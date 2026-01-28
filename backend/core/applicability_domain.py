"""
Applicability Domain Checking (Phase 4)

This module determines whether predictions are within the reliable domain
of each predictor. It detects when models are extrapolating and adjusts
confidence accordingly.

Design Philosophy:
- Know the limits of our knowledge
- Detect when predictions are unreliable
- Communicate uncertainty honestly
- Prevent false confidence in extrapolations
"""

import logging
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field
from enum import Enum

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, DataStructs, rdMolDescriptors
import numpy as np

logger = logging.getLogger(__name__)


class DomainStatus(Enum):
    """Status of a molecule relative to applicability domain."""
    INSIDE = "inside"          # Well within domain
    BORDERLINE = "borderline"  # Near the edge of domain
    OUTSIDE = "outside"        # Outside domain
    UNKNOWN = "unknown"        # Cannot determine


@dataclass
class DomainCheckResult:
    """Result of an applicability domain check."""
    status: DomainStatus
    confidence_modifier: float  # Multiply confidence by this (0.0 - 1.0)
    reasons: List[str] = field(default_factory=list)
    details: Dict[str, Any] = field(default_factory=dict)

    @property
    def is_reliable(self) -> bool:
        """Check if prediction is reliable."""
        return self.status == DomainStatus.INSIDE

    @property
    def needs_caution(self) -> bool:
        """Check if prediction needs extra caution."""
        return self.status in (DomainStatus.BORDERLINE, DomainStatus.OUTSIDE)


class ApplicabilityDomainChecker:
    """
    Checks if molecules are within the applicability domain of predictors.

    This is crucial for honest predictions. A model trained on drug-like
    molecules shouldn't confidently predict properties of molecules
    outside that space.

    Methods:
    1. Property-based: Check if properties are in expected ranges
    2. Similarity-based: Check similarity to reference compounds
    3. Structural-based: Check for unusual substructures
    """

    # Standard drug-like property ranges (from large drug databases)
    DRUG_LIKE_RANGES = {
        "molecular_weight": (100, 900),
        "logp": (-3, 8),
        "tpsa": (0, 200),
        "hbd": (0, 8),
        "hba": (0, 15),
        "rotatable_bonds": (0, 15),
        "heavy_atoms": (5, 70),
        "aromatic_rings": (0, 5),
        "fraction_sp3": (0, 1),
    }

    # Stricter ranges for oral drugs
    ORAL_DRUG_RANGES = {
        "molecular_weight": (150, 600),
        "logp": (-1, 6),
        "tpsa": (0, 150),
        "hbd": (0, 5),
        "hba": (0, 10),
        "rotatable_bonds": (0, 10),
    }

    # Unusual elements that suggest non-drug-like chemistry
    UNUSUAL_ELEMENTS = {
        "Si", "B", "As", "Se", "Te", "metal"  # Metals are special case
    }

    def __init__(self):
        self._reference_fps = {}  # Cache for reference fingerprints

    def check_domain(
        self,
        smiles: str,
        domain_type: str = "drug_like",
        reference_smiles: Optional[List[str]] = None
    ) -> DomainCheckResult:
        """
        Check if a molecule is within the applicability domain.

        Args:
            smiles: SMILES string to check
            domain_type: "drug_like", "oral_drug", or "custom"
            reference_smiles: Optional reference compounds for similarity check

        Returns:
            DomainCheckResult with status and details
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return DomainCheckResult(
                status=DomainStatus.UNKNOWN,
                confidence_modifier=0.0,
                reasons=["Invalid SMILES - cannot parse molecule"],
            )

        # Run all checks
        property_result = self._check_properties(mol, domain_type)
        structural_result = self._check_structure(mol)
        similarity_result = self._check_similarity(mol, reference_smiles) if reference_smiles else None

        # Combine results
        return self._combine_results(property_result, structural_result, similarity_result)

    def _check_properties(self, mol: Chem.Mol, domain_type: str) -> DomainCheckResult:
        """Check if molecular properties are in expected ranges."""
        if domain_type == "oral_drug":
            ranges = self.ORAL_DRUG_RANGES
        else:
            ranges = self.DRUG_LIKE_RANGES

        violations = []
        borderline = []
        details = {}

        # Calculate properties
        props = {
            "molecular_weight": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "tpsa": Descriptors.TPSA(mol),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
            "heavy_atoms": mol.GetNumHeavyAtoms(),
            "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        }

        details["properties"] = props

        # Check each property
        for prop, (min_val, max_val) in ranges.items():
            if prop not in props:
                continue

            value = props[prop]

            # Calculate how far outside the range
            if value < min_val:
                deviation = (min_val - value) / (max_val - min_val) if max_val != min_val else 1
                if deviation > 0.2:
                    violations.append(f"{prop} ({value:.1f}) below range ({min_val}-{max_val})")
                else:
                    borderline.append(f"{prop} ({value:.1f}) near lower limit")
            elif value > max_val:
                deviation = (value - max_val) / (max_val - min_val) if max_val != min_val else 1
                if deviation > 0.2:
                    violations.append(f"{prop} ({value:.1f}) above range ({min_val}-{max_val})")
                else:
                    borderline.append(f"{prop} ({value:.1f}) near upper limit")

        # Determine status
        if violations:
            return DomainCheckResult(
                status=DomainStatus.OUTSIDE,
                confidence_modifier=0.5,
                reasons=violations,
                details=details,
            )
        elif borderline:
            return DomainCheckResult(
                status=DomainStatus.BORDERLINE,
                confidence_modifier=0.8,
                reasons=borderline,
                details=details,
            )
        else:
            return DomainCheckResult(
                status=DomainStatus.INSIDE,
                confidence_modifier=1.0,
                reasons=[],
                details=details,
            )

    def _check_structure(self, mol: Chem.Mol) -> DomainCheckResult:
        """Check for unusual structural features."""
        issues = []
        details = {}

        # Check for unusual elements
        elements = set()
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            elements.add(symbol)

            # Check for metals
            if atom.GetAtomicNum() in [
                3, 4, 11, 12, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
                37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
                55, 56, 72, 73, 74, 75, 76, 77, 78, 79, 80,
            ]:
                issues.append(f"Contains metal atom: {symbol}")

        details["elements"] = list(elements)

        unusual = elements & self.UNUSUAL_ELEMENTS
        if unusual:
            issues.extend([f"Unusual element: {e}" for e in unusual])

        # Check for very large ring systems
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if len(ring) > 8:
                issues.append(f"Large ring system ({len(ring)} atoms)")

        # Check for many disconnected fragments
        frags = Chem.GetMolFrags(mol)
        if len(frags) > 1:
            issues.append(f"Multiple disconnected fragments ({len(frags)})")

        # Check for unusual formal charges
        total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        if abs(total_charge) > 2:
            issues.append(f"High formal charge: {total_charge:+d}")

        # Determine status
        if len(issues) >= 2:
            return DomainCheckResult(
                status=DomainStatus.OUTSIDE,
                confidence_modifier=0.4,
                reasons=issues,
                details=details,
            )
        elif issues:
            return DomainCheckResult(
                status=DomainStatus.BORDERLINE,
                confidence_modifier=0.7,
                reasons=issues,
                details=details,
            )
        else:
            return DomainCheckResult(
                status=DomainStatus.INSIDE,
                confidence_modifier=1.0,
                reasons=[],
                details=details,
            )

    def _check_similarity(
        self,
        mol: Chem.Mol,
        reference_smiles: List[str]
    ) -> DomainCheckResult:
        """Check similarity to reference compounds."""
        if not reference_smiles:
            return DomainCheckResult(
                status=DomainStatus.UNKNOWN,
                confidence_modifier=1.0,
                reasons=["No reference compounds provided"],
            )

        # Calculate fingerprint for query
        query_fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)

        # Calculate similarities
        similarities = []
        for ref_smiles in reference_smiles:
            # Use cached fingerprint if available
            if ref_smiles not in self._reference_fps:
                ref_mol = Chem.MolFromSmiles(ref_smiles)
                if ref_mol:
                    self._reference_fps[ref_smiles] = AllChem.GetMorganFingerprintAsBitVect(
                        ref_mol, 2, nBits=2048
                    )

            if ref_smiles in self._reference_fps:
                sim = DataStructs.TanimotoSimilarity(query_fp, self._reference_fps[ref_smiles])
                similarities.append(sim)

        if not similarities:
            return DomainCheckResult(
                status=DomainStatus.UNKNOWN,
                confidence_modifier=0.8,
                reasons=["Could not compare to reference compounds"],
            )

        max_sim = max(similarities)
        avg_sim = sum(similarities) / len(similarities)

        details = {
            "max_similarity": round(max_sim, 3),
            "avg_similarity": round(avg_sim, 3),
            "num_references": len(similarities),
        }

        # Determine status based on similarity
        if max_sim >= 0.5:
            return DomainCheckResult(
                status=DomainStatus.INSIDE,
                confidence_modifier=1.0,
                reasons=[f"Similar to known compounds (max Tanimoto: {max_sim:.2f})"],
                details=details,
            )
        elif max_sim >= 0.3:
            return DomainCheckResult(
                status=DomainStatus.BORDERLINE,
                confidence_modifier=0.7,
                reasons=[f"Moderate similarity to references (max: {max_sim:.2f})"],
                details=details,
            )
        else:
            return DomainCheckResult(
                status=DomainStatus.OUTSIDE,
                confidence_modifier=0.4,
                reasons=[f"Low similarity to known compounds (max: {max_sim:.2f})"],
                details=details,
            )

    def _combine_results(
        self,
        property_result: DomainCheckResult,
        structural_result: DomainCheckResult,
        similarity_result: Optional[DomainCheckResult]
    ) -> DomainCheckResult:
        """Combine multiple domain check results."""
        all_results = [property_result, structural_result]
        if similarity_result:
            all_results.append(similarity_result)

        # Take most severe status
        status_priority = {
            DomainStatus.OUTSIDE: 0,
            DomainStatus.BORDERLINE: 1,
            DomainStatus.UNKNOWN: 2,
            DomainStatus.INSIDE: 3,
        }

        worst_status = min(all_results, key=lambda r: status_priority[r.status]).status

        # Combine confidence modifiers (multiply)
        combined_modifier = 1.0
        for result in all_results:
            combined_modifier *= result.confidence_modifier

        # Collect all reasons
        all_reasons = []
        for result in all_results:
            all_reasons.extend(result.reasons)

        # Combine details
        combined_details = {}
        for result in all_results:
            combined_details.update(result.details)

        return DomainCheckResult(
            status=worst_status,
            confidence_modifier=combined_modifier,
            reasons=all_reasons,
            details=combined_details,
        )

    def adjust_prediction_confidence(
        self,
        smiles: str,
        original_confidence_score: float,
        reference_smiles: Optional[List[str]] = None
    ) -> Tuple[float, List[str]]:
        """
        Adjust prediction confidence based on applicability domain.

        Args:
            smiles: SMILES of molecule
            original_confidence_score: Original confidence (0-1)
            reference_smiles: Optional reference compounds

        Returns:
            Tuple of (adjusted_confidence, list of reasons for adjustment)
        """
        domain_check = self.check_domain(smiles, reference_smiles=reference_smiles)

        adjusted_confidence = original_confidence_score * domain_check.confidence_modifier
        adjusted_confidence = max(0.1, min(1.0, adjusted_confidence))  # Keep in bounds

        reasons = []
        if domain_check.status == DomainStatus.OUTSIDE:
            reasons.append("Confidence reduced - molecule outside typical drug-like space")
            reasons.extend(domain_check.reasons)
        elif domain_check.status == DomainStatus.BORDERLINE:
            reasons.append("Confidence slightly reduced - molecule near domain boundary")

        return adjusted_confidence, reasons


# Global instance
applicability_checker = ApplicabilityDomainChecker()
