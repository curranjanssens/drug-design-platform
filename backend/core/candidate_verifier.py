"""
Candidate Verification Layer

The key insight: LLMs find generation easy but rigorous verification hard.
This module provides the "senior medicinal chemist intuition" that catches:
- SMILES/name mismatches
- Pharmacophore inconsistencies
- Known drug duplicates
- Structural alert violations
- Claims that don't match the actual structure

This is where the ROI is - a strong verifier makes a decent generator excellent.
"""

import logging
import re
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field
from enum import Enum

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import DataStructs

logger = logging.getLogger(__name__)


class VerificationSeverity(Enum):
    """Severity levels for verification issues."""
    CRITICAL = "critical"      # Reject candidate
    WARNING = "warning"        # Flag but allow
    INFO = "info"              # Note for transparency


@dataclass
class VerificationIssue:
    """A single verification issue found."""
    category: str              # e.g., "structure_mismatch", "known_drug", "pharmacophore"
    severity: VerificationSeverity
    message: str
    details: Dict = field(default_factory=dict)


@dataclass
class VerificationResult:
    """Complete verification result for a candidate."""
    smiles: str
    name: str
    is_valid: bool
    issues: List[VerificationIssue] = field(default_factory=list)
    verified_properties: Dict = field(default_factory=dict)

    @property
    def critical_issues(self) -> List[VerificationIssue]:
        return [i for i in self.issues if i.severity == VerificationSeverity.CRITICAL]

    @property
    def warnings(self) -> List[VerificationIssue]:
        return [i for i in self.issues if i.severity == VerificationSeverity.WARNING]


# ============================================================================
# KNOWN COMPOUNDS DATABASE
# ============================================================================

# Known recreational drugs (high similarity = reject)
KNOWN_RECREATIONAL_DRUGS = [
    # Phenethylamines/Amphetamines
    ("CC(N)Cc1ccccc1", "Amphetamine", "stimulant"),
    ("CC(NC)Cc1ccccc1", "Methamphetamine", "stimulant"),
    ("CNC(C)Cc1ccc2c(c1)OCO2", "MDMA", "entactogen"),
    ("CC(N)Cc1ccc2c(c1)OCO2", "MDA", "entactogen"),
    ("CNC(C)Cc1ccc(OC)c(OC)c1", "MMDA", "entactogen"),
    # Benzofurans - commonly missed!
    ("CC(N)Cc1ccc2occc2c1", "6-APB", "entactogen"),
    ("CNC(C)Cc1ccc2occc2c1", "6-MAPB", "entactogen"),
    ("CC(N)Cc1ccc2ccoc2c1", "5-APB", "entactogen"),
    ("CNC(C)Cc1ccc2ccoc2c1", "5-MAPB", "entactogen"),
    ("CC(N)Cc1cccc2occc12", "4-APB", "entactogen"),
    # Cathinones
    ("CNC(C)C(=O)c1ccc(C)cc1", "Mephedrone", "stimulant"),
    ("CCNC(C)C(=O)c1ccc(C)cc1", "4-MEC", "stimulant"),
    ("CNC(C)C(=O)c1ccc2c(c1)OCO2", "Methylone", "stimulant"),
    ("CNC(C)C(=O)c1ccccc1", "Methcathinone", "stimulant"),
    ("CC(NC(C)C)C(=O)c1ccc2c(c1)OCO2", "MDPV", "stimulant"),
    # Tryptamines
    ("CN(C)CCc1c[nH]c2ccccc12", "DMT", "psychedelic"),
    ("CCN(CC)CCc1c[nH]c2ccccc12", "DET", "psychedelic"),
    ("CCCN(CCC)CCc1c[nH]c2ccccc12", "DPT", "psychedelic"),
    ("COc1ccc2[nH]cc(CCN(C)C)c2c1", "5-MeO-DMT", "psychedelic"),
    ("CCCc1c[nH]c2cccc(O)c12", "Psilocin", "psychedelic"),
    # 2C series
    ("COc1cc(CCN)c(OC)cc1Br", "2C-B", "psychedelic"),
    ("COc1cc(CCN)c(OC)cc1I", "2C-I", "psychedelic"),
    ("COc1cc(CCN)c(OC)cc1C", "2C-D", "psychedelic"),
    # NBOMes
    ("COc1ccc(CCNCc2ccccc2OC)c(OC)c1Br", "25B-NBOMe", "psychedelic"),
    ("COc1ccc(CCNCc2ccccc2OC)c(OC)c1I", "25I-NBOMe", "psychedelic"),
    # Opioids
    ("CCC(=O)N(c1ccccc1)C1CCN(CCc2ccccc2)CC1", "Fentanyl", "opioid"),
    ("CCC(=O)N(c1ccccc1)C1CCN(CC(C)c2ccccc2)CC1", "α-Methylfentanyl", "opioid"),
    # GHB
    ("OCC1CCCO1", "GBL", "depressant"),
]

# Known approved drugs (flag similarity but context matters)
KNOWN_APPROVED_DRUGS = [
    ("CN1C(=O)N(C)c2nc[nH]c2C1=O", "Theophylline", "bronchodilator"),
    ("Cn1cnc2c1c(=O)n(C)c(=O)n2C", "Caffeine", "stimulant"),
    # Add more as needed
]


# ============================================================================
# NAME-STRUCTURE CONSISTENCY PATTERNS
# ============================================================================

# Substructure indicators that should be present if name contains certain terms
NAME_STRUCTURE_PATTERNS = {
    # Heterocycles
    "benzofuran": ("c1ccc2occc2c1", "Benzofuran ring"),
    "benzothiophene": ("c1ccc2sccc2c1", "Benzothiophene ring"),
    "indole": ("c1ccc2[nH]ccc2c1", "Indole ring"),
    "indazole": ("c1ccc2[nH]ncc2c1", "Indazole ring"),
    "quinoline": ("c1ccc2ncccc2c1", "Quinoline ring"),
    "isoquinoline": ("c1ccc2cnccc2c1", "Isoquinoline ring"),
    "pyridine": ("c1ccncc1", "Pyridine ring"),
    "pyrimidine": ("c1cncnc1", "Pyrimidine ring"),
    "pyrazine": ("c1cnccn1", "Pyrazine ring"),
    "thiophene": ("c1ccsc1", "Thiophene ring"),
    "furan": ("c1ccoc1", "Furan ring"),
    "pyrrole": ("c1cc[nH]c1", "Pyrrole ring"),
    "imidazole": ("c1c[nH]cn1", "Imidazole ring"),
    "oxazole": ("c1cocn1", "Oxazole ring"),
    "thiazole": ("c1cscn1", "Thiazole ring"),
    "triazole": ("c1nc[nH]n1", "Triazole ring"),
    "tetrazole": ("c1nn[nH]n1", "Tetrazole ring"),

    # Functional groups
    "cyano": ("C#N", "Cyano/nitrile group"),
    "nitrile": ("C#N", "Cyano/nitrile group"),
    "methoxy": ("COc", "Methoxy group"),
    "fluoro": ("F", "Fluorine"),
    "chloro": ("Cl", "Chlorine"),
    "bromo": ("Br", "Bromine"),
    "trifluoromethyl": ("C(F)(F)F", "Trifluoromethyl group"),
    "cf3": ("C(F)(F)F", "Trifluoromethyl group"),
    "hydroxyl": ("O", "Hydroxyl group"),  # Too generic - use carefully
    "amino": ("N", "Amino group"),  # Too generic - use carefully
    "nitro": ("[N+](=O)[O-]", "Nitro group"),
    "sulfonyl": ("S(=O)(=O)", "Sulfonyl group"),
    "carbonyl": ("C=O", "Carbonyl group"),
    "carboxyl": ("C(=O)O", "Carboxyl group"),
    "amide": ("C(=O)N", "Amide group"),
    "ester": ("C(=O)O", "Ester group"),
    "ether": ("COC", "Ether linkage"),

    # Scaffolds
    "phenyl": ("c1ccccc1", "Phenyl ring"),
    "benzyl": ("Cc1ccccc1", "Benzyl group"),
    "naphthyl": ("c1ccc2ccccc2c1", "Naphthalene"),
    "cyclopropyl": ("C1CC1", "Cyclopropyl ring"),
    "cyclobutyl": ("C1CCC1", "Cyclobutyl ring"),
    "cyclopentyl": ("C1CCCC1", "Cyclopentyl ring"),
    "cyclohexyl": ("C1CCCCC1", "Cyclohexyl ring"),
    "piperidine": ("C1CCNCC1", "Piperidine ring"),
    "piperazine": ("C1CNCCN1", "Piperazine ring"),
    "morpholine": ("C1COCCN1", "Morpholine ring"),
    "pyrrolidine": ("C1CCNC1", "Pyrrolidine ring"),

    # Cathinone/amphetamine features
    "cathinone": ("C(=O)c1ccccc1", "Cathinone core (beta-keto)"),
    "amphetamine": ("CC(N)Cc1ccccc1", "Amphetamine core"),
    "methamphetamine": ("CC(NC)Cc1ccccc1", "Methamphetamine core"),
}


# ============================================================================
# PHARMACOPHORE PATTERNS
# ============================================================================

# Common pharmacophore requirements by target class
PHARMACOPHORE_REQUIREMENTS = {
    "sert": {
        "required": [
            ("N", "Basic nitrogen for SERT binding"),
        ],
        "common": [
            ("c1ccccc1", "Aromatic ring for hydrophobic binding"),
        ]
    },
    "dat": {
        "required": [
            ("N", "Basic nitrogen for DAT binding"),
        ],
        "common": [
            ("c1ccccc1", "Aromatic ring"),
        ]
    },
    "5ht": {
        "required": [
            ("N", "Basic nitrogen for 5-HT binding"),
        ],
        "common": [
            ("c1c[nH]c2ccccc12", "Indole-like system"),
        ]
    },
    "kinase": {
        "required": [
            ("c1cncnc1", "Hinge-binding heterocycle (pyrimidine/purine-like)"),
        ],
        "common": [
            ("NC(=O)", "H-bond donor/acceptor"),
        ]
    },
    "gpcr": {
        "required": [
            ("N", "Basic nitrogen"),
        ],
        "common": [
            ("c1ccccc1", "Aromatic system"),
        ]
    },
}


class CandidateVerifier:
    """
    Comprehensive verification layer for drug candidates.

    Catches issues that LLM generation commonly produces:
    1. SMILES/name mismatches
    2. Known drug duplicates
    3. Pharmacophore inconsistencies
    4. Structural impossibilities
    5. Claims contradicting structure
    """

    def __init__(self):
        # Pre-compute fingerprints for known drugs
        self._recreational_fps = []
        for smiles, name, category in KNOWN_RECREATIONAL_DRUGS:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)
                self._recreational_fps.append((fp, name, category, smiles))

        logger.info(f"CandidateVerifier initialized with {len(self._recreational_fps)} known recreational drugs")

    def verify(
        self,
        smiles: str,
        name: str,
        rationale: str = "",
        target_class: str = "",
        claimed_features: List[str] = None
    ) -> VerificationResult:
        """
        Comprehensive verification of a drug candidate.

        Args:
            smiles: SMILES string
            name: Candidate name
            rationale: Design rationale/description
            target_class: Target class (e.g., "sert", "kinase")
            claimed_features: Features claimed in the design

        Returns:
            VerificationResult with all issues found
        """
        result = VerificationResult(smiles=smiles, name=name, is_valid=True, issues=[])

        # Step 1: Basic SMILES validity
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            result.issues.append(VerificationIssue(
                category="invalid_smiles",
                severity=VerificationSeverity.CRITICAL,
                message=f"Invalid SMILES: {smiles}"
            ))
            result.is_valid = False
            return result

        # Compute basic properties
        result.verified_properties = self._compute_properties(mol)

        # Step 2: Check against known recreational drugs
        rec_issues = self._check_known_recreational(mol, smiles)
        result.issues.extend(rec_issues)

        # Step 3: Name-structure consistency
        name_issues = self._check_name_structure_consistency(mol, name)
        result.issues.extend(name_issues)

        # Step 4: Rationale-structure consistency
        if rationale:
            rationale_issues = self._check_rationale_consistency(mol, rationale, name)
            result.issues.extend(rationale_issues)

        # Step 5: Pharmacophore requirements
        if target_class:
            pharma_issues = self._check_pharmacophore(mol, target_class)
            result.issues.extend(pharma_issues)

        # Step 6: Claimed features verification
        if claimed_features:
            claim_issues = self._verify_claimed_features(mol, claimed_features)
            result.issues.extend(claim_issues)

        # Step 7: Structural alerts (stability)
        stability_issues = self._check_stability_alerts(mol)
        result.issues.extend(stability_issues)

        # Determine overall validity
        if result.critical_issues:
            result.is_valid = False

        return result

    def _compute_properties(self, mol: Chem.Mol) -> Dict:
        """Compute verified molecular properties."""
        return {
            "mw": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "hbd": rdMolDescriptors.CalcNumHBD(mol),
            "hba": rdMolDescriptors.CalcNumHBA(mol),
            "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
            "tpsa": rdMolDescriptors.CalcTPSA(mol),
            "rings": rdMolDescriptors.CalcNumRings(mol),
            "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
            "heavy_atoms": mol.GetNumHeavyAtoms(),
            "has_nitrogen": mol.HasSubstructMatch(Chem.MolFromSmarts("[N]")),
            "has_basic_nitrogen": mol.HasSubstructMatch(Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(NC=O)]")),
        }

    def _check_known_recreational(self, mol: Chem.Mol, smiles: str) -> List[VerificationIssue]:
        """Check if candidate is too similar to known recreational drugs."""
        issues = []

        mol_fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)

        for ref_fp, ref_name, category, ref_smiles in self._recreational_fps:
            similarity = DataStructs.TanimotoSimilarity(mol_fp, ref_fp)

            if similarity >= 0.95:
                # Very high similarity - essentially the same compound
                issues.append(VerificationIssue(
                    category="known_drug_duplicate",
                    severity=VerificationSeverity.CRITICAL,
                    message=f"Structure is {similarity:.0%} similar to {ref_name} - essentially identical",
                    details={
                        "known_drug": ref_name,
                        "category": category,
                        "similarity": similarity,
                        "known_smiles": ref_smiles
                    }
                ))
            elif similarity >= 0.85:
                # High similarity - very close analog
                issues.append(VerificationIssue(
                    category="known_drug_analog",
                    severity=VerificationSeverity.CRITICAL,
                    message=f"Structure is {similarity:.0%} similar to {ref_name} ({category}) - too close for novelty",
                    details={
                        "known_drug": ref_name,
                        "category": category,
                        "similarity": similarity,
                        "known_smiles": ref_smiles
                    }
                ))
            elif similarity >= 0.70:
                # Moderate similarity - flag as warning
                issues.append(VerificationIssue(
                    category="known_drug_related",
                    severity=VerificationSeverity.WARNING,
                    message=f"Structure has {similarity:.0%} similarity to {ref_name} ({category})",
                    details={
                        "known_drug": ref_name,
                        "category": category,
                        "similarity": similarity
                    }
                ))

        return issues

    def _check_name_structure_consistency(self, mol: Chem.Mol, name: str) -> List[VerificationIssue]:
        """Verify that the name is consistent with the structure."""
        issues = []
        name_lower = name.lower()

        for term, (smarts, description) in NAME_STRUCTURE_PATTERNS.items():
            if term in name_lower:
                pattern = Chem.MolFromSmarts(smarts)
                if pattern and not mol.HasSubstructMatch(pattern):
                    issues.append(VerificationIssue(
                        category="name_structure_mismatch",
                        severity=VerificationSeverity.CRITICAL,
                        message=f"Name contains '{term}' but structure lacks {description}",
                        details={
                            "claimed_feature": term,
                            "expected_smarts": smarts,
                            "description": description
                        }
                    ))

        return issues

    def _check_rationale_consistency(self, mol: Chem.Mol, rationale: str, name: str) -> List[VerificationIssue]:
        """Verify that claims in the rationale match the actual structure."""
        issues = []
        rationale_lower = rationale.lower()

        # Check for specific claims about functional groups
        claim_patterns = [
            ("cyano", "C#N", "cyano/nitrile group"),
            ("nitrile", "C#N", "cyano/nitrile group"),
            ("fluorine", "F", "fluorine"),
            ("fluoro", "F", "fluorine"),
            ("chlorine", "Cl", "chlorine"),
            ("chloro", "Cl", "chlorine"),
            ("methoxy", "COc", "methoxy group"),
            ("trifluoromethyl", "C(F)(F)F", "CF3 group"),
            ("-cf3", "C(F)(F)F", "CF3 group"),
            ("hydroxyl", "[OH]", "hydroxyl group"),
            ("carboxylic", "C(=O)O", "carboxylic acid"),
            ("amide", "C(=O)N", "amide"),
            ("sulfonamide", "S(=O)(=O)N", "sulfonamide"),
        ]

        for term, smarts, description in claim_patterns:
            if term in rationale_lower:
                pattern = Chem.MolFromSmarts(smarts)
                if pattern and not mol.HasSubstructMatch(pattern):
                    issues.append(VerificationIssue(
                        category="rationale_structure_mismatch",
                        severity=VerificationSeverity.WARNING,
                        message=f"Rationale mentions '{term}' but structure lacks {description}",
                        details={
                            "claimed": term,
                            "expected_smarts": smarts
                        }
                    ))

        return issues

    def _check_pharmacophore(self, mol: Chem.Mol, target_class: str) -> List[VerificationIssue]:
        """Check if structure has required pharmacophore features for target class."""
        issues = []
        target_lower = target_class.lower()

        # Find matching pharmacophore requirements
        for target_key, requirements in PHARMACOPHORE_REQUIREMENTS.items():
            if target_key in target_lower:
                # Check required features
                for smarts, description in requirements.get("required", []):
                    pattern = Chem.MolFromSmarts(smarts)
                    if pattern and not mol.HasSubstructMatch(pattern):
                        issues.append(VerificationIssue(
                            category="missing_pharmacophore",
                            severity=VerificationSeverity.WARNING,
                            message=f"Missing required pharmacophore for {target_key}: {description}",
                            details={
                                "target": target_key,
                                "required_feature": description,
                                "smarts": smarts
                            }
                        ))
                break

        return issues

    def _verify_claimed_features(self, mol: Chem.Mol, claimed_features: List[str]) -> List[VerificationIssue]:
        """Verify that claimed features are actually present."""
        issues = []

        for feature in claimed_features:
            feature_lower = feature.lower()

            # Map feature claims to SMARTS
            for term, (smarts, description) in NAME_STRUCTURE_PATTERNS.items():
                if term in feature_lower:
                    pattern = Chem.MolFromSmarts(smarts)
                    if pattern and not mol.HasSubstructMatch(pattern):
                        issues.append(VerificationIssue(
                            category="missing_claimed_feature",
                            severity=VerificationSeverity.WARNING,
                            message=f"Claimed feature '{feature}' not found in structure",
                            details={
                                "claimed": feature,
                                "expected": description
                            }
                        ))
                    break

        return issues

    def _check_stability_alerts(self, mol: Chem.Mol) -> List[VerificationIssue]:
        """Check for stability issues (reactive/unstable groups)."""
        issues = []

        stability_patterns = [
            ("[SH]", "Free thiol - oxidizes in plasma, binds albumin", VerificationSeverity.CRITICAL),
            ("[NX3]C(F)(F)F", "N-CF3 - chemically unstable", VerificationSeverity.CRITICAL),
            ("OO", "Peroxide - unstable/explosive", VerificationSeverity.CRITICAL),
            ("C(=O)Cl", "Acid chloride - too reactive", VerificationSeverity.CRITICAL),
            ("C(=O)F", "Acyl fluoride - too reactive", VerificationSeverity.CRITICAL),
            ("N=NN", "Triazene - photolabile", VerificationSeverity.CRITICAL),
            ("[CH]=O", "Aldehyde - reactive", VerificationSeverity.WARNING),
            ("C1OC1", "Epoxide - reactive/mutagenic risk", VerificationSeverity.WARNING),
            ("SS", "Disulfide - may scramble", VerificationSeverity.WARNING),
            ("c[NH2]", "Primary aromatic amine - metabolic liability", VerificationSeverity.WARNING),
            ("c[N+](=O)[O-]", "Nitroaromatic - mutagenicity concern", VerificationSeverity.WARNING),
        ]

        for smarts, message, severity in stability_patterns:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                issues.append(VerificationIssue(
                    category="stability_alert",
                    severity=severity,
                    message=message,
                    details={"pattern": smarts}
                ))

        return issues

    def verify_batch(
        self,
        candidates: List[Dict],
        target_class: str = ""
    ) -> List[VerificationResult]:
        """
        Verify a batch of candidates.

        Args:
            candidates: List of dicts with 'smiles', 'name', optionally 'rationale'
            target_class: Target class for pharmacophore checking

        Returns:
            List of VerificationResults
        """
        results = []

        for candidate in candidates:
            result = self.verify(
                smiles=candidate.get("smiles", ""),
                name=candidate.get("name", ""),
                rationale=candidate.get("rationale", ""),
                target_class=target_class,
                claimed_features=candidate.get("features", [])
            )
            results.append(result)

        return results

    def filter_valid(
        self,
        candidates: List[Dict],
        target_class: str = ""
    ) -> Tuple[List[Dict], List[Dict]]:
        """
        Filter candidates into valid and rejected.

        Returns:
            (valid_candidates, rejected_candidates)
        """
        results = self.verify_batch(candidates, target_class)

        valid = []
        rejected = []

        for candidate, result in zip(candidates, results):
            if result.is_valid:
                # Add verification info
                candidate["verification"] = {
                    "verified": True,
                    "warnings": [i.message for i in result.warnings],
                    "properties": result.verified_properties
                }
                valid.append(candidate)
            else:
                candidate["verification"] = {
                    "verified": False,
                    "rejection_reasons": [i.message for i in result.critical_issues],
                    "warnings": [i.message for i in result.warnings]
                }
                rejected.append(candidate)

        return valid, rejected


# Global instance
candidate_verifier = CandidateVerifier()
