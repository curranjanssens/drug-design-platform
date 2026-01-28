"""
Property Prediction Service using RDKit and additional models
"""
from typing import Dict, Any, Optional, List
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, Lipinski, QED
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
import numpy as np
from loguru import logger

from models import PropertyPrediction


class PropertyPredictor:
    """Service for predicting molecular properties"""

    def __init__(self):
        # Initialize PAINS filter
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        self.pains_catalog = FilterCatalog(params)

    def predict_properties(self, smiles: str) -> Optional[PropertyPrediction]:
        """Predict all relevant properties for a compound"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.error(f"Could not parse SMILES: {smiles}")
                return None

            # Basic descriptors
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            tpsa = Descriptors.TPSA(mol)
            rotatable = Descriptors.NumRotatableBonds(mol)

            # QED (Quantitative Estimate of Drug-likeness)
            qed_score = QED.qed(mol)

            # Synthetic accessibility
            sa_score = self._calculate_sa_score(mol)

            # Lipinski violations
            lipinski_violations = self._count_lipinski_violations(mol)

            # Additional predictions
            solubility = self._predict_solubility(mol)
            half_life = self._predict_peptide_half_life(mol)

            return PropertyPrediction(
                molecular_weight=round(mw, 2),
                logp=round(logp, 2),
                hbd=hbd,
                hba=hba,
                tpsa=round(tpsa, 2),
                rotatable_bonds=rotatable,
                qed=round(qed_score, 3),
                synthetic_accessibility=round(sa_score, 2),
                lipinski_violations=lipinski_violations,
                predicted_solubility=solubility,
                predicted_half_life=half_life
            )

        except Exception as e:
            logger.error(f"Error predicting properties: {e}")
            return None

    def _calculate_sa_score(self, mol: Chem.Mol) -> float:
        """
        Calculate Synthetic Accessibility Score (1-10)
        Based on Ertl & Schuffenhauer method
        """
        try:
            from rdkit.Chem import RDConfig
            import os
            import sys

            # Simplified SA score calculation
            # Full implementation would use the SA_Score module

            # Heuristic based on molecular complexity
            num_rings = rdMolDescriptors.CalcNumRings(mol)
            num_stereo = len(Chem.FindMolChiralCenters(mol))
            num_atoms = mol.GetNumAtoms()
            num_rotatable = Descriptors.NumRotatableBonds(mol)

            # Peptide-specific adjustments
            num_amide = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[NX3][CX3](=[OX1])')))

            # Base score
            score = 3.0

            # Ring complexity
            score += num_rings * 0.3

            # Stereocenters
            score += num_stereo * 0.2

            # Size penalty
            if num_atoms > 50:
                score += (num_atoms - 50) * 0.05

            # Peptide bonus (peptides are generally synthesizable)
            if num_amide >= 3:
                score -= 1.0

            return max(1.0, min(10.0, score))

        except Exception as e:
            logger.warning(f"SA score calculation error: {e}")
            return 5.0  # Default moderate difficulty

    def _count_lipinski_violations(self, mol: Chem.Mol) -> int:
        """Count Lipinski Rule of Five violations"""
        violations = 0

        if Descriptors.MolWt(mol) > 500:
            violations += 1
        if Descriptors.MolLogP(mol) > 5:
            violations += 1
        if Descriptors.NumHDonors(mol) > 5:
            violations += 1
        if Descriptors.NumHAcceptors(mol) > 10:
            violations += 1

        return violations

    def _predict_solubility(self, mol: Chem.Mol) -> str:
        """Predict aqueous solubility category"""
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        mw = Descriptors.MolWt(mol)

        # Simple heuristic model
        if logp < 0 and tpsa > 100:
            return "High"
        elif logp < 2 and tpsa > 60:
            return "Moderate"
        elif logp < 4:
            return "Low"
        else:
            return "Very Low"

    def _predict_peptide_half_life(self, mol: Chem.Mol) -> str:
        """Predict relative plasma half-life for peptide-like molecules"""
        # Count protective modifications
        smiles = Chem.MolToSmiles(mol)

        # Check for N-terminal protection
        has_n_protection = any([
            "CC(C)(C)C(=O)N" in smiles,  # pivaloyl
            "C1CC1C(=O)N" in smiles,  # cyclopropyl
        ])

        # Check for D-amino acids (simplified check)
        has_d_amino = "[C@H]" in smiles or "[C@@H]" in smiles

        # Check for N-methylation
        has_n_methyl = "N(C)C" in smiles or "N(C)CC" in smiles

        # Count amide bonds
        num_amides = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[NX3][CX3](=[OX1])')))

        # Estimate half-life category
        protection_score = sum([
            has_n_protection * 2,
            has_d_amino * 1,
            has_n_methyl * 1,
        ])

        if num_amides < 3:
            return "Short (minutes)"
        elif protection_score >= 2:
            return "Extended (hours)"
        elif protection_score >= 1:
            return "Moderate (30-60 min)"
        else:
            return "Short (minutes)"

    def check_pains(self, mol: Chem.Mol) -> List[str]:
        """Check for PAINS (Pan Assay Interference Compounds) alerts"""
        alerts = []
        entries = self.pains_catalog.GetMatches(mol)
        for entry in entries:
            alerts.append(entry.GetDescription())
        return alerts

    def calculate_druglikeness_profile(self, smiles: str) -> Dict[str, Any]:
        """Calculate comprehensive drug-likeness profile"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {}

        profile = {
            "lipinski": self._check_lipinski(mol),
            "veber": self._check_veber(mol),
            "ghose": self._check_ghose(mol),
            "egan": self._check_egan(mol),
            "muegge": self._check_muegge(mol),
            "pains_alerts": self.check_pains(mol),
            "brenk_alerts": self._check_brenk(mol),
        }

        # Overall drug-likeness
        passes = sum([
            profile["lipinski"]["pass"],
            profile["veber"]["pass"],
            profile["ghose"]["pass"],
        ])
        profile["overall_score"] = passes / 3.0

        return profile

    def _check_lipinski(self, mol: Chem.Mol) -> Dict[str, Any]:
        """Lipinski Rule of Five"""
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)

        violations = []
        if mw > 500:
            violations.append(f"MW={mw:.1f} > 500")
        if logp > 5:
            violations.append(f"LogP={logp:.2f} > 5")
        if hbd > 5:
            violations.append(f"HBD={hbd} > 5")
        if hba > 10:
            violations.append(f"HBA={hba} > 10")

        return {
            "pass": len(violations) <= 1,
            "violations": violations,
            "values": {"MW": mw, "LogP": logp, "HBD": hbd, "HBA": hba}
        }

    def _check_veber(self, mol: Chem.Mol) -> Dict[str, Any]:
        """Veber rules for oral bioavailability"""
        tpsa = Descriptors.TPSA(mol)
        rotatable = Descriptors.NumRotatableBonds(mol)

        violations = []
        if tpsa > 140:
            violations.append(f"TPSA={tpsa:.1f} > 140")
        if rotatable > 10:
            violations.append(f"RotBonds={rotatable} > 10")

        return {
            "pass": len(violations) == 0,
            "violations": violations,
            "values": {"TPSA": tpsa, "RotatableBonds": rotatable}
        }

    def _check_ghose(self, mol: Chem.Mol) -> Dict[str, Any]:
        """Ghose filter"""
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        mr = Descriptors.MolMR(mol)
        num_atoms = mol.GetNumAtoms()

        violations = []
        if not (160 <= mw <= 480):
            violations.append(f"MW={mw:.1f} not in [160,480]")
        if not (-0.4 <= logp <= 5.6):
            violations.append(f"LogP={logp:.2f} not in [-0.4,5.6]")
        if not (40 <= mr <= 130):
            violations.append(f"MR={mr:.1f} not in [40,130]")
        if not (20 <= num_atoms <= 70):
            violations.append(f"Atoms={num_atoms} not in [20,70]")

        return {
            "pass": len(violations) == 0,
            "violations": violations
        }

    def _check_egan(self, mol: Chem.Mol) -> Dict[str, Any]:
        """Egan egg (absorption)"""
        tpsa = Descriptors.TPSA(mol)
        logp = Descriptors.MolLogP(mol)

        violations = []
        if tpsa > 131.6:
            violations.append(f"TPSA={tpsa:.1f} > 131.6")
        if logp > 5.88:
            violations.append(f"LogP={logp:.2f} > 5.88")

        return {
            "pass": len(violations) == 0,
            "violations": violations
        }

    def _check_muegge(self, mol: Chem.Mol) -> Dict[str, Any]:
        """Muegge filter"""
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rotatable = Descriptors.NumRotatableBonds(mol)
        rings = rdMolDescriptors.CalcNumRings(mol)
        num_atoms = mol.GetNumAtoms()

        violations = []
        if not (200 <= mw <= 600):
            violations.append(f"MW={mw:.1f}")
        if not (-2 <= logp <= 5):
            violations.append(f"LogP={logp:.2f}")
        if tpsa > 150:
            violations.append(f"TPSA={tpsa:.1f}")
        if hbd > 5:
            violations.append(f"HBD={hbd}")
        if hba > 10:
            violations.append(f"HBA={hba}")
        if rotatable > 15:
            violations.append(f"RotBonds={rotatable}")
        if rings > 7:
            violations.append(f"Rings={rings}")
        if num_atoms < 5:
            violations.append(f"Atoms={num_atoms}")

        return {
            "pass": len(violations) == 0,
            "violations": violations
        }

    def _check_brenk(self, mol: Chem.Mol) -> List[str]:
        """Check for Brenk structural alerts"""
        # Simplified - full implementation would use complete Brenk filter set
        alerts = []

        # Check for reactive groups
        reactive_smarts = [
            ("[N+]([O-])=O", "Nitro group"),
            ("C#N", "Nitrile"),
            ("SC#N", "Thiocyanate"),
            ("[N-]=[N+]=[N-]", "Azide"),
        ]

        for smarts, name in reactive_smarts:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                alerts.append(name)

        return alerts

    def predict_admet(self, smiles: str) -> Dict[str, Any]:
        """Predict ADMET properties"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {}

        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = Descriptors.NumHDonors(mol)

        # Absorption predictions (simplified models)
        caco2_permeable = tpsa < 140 and mw < 500
        intestinal_absorption = "High" if (tpsa < 100 and logp > 0) else "Moderate" if tpsa < 140 else "Low"

        # Distribution
        bbb_permeable = tpsa < 90 and mw < 400 and logp > 0

        # Metabolism
        cyp_substrate_risk = "High" if logp > 3 else "Moderate" if logp > 1 else "Low"

        # Excretion (renal)
        renal_clearance = "High" if (mw < 300 and logp < 0) else "Moderate" if mw < 500 else "Low"

        # Toxicity alerts
        toxicity_alerts = self.check_pains(mol) + self._check_brenk(mol)

        return {
            "absorption": {
                "caco2_permeable": caco2_permeable,
                "intestinal_absorption": intestinal_absorption,
                "oral_bioavailability": "Likely" if caco2_permeable else "Unlikely"
            },
            "distribution": {
                "bbb_permeable": bbb_permeable,
                "volume_distribution": "Moderate"  # Would need ML model
            },
            "metabolism": {
                "cyp_substrate_risk": cyp_substrate_risk,
                "likely_cyp_enzymes": ["CYP3A4"] if logp > 2 else ["CYP2D6"]
            },
            "excretion": {
                "renal_clearance": renal_clearance,
                "half_life_estimate": self._predict_peptide_half_life(mol)
            },
            "toxicity": {
                "alerts": toxicity_alerts,
                "herg_risk": "Low" if logp < 3 else "Moderate",
                "hepatotoxicity_risk": "Low" if len(toxicity_alerts) == 0 else "Moderate"
            }
        }
