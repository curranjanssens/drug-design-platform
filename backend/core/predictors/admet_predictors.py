"""
ADMET Predictors

Absorption, Distribution, Metabolism, Excretion, Toxicity predictions.

Currently implements rule-based predictions with appropriate uncertainty.
TODO: Add ML model-based predictions when models are available.
"""

from typing import Dict, Any, List, Optional
from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski

from ..prediction import Prediction, ConfidenceLevel, PredictionBasis
from ..predictor_registry import Predictor, PredictorCapability, registry


@dataclass
class ADMETPrediction:
    """Structured ADMET prediction result."""
    property_name: str
    value: Any
    unit: str
    threshold: Optional[float]
    passes_threshold: Optional[bool]
    interpretation: str


@registry.register
class RuleBasedADMETPredictor(Predictor):
    """
    Rule-based ADMET predictions using molecular descriptors.

    These are HEURISTIC predictions based on empirical rules (Lipinski, Veber, etc.)
    They should be treated as rough estimates, not validated predictions.

    Confidence is LOW to MODERATE because these are rules, not validated models.
    """

    capabilities = [
        # Absorption
        PredictorCapability(
            predicts="oral_absorption_estimate",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.LOW,
            description="Estimated oral absorption (rule-based)",
            tags=["admet", "absorption"],
        ),
        PredictorCapability(
            predicts="caco2_estimate",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.LOW,
            description="Estimated Caco-2 permeability category",
            tags=["admet", "absorption", "permeability"],
        ),

        # Distribution
        PredictorCapability(
            predicts="bbb_penetration_estimate",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.LOW,
            description="Estimated BBB penetration (rule-based)",
            tags=["admet", "distribution", "cns"],
        ),
        PredictorCapability(
            predicts="pgp_substrate_estimate",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.LOW,
            description="Estimated P-gp substrate likelihood",
            tags=["admet", "distribution", "transporter"],
        ),

        # Metabolism
        PredictorCapability(
            predicts="cyp_liability_estimate",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.LOW,
            description="Estimated CYP inhibition liability",
            tags=["admet", "metabolism", "cyp"],
        ),

        # Toxicity
        PredictorCapability(
            predicts="herg_liability_estimate",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.LOW,
            description="Estimated hERG liability (rule-based)",
            tags=["admet", "toxicity", "cardiac"],
        ),
        PredictorCapability(
            predicts="ames_estimate",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.LOW,
            description="Estimated Ames mutagenicity",
            tags=["admet", "toxicity", "mutagenicity"],
        ),

        # Overall
        PredictorCapability(
            predicts="admet_safety_score",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.LOW,
            description="Overall ADMET safety score (0-1)",
            tags=["admet", "safety", "score"],
        ),
    ]

    def _get_mol(self, smiles: str) -> Chem.Mol:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        return mol

    async def predict(
        self,
        inputs: Dict[str, Any],
        capability: PredictorCapability
    ) -> Prediction:
        smiles = inputs["smiles"]
        mol = self._get_mol(smiles)
        prop = capability.predicts

        # Calculate common properties once
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rotatable = Descriptors.NumRotatableBonds(mol)

        # Standard caveat for all rule-based predictions
        base_caveat = "Rule-based estimate - experimental validation required"

        if prop == "oral_absorption_estimate":
            # Lipinski + Veber rules
            violations = 0
            details = []

            if mw > 500:
                violations += 1
                details.append(f"MW ({mw:.0f}) > 500")
            if logp > 5:
                violations += 1
                details.append(f"LogP ({logp:.1f}) > 5")
            if hbd > 5:
                violations += 1
                details.append(f"HBD ({hbd}) > 5")
            if hba > 10:
                violations += 1
                details.append(f"HBA ({hba}) > 10")
            if tpsa > 140:
                violations += 1
                details.append(f"TPSA ({tpsa:.0f}) > 140")
            if rotatable > 10:
                violations += 1
                details.append(f"Rotatable bonds ({rotatable}) > 10")

            if violations == 0:
                estimate = "high"
                score = 0.85
            elif violations <= 2:
                estimate = "moderate"
                score = 0.6
            else:
                estimate = "low"
                score = 0.3

            caveats = [base_caveat]
            if details:
                caveats.append(f"Violations: {'; '.join(details)}")

            return Prediction(
                value=estimate,
                confidence=ConfidenceLevel.LOW,
                confidence_score=0.4,
                basis=PredictionBasis.HEURISTIC,
                caveats=caveats,
                metadata={
                    "violations": violations,
                    "details": details,
                    "score": score,
                },
            )

        elif prop == "caco2_estimate":
            # Based on TPSA and LogP
            if tpsa < 60 and 0 < logp < 3:
                estimate = "high"
                score = 0.8
            elif tpsa < 100 and -1 < logp < 5:
                estimate = "moderate"
                score = 0.5
            else:
                estimate = "low"
                score = 0.2

            caveats = [base_caveat]
            if tpsa > 100:
                caveats.append(f"High TPSA ({tpsa:.0f} Å²) suggests poor permeability")

            return Prediction(
                value=estimate,
                confidence=ConfidenceLevel.LOW,
                confidence_score=0.35,
                basis=PredictionBasis.HEURISTIC,
                caveats=caveats,
                metadata={"tpsa": tpsa, "logp": logp, "score": score},
            )

        elif prop == "bbb_penetration_estimate":
            # Simple BBB rules: MW < 450, TPSA < 90, HBD < 3
            passes = mw < 450 and tpsa < 90 and hbd < 3 and 1 < logp < 4

            caveats = [base_caveat]
            reasons = []

            if mw >= 450:
                reasons.append(f"MW ({mw:.0f}) >= 450")
            if tpsa >= 90:
                reasons.append(f"TPSA ({tpsa:.0f}) >= 90")
            if hbd >= 3:
                reasons.append(f"HBD ({hbd}) >= 3")
            if not (1 < logp < 4):
                reasons.append(f"LogP ({logp:.1f}) outside 1-4 range")

            if reasons:
                caveats.append(f"BBB concerns: {'; '.join(reasons)}")

            return Prediction(
                value=passes,
                confidence=ConfidenceLevel.LOW,
                confidence_score=0.35,
                basis=PredictionBasis.HEURISTIC,
                caveats=caveats,
                metadata={"reasons": reasons},
            )

        elif prop == "pgp_substrate_estimate":
            # Rough heuristic: high MW, high TPSA, many HBD/HBA
            pgp_risk = mw > 400 or tpsa > 100 or (hbd + hba) > 10

            caveats = [base_caveat]
            if pgp_risk:
                caveats.append("Properties suggest possible P-gp substrate")

            return Prediction(
                value=pgp_risk,
                confidence=ConfidenceLevel.LOW,
                confidence_score=0.3,
                basis=PredictionBasis.HEURISTIC,
                caveats=caveats,
            )

        elif prop == "cyp_liability_estimate":
            # Very rough: lipophilic compounds more likely to inhibit CYPs
            liability = "high" if logp > 4 else "moderate" if logp > 2 else "low"

            return Prediction(
                value=liability,
                confidence=ConfidenceLevel.LOW,
                confidence_score=0.25,
                basis=PredictionBasis.HEURISTIC,
                caveats=[base_caveat, "CYP inhibition is highly dependent on specific structure"],
            )

        elif prop == "herg_liability_estimate":
            # Basic positive nitrogen and lipophilicity check
            # hERG blockers often have: LogP > 3, basic nitrogen, aromatic
            num_basic_n = sum(1 for atom in mol.GetAtoms()
                             if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() >= 0)
            aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)

            risk_score = 0
            reasons = []

            if logp > 3:
                risk_score += 1
                reasons.append(f"LogP ({logp:.1f}) > 3")
            if num_basic_n > 0:
                risk_score += 1
                reasons.append(f"Contains basic nitrogen")
            if aromatic_rings >= 2:
                risk_score += 1
                reasons.append(f"{aromatic_rings} aromatic rings")

            if risk_score >= 2:
                estimate = "high"
            elif risk_score == 1:
                estimate = "moderate"
            else:
                estimate = "low"

            caveats = [base_caveat]
            if reasons:
                caveats.append(f"hERG risk factors: {'; '.join(reasons)}")

            return Prediction(
                value=estimate,
                confidence=ConfidenceLevel.LOW,
                confidence_score=0.3,
                basis=PredictionBasis.HEURISTIC,
                caveats=caveats,
                metadata={"risk_factors": reasons, "risk_score": risk_score},
            )

        elif prop == "ames_estimate":
            # Check for common mutagenic structural alerts
            # This is very simplified - real Ames prediction needs proper alerts
            mutagenic_smarts = [
                "[N+](=O)[O-]",  # Nitro
                "[N]=[N]",       # Azo
                "N(=O)O",        # Nitroso-like
            ]

            alerts_found = []
            for smarts in mutagenic_smarts:
                pattern = Chem.MolFromSmarts(smarts)
                if pattern and mol.HasSubstructMatch(pattern):
                    alerts_found.append(smarts)

            has_alerts = len(alerts_found) > 0

            caveats = [base_caveat]
            if alerts_found:
                caveats.append(f"Mutagenic alerts: {len(alerts_found)} patterns found")

            return Prediction(
                value=has_alerts,
                confidence=ConfidenceLevel.LOW,
                confidence_score=0.35,
                basis=PredictionBasis.HEURISTIC,
                caveats=caveats,
                metadata={"alerts": alerts_found},
            )

        elif prop == "admet_safety_score":
            # Aggregate safety score based on multiple factors
            score = 1.0
            issues = []

            # Lipinski compliance
            lipinski_violations = sum([
                mw > 500, logp > 5, hbd > 5, hba > 10
            ])
            if lipinski_violations > 0:
                score -= 0.1 * lipinski_violations
                issues.append(f"{lipinski_violations} Lipinski violations")

            # Permeability concerns
            if tpsa > 140:
                score -= 0.15
                issues.append("High TPSA")

            # Flexibility
            if rotatable > 10:
                score -= 0.1
                issues.append("High flexibility")

            # hERG risk (simplified)
            if logp > 3:
                score -= 0.1
                issues.append("Lipophilic (hERG risk)")

            score = max(0.0, min(1.0, score))

            caveats = [base_caveat]
            if issues:
                caveats.append(f"Issues: {'; '.join(issues)}")

            return Prediction(
                value=round(score, 2),
                confidence=ConfidenceLevel.LOW,
                confidence_score=0.35,
                basis=PredictionBasis.HEURISTIC,
                caveats=caveats,
                metadata={"issues": issues},
            )

        raise ValueError(f"Unknown property: {prop}")


@registry.register
class ChemicalStabilityPredictor(Predictor):
    """
    Predicts chemical stability based on known unstable functional group patterns.

    This catches issues like:
    - N-CF3 groups (unstable N-C bond with electron-withdrawing CF3)
    - Peroxides, N-oxides in sensitive contexts
    - Labile esters, hydrazines, etc.
    """

    capabilities = [
        PredictorCapability(
            predicts="chemical_stability",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.MODERATE,
            description="Chemical stability assessment based on structural alerts",
            tags=["stability", "chemistry", "alerts"],
        ),
        PredictorCapability(
            predicts="stability_issues",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.MODERATE,
            description="List of potential stability issues",
            tags=["stability", "chemistry", "alerts"],
        ),
    ]

    # SMARTS patterns for unstable functional groups
    # Each entry: (SMARTS, issue_description, severity: critical/warning)
    STABILITY_ALERTS = [
        # N-CF3 is chemically unstable - the N-C bond is weak
        ("[NX3]C(F)(F)F", "N-CF3 group is chemically unstable (weak N-C bond)", "critical"),
        ("[nX2]C(F)(F)F", "N-CF3 on aromatic N is unstable", "critical"),

        # Acyl fluorides - very reactive
        ("C(=O)F", "Acyl fluoride is highly reactive", "critical"),

        # Peroxides
        ("OO", "Peroxide group is unstable and potentially explosive", "critical"),

        # Triazene
        ("N=NN", "Triazene is photolabile and potentially mutagenic", "critical"),

        # Acid chlorides in final compounds
        ("C(=O)Cl", "Acid chloride is too reactive for drug use", "critical"),

        # Vinyl halides (some are carcinogenic)
        ("C=CCl", "Vinyl chloride is potentially carcinogenic", "warning"),
        ("C=CBr", "Vinyl bromide is potentially carcinogenic", "warning"),

        # Epoxides (reactive, potential mutagen)
        ("C1OC1", "Epoxide is reactive and potentially mutagenic", "warning"),

        # Hydrazines (oxidatively labile, potential mutagen)
        ("[NH2][NH2]", "Hydrazine is oxidatively labile", "warning"),
        ("N(N)C=O", "Acyl hydrazine is oxidatively labile", "warning"),

        # Michael acceptors (reactive, electrophilic)
        # Only flag if not intentional (like in covalent drugs)
        # ("[CX3]=[CX3][CX3]=O", "Michael acceptor - reactive electrophile", "warning"),

        # Thiols - oxidize in plasma, bind albumin Cys34, poor oral bioavailability
        ("[SH]", "Free thiol is reactive/oxidizable - not suitable for drugs", "critical"),

        # Disulfides - can scramble
        ("SS", "Disulfide may scramble or reduce", "warning"),

        # Aldehydes - reactive
        ("[CH]=O", "Aldehyde is reactive (consider masked form)", "warning"),

        # Primary aromatic amines - potential metabolic activation
        ("c[NH2]", "Primary aromatic amine - metabolic liability", "warning"),

        # Nitro aromatics - potential mutagenicity
        ("c[N+](=O)[O-]", "Aromatic nitro - mutagenicity concern", "warning"),

        # Phenolic ester - hydrolytically labile
        ("cOC(=O)", "Phenolic ester is hydrolytically labile", "warning"),

        # Anhydrides
        ("C(=O)OC(=O)", "Anhydride is hydrolytically unstable", "critical"),

        # Carbamic acid (decomposes to CO2 + amine)
        ("[NH]C(=O)O", "Carbamic acid derivative may be unstable", "warning"),
    ]

    def _get_mol(self, smiles: str) -> Chem.Mol:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        return mol

    async def predict(
        self,
        inputs: Dict[str, Any],
        capability: PredictorCapability
    ) -> Prediction:
        smiles = inputs["smiles"]
        mol = self._get_mol(smiles)
        prop = capability.predicts

        critical_issues = []
        warnings = []

        for smarts, description, severity in self.STABILITY_ALERTS:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                matches = mol.GetSubstructMatches(pattern)
                issue = {
                    "pattern": smarts,
                    "description": description,
                    "count": len(matches),
                }
                if severity == "critical":
                    critical_issues.append(issue)
                else:
                    warnings.append(issue)

        if prop == "stability_issues":
            all_issues = critical_issues + warnings
            return Prediction(
                value=all_issues,
                confidence=ConfidenceLevel.MODERATE,
                confidence_score=0.7,
                basis=PredictionBasis.HEURISTIC,
                caveats=["Based on known structural alerts - context matters"],
                metadata={
                    "critical_count": len(critical_issues),
                    "warning_count": len(warnings),
                },
            )

        elif prop == "chemical_stability":
            # Score from 0 (unstable) to 1 (stable)
            score = 1.0
            score -= 0.3 * len(critical_issues)
            score -= 0.1 * len(warnings)
            score = max(0.0, score)

            if critical_issues:
                category = "unstable"
            elif warnings:
                category = "caution"
            else:
                category = "stable"

            caveats = []
            if critical_issues:
                caveats.append(f"Critical issues: {[i['description'] for i in critical_issues]}")
            if warnings:
                caveats.append(f"Warnings: {[i['description'] for i in warnings]}")

            return Prediction(
                value=category,
                confidence=ConfidenceLevel.MODERATE,
                confidence_score=0.7,
                basis=PredictionBasis.HEURISTIC,
                caveats=caveats if caveats else ["No stability alerts found"],
                metadata={
                    "score": score,
                    "critical_issues": critical_issues,
                    "warnings": warnings,
                },
            )

        raise ValueError(f"Unknown property: {prop}")


# TODO: Add ModelBasedADMETPredictor when ML models are integrated
# This would use actual trained models with proper applicability domain checking
