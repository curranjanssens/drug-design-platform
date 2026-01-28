"""
Molecular Property Predictors

High-confidence predictors using RDKit for fundamental molecular properties.
These are validated, well-understood calculations.
"""

from typing import Dict, Any, List
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski, FilterCatalog

from ..prediction import Prediction, ConfidenceLevel, PredictionBasis
from ..predictor_registry import Predictor, PredictorCapability, registry


@registry.register
class RDKitPropertyPredictor(Predictor):
    """
    Calculates fundamental molecular properties using RDKit.

    These are exact calculations (not predictions), so confidence is HIGH.
    """

    capabilities = [
        PredictorCapability(
            predicts="molecular_weight",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.HIGH,
            description="Molecular weight (Da)",
            tags=["property", "basic"],
        ),
        PredictorCapability(
            predicts="logp",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.HIGH,
            description="Calculated LogP (Wildman-Crippen)",
            tags=["property", "lipophilicity"],
        ),
        PredictorCapability(
            predicts="tpsa",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.HIGH,
            description="Topological Polar Surface Area (Å²)",
            tags=["property", "polarity"],
        ),
        PredictorCapability(
            predicts="hbd",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.HIGH,
            description="Hydrogen Bond Donors",
            tags=["property", "hbond"],
        ),
        PredictorCapability(
            predicts="hba",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.HIGH,
            description="Hydrogen Bond Acceptors",
            tags=["property", "hbond"],
        ),
        PredictorCapability(
            predicts="rotatable_bonds",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.HIGH,
            description="Number of rotatable bonds",
            tags=["property", "flexibility"],
        ),
        PredictorCapability(
            predicts="aromatic_rings",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.HIGH,
            description="Number of aromatic rings",
            tags=["property", "aromaticity"],
        ),
        PredictorCapability(
            predicts="heavy_atoms",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.HIGH,
            description="Number of heavy atoms",
            tags=["property", "size"],
        ),
        PredictorCapability(
            predicts="fraction_sp3",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.HIGH,
            description="Fraction of sp3 carbons",
            tags=["property", "saturation"],
        ),
        PredictorCapability(
            predicts="lipinski_violations",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.HIGH,
            description="Number of Lipinski Rule of 5 violations",
            tags=["property", "druglikeness"],
        ),
    ]

    def _get_mol(self, smiles: str) -> Chem.Mol:
        """Parse SMILES to RDKit mol, raise if invalid."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        return mol

    async def predict(
        self,
        inputs: Dict[str, Any],
        capability: PredictorCapability
    ) -> Prediction:
        smiles = inputs.get("smiles", "")
        if not smiles or not smiles.strip():
            return Prediction.not_assessed(
                "No SMILES provided",
                predictor_id=self.predictor_id
            )
        mol = self._get_mol(smiles)

        prop = capability.predicts

        # Calculate the requested property
        if prop == "molecular_weight":
            value = Descriptors.MolWt(mol)
        elif prop == "logp":
            value = Descriptors.MolLogP(mol)
        elif prop == "tpsa":
            value = Descriptors.TPSA(mol)
        elif prop == "hbd":
            value = Descriptors.NumHDonors(mol)
        elif prop == "hba":
            value = Descriptors.NumHAcceptors(mol)
        elif prop == "rotatable_bonds":
            value = Descriptors.NumRotatableBonds(mol)
        elif prop == "aromatic_rings":
            value = rdMolDescriptors.CalcNumAromaticRings(mol)
        elif prop == "heavy_atoms":
            value = mol.GetNumHeavyAtoms()
        elif prop == "fraction_sp3":
            value = rdMolDescriptors.CalcFractionCSP3(mol)
        elif prop == "lipinski_violations":
            violations = 0
            if Descriptors.MolWt(mol) > 500:
                violations += 1
            if Descriptors.MolLogP(mol) > 5:
                violations += 1
            if Descriptors.NumHDonors(mol) > 5:
                violations += 1
            if Descriptors.NumHAcceptors(mol) > 10:
                violations += 1
            value = violations
        else:
            raise ValueError(f"Unknown property: {prop}")

        return Prediction(
            value=round(value, 3) if isinstance(value, float) else value,
            confidence=ConfidenceLevel.HIGH,
            confidence_score=0.99,
            basis=PredictionBasis.HEURISTIC,  # These are calculations, not predictions
            caveats=[],
        )


@registry.register
class DrugLikenessPredictor(Predictor):
    """
    Drug-likeness scores and filters.
    """

    capabilities = [
        PredictorCapability(
            predicts="qed",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.HIGH,
            description="Quantitative Estimate of Drug-likeness (0-1)",
            tags=["druglikeness", "score"],
        ),
        PredictorCapability(
            predicts="pains_alerts",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.HIGH,
            description="PAINS (Pan Assay Interference) alerts",
            tags=["filter", "pains"],
        ),
        PredictorCapability(
            predicts="brenk_alerts",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.HIGH,
            description="Brenk structural alerts",
            tags=["filter", "structural_alert"],
        ),
    ]

    def __init__(self):
        super().__init__()
        self._pains_catalog = None
        self._brenk_catalog = None

    async def initialize(self):
        """Load filter catalogs."""
        # PAINS filters
        pains_params = FilterCatalog.FilterCatalogParams()
        pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        self._pains_catalog = FilterCatalog.FilterCatalog(pains_params)

        # Brenk filters
        brenk_params = FilterCatalog.FilterCatalogParams()
        brenk_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
        self._brenk_catalog = FilterCatalog.FilterCatalog(brenk_params)

    async def predict(
        self,
        inputs: Dict[str, Any],
        capability: PredictorCapability
    ) -> Prediction:
        smiles = inputs.get("smiles", "")
        if not smiles or not smiles.strip():
            return Prediction.not_assessed(
                "No SMILES provided",
                predictor_id=self.predictor_id
            )
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")

        prop = capability.predicts

        if prop == "qed":
            try:
                value = Descriptors.qed(mol)
                return Prediction(
                    value=round(value, 3),
                    confidence=ConfidenceLevel.HIGH,
                    confidence_score=0.95,
                    basis=PredictionBasis.VALIDATED_MODEL,
                    caveats=[],
                )
            except Exception:
                return Prediction(
                    value=0.5,
                    confidence=ConfidenceLevel.LOW,
                    confidence_score=0.3,
                    basis=PredictionBasis.HEURISTIC,
                    caveats=["QED calculation failed, using default"],
                )

        elif prop == "pains_alerts":
            entries = self._pains_catalog.GetMatches(mol)
            alerts = [entry.GetDescription() for entry in entries]

            caveats = []
            if alerts:
                caveats = [f"PAINS alert: {a}" for a in alerts[:3]]  # Limit to 3

            return Prediction(
                value=len(alerts),
                confidence=ConfidenceLevel.HIGH,
                confidence_score=0.95,
                basis=PredictionBasis.HEURISTIC,
                caveats=caveats,
                metadata={"alerts": alerts},
            )

        elif prop == "brenk_alerts":
            entries = self._brenk_catalog.GetMatches(mol)
            alerts = [entry.GetDescription() for entry in entries]

            caveats = []
            if alerts:
                caveats = [f"Structural alert: {a}" for a in alerts[:3]]

            return Prediction(
                value=len(alerts),
                confidence=ConfidenceLevel.HIGH,
                confidence_score=0.95,
                basis=PredictionBasis.HEURISTIC,
                caveats=caveats,
                metadata={"alerts": alerts},
            )

        raise ValueError(f"Unknown property: {prop}")
