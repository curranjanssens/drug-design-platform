"""
Synthesis Feasibility Predictors

Predictions for how difficult a molecule is to synthesize.

Implements:
- SA Score (Synthetic Accessibility Score)
- TODO: Retrosynthesis-based assessment
"""

from typing import Dict, Any, List, Optional
from dataclasses import dataclass
import logging
import math

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem

# SA Score module is in Contrib which may not always be available
try:
    from rdkit.Contrib.SA_Score import sascorer
    SASCORER_AVAILABLE = True
except ImportError:
    sascorer = None
    SASCORER_AVAILABLE = False

from ..prediction import Prediction, ConfidenceLevel, PredictionBasis
from ..predictor_registry import Predictor, PredictorCapability, registry

logger = logging.getLogger(__name__)


@registry.register
class SASynthesisPredictor(Predictor):
    """
    Predict synthetic accessibility using the SA Score.

    The SA Score (Ertl and Schuffenhauer, 2009) estimates how easy
    it is to synthesize a molecule based on:
    - Fragment contributions (common vs rare substructures)
    - Complexity penalties (stereocenters, rings, etc.)

    Score ranges from 1 (easy) to 10 (hard).

    CAVEATS:
    - Based on fragment frequencies, not actual synthesis routes
    - May not capture all synthetic challenges
    - Should be combined with retrosynthesis for best results
    """

    capabilities = [
        PredictorCapability(
            predicts="synthetic_accessibility",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.MODERATE,
            description="Synthetic Accessibility Score (1=easy, 10=hard)",
            tags=["synthesis", "accessibility", "score"],
        ),
        PredictorCapability(
            predicts="synthesis_category",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.MODERATE,
            description="Synthesis difficulty category (easy/moderate/hard)",
            tags=["synthesis", "category"],
        ),
        PredictorCapability(
            predicts="synthesis_feasibility_score",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.MODERATE,
            description="Normalized synthesis feasibility (0=hard, 1=easy)",
            tags=["synthesis", "feasibility", "normalized"],
        ),
    ]

    def _get_mol(self, smiles: str) -> Chem.Mol:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        return mol

    def _analyze_complexity(self, mol: Chem.Mol) -> Dict[str, Any]:
        """Analyze structural complexity factors."""
        complexity = {
            "ring_count": rdMolDescriptors.CalcNumRings(mol),
            "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
            "stereocenters": len(Chem.FindMolChiralCenters(mol, includeUnassigned=True)),
            "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
            "heavy_atoms": mol.GetNumHeavyAtoms(),
            "heteroatoms": sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [6, 1]),
            "spiro_centers": rdMolDescriptors.CalcNumSpiroAtoms(mol),
            "bridgehead_atoms": rdMolDescriptors.CalcNumBridgeheadAtoms(mol),
        }

        # Count fused ring systems (complex fused systems are harder)
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() > 0:
            complexity["fused_ring_systems"] = self._count_fused_systems(mol)
        else:
            complexity["fused_ring_systems"] = 0

        return complexity

    def _count_fused_systems(self, mol: Chem.Mol) -> int:
        """Count the number of fused ring systems."""
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() <= 1:
            return 0

        # Simple heuristic: count shared atoms between rings
        atom_rings = ring_info.AtomRings()
        shared_count = 0
        for i, ring1 in enumerate(atom_rings):
            for ring2 in atom_rings[i+1:]:
                if set(ring1) & set(ring2):
                    shared_count += 1

        return shared_count

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

        # Calculate SA Score
        complexity = self._analyze_complexity(mol)
        if SASCORER_AVAILABLE and sascorer is not None:
            try:
                sa_score = sascorer.calculateScore(mol)
            except Exception as e:
                logger.warning(f"SA Score calculation failed: {e}")
                # Fall back to simple complexity-based estimate
                sa_score = 1 + (complexity["heavy_atoms"] / 10) + \
                          (complexity["stereocenters"] * 0.5) + \
                          (complexity["ring_count"] * 0.3)
                sa_score = min(10, max(1, sa_score))
        else:
            # SA Score module not available, use simple complexity-based estimate
            sa_score = 1 + (complexity["heavy_atoms"] / 10) + \
                      (complexity["stereocenters"] * 0.5) + \
                      (complexity["ring_count"] * 0.3)
            sa_score = min(10, max(1, sa_score))

        # Analyze complexity for detailed feedback
        complexity = self._analyze_complexity(mol)

        # Generate caveats based on complexity
        caveats = ["SA Score is a statistical estimate, not a synthesis plan"]
        synthesis_challenges = []

        if complexity["stereocenters"] > 2:
            synthesis_challenges.append(
                f"Multiple stereocenters ({complexity['stereocenters']}) - may require chiral synthesis"
            )
        if complexity["spiro_centers"] > 0:
            synthesis_challenges.append(
                f"Spiro centers ({complexity['spiro_centers']}) - specialized chemistry needed"
            )
        if complexity["bridgehead_atoms"] > 0:
            synthesis_challenges.append(
                f"Bridgehead atoms ({complexity['bridgehead_atoms']}) - complex bicyclic systems"
            )
        if complexity["fused_ring_systems"] > 2:
            synthesis_challenges.append(
                f"Multiple fused ring systems - may require multi-step synthesis"
            )
        if complexity["heavy_atoms"] > 40:
            synthesis_challenges.append(
                "Large molecule - synthesis may involve many steps"
            )

        if synthesis_challenges:
            caveats.extend(synthesis_challenges)

        if prop == "synthetic_accessibility":
            return Prediction(
                value=round(sa_score, 2),
                confidence=ConfidenceLevel.MODERATE,
                confidence_score=0.7,
                basis=PredictionBasis.HEURISTIC,
                caveats=caveats,
                metadata={
                    "complexity": complexity,
                    "interpretation": self._interpret_sa_score(sa_score),
                },
            )

        elif prop == "synthesis_category":
            if sa_score <= 4:
                category = "easy"
            elif sa_score <= 6:
                category = "moderate"
            else:
                category = "hard"

            return Prediction(
                value=category,
                confidence=ConfidenceLevel.MODERATE,
                confidence_score=0.65,
                basis=PredictionBasis.HEURISTIC,
                caveats=caveats,
                metadata={
                    "sa_score": round(sa_score, 2),
                    "complexity": complexity,
                },
            )

        elif prop == "synthesis_feasibility_score":
            # Normalize to 0-1 where 1 is easy
            # SA Score 1 -> 1.0, SA Score 10 -> 0.0
            normalized = max(0, (10 - sa_score) / 9)

            return Prediction(
                value=round(normalized, 3),
                confidence=ConfidenceLevel.MODERATE,
                confidence_score=0.65,
                basis=PredictionBasis.HEURISTIC,
                caveats=caveats,
                metadata={
                    "sa_score": round(sa_score, 2),
                    "complexity": complexity,
                },
            )

        raise ValueError(f"Unknown property: {prop}")

    def _interpret_sa_score(self, score: float) -> str:
        """Provide human-readable interpretation of SA Score."""
        if score <= 2:
            return "Very easy to synthesize - simple, common building blocks"
        elif score <= 3:
            return "Easy to synthesize - standard chemistry"
        elif score <= 4:
            return "Moderately easy - routine synthesis"
        elif score <= 5:
            return "Moderate difficulty - may require some optimization"
        elif score <= 6:
            return "Moderately difficult - multiple synthetic steps"
        elif score <= 7:
            return "Difficult - requires careful planning"
        elif score <= 8:
            return "Very difficult - specialized chemistry likely needed"
        elif score <= 9:
            return "Extremely difficult - significant synthetic challenge"
        else:
            return "Near impossible - may require novel synthetic methods"


# TODO: Add RetrosynthesisPredictor when integration with retrosynthesis tools is ready
"""
@registry.register
class RetrosynthesisPredictor(Predictor):
    '''
    Predict synthesis feasibility using retrosynthetic analysis.

    Uses retrosynthesis engines (ASKCOS, IBM RXN, etc.) to find
    actual synthesis routes and assess their feasibility.

    This is more reliable than SA Score but requires external API calls.
    '''

    capabilities = [
        PredictorCapability(
            predicts="retrosynthesis_routes",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.HIGH,
            description="Retrosynthetic routes to target molecule",
            tags=["synthesis", "retrosynthesis", "routes"],
        ),
        PredictorCapability(
            predicts="synthesis_steps",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.MODERATE,
            description="Estimated number of synthesis steps",
            tags=["synthesis", "steps"],
        ),
    ]

    async def predict(self, inputs: Dict, capability: PredictorCapability) -> Prediction:
        # TODO: Integrate with retrosynthesis API
        pass
"""
