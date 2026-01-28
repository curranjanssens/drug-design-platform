"""
Binding Affinity Predictors

Predictions for how strongly a molecule binds to a target.

Currently implements:
- Similarity-based predictions (using known active compounds)
- TODO: Structure-based docking when 3D structures available
"""

from typing import Dict, Any, List, Optional
from dataclasses import dataclass
import logging

from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit import RDLogger

from ..prediction import Prediction, ConfidenceLevel, PredictionBasis
from ..predictor_registry import Predictor, PredictorCapability, registry

logger = logging.getLogger(__name__)

# Suppress RDKit warnings for cleaner output
RDLogger.DisableLog('rdApp.*')


@dataclass
class ReferenceCompound:
    """A known active compound for similarity comparison."""
    smiles: str
    name: str
    activity_value: float  # e.g., pIC50, pKi
    activity_type: str  # "pIC50", "pKi", "IC50_nM", etc.
    target: str
    source: str  # Where data came from


# Known reference compounds for common targets
# In production, these would be loaded from ChEMBL or a database
REFERENCE_COMPOUNDS: Dict[str, List[ReferenceCompound]] = {
    "FAAH": [
        ReferenceCompound(
            smiles="O=C(Nc1ccc(OC(F)(F)F)cc1)N1CCC(c2ccnc3ccccc23)CC1",
            name="PF-04457845",
            activity_value=7.5,  # pIC50
            activity_type="pIC50",
            target="FAAH",
            source="ChEMBL",
        ),
        ReferenceCompound(
            smiles="CC(C)(C)c1ccc(C(=O)Nc2ccc(C(=O)c3ccccc3)cc2)cc1",
            name="URB597",
            activity_value=6.5,
            activity_type="pIC50",
            target="FAAH",
            source="ChEMBL",
        ),
    ],
    "COX2": [
        ReferenceCompound(
            smiles="CS(=O)(=O)c1ccc(C2=C(c3ccccc3)C(=O)OC2)cc1",
            name="Rofecoxib",
            activity_value=8.0,
            activity_type="pIC50",
            target="COX2",
            source="ChEMBL",
        ),
        ReferenceCompound(
            smiles="Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1",
            name="Celecoxib",
            activity_value=7.8,
            activity_type="pIC50",
            target="COX2",
            source="ChEMBL",
        ),
    ],
    # Add more targets as needed
}


@registry.register
class SimilarityBindingPredictor(Predictor):
    """
    Predict binding affinity based on similarity to known actives.

    This uses a similarity-based approach where:
    1. Find most similar known active compounds
    2. Estimate activity based on similarity and known activity
    3. Confidence depends on how similar the closest match is

    IMPORTANT CAVEATS:
    - This assumes similar structures have similar activity (SAR)
    - Cliff compounds (similar but different activity) will be missed
    - Should be used as rough guidance, not definitive prediction
    """

    capabilities = [
        PredictorCapability(
            predicts="binding_affinity_estimate",
            requires=["smiles", "target"],
            confidence_tier=ConfidenceLevel.LOW,
            description="Estimated binding affinity based on similarity to known actives",
            tags=["binding", "affinity", "similarity"],
        ),
        PredictorCapability(
            predicts="similarity_to_known_actives",
            requires=["smiles", "target"],
            confidence_tier=ConfidenceLevel.MODERATE,
            description="Tanimoto similarity to known active compounds",
            tags=["similarity", "novelty"],
        ),
    ]

    def __init__(self):
        super().__init__()
        self._fingerprint_cache: Dict[str, Any] = {}

    def _get_fingerprint(self, smiles: str):
        """Get Morgan fingerprint for a SMILES string."""
        if smiles in self._fingerprint_cache:
            return self._fingerprint_cache[smiles]

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        self._fingerprint_cache[smiles] = fp
        return fp

    def _calculate_similarity(self, smiles1: str, smiles2: str) -> Optional[float]:
        """Calculate Tanimoto similarity between two molecules."""
        fp1 = self._get_fingerprint(smiles1)
        fp2 = self._get_fingerprint(smiles2)

        if fp1 is None or fp2 is None:
            return None

        return DataStructs.TanimotoSimilarity(fp1, fp2)

    async def predict(
        self,
        inputs: Dict[str, Any],
        capability: PredictorCapability
    ) -> Prediction:
        smiles = inputs["smiles"]
        target = inputs["target"].upper()

        # Validate SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")

        # Get reference compounds for this target
        references = REFERENCE_COMPOUNDS.get(target, [])

        if not references:
            return Prediction.not_assessed(
                f"No reference compounds available for target: {target}",
                predictor_id=self.predictor_id
            )

        # Calculate similarities to all references
        similarities = []
        for ref in references:
            sim = self._calculate_similarity(smiles, ref.smiles)
            if sim is not None:
                similarities.append((ref, sim))

        if not similarities:
            return Prediction.not_assessed(
                "Could not calculate similarity to reference compounds",
                predictor_id=self.predictor_id
            )

        # Sort by similarity (highest first)
        similarities.sort(key=lambda x: x[1], reverse=True)
        best_match, best_sim = similarities[0]

        prop = capability.predicts

        if prop == "similarity_to_known_actives":
            # Return the highest similarity
            caveats = []
            if best_sim < 0.3:
                caveats.append("Low similarity to known actives - novel scaffold")
            elif best_sim < 0.5:
                caveats.append("Moderate similarity - may have different SAR")

            return Prediction(
                value=round(best_sim, 3),
                confidence=ConfidenceLevel.MODERATE,
                confidence_score=0.8,
                basis=PredictionBasis.SIMILARITY,
                caveats=caveats,
                metadata={
                    "best_match_name": best_match.name,
                    "best_match_smiles": best_match.smiles,
                    "top_similarities": [
                        {"name": ref.name, "similarity": round(sim, 3)}
                        for ref, sim in similarities[:3]
                    ],
                },
            )

        elif prop == "binding_affinity_estimate":
            # Estimate binding affinity based on similarity
            # Simple approach: interpolate from known activity

            # Determine confidence based on similarity
            if best_sim >= 0.7:
                # High similarity - fairly confident in prediction
                confidence = ConfidenceLevel.MODERATE
                conf_score = 0.6
                estimated_activity = best_match.activity_value
                caveats = [
                    f"Estimated from {best_match.name} (similarity: {best_sim:.2f})",
                    "High similarity suggests comparable activity"
                ]
            elif best_sim >= 0.4:
                # Moderate similarity - less confident
                confidence = ConfidenceLevel.LOW
                conf_score = 0.35
                # Average of top matches, weighted by similarity
                weighted_sum = sum(ref.activity_value * sim for ref, sim in similarities[:3])
                total_sim = sum(sim for _, sim in similarities[:3])
                estimated_activity = weighted_sum / total_sim if total_sim > 0 else 5.0
                caveats = [
                    "Moderate similarity to known actives",
                    "Activity cliff possible - verify experimentally"
                ]
            else:
                # Low similarity - prediction is unreliable
                confidence = ConfidenceLevel.OUTSIDE_DOMAIN
                conf_score = 0.15
                # Use average of known activities as rough guess
                estimated_activity = sum(ref.activity_value for ref in references) / len(references)
                caveats = [
                    f"Low similarity to known actives (max: {best_sim:.2f})",
                    "Novel scaffold - prediction highly uncertain",
                    "Consider this a rough estimate only"
                ]

            # Add general caveat about similarity-based prediction
            caveats.append("Similarity-based estimate - does not account for SAR cliffs")

            return Prediction(
                value=round(estimated_activity, 2),
                confidence=confidence,
                confidence_score=conf_score,
                basis=PredictionBasis.SIMILARITY,
                caveats=caveats,
                metadata={
                    "activity_type": best_match.activity_type,
                    "best_match": best_match.name,
                    "similarity": best_sim,
                    "unit": "pIC50" if "pIC50" in best_match.activity_type else best_match.activity_type,
                },
            )

        raise ValueError(f"Unknown property: {prop}")


# TODO: Add DockingPredictor when structure-based docking is available
# This would use AutoDock Vina, Glide, or similar docking engines
# and requires target 3D structure
"""
@registry.register
class DockingPredictor(Predictor):
    '''
    Structure-based binding prediction using molecular docking.

    Requires:
    - 3D structure of target (PDB file or AlphaFold model)
    - Binding site definition
    - Docking engine (AutoDock Vina, etc.)
    '''

    capabilities = [
        PredictorCapability(
            predicts="docking_score",
            requires=["smiles", "target_structure"],
            confidence_tier=ConfidenceLevel.MODERATE,
            description="Docking score from structure-based docking",
            tags=["binding", "docking", "structure"],
        ),
    ]

    async def predict(self, inputs: Dict, capability: PredictorCapability) -> Prediction:
        # TODO: Implement docking
        pass
"""
