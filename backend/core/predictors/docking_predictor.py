"""
Docking Predictor Stub

This module provides a placeholder for molecular docking predictions.
Currently returns NOT_ASSESSED for all docking requests.

Future integration options:
- AutoDock Vina (local, open-source, good accuracy)
- GNINA (local, CNN-based, better for flexible docking)
- Glide (commercial, highest accuracy for rigid docking)

The PDB fetcher is still available for target structure retrieval.
"""

import logging
from typing import Dict, Any, Optional

from ..prediction import Prediction, ConfidenceLevel, PredictionBasis
from ..predictor_registry import Predictor, PredictorCapability, registry

logger = logging.getLogger(__name__)


@registry.register
class DiffDockPredictor(Predictor):
    """
    Placeholder for molecular docking predictions.

    This predictor currently returns NOT_ASSESSED for all requests.
    The interface is preserved for future integration with local
    docking tools like AutoDock Vina or GNINA.

    Capabilities are registered but will return "not assessed" status
    to indicate docking is not currently available.
    """

    capabilities = [
        PredictorCapability(
            predicts="docking_score",
            requires=["smiles", "target"],
            confidence_tier=ConfidenceLevel.MODERATE,
            description="Molecular docking score (not currently available)",
            tags=["binding", "docking", "structure"],
        ),
        PredictorCapability(
            predicts="binding_pose_confidence",
            requires=["smiles", "target"],
            confidence_tier=ConfidenceLevel.MODERATE,
            description="Confidence in predicted binding pose (not currently available)",
            tags=["binding", "pose", "confidence"],
        ),
        PredictorCapability(
            predicts="binding_affinity_estimate",
            requires=["smiles", "target"],
            confidence_tier=ConfidenceLevel.LOW,
            description="Estimated binding affinity (not currently available)",
            tags=["binding", "affinity", "energy"],
        ),
        PredictorCapability(
            predicts="docking_validated",
            requires=["smiles", "target"],
            confidence_tier=ConfidenceLevel.MODERATE,
            description="Whether docking suggests valid binding (not currently available)",
            tags=["binding", "validation", "boolean"],
        ),
    ]

    async def predict(
        self,
        inputs: Dict[str, Any],
        capability: PredictorCapability,
    ) -> Prediction:
        """
        Returns NOT_ASSESSED for all docking predictions.

        Docking requires integration with a local docking tool.
        This stub preserves the interface for future implementation.
        """
        return Prediction.not_assessed(
            "Molecular docking not currently configured. "
            "Consider integrating AutoDock Vina or GNINA for local docking.",
            predictor_id=self.predictor_id
        )


# Utility function stub for standalone docking
async def dock_molecule(
    smiles: str,
    target: str,
    num_samples: int = 10,
) -> Optional[Dict[str, Any]]:
    """
    Convenience function to dock a molecule against a target.

    Currently returns None as docking is not configured.

    Future implementation should return:
        - confidence: Best pose confidence (0-1)
        - affinity: Estimated binding affinity (kcal/mol)
        - validated: True if confidence > threshold
        - runtime: Seconds taken
    """
    logger.info(
        f"dock_molecule called but docking not configured. "
        f"SMILES: {smiles[:30]}..., target: {target}"
    )
    return None
