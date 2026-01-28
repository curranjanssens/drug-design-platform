"""
Prediction Service

High-level service that integrates the predictor registry into the drug design pipeline.
This is the main interface for making predictions throughout the application.

Usage:
    from core.prediction_service import prediction_service

    # Score a molecule comprehensively
    result = await prediction_service.score_molecule(
        smiles="CCO",
        target="FAAH",
        reference_compounds=["CC(=O)Nc1ccc(O)cc1"]
    )

    # Get specific predictions
    predictions = await prediction_service.predict(
        smiles="CCO",
        properties=["molecular_weight", "logp", "qed"]
    )
"""

import asyncio
import logging
from typing import Dict, List, Optional, Any, Set
from dataclasses import dataclass, field

from .prediction import Prediction, ConfidenceLevel, PredictionBasis, aggregate_predictions
from .predictor_registry import registry

# Import all predictors to register them with the registry
# This ensures predictors are available when the service is imported
from .predictors import (
    property_predictors,
    admet_predictors,
    synthesis_predictors,
    binding_predictors,
    novelty_predictors,
)

logger = logging.getLogger(__name__)


@dataclass
class MoleculeScore:
    """
    Comprehensive scoring result for a molecule.

    This aggregates multiple predictions into a single actionable result,
    while preserving uncertainty information.
    """
    smiles: str

    # Overall scores (0-1 where 1 is best)
    overall_score: float
    overall_confidence: ConfidenceLevel

    # Component scores
    drug_likeness_score: float
    safety_score: float
    synthesis_score: float
    novelty_score: float
    binding_score: Optional[float] = None  # Only if target provided

    # Detailed predictions (property -> Prediction)
    predictions: Dict[str, Prediction] = field(default_factory=dict)

    # Issues and warnings
    critical_issues: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    # Metadata
    target: Optional[str] = None
    assessed_properties: List[str] = field(default_factory=list)
    unassessed_properties: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "smiles": self.smiles,
            "overall_score": round(self.overall_score, 3),
            "overall_confidence": self.overall_confidence.value,
            "component_scores": {
                "drug_likeness": round(self.drug_likeness_score, 3),
                "safety": round(self.safety_score, 3),
                "synthesis": round(self.synthesis_score, 3),
                "novelty": round(self.novelty_score, 3),
                "binding": round(self.binding_score, 3) if self.binding_score else None,
            },
            "critical_issues": self.critical_issues,
            "warnings": self.warnings,
            "target": self.target,
            "predictions": {k: v.to_dict() for k, v in self.predictions.items()},
            "assessed_properties": self.assessed_properties,
            "unassessed_properties": self.unassessed_properties,
        }


class PredictionService:
    """
    High-level prediction service for drug design.

    This wraps the predictor registry and provides:
    - Comprehensive molecule scoring
    - Property-specific predictions
    - Batch processing
    - Caching and optimization
    """

    def __init__(self):
        self._registry = registry
        self._initialized = False

    async def initialize(self):
        """Initialize all predictors."""
        if self._initialized:
            return

        # Import predictors to register them
        from . import predictors  # noqa: F401

        self._initialized = True
        logger.info(f"PredictionService initialized with {len(self._registry.list_predictors())} predictors")

    async def predict(
        self,
        smiles: str,
        properties: Optional[List[str]] = None,
        target: Optional[str] = None,
        reference_smiles: Optional[List[str]] = None,
    ) -> Dict[str, Prediction]:
        """
        Predict specific properties for a molecule.

        Args:
            smiles: SMILES string of molecule
            properties: List of properties to predict (None = all available)
            target: Target protein (for binding predictions)
            reference_smiles: Reference compounds (for novelty predictions)

        Returns:
            Dictionary mapping property name to Prediction
        """
        await self.initialize()

        # Build input dictionary
        inputs = {"smiles": smiles}
        if target:
            inputs["target"] = target
        if reference_smiles:
            inputs["reference_smiles"] = reference_smiles

        # Get available properties
        available = self._registry.what_can_i_predict(inputs)

        if properties:
            # Filter to requested properties that are available
            to_predict = [p for p in properties if p in available]
            unavailable = [p for p in properties if p not in available]

            if unavailable:
                logger.warning(f"Cannot predict: {unavailable}")
        else:
            to_predict = available

        # Make predictions
        predictions = await self._registry.predict_all(inputs, to_predict)

        return predictions

    async def score_molecule(
        self,
        smiles: str,
        target: Optional[str] = None,
        reference_compounds: Optional[List[str]] = None,
        scoring_weights: Optional[Dict[str, float]] = None,
    ) -> MoleculeScore:
        """
        Comprehensively score a molecule for drug design.

        This is the main scoring function that evaluates a molecule across
        all relevant dimensions and returns an actionable score.

        Args:
            smiles: SMILES string of molecule
            target: Target protein (optional, for binding predictions)
            reference_compounds: Known actives/leads to compare against
            scoring_weights: Custom weights for score components
                Default: {"drug_likeness": 0.25, "safety": 0.25,
                         "synthesis": 0.20, "novelty": 0.15, "binding": 0.15}

        Returns:
            MoleculeScore with comprehensive evaluation
        """
        await self.initialize()

        # Default weights
        weights = scoring_weights or {
            "drug_likeness": 0.25,
            "safety": 0.25,
            "synthesis": 0.20,
            "novelty": 0.15,
            "binding": 0.15,
        }

        # Build inputs
        inputs = {"smiles": smiles}
        if target:
            inputs["target"] = target
        if reference_compounds:
            inputs["reference_smiles"] = reference_compounds

        # Define property groups
        drug_likeness_props = [
            "molecular_weight", "logp", "tpsa", "hbd", "hba",
            "rotatable_bonds", "lipinski_violations", "qed"
        ]

        safety_props = [
            "pains_alerts", "brenk_alerts", "herg_liability_estimate",
            "ames_estimate", "admet_safety_score"
        ]

        synthesis_props = [
            "synthetic_accessibility", "synthesis_feasibility_score"
        ]

        novelty_props = [
            "novelty_score", "reference_similarity",
            "pubchem_novelty", "database_novelty"
        ]

        binding_props = [
            "binding_affinity_estimate", "similarity_to_known_actives"
        ]

        # Collect all properties to predict
        all_props = (
            drug_likeness_props + safety_props + synthesis_props +
            novelty_props + (binding_props if target else [])
        )

        # Get available predictions
        available = self._registry.what_can_i_predict(inputs)
        to_predict = [p for p in all_props if p in available]

        # Make predictions
        predictions = await self._registry.predict_all(inputs, to_predict)

        # Calculate component scores
        drug_likeness_score = self._calculate_drug_likeness(predictions)
        safety_score = self._calculate_safety_score(predictions)
        synthesis_score = self._calculate_synthesis_score(predictions)
        novelty_score = self._calculate_novelty_score(predictions)
        binding_score = self._calculate_binding_score(predictions) if target else None

        # Calculate overall score
        if target and binding_score is not None:
            overall_score = (
                weights["drug_likeness"] * drug_likeness_score +
                weights["safety"] * safety_score +
                weights["synthesis"] * synthesis_score +
                weights["novelty"] * novelty_score +
                weights["binding"] * binding_score
            )
        else:
            # Redistribute binding weight to others
            adjusted_weights = {
                "drug_likeness": weights["drug_likeness"] / (1 - weights["binding"]),
                "safety": weights["safety"] / (1 - weights["binding"]),
                "synthesis": weights["synthesis"] / (1 - weights["binding"]),
                "novelty": weights["novelty"] / (1 - weights["binding"]),
            }
            overall_score = (
                adjusted_weights["drug_likeness"] * drug_likeness_score +
                adjusted_weights["safety"] * safety_score +
                adjusted_weights["synthesis"] * synthesis_score +
                adjusted_weights["novelty"] * novelty_score
            )

        # Determine overall confidence
        assessed_preds = [p for p in predictions.values() if p.assessed]
        if not assessed_preds:
            overall_confidence = ConfidenceLevel.NOT_ASSESSED
        else:
            min_conf = min(p.confidence for p in assessed_preds)
            avg_conf_score = sum(p.confidence_score or 0.5 for p in assessed_preds) / len(assessed_preds)

            if avg_conf_score >= 0.7:
                overall_confidence = ConfidenceLevel.MODERATE
            elif avg_conf_score >= 0.4:
                overall_confidence = ConfidenceLevel.LOW
            else:
                overall_confidence = min_conf

        # Collect issues and warnings
        critical_issues = []
        warnings = []

        for prop, pred in predictions.items():
            if not pred.assessed:
                continue

            for caveat in pred.caveats:
                if any(kw in caveat.lower() for kw in ["alert", "violation", "fail", "high risk"]):
                    critical_issues.append(caveat)
                elif any(kw in caveat.lower() for kw in ["concern", "warning", "caution"]):
                    warnings.append(caveat)

        # Check specific critical conditions
        if predictions.get("pains_alerts"):
            pains = predictions["pains_alerts"]
            if pains.assessed and pains.value and pains.value > 0:
                critical_issues.append(f"PAINS alerts detected: {pains.value}")

        if predictions.get("lipinski_violations"):
            lip = predictions["lipinski_violations"]
            if lip.assessed and lip.value and lip.value > 2:
                warnings.append(f"Multiple Lipinski violations: {lip.value}")

        # Categorize properties
        assessed = [p for p, pred in predictions.items() if pred.assessed]
        unassessed = [p for p, pred in predictions.items() if not pred.assessed]

        return MoleculeScore(
            smiles=smiles,
            overall_score=overall_score,
            overall_confidence=overall_confidence,
            drug_likeness_score=drug_likeness_score,
            safety_score=safety_score,
            synthesis_score=synthesis_score,
            novelty_score=novelty_score,
            binding_score=binding_score,
            predictions=predictions,
            critical_issues=critical_issues,
            warnings=warnings,
            target=target,
            assessed_properties=assessed,
            unassessed_properties=unassessed,
        )

    def _calculate_drug_likeness(self, predictions: Dict[str, Prediction]) -> float:
        """Calculate drug-likeness component score."""
        score = 1.0

        # QED is the primary drug-likeness metric
        qed = predictions.get("qed")
        if qed and qed.assessed and qed.value is not None:
            return float(qed.value)  # QED is already 0-1

        # Fallback to Lipinski-based score
        lip = predictions.get("lipinski_violations")
        if lip and lip.assessed and lip.value is not None:
            violations = int(lip.value)
            score -= violations * 0.15

        # Check individual properties
        mw = predictions.get("molecular_weight")
        if mw and mw.assessed and mw.value:
            if mw.value > 500:
                score -= 0.1
            if mw.value > 700:
                score -= 0.1

        logp = predictions.get("logp")
        if logp and logp.assessed and logp.value is not None:
            if logp.value > 5 or logp.value < -1:
                score -= 0.1

        return max(0.0, min(1.0, score))

    def _calculate_safety_score(self, predictions: Dict[str, Prediction]) -> float:
        """Calculate safety component score."""
        score = 1.0

        # ADMET safety score if available
        admet = predictions.get("admet_safety_score")
        if admet and admet.assessed and admet.value is not None:
            return float(admet.value)

        # Fallback to individual safety checks
        pains = predictions.get("pains_alerts")
        if pains and pains.assessed and pains.value:
            score -= 0.3 * min(pains.value, 3)  # Max penalty 0.9

        brenk = predictions.get("brenk_alerts")
        if brenk and brenk.assessed and brenk.value:
            score -= 0.15 * min(brenk.value, 4)

        herg = predictions.get("herg_liability_estimate")
        if herg and herg.assessed:
            if herg.value == "high":
                score -= 0.3
            elif herg.value == "moderate":
                score -= 0.15

        ames = predictions.get("ames_estimate")
        if ames and ames.assessed and ames.value:
            score -= 0.25

        return max(0.0, min(1.0, score))

    def _calculate_synthesis_score(self, predictions: Dict[str, Prediction]) -> float:
        """Calculate synthesis feasibility component score."""
        # Use normalized synthesis score if available
        synth = predictions.get("synthesis_feasibility_score")
        if synth and synth.assessed and synth.value is not None:
            return float(synth.value)

        # Fallback to SA score
        sa = predictions.get("synthetic_accessibility")
        if sa and sa.assessed and sa.value is not None:
            # SA Score 1-10, convert to 0-1 where 1 is easy
            return max(0, (10 - sa.value) / 9)

        return 0.5  # Default moderate if no synthesis prediction

    def _calculate_novelty_score(self, predictions: Dict[str, Prediction]) -> float:
        """Calculate novelty component score."""
        # Prefer explicit novelty score
        novelty = predictions.get("novelty_score")
        if novelty and novelty.assessed and novelty.value is not None:
            return float(novelty.value)

        # Use similarity as inverse
        sim = predictions.get("reference_similarity")
        if sim and sim.assessed and sim.value is not None:
            return 1.0 - float(sim.value)

        # Database novelty (boolean)
        db_novel = predictions.get("database_novelty")
        if db_novel and db_novel.assessed and db_novel.value is not None:
            return 0.8 if db_novel.value else 0.3

        return 0.5  # Default moderate

    def _calculate_binding_score(self, predictions: Dict[str, Prediction]) -> float:
        """Calculate binding affinity component score."""
        # Use similarity to known actives as proxy
        sim = predictions.get("similarity_to_known_actives")
        if sim and sim.assessed and sim.value is not None:
            # Higher similarity to actives = better binding estimate
            return float(sim.value)

        # Use binding estimate (usually pIC50-like)
        binding = predictions.get("binding_affinity_estimate")
        if binding and binding.assessed and binding.value is not None:
            # pIC50 5-10 is typical range, normalize
            normalized = (float(binding.value) - 5) / 5
            return max(0.0, min(1.0, normalized))

        return 0.5  # Default moderate

    async def batch_score(
        self,
        molecules: List[str],
        target: Optional[str] = None,
        reference_compounds: Optional[List[str]] = None,
        max_concurrent: int = 10,
    ) -> List[MoleculeScore]:
        """
        Score multiple molecules concurrently.

        Args:
            molecules: List of SMILES strings
            target: Target protein
            reference_compounds: Reference compounds
            max_concurrent: Maximum concurrent scoring operations

        Returns:
            List of MoleculeScore objects
        """
        semaphore = asyncio.Semaphore(max_concurrent)

        async def score_with_semaphore(smiles: str) -> MoleculeScore:
            async with semaphore:
                try:
                    return await self.score_molecule(
                        smiles=smiles,
                        target=target,
                        reference_compounds=reference_compounds,
                    )
                except Exception as e:
                    logger.error(f"Error scoring {smiles}: {e}")
                    # Return a minimal error score
                    return MoleculeScore(
                        smiles=smiles,
                        overall_score=0.0,
                        overall_confidence=ConfidenceLevel.NOT_ASSESSED,
                        drug_likeness_score=0.0,
                        safety_score=0.0,
                        synthesis_score=0.0,
                        novelty_score=0.0,
                        critical_issues=[f"Scoring failed: {str(e)}"],
                    )

        scores = await asyncio.gather(*[score_with_semaphore(s) for s in molecules])
        return list(scores)

    def list_available_predictions(self) -> List[str]:
        """List all properties that can be predicted."""
        return self._registry.list_predictable_properties()

    def list_predictors(self) -> List[Dict[str, Any]]:
        """List all registered predictors and their capabilities."""
        return self._registry.list_predictors()


# Global singleton instance
prediction_service = PredictionService()
