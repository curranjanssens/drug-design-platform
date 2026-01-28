"""
Prediction Abstraction Layer (Task 1.1)

Unified interface for all predictions that carries uncertainty metadata.

Design Philosophy:
- Every prediction carries information about its reliability
- The system distinguishes between "confident," "uncertain," and "outside capability"
- New prediction sources can be added without changing downstream code
- Aggregation of multiple predictions is explicit and principled

Example usage:
    # A confident prediction
    pred = Prediction(
        value=0.85,
        confidence=ConfidenceLevel.HIGH,
        confidence_score=0.92,
        basis=PredictionBasis.VALIDATED_MODEL,
        caveats=[]
    )

    # An uncertain prediction
    pred = Prediction(
        value=0.5,
        confidence=ConfidenceLevel.LOW,
        confidence_score=0.3,
        basis=PredictionBasis.HEURISTIC,
        caveats=["No similar compounds in training data", "Extrapolating beyond known range"]
    )

    # Prediction that couldn't be made
    pred = Prediction.not_assessed("No ADMET predictor available for this property")
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import TypeVar, Generic, Optional, List, Dict, Any, Tuple, Union
from datetime import datetime


class ConfidenceLevel(Enum):
    """
    Categorical confidence levels for predictions.

    These map to semantic meanings that guide decision-making:
    - HIGH: Safe to use for primary decision-making
    - MODERATE: Use with awareness of uncertainty
    - LOW: Treat as rough estimate, verify experimentally
    - OUTSIDE_DOMAIN: Model extrapolating, prediction may be unreliable
    - NOT_ASSESSED: No prediction was made (missing capability)
    """
    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"
    OUTSIDE_DOMAIN = "outside_domain"
    NOT_ASSESSED = "not_assessed"

    def __lt__(self, other):
        """Enable sorting by confidence level."""
        order = [
            ConfidenceLevel.NOT_ASSESSED,
            ConfidenceLevel.OUTSIDE_DOMAIN,
            ConfidenceLevel.LOW,
            ConfidenceLevel.MODERATE,
            ConfidenceLevel.HIGH,
        ]
        return order.index(self) < order.index(other)


class PredictionBasis(Enum):
    """
    What the prediction is based on.

    This helps users understand WHERE the prediction came from,
    which informs how much to trust it.
    """
    # Model-based predictions
    VALIDATED_MODEL = "validated_model"      # Prospectively validated ML model
    ML_MODEL = "ml_model"                    # ML model without prospective validation
    DOCKING = "docking"                      # Structure-based docking
    PHARMACOPHORE = "pharmacophore"          # Pharmacophore matching

    # Rule/heuristic-based
    HEURISTIC = "heuristic"                  # Rule-based (Lipinski, PAINS, etc.)
    SIMILARITY = "similarity"               # Based on similar known compounds
    LITERATURE = "literature"               # From literature/database lookup

    # LLM-based
    LLM_REASONING = "llm_reasoning"          # Claude's chemical reasoning
    LLM_EXTRACTION = "llm_extraction"        # Extracted from LLM knowledge

    # Meta
    AGGREGATED = "aggregated"                # Combined from multiple sources
    UNKNOWN = "unknown"                      # Basis not specified
    NOT_APPLICABLE = "not_applicable"        # For NOT_ASSESSED predictions


T = TypeVar('T')


@dataclass
class Prediction(Generic[T]):
    """
    A prediction with uncertainty metadata.

    This is the core abstraction - every prediction in the system
    should be wrapped in this type.

    Attributes:
        value: The predicted value (any type)
        confidence: Categorical confidence level
        confidence_score: Continuous confidence (0.0-1.0), optional
        bounds: Lower and upper bounds if applicable
        basis: What method produced this prediction
        caveats: List of warnings/limitations
        assessed: Whether a prediction was actually made
        timestamp: When the prediction was made
        predictor_id: Which predictor made this (for tracing)
        raw_data: Original predictor output (for debugging)
        metadata: Additional predictor-specific data
    """
    value: T
    confidence: ConfidenceLevel = ConfidenceLevel.MODERATE
    confidence_score: Optional[float] = None
    bounds: Optional[Tuple[T, T]] = None
    basis: PredictionBasis = PredictionBasis.UNKNOWN
    caveats: List[str] = field(default_factory=list)
    assessed: bool = True
    timestamp: datetime = field(default_factory=datetime.utcnow)
    predictor_id: Optional[str] = None
    raw_data: Optional[Dict[str, Any]] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Validate and normalize prediction."""
        # Ensure confidence_score is in valid range
        if self.confidence_score is not None:
            self.confidence_score = max(0.0, min(1.0, self.confidence_score))

        # NOT_ASSESSED predictions should have assessed=False
        if self.confidence == ConfidenceLevel.NOT_ASSESSED:
            self.assessed = False

    @classmethod
    def not_assessed(cls, reason: str, predictor_id: Optional[str] = None) -> "Prediction[None]":
        """
        Create a prediction indicating the property wasn't assessed.

        Use this when:
        - No predictor is available for this property
        - The predictor couldn't run (missing data, error, etc.)
        - The prediction was intentionally skipped
        """
        return cls(
            value=None,
            confidence=ConfidenceLevel.NOT_ASSESSED,
            confidence_score=0.0,
            basis=PredictionBasis.NOT_APPLICABLE,
            caveats=[reason],
            assessed=False,
            predictor_id=predictor_id,
        )

    @classmethod
    def outside_domain(
        cls,
        value: T,
        reason: str,
        basis: PredictionBasis = PredictionBasis.UNKNOWN,
        predictor_id: Optional[str] = None
    ) -> "Prediction[T]":
        """
        Create a prediction for queries outside the model's domain.

        Use this when:
        - The query compound is dissimilar to training data
        - The prediction involves extrapolation
        - The model explicitly reports low confidence
        """
        return cls(
            value=value,
            confidence=ConfidenceLevel.OUTSIDE_DOMAIN,
            confidence_score=0.2,  # Low but not zero
            basis=basis,
            caveats=[f"Outside applicability domain: {reason}"],
            assessed=True,
            predictor_id=predictor_id,
        )

    @property
    def is_reliable(self) -> bool:
        """Check if prediction is reliable enough for primary decision-making."""
        return self.confidence in (ConfidenceLevel.HIGH, ConfidenceLevel.MODERATE)

    @property
    def needs_verification(self) -> bool:
        """Check if prediction should be experimentally verified."""
        return self.confidence in (ConfidenceLevel.LOW, ConfidenceLevel.OUTSIDE_DOMAIN)

    def with_caveat(self, caveat: str) -> "Prediction[T]":
        """Return a copy with an additional caveat."""
        new_caveats = self.caveats + [caveat]
        return Prediction(
            value=self.value,
            confidence=self.confidence,
            confidence_score=self.confidence_score,
            bounds=self.bounds,
            basis=self.basis,
            caveats=new_caveats,
            assessed=self.assessed,
            timestamp=self.timestamp,
            predictor_id=self.predictor_id,
            raw_data=self.raw_data,
            metadata=self.metadata,
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "value": self.value,
            "confidence": self.confidence.value,
            "confidence_score": self.confidence_score,
            "bounds": self.bounds,
            "basis": self.basis.value,
            "caveats": self.caveats,
            "assessed": self.assessed,
            "timestamp": self.timestamp.isoformat(),
            "predictor_id": self.predictor_id,
            "metadata": self.metadata,
        }

    def __repr__(self) -> str:
        if not self.assessed:
            return f"Prediction(NOT_ASSESSED: {self.caveats[0] if self.caveats else 'unknown reason'})"

        conf_str = f"{self.confidence.value}"
        if self.confidence_score is not None:
            conf_str += f" ({self.confidence_score:.0%})"

        caveat_str = f", {len(self.caveats)} caveats" if self.caveats else ""
        return f"Prediction({self.value}, {conf_str}, {self.basis.value}{caveat_str})"


@dataclass
class AggregatePrediction(Generic[T]):
    """
    A prediction aggregated from multiple sources.

    When multiple predictors assess the same property, this type
    captures the aggregate result along with information about
    agreement/disagreement.
    """
    value: T
    confidence: ConfidenceLevel
    confidence_score: float
    sources: List[Prediction[T]]
    agreement_score: float  # 0.0 = complete disagreement, 1.0 = complete agreement
    caveats: List[str] = field(default_factory=list)

    @property
    def source_count(self) -> int:
        """Number of sources that contributed."""
        return len([s for s in self.sources if s.assessed])

    @property
    def has_disagreement(self) -> bool:
        """Check if sources significantly disagree."""
        return self.agreement_score < 0.7

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "value": self.value,
            "confidence": self.confidence.value,
            "confidence_score": self.confidence_score,
            "source_count": self.source_count,
            "agreement_score": self.agreement_score,
            "caveats": self.caveats,
            "sources": [s.to_dict() for s in self.sources],
        }


def aggregate_predictions(
    predictions: List[Prediction[float]],
    method: str = "confidence_weighted"
) -> AggregatePrediction[float]:
    """
    Aggregate multiple numeric predictions into one.

    Methods:
    - "confidence_weighted": Weight by confidence score
    - "median": Take median value
    - "conservative": Take most pessimistic value

    Args:
        predictions: List of predictions to aggregate
        method: Aggregation method

    Returns:
        AggregatePrediction with combined result
    """
    # Filter to assessed predictions
    assessed = [p for p in predictions if p.assessed and p.value is not None]

    if not assessed:
        return AggregatePrediction(
            value=None,
            confidence=ConfidenceLevel.NOT_ASSESSED,
            confidence_score=0.0,
            sources=predictions,
            agreement_score=0.0,
            caveats=["No predictions available to aggregate"],
        )

    if len(assessed) == 1:
        p = assessed[0]
        return AggregatePrediction(
            value=p.value,
            confidence=p.confidence,
            confidence_score=p.confidence_score or 0.5,
            sources=predictions,
            agreement_score=1.0,
            caveats=p.caveats + ["Single source prediction"],
        )

    # Get values and weights
    values = [p.value for p in assessed]
    weights = [p.confidence_score or 0.5 for p in assessed]

    # Calculate aggregated value
    if method == "confidence_weighted":
        total_weight = sum(weights)
        agg_value = sum(v * w for v, w in zip(values, weights)) / total_weight
    elif method == "median":
        sorted_values = sorted(values)
        mid = len(sorted_values) // 2
        agg_value = sorted_values[mid] if len(sorted_values) % 2 else (sorted_values[mid-1] + sorted_values[mid]) / 2
    elif method == "conservative":
        # For scores where higher is better, take minimum
        # For risks where higher is worse, take maximum
        # Default to minimum (conservative for "good" predictions)
        agg_value = min(values)
    else:
        raise ValueError(f"Unknown aggregation method: {method}")

    # Calculate agreement (using coefficient of variation)
    mean_val = sum(values) / len(values)
    if mean_val != 0:
        variance = sum((v - mean_val) ** 2 for v in values) / len(values)
        std_dev = variance ** 0.5
        cv = std_dev / abs(mean_val)
        agreement = max(0.0, 1.0 - cv)  # Lower CV = higher agreement
    else:
        agreement = 1.0 if all(v == 0 for v in values) else 0.0

    # Determine aggregate confidence
    min_confidence = min(p.confidence for p in assessed)
    if agreement < 0.5:
        # Significant disagreement reduces confidence
        agg_confidence = ConfidenceLevel.LOW
    elif agreement < 0.7:
        agg_confidence = min(min_confidence, ConfidenceLevel.MODERATE)
    else:
        agg_confidence = min_confidence

    # Aggregate confidence score
    agg_conf_score = sum(weights) / len(weights) * agreement

    # Collect caveats
    all_caveats = []
    for p in assessed:
        all_caveats.extend(p.caveats)

    if agreement < 0.7:
        all_caveats.append(f"Predictors disagree (agreement={agreement:.0%})")

    return AggregatePrediction(
        value=agg_value,
        confidence=agg_confidence,
        confidence_score=agg_conf_score,
        sources=predictions,
        agreement_score=agreement,
        caveats=list(set(all_caveats)),  # Deduplicate
    )


# Type aliases for common prediction types
BindingPrediction = Prediction[float]  # Binding affinity (e.g., pKi, docking score)
PropertyPrediction = Prediction[float]  # Molecular property (e.g., LogP, MW)
BooleanPrediction = Prediction[bool]   # Yes/no prediction (e.g., BBB penetration)
CategoryPrediction = Prediction[str]   # Categorical prediction (e.g., toxicity class)
