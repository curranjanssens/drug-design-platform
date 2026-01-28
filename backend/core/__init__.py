"""
Core abstractions for the Drug Design Platform.

This module provides foundational types and patterns used throughout the system,
including the unified Prediction interface and Predictor registry.
"""

from .prediction import (
    Prediction,
    ConfidenceLevel,
    PredictionBasis,
    AggregatePrediction,
    aggregate_predictions,
)

from .predictor_registry import (
    Predictor,
    PredictorCapability,
    PredictorRegistry,
    registry,  # Global singleton
)

from .prediction_service import (
    PredictionService,
    MoleculeScore,
    prediction_service,  # Global singleton
)

from .target_knowledge import (
    TargetKnowledgeBase,
    TargetKnowledgeExtractor,
    TargetType,
    BindingMechanism,
    target_knowledge_extractor,  # Global singleton
)

from .honest_output import (
    HonestOutputFormatter,
    FormattedPrediction,
    FormattedScoreCard,
    OutputTone,
    honest_formatter,  # Global singleton
)

from .applicability_domain import (
    ApplicabilityDomainChecker,
    DomainCheckResult,
    DomainStatus,
    applicability_checker,  # Global singleton
)

from .design_report import (
    DesignReport,
    DesignInsight,
    DesignReportGenerator,
    Recommendation,
    report_generator,  # Global singleton
)

__all__ = [
    # Prediction types
    "Prediction",
    "ConfidenceLevel",
    "PredictionBasis",
    "AggregatePrediction",
    "aggregate_predictions",
    # Predictor types
    "Predictor",
    "PredictorCapability",
    "PredictorRegistry",
    "registry",
    # Prediction service
    "PredictionService",
    "MoleculeScore",
    "prediction_service",
    # Target knowledge
    "TargetKnowledgeBase",
    "TargetKnowledgeExtractor",
    "TargetType",
    "BindingMechanism",
    "target_knowledge_extractor",
    # Honest output
    "HonestOutputFormatter",
    "FormattedPrediction",
    "FormattedScoreCard",
    "OutputTone",
    "honest_formatter",
    # Applicability domain
    "ApplicabilityDomainChecker",
    "DomainCheckResult",
    "DomainStatus",
    "applicability_checker",
    # Design reports
    "DesignReport",
    "DesignInsight",
    "DesignReportGenerator",
    "Recommendation",
    "report_generator",
]
