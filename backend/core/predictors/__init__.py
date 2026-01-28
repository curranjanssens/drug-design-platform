"""
Concrete predictor implementations.

This module contains all predictor implementations organized by category:
- properties: Molecular property predictors (LogP, MW, TPSA, etc.)
- admet: ADMET predictions (absorption, toxicity, etc.)
- binding: Binding affinity predictions (docking, similarity-based)
- novelty: Novelty and patentability assessment
- synthesis: Synthesis feasibility predictions

All predictors are automatically registered with the global registry
when this module is imported.
"""

from .property_predictors import (
    RDKitPropertyPredictor,
    DrugLikenessPredictor,
)

from .novelty_predictors import (
    FingerprintNoveltyPredictor,
    DatabaseNoveltyPredictor,
)

from .admet_predictors import (
    RuleBasedADMETPredictor,
    # ModelBasedADMETPredictor,  # TODO: Add when we have real models
)

from .binding_predictors import (
    SimilarityBindingPredictor,
)

from .docking_predictor import (
    DiffDockPredictor,
    dock_molecule,
)

from .synthesis_predictors import (
    SASynthesisPredictor,
)

__all__ = [
    "RDKitPropertyPredictor",
    "DrugLikenessPredictor",
    "FingerprintNoveltyPredictor",
    "DatabaseNoveltyPredictor",
    "RuleBasedADMETPredictor",
    "SimilarityBindingPredictor",
    "DiffDockPredictor",
    "dock_molecule",
    "SASynthesisPredictor",
]
