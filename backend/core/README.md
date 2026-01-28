# Drug Design Platform - Core Prediction Infrastructure

This module implements the "Honest Scoring Infrastructure" from the Implementation Roadmap.

## Philosophy

- **Honest Uncertainty**: Every prediction carries information about its reliability
- **Target-Agnostic Design**: The system works for any drug target without hardcoded classes
- **Plugin Architecture**: New predictors can be added without changing core code
- **Graceful Degradation**: Missing predictors degrade to "not assessed" instead of crashing

## Components

### 1. Prediction Abstraction (`prediction.py`)

The core `Prediction[T]` type that wraps all predictions with uncertainty metadata:

```python
@dataclass
class Prediction(Generic[T]):
    value: T                              # The predicted value
    confidence: ConfidenceLevel           # HIGH, MODERATE, LOW, OUTSIDE_DOMAIN, NOT_ASSESSED
    confidence_score: Optional[float]     # 0.0 - 1.0
    basis: PredictionBasis               # How prediction was made
    caveats: List[str]                   # Warnings and limitations
    assessed: bool                       # Whether prediction was made
```

### 2. Predictor Registry (`predictor_registry.py`)

Plugin system for prediction capabilities:

```python
# Register a predictor
@registry.register
class MyPredictor(Predictor):
    capabilities = [
        PredictorCapability(
            predicts="my_property",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.MODERATE,
        )
    ]

    async def predict(self, inputs, capability) -> Prediction:
        # ... make prediction
        return Prediction(value=result, confidence=ConfidenceLevel.MODERATE)

# Discover what's possible
available = registry.what_can_i_predict({"smiles": "CCO"})

# Make predictions
predictions = await registry.predict_all(inputs, properties)
```

### 3. Prediction Service (`prediction_service.py`)

High-level interface for comprehensive molecule scoring:

```python
# Score a single molecule
score = await prediction_service.score_molecule(
    smiles="CCO",
    target="FAAH",
    reference_compounds=["..."]
)

print(f"Overall: {score.overall_score}")
print(f"Confidence: {score.overall_confidence.value}")
print(f"Issues: {score.critical_issues}")
```

### 4. Predictors (`predictors/`)

Concrete predictor implementations:

| File | Predictors | Properties |
|------|-----------|------------|
| `property_predictors.py` | RDKitPropertyPredictor, DrugLikenessPredictor | MW, LogP, TPSA, HBD, HBA, QED, PAINS, etc. |
| `admet_predictors.py` | RuleBasedADMETPredictor | Oral absorption, BBB, hERG, Ames, etc. |
| `synthesis_predictors.py` | SASynthesisPredictor | Synthetic accessibility, feasibility |
| `binding_predictors.py` | SimilarityBindingPredictor | Binding affinity, similarity to actives |
| `novelty_predictors.py` | FingerprintNoveltyPredictor, DatabaseNoveltyPredictor | Novelty, database checks |

### 5. Target Knowledge (`target_knowledge.py`)

Dynamic extraction of target-specific knowledge:

```python
knowledge = await target_knowledge_extractor.extract_knowledge(
    "Design FAAH inhibitors for pain treatment"
)

print(f"Target: {knowledge.target_name}")
print(f"Mechanism: {knowledge.binding_mechanism}")
print(f"Is covalent: {knowledge.is_covalent}")
```

### 6. Honest Output (`honest_output.py`)

Formats predictions for clear, honest communication:

```python
score_card = honest_formatter.format_score_card(score)

print(f"Grade: {score_card.overall_grade}")
print(f"Summary: {score_card.overall_summary}")
print(f"Recommendations: {score_card.recommendations}")
```

### 7. Applicability Domain (`applicability_domain.py`)

Checks if molecules are within reliable prediction space:

```python
result = applicability_checker.check_domain(smiles)

if result.status == DomainStatus.OUTSIDE:
    print("Predictions may be unreliable")
    print(f"Confidence modifier: {result.confidence_modifier}")
```

## API Endpoints (v2)

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/v2/score` | POST | Comprehensive molecule scoring |
| `/api/v2/score/batch` | POST | Batch scoring with summary |
| `/api/v2/predict` | POST | Specific property predictions |
| `/api/v2/capabilities` | GET | List available predictions |
| `/api/v2/domain-check` | POST | Check applicability domain |

## Confidence Levels

| Level | Meaning | Action |
|-------|---------|--------|
| HIGH | Reliable calculation | Use for decisions |
| MODERATE | Reasonable estimate | Consider verification |
| LOW | Rough heuristic | Verify experimentally |
| OUTSIDE_DOMAIN | Model extrapolating | Treat with caution |
| NOT_ASSESSED | No prediction made | Cannot evaluate |

## Prediction Basis

| Basis | Description |
|-------|-------------|
| HEURISTIC | Rule-based calculation (Lipinski, etc.) |
| VALIDATED_MODEL | Prospectively validated ML model |
| ML_MODEL | ML model without prospective validation |
| SIMILARITY | Based on similar known compounds |
| LITERATURE | From database/literature |
| LLM_REASONING | Claude's chemical reasoning |
| DOCKING | Structure-based docking |

## Usage Example

```python
import asyncio
from core import prediction_service, honest_formatter, applicability_checker

async def evaluate_molecule(smiles: str, target: str = None):
    # Check if predictions will be reliable
    domain = applicability_checker.check_domain(smiles)

    # Get comprehensive score
    score = await prediction_service.score_molecule(
        smiles=smiles,
        target=target
    )

    # Format for display
    card = honest_formatter.format_score_card(score)

    return {
        "grade": card.overall_grade,
        "summary": card.overall_summary,
        "confidence": card.overall_confidence,
        "issues": card.critical_issues,
        "recommendations": card.recommendations,
        "domain_status": domain.status.value,
    }

# Run
result = asyncio.run(evaluate_molecule("CC(=O)OC1=CC=CC=C1C(=O)O"))
```

## Adding New Predictors

1. Create predictor class extending `Predictor`
2. Define `capabilities` list
3. Implement `async def predict()` method
4. Decorate with `@registry.register`

```python
from core import Predictor, PredictorCapability, ConfidenceLevel, registry

@registry.register
class MyNewPredictor(Predictor):
    capabilities = [
        PredictorCapability(
            predicts="my_new_property",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.MODERATE,
            description="My new prediction",
            tags=["custom"],
        )
    ]

    async def predict(self, inputs, capability):
        smiles = inputs["smiles"]
        # ... calculate prediction
        return Prediction(
            value=result,
            confidence=ConfidenceLevel.MODERATE,
            confidence_score=0.7,
            basis=PredictionBasis.HEURISTIC,
            caveats=["Add relevant caveats"],
        )
```

The predictor is automatically discovered and available through the registry and prediction service.
