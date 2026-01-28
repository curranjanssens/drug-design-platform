"""
Pluggable Predictor Registry (Task 1.2)

Allows prediction capabilities to be registered, discovered, and invoked dynamically.

Design Philosophy:
- Predictors register themselves with capability declarations
- The system can answer "what can I predict for this target?" at runtime
- Missing predictors degrade gracefully with explicit "not assessed" outputs
- Supports both sync and async predictors
- Handles expensive setup (model loading) via lazy initialization

Example usage:
    # Register a predictor
    @registry.register
    class LogPPredictor(Predictor):
        capabilities = [
            PredictorCapability(
                predicts="logp",
                requires=["smiles"],
                confidence_tier=ConfidenceLevel.HIGH,
            )
        ]

        async def predict(self, inputs: Dict) -> Prediction:
            smiles = inputs["smiles"]
            # ... calculate LogP ...
            return Prediction(value=logp, confidence=ConfidenceLevel.HIGH, ...)

    # Discover available predictions
    available = registry.what_can_i_predict({"smiles": "CCO", "target": "FAAH"})
    # Returns: ["logp", "molecular_weight", "binding_affinity", ...]

    # Make predictions
    predictions = await registry.predict_all(
        inputs={"smiles": "CCO"},
        properties=["logp", "tpsa"]
    )
"""

import asyncio
import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any, Set, Type, Callable, Union
from enum import Enum

from .prediction import Prediction, ConfidenceLevel, PredictionBasis

logger = logging.getLogger(__name__)


@dataclass
class PredictorCapability:
    """
    Declaration of what a predictor can predict.

    Attributes:
        predicts: What property this predicts (e.g., "logp", "binding_affinity")
        requires: What inputs are needed (e.g., ["smiles"], ["smiles", "target_structure"])
        confidence_tier: Expected confidence level for predictions
        description: Human-readable description
        tags: Searchable tags (e.g., ["admet", "toxicity"])
    """
    predicts: str
    requires: List[str]
    confidence_tier: ConfidenceLevel = ConfidenceLevel.MODERATE
    description: str = ""
    tags: List[str] = field(default_factory=list)

    def can_predict_with(self, available_inputs: Set[str]) -> bool:
        """Check if this capability can be used given available inputs."""
        return all(req in available_inputs for req in self.requires)


class Predictor(ABC):
    """
    Base class for all predictors.

    Subclasses must:
    1. Define `capabilities` class attribute
    2. Implement `predict()` method
    3. Optionally implement `initialize()` for expensive setup

    Example:
        class MyPredictor(Predictor):
            capabilities = [
                PredictorCapability(predicts="my_property", requires=["smiles"])
            ]

            async def predict(self, inputs: Dict, capability: PredictorCapability) -> Prediction:
                smiles = inputs["smiles"]
                value = calculate_something(smiles)
                return Prediction(value=value, confidence=ConfidenceLevel.MODERATE)
    """

    # Subclasses must define this
    capabilities: List[PredictorCapability] = []

    def __init__(self):
        self._initialized = False
        self._initialization_error: Optional[str] = None

    @property
    def predictor_id(self) -> str:
        """Unique identifier for this predictor."""
        return self.__class__.__name__

    async def initialize(self) -> None:
        """
        Expensive initialization (model loading, API auth, etc.)

        Called lazily on first prediction. Override in subclasses
        that need setup.
        """
        pass

    async def ensure_initialized(self) -> bool:
        """Ensure predictor is initialized. Returns False if initialization failed."""
        if self._initialized:
            return self._initialization_error is None

        try:
            await self.initialize()
            self._initialized = True
            return True
        except Exception as e:
            self._initialization_error = str(e)
            self._initialized = True
            logger.error(f"Predictor {self.predictor_id} initialization failed: {e}")
            return False

    @abstractmethod
    async def predict(
        self,
        inputs: Dict[str, Any],
        capability: PredictorCapability
    ) -> Prediction:
        """
        Make a prediction.

        Args:
            inputs: Dictionary of input values (e.g., {"smiles": "CCO"})
            capability: Which capability is being invoked

        Returns:
            Prediction with value and uncertainty metadata
        """
        pass

    async def safe_predict(
        self,
        inputs: Dict[str, Any],
        capability: PredictorCapability
    ) -> Prediction:
        """
        Safely make a prediction, handling errors gracefully.

        This is what the registry calls - it wraps predict() with
        error handling and initialization.
        """
        # Ensure initialized
        if not await self.ensure_initialized():
            return Prediction.not_assessed(
                f"Predictor initialization failed: {self._initialization_error}",
                predictor_id=self.predictor_id
            )

        # Check inputs
        missing = [r for r in capability.requires if r not in inputs]
        if missing:
            return Prediction.not_assessed(
                f"Missing required inputs: {missing}",
                predictor_id=self.predictor_id
            )

        # Make prediction
        try:
            prediction = await self.predict(inputs, capability)
            prediction.predictor_id = self.predictor_id
            return prediction
        except Exception as e:
            logger.error(f"Predictor {self.predictor_id} failed: {e}")
            return Prediction.not_assessed(
                f"Prediction failed: {str(e)}",
                predictor_id=self.predictor_id
            )


class PredictorRegistry:
    """
    Central registry for all predictors.

    Provides:
    - Registration of predictors
    - Discovery of available predictions
    - Invocation of predictions
    - Caching of expensive predictions
    """

    def __init__(self):
        self._predictors: Dict[str, Predictor] = {}
        self._capability_index: Dict[str, List[tuple[str, PredictorCapability]]] = {}
        self._cache: Dict[str, Prediction] = {}
        self._cache_enabled = True

    def register(self, predictor_class: Type[Predictor]) -> Type[Predictor]:
        """
        Register a predictor class.

        Can be used as a decorator:
            @registry.register
            class MyPredictor(Predictor):
                ...
        """
        instance = predictor_class()
        predictor_id = instance.predictor_id

        if predictor_id in self._predictors:
            logger.warning(f"Overwriting existing predictor: {predictor_id}")

        self._predictors[predictor_id] = instance

        # Index capabilities
        for cap in instance.capabilities:
            if cap.predicts not in self._capability_index:
                self._capability_index[cap.predicts] = []
            self._capability_index[cap.predicts].append((predictor_id, cap))

        logger.info(f"Registered predictor {predictor_id} with capabilities: {[c.predicts for c in instance.capabilities]}")

        return predictor_class

    def register_instance(self, predictor: Predictor) -> None:
        """Register an already-instantiated predictor."""
        predictor_id = predictor.predictor_id

        if predictor_id in self._predictors:
            logger.warning(f"Overwriting existing predictor: {predictor_id}")

        self._predictors[predictor_id] = predictor

        for cap in predictor.capabilities:
            if cap.predicts not in self._capability_index:
                self._capability_index[cap.predicts] = []
            self._capability_index[cap.predicts].append((predictor_id, cap))

    def what_can_i_predict(self, available_inputs: Dict[str, Any]) -> List[str]:
        """
        Given available inputs, what properties can be predicted?

        Args:
            available_inputs: Dict of available input values

        Returns:
            List of property names that can be predicted
        """
        input_keys = set(available_inputs.keys())
        predictable = []

        for property_name, predictors in self._capability_index.items():
            for predictor_id, cap in predictors:
                if cap.can_predict_with(input_keys):
                    predictable.append(property_name)
                    break  # Only need one predictor to make it predictable

        return predictable

    def get_predictors_for(
        self,
        property_name: str,
        available_inputs: Dict[str, Any]
    ) -> List[tuple[Predictor, PredictorCapability]]:
        """
        Get all predictors that can predict a property given inputs.

        Returns list of (predictor, capability) tuples, sorted by confidence tier.
        """
        if property_name not in self._capability_index:
            return []

        input_keys = set(available_inputs.keys())
        matches = []

        for predictor_id, cap in self._capability_index[property_name]:
            if cap.can_predict_with(input_keys):
                predictor = self._predictors[predictor_id]
                matches.append((predictor, cap))

        # Sort by confidence tier (highest first)
        matches.sort(key=lambda x: x[1].confidence_tier, reverse=True)

        return matches

    def _cache_key(self, predictor_id: str, property_name: str, inputs: Dict) -> str:
        """Generate cache key for a prediction."""
        # Use sorted items for consistent key
        input_str = str(sorted((k, str(v)[:100]) for k, v in inputs.items()))
        return f"{predictor_id}:{property_name}:{hash(input_str)}"

    async def predict(
        self,
        inputs: Dict[str, Any],
        property_name: str,
        predictor_id: Optional[str] = None,
        use_cache: bool = True
    ) -> Prediction:
        """
        Make a single prediction.

        Args:
            inputs: Input values (e.g., {"smiles": "CCO"})
            property_name: What to predict (e.g., "logp")
            predictor_id: Optionally specify which predictor to use
            use_cache: Whether to use cached results

        Returns:
            Prediction result
        """
        # Find predictor
        if predictor_id:
            if predictor_id not in self._predictors:
                return Prediction.not_assessed(f"Unknown predictor: {predictor_id}")
            predictor = self._predictors[predictor_id]
            caps = [c for c in predictor.capabilities if c.predicts == property_name]
            if not caps:
                return Prediction.not_assessed(
                    f"Predictor {predictor_id} does not predict {property_name}"
                )
            capability = caps[0]
        else:
            matches = self.get_predictors_for(property_name, inputs)
            if not matches:
                return Prediction.not_assessed(
                    f"No predictor available for {property_name} with given inputs"
                )
            predictor, capability = matches[0]  # Use highest confidence predictor

        # Check cache
        cache_key = self._cache_key(predictor.predictor_id, property_name, inputs)
        if use_cache and self._cache_enabled and cache_key in self._cache:
            cached = self._cache[cache_key]
            cached.metadata["cached"] = True
            return cached

        # Make prediction
        prediction = await predictor.safe_predict(inputs, capability)

        # Cache result
        if use_cache and self._cache_enabled and prediction.assessed:
            self._cache[cache_key] = prediction

        return prediction

    async def predict_all(
        self,
        inputs: Dict[str, Any],
        properties: Optional[List[str]] = None,
        use_cache: bool = True
    ) -> Dict[str, Prediction]:
        """
        Predict multiple properties in parallel.

        Args:
            inputs: Input values
            properties: List of properties to predict (None = all possible)
            use_cache: Whether to use cached results

        Returns:
            Dict mapping property name to Prediction
        """
        if properties is None:
            properties = self.what_can_i_predict(inputs)

        # Run predictions in parallel
        tasks = [
            self.predict(inputs, prop, use_cache=use_cache)
            for prop in properties
        ]

        results = await asyncio.gather(*tasks, return_exceptions=True)

        output = {}
        for prop, result in zip(properties, results):
            if isinstance(result, Exception):
                output[prop] = Prediction.not_assessed(f"Prediction error: {str(result)}")
            else:
                output[prop] = result

        return output

    async def predict_with_all_predictors(
        self,
        inputs: Dict[str, Any],
        property_name: str
    ) -> List[Prediction]:
        """
        Get predictions from ALL available predictors for a property.

        Useful when you want to aggregate or compare predictions
        from multiple sources.
        """
        matches = self.get_predictors_for(property_name, inputs)

        if not matches:
            return [Prediction.not_assessed(f"No predictor for {property_name}")]

        tasks = [
            predictor.safe_predict(inputs, cap)
            for predictor, cap in matches
        ]

        return await asyncio.gather(*tasks)

    def clear_cache(self) -> None:
        """Clear the prediction cache."""
        self._cache.clear()

    def list_predictors(self) -> List[Dict[str, Any]]:
        """List all registered predictors and their capabilities."""
        return [
            {
                "id": predictor_id,
                "capabilities": [
                    {
                        "predicts": cap.predicts,
                        "requires": cap.requires,
                        "confidence_tier": cap.confidence_tier.value,
                        "description": cap.description,
                        "tags": cap.tags,
                    }
                    for cap in predictor.capabilities
                ]
            }
            for predictor_id, predictor in self._predictors.items()
        ]

    def list_predictable_properties(self) -> List[str]:
        """List all properties that can be predicted."""
        return list(self._capability_index.keys())


# Global singleton registry
registry = PredictorRegistry()
