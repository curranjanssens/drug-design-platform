"""
Retrosynthesis planning service.
Integrates with IBM RXN, AiZynthFinder, and ASKCOS.
"""
import asyncio
import os
from typing import Dict, List, Optional, Any
from dataclasses import dataclass
import json
import logging

from models import SynthesisRoute, ReactionStep

logger = logging.getLogger(__name__)


class IBMRXNClient:
    """Client for IBM RXN for Chemistry API."""

    def __init__(self, api_key: str = None):
        self.api_key = api_key or os.getenv("IBM_RXN_API_KEY")
        self.base_url = "https://rxn.res.ibm.com/rxn/api/api/v1"
        self.project_id = None

    async def _ensure_project(self):
        """Ensure we have a project created."""
        if self.project_id:
            return

        # In real implementation, create project via API
        # For now, assume project exists
        self.project_id = "drug_design_platform"

    async def predict_retrosynthesis(
        self,
        smiles: str,
        max_steps: int = 6,
        availability_pricing: bool = True
    ) -> List[SynthesisRoute]:
        """
        Predict retrosynthesis routes using IBM RXN.

        Args:
            smiles: Target molecule SMILES
            max_steps: Maximum synthesis steps
            availability_pricing: Whether to include pricing info

        Returns:
            List of synthesis routes
        """
        if not self.api_key:
            logger.warning("IBM RXN API key not set, returning mock data")
            return self._mock_routes(smiles)

        await self._ensure_project()

        try:
            # This would be the actual API call
            # response = await self._call_api(
            #     "POST",
            #     f"/retrosynthesis/auto/{self.project_id}",
            #     json={"product": smiles, "maxSteps": max_steps}
            # )

            # For now, return mock data
            return self._mock_routes(smiles)

        except Exception as e:
            logger.error(f"IBM RXN API error: {e}")
            return []

    def _mock_routes(self, smiles: str) -> List[SynthesisRoute]:
        """Generate mock synthesis routes for development."""
        route = SynthesisRoute(
            steps=[
                ReactionStep(
                    step_number=1,
                    reaction_smiles=f"A.B>>{smiles[:10]}...",
                    reaction_name="Suzuki Coupling",
                    reactants=["Cc1ccc(B(O)O)cc1", "Brc1ccccc1"],
                    reagents=["Pd(PPh3)4", "K2CO3"],
                    solvent="THF/H2O",
                    temperature="80°C",
                    time="12h",
                    atmosphere="N2",
                    expected_yield=0.85,
                    confidence=0.9,
                    reference="Mock reference"
                )
            ],
            starting_materials=[
                {
                    "smiles": "Cc1ccc(B(O)O)cc1",
                    "name": "4-Methylphenylboronic acid",
                    "vendor": "Sigma-Aldrich",
                    "cas": "5720-05-8",
                    "price_per_g": 35.0,
                    "available": True
                }
            ],
            source="ibm_rxn"
        )
        route.calculate_metrics()
        return [route]


class AiZynthFinderWrapper:
    """Wrapper for AiZynthFinder retrosynthesis tool."""

    def __init__(self, config_path: str = None):
        self.config_path = config_path
        self.finder = None

    def _init_finder(self):
        """Initialize AiZynthFinder (lazy loading)."""
        if self.finder is not None:
            return

        try:
            from aizynthfinder.aizynthfinder import AiZynthFinder
            self.finder = AiZynthFinder(configfile=self.config_path)
        except ImportError:
            logger.warning("AiZynthFinder not installed")
            self.finder = None

    async def predict_retrosynthesis(
        self,
        smiles: str,
        time_limit: int = 120
    ) -> List[SynthesisRoute]:
        """
        Predict retrosynthesis using AiZynthFinder.

        Args:
            smiles: Target molecule SMILES
            time_limit: Time limit in seconds

        Returns:
            List of synthesis routes
        """
        self._init_finder()

        if self.finder is None:
            logger.warning("AiZynthFinder not available, returning empty")
            return []

        try:
            self.finder.config.search.time_limit = time_limit
            self.finder.target_smiles = smiles

            # Run search in thread pool to not block
            loop = asyncio.get_event_loop()
            await loop.run_in_executor(None, self.finder.tree_search)
            self.finder.build_routes()

            # Convert to our format
            routes = []
            for route_dict in self.finder.routes.dicts:
                route = self._convert_route(route_dict)
                if route:
                    routes.append(route)

            return routes

        except Exception as e:
            logger.error(f"AiZynthFinder error: {e}")
            return []

    def _convert_route(self, route_dict: Dict) -> Optional[SynthesisRoute]:
        """Convert AiZynthFinder route to our format."""
        try:
            steps = []
            for i, reaction in enumerate(route_dict.get("reactions", [])):
                step = ReactionStep(
                    step_number=i + 1,
                    reaction_smiles=reaction.get("reaction_smiles", ""),
                    reaction_name=reaction.get("template_name", "Unknown"),
                    reactants=reaction.get("reactants", []),
                    reagents=[],
                    expected_yield=reaction.get("score", 0.7),
                    confidence=reaction.get("confidence", 0.5)
                )
                steps.append(step)

            route = SynthesisRoute(
                steps=steps,
                starting_materials=[
                    {"smiles": sm, "available": True}
                    for sm in route_dict.get("starting_materials", [])
                ],
                source="aizynthfinder"
            )
            route.calculate_metrics()
            return route

        except Exception as e:
            logger.error(f"Error converting route: {e}")
            return None


class RetrosynthesisService:
    """
    Main retrosynthesis service that orchestrates multiple backends.
    """

    def __init__(self):
        self.ibm_rxn = IBMRXNClient()
        self.aizynthfinder = AiZynthFinderWrapper()

    async def plan_synthesis(
        self,
        smiles: str,
        max_steps: int = 6,
        use_ibm_rxn: bool = True,
        use_aizynthfinder: bool = True,
        time_limit: int = 120
    ) -> List[SynthesisRoute]:
        """
        Plan synthesis routes using multiple backends.

        Args:
            smiles: Target molecule SMILES
            max_steps: Maximum synthesis steps
            use_ibm_rxn: Use IBM RXN backend
            use_aizynthfinder: Use AiZynthFinder backend
            time_limit: Time limit per backend

        Returns:
            List of synthesis routes from all backends, deduplicated and ranked
        """
        all_routes = []
        tasks = []

        if use_ibm_rxn:
            tasks.append(self.ibm_rxn.predict_retrosynthesis(smiles, max_steps))

        if use_aizynthfinder:
            tasks.append(self.aizynthfinder.predict_retrosynthesis(smiles, time_limit))

        if tasks:
            results = await asyncio.gather(*tasks, return_exceptions=True)

            for result in results:
                if isinstance(result, Exception):
                    logger.error(f"Retrosynthesis backend error: {result}")
                elif isinstance(result, list):
                    all_routes.extend(result)

        # Rank routes
        ranked_routes = self._rank_routes(all_routes)

        return ranked_routes

    def _rank_routes(self, routes: List[SynthesisRoute]) -> List[SynthesisRoute]:
        """Rank synthesis routes by quality metrics."""
        if not routes:
            return []

        # Score each route
        scored = []
        for route in routes:
            score = self._calculate_route_score(route)
            scored.append((route, score))

        # Sort by score descending
        scored.sort(key=lambda x: x[1], reverse=True)

        return [route for route, _ in scored]

    def _calculate_route_score(self, route: SynthesisRoute) -> float:
        """Calculate composite score for a route."""
        score = 0.0

        # Confidence (0-1) * 40%
        score += route.confidence * 0.4

        # Yield (0-1) * 30%
        score += route.overall_yield * 0.3

        # Fewer steps is better (max 20%)
        steps_score = max(0, 1 - (route.total_steps - 1) * 0.1)
        score += steps_score * 0.2

        # Starting material availability (10%)
        available_count = sum(
            1 for sm in route.starting_materials
            if sm.get("available", False)
        )
        if route.starting_materials:
            availability = available_count / len(route.starting_materials)
            score += availability * 0.1

        return score

    async def get_starting_material_info(
        self,
        smiles: str
    ) -> Dict[str, Any]:
        """
        Get pricing and availability for a starting material.

        Would integrate with:
        - ZINC database
        - eMolecules
        - Sigma-Aldrich/MilliporeSigma API
        - ChemSpace
        """
        # Mock implementation
        return {
            "smiles": smiles,
            "available": True,
            "vendors": [
                {
                    "name": "Sigma-Aldrich",
                    "price_per_g": 45.0,
                    "currency": "USD",
                    "catalog_number": "MOCK-12345",
                    "purity": "95%"
                }
            ],
            "estimated_delivery": "3-5 business days"
        }

    def format_synthesis_report(self, route: SynthesisRoute) -> str:
        """Generate a human-readable synthesis report."""
        lines = [
            "=" * 60,
            "SYNTHESIS ROUTE REPORT",
            "=" * 60,
            f"Total Steps: {route.total_steps}",
            f"Estimated Overall Yield: {route.overall_yield * 100:.1f}%",
            f"Route Confidence: {route.confidence * 100:.1f}%",
            f"Source: {route.source}",
            "",
            "STARTING MATERIALS:",
            "-" * 40
        ]

        for sm in route.starting_materials:
            lines.append(f"  • {sm.get('name', 'Unknown')} ({sm['smiles']})")
            if sm.get("vendor"):
                lines.append(f"    Vendor: {sm['vendor']} @ ${sm.get('price_per_g', 'N/A')}/g")

        lines.append("")
        lines.append("REACTION STEPS:")
        lines.append("-" * 40)

        for step in route.steps:
            lines.append(f"\nStep {step.step_number}: {step.reaction_name}")
            lines.append(f"  Reaction: {step.reaction_smiles}")
            if step.reagents:
                lines.append(f"  Reagents: {', '.join(step.reagents)}")
            if step.solvent:
                lines.append(f"  Solvent: {step.solvent}")
            if step.temperature:
                lines.append(f"  Temperature: {step.temperature}")
            if step.time:
                lines.append(f"  Time: {step.time}")
            lines.append(f"  Expected Yield: {step.expected_yield * 100:.0f}%")
            lines.append(f"  Confidence: {step.confidence * 100:.0f}%")

        lines.append("")
        lines.append("=" * 60)

        return "\n".join(lines)
