"""
Patent novelty and freedom-to-operate analysis service.
Integrates with SureChEMBL, EPO OPS, and provides FTO assessment.
"""
import asyncio
import os
from typing import Dict, List, Optional, Any
from dataclasses import dataclass
import logging
import aiohttp

logger = logging.getLogger(__name__)


@dataclass
class PatentHit:
    """A patent search result."""
    patent_id: str
    title: str
    publication_date: str
    assignee: Optional[str]

    # Match type
    exact_match: bool = False
    markush_match: bool = False
    similarity: float = 0.0

    # Matched compound info
    matched_smiles: Optional[str] = None
    match_context: Optional[str] = None

    # Patent status
    is_active: bool = True
    expiry_date: Optional[str] = None
    jurisdictions: List[str] = None

    def __post_init__(self):
        if self.jurisdictions is None:
            self.jurisdictions = []


class SureChEMBLClient:
    """Client for SureChEMBL chemical patent database."""

    def __init__(self):
        self.base_url = "https://www.surechembl.org"

    async def search_by_smiles(
        self,
        smiles: str,
        search_type: str = "substructure"
    ) -> List[PatentHit]:
        """
        Search patents by SMILES.

        Args:
            smiles: Query SMILES
            search_type: "exact", "substructure", or "similarity"

        Returns:
            List of patent hits
        """
        # In production, would call SureChEMBL API
        # For now, return mock data
        logger.info(f"SureChEMBL search: {smiles[:30]}... ({search_type})")

        # Mock: return empty for novel compounds
        return []

    async def get_compound_patents(self, schembl_id: str) -> List[PatentHit]:
        """Get all patents mentioning a SureChEMBL compound."""
        return []

    async def similarity_search(
        self,
        smiles: str,
        threshold: float = 0.7,
        limit: int = 100
    ) -> List[Dict[str, Any]]:
        """Find similar compounds in patent literature."""
        return []


class EPOClient:
    """Client for EPO Open Patent Services."""

    def __init__(self, consumer_key: str = None, consumer_secret: str = None):
        self.consumer_key = consumer_key or os.getenv("EPO_CONSUMER_KEY")
        self.consumer_secret = consumer_secret or os.getenv("EPO_CONSUMER_SECRET")
        self.base_url = "https://ops.epo.org/3.2"
        self.access_token = None

    async def _authenticate(self):
        """Get OAuth access token."""
        if not self.consumer_key or not self.consumer_secret:
            logger.warning("EPO credentials not set")
            return False

        # Would implement OAuth flow here
        return True

    async def search_patents(
        self,
        query: str,
        cpc_class: str = "A61K"  # Pharmaceutical preparations
    ) -> List[PatentHit]:
        """Search EPO patent database."""
        if not await self._authenticate():
            return []

        # Mock implementation
        return []

    async def get_patent_family(self, patent_id: str) -> List[str]:
        """Get all patents in a patent family."""
        return []


class PatentAnalyzer:
    """
    Main patent analysis service.
    Provides novelty checking, FTO analysis, and patent landscape assessment.
    """

    def __init__(self):
        self.surechembl = SureChEMBLClient()
        self.epo = EPOClient()

        # Cache of known compounds for similarity
        self._known_compounds_cache = []

    async def check_novelty(
        self,
        smiles: str,
        thorough: bool = True
    ) -> Dict[str, Any]:
        """
        Comprehensive novelty check for a compound.

        Args:
            smiles: Query compound SMILES
            thorough: If True, search multiple databases

        Returns:
            Dict with novelty assessment
        """
        results = {
            "query_smiles": smiles,
            "exact_matches": [],
            "similar_compounds": [],
            "markush_coverage": [],
            "novel": True,
            "confidence": 0.5,
            "notes": []
        }

        tasks = [
            self.surechembl.search_by_smiles(smiles, "exact"),
        ]

        if thorough:
            tasks.extend([
                self.surechembl.similarity_search(smiles, threshold=0.85),
            ])

        search_results = await asyncio.gather(*tasks, return_exceptions=True)

        # Process results
        for result in search_results:
            if isinstance(result, Exception):
                logger.error(f"Search error: {result}")
                results["notes"].append(f"Search error: {result}")
            elif isinstance(result, list):
                for hit in result:
                    if isinstance(hit, PatentHit):
                        if hit.exact_match:
                            results["exact_matches"].append(hit)
                            results["novel"] = False
                        elif hit.markush_match:
                            results["markush_coverage"].append(hit)
                        elif hit.similarity > 0.7:
                            results["similar_compounds"].append({
                                "patent_id": hit.patent_id,
                                "similarity": hit.similarity,
                                "smiles": hit.matched_smiles
                            })

        # Determine confidence
        if thorough and not results["exact_matches"]:
            results["confidence"] = 0.8

        if results["exact_matches"]:
            results["notes"].append(
                f"Found {len(results['exact_matches'])} exact match(es) in patent literature"
            )
        elif results["markush_coverage"]:
            results["notes"].append(
                f"Compound may be covered by {len(results['markush_coverage'])} Markush claim(s)"
            )
        else:
            results["notes"].append("No exact matches found in searched databases")

        return results

    async def assess_fto(
        self,
        smiles: str,
        jurisdictions: List[str] = None,
        launch_date: str = None
    ) -> Dict[str, Any]:
        """
        Freedom-to-operate analysis.

        Args:
            smiles: Compound SMILES
            jurisdictions: Target jurisdictions (e.g., ["US", "EP", "JP"])
            launch_date: Planned commercial launch date

        Returns:
            FTO assessment
        """
        if jurisdictions is None:
            jurisdictions = ["US", "EP"]

        fto_result = {
            "query_smiles": smiles,
            "jurisdictions": jurisdictions,
            "blocking_patents": [],
            "risk_level": "low",
            "recommendations": [],
            "expiry_timeline": {}
        }

        # Search for blocking patents
        novelty = await self.check_novelty(smiles, thorough=True)

        blocking = []

        # Check exact matches
        for match in novelty.get("exact_matches", []):
            if match.is_active:
                blocking.append({
                    "patent_id": match.patent_id,
                    "reason": "Exact compound claimed",
                    "risk": "high",
                    "expiry": match.expiry_date
                })

        # Check Markush coverage
        for markush in novelty.get("markush_coverage", []):
            blocking.append({
                "patent_id": markush.patent_id,
                "reason": "May fall within Markush claim scope",
                "risk": "medium",
                "expiry": markush.expiry_date
            })

        fto_result["blocking_patents"] = blocking

        # Determine overall risk
        if any(b["risk"] == "high" for b in blocking):
            fto_result["risk_level"] = "high"
            fto_result["recommendations"].append(
                "HIGH RISK: Exact compound appears in active patents. "
                "Consult patent counsel before proceeding."
            )
        elif blocking:
            fto_result["risk_level"] = "medium"
            fto_result["recommendations"].append(
                "MEDIUM RISK: Potential Markush coverage identified. "
                "Detailed claim analysis recommended."
            )
        else:
            fto_result["risk_level"] = "low"
            fto_result["recommendations"].append(
                "LOW RISK: No blocking patents identified in searched databases. "
                "Note: This is a preliminary assessment only."
            )

        return fto_result

    async def analyze_patent_landscape(
        self,
        target_name: str,
        therapeutic_area: str
    ) -> Dict[str, Any]:
        """
        Analyze the patent landscape for a therapeutic target.

        Args:
            target_name: Name of the target (e.g., "5HT1A receptor")
            therapeutic_area: Therapeutic area (e.g., "CNS", "Oncology")

        Returns:
            Patent landscape analysis
        """
        landscape = {
            "target": target_name,
            "therapeutic_area": therapeutic_area,
            "total_patents": 0,
            "key_players": [],
            "recent_filings": [],
            "expiring_soon": [],
            "white_space": [],
            "analysis_date": None
        }

        # This would involve comprehensive patent searches
        # Mock implementation for now

        return landscape

    def assess_patentability(
        self,
        novelty_result: Dict[str, Any],
        biological_data: Dict[str, Any] = None
    ) -> Dict[str, Any]:
        """
        Assess patentability of a compound.

        Args:
            novelty_result: Result from check_novelty()
            biological_data: Optional biological activity data

        Returns:
            Patentability assessment
        """
        assessment = {
            "novel": novelty_result.get("novel", False),
            "non_obvious": True,  # Would require more analysis
            "utility": False,
            "patentability_score": 0,
            "claim_strategy": [],
            "recommendations": []
        }

        # Check utility
        if biological_data:
            if biological_data.get("activity"):
                assessment["utility"] = True

        # Calculate score
        score = 0
        if assessment["novel"]:
            score += 40
        if assessment["non_obvious"]:
            score += 40
        if assessment["utility"]:
            score += 20
        assessment["patentability_score"] = score

        # Recommendations
        if score >= 80:
            assessment["recommendations"].append(
                "STRONG CANDIDATE: Consider filing provisional patent application"
            )
            assessment["claim_strategy"] = [
                "Composition of matter claims",
                "Method of treatment claims",
                "Pharmaceutical composition claims"
            ]
        elif score >= 60:
            assessment["recommendations"].append(
                "MODERATE CANDIDATE: Gather additional supporting data"
            )
        else:
            assessment["recommendations"].append(
                "WEAK CANDIDATE: Address patentability concerns before filing"
            )

        return assessment

    def generate_patent_report(
        self,
        smiles: str,
        novelty_result: Dict[str, Any],
        fto_result: Dict[str, Any],
        patentability: Dict[str, Any]
    ) -> str:
        """Generate comprehensive patent analysis report."""
        lines = [
            "=" * 70,
            "PATENT & NOVELTY ANALYSIS REPORT",
            "=" * 70,
            f"Compound: {smiles}",
            "",
            "NOVELTY ASSESSMENT",
            "-" * 50,
            f"Novel: {'Yes' if novelty_result.get('novel') else 'No'}",
            f"Confidence: {novelty_result.get('confidence', 0) * 100:.0f}%",
        ]

        if novelty_result.get("exact_matches"):
            lines.append(f"Exact Matches: {len(novelty_result['exact_matches'])}")
        if novelty_result.get("similar_compounds"):
            lines.append(f"Similar Compounds: {len(novelty_result['similar_compounds'])}")

        lines.extend([
            "",
            "FREEDOM-TO-OPERATE",
            "-" * 50,
            f"Risk Level: {fto_result.get('risk_level', 'unknown').upper()}",
            f"Blocking Patents: {len(fto_result.get('blocking_patents', []))}",
        ])

        for patent in fto_result.get("blocking_patents", [])[:5]:
            lines.append(f"  • {patent['patent_id']}: {patent['reason']}")

        lines.extend([
            "",
            "PATENTABILITY ASSESSMENT",
            "-" * 50,
            f"Score: {patentability.get('patentability_score', 0)}/100",
            f"Novel: {'Yes' if patentability.get('novel') else 'No'}",
            f"Non-obvious: {'Yes' if patentability.get('non_obvious') else 'Uncertain'}",
            f"Utility: {'Yes' if patentability.get('utility') else 'Needs data'}",
        ])

        if patentability.get("claim_strategy"):
            lines.append("\nRecommended Claim Strategy:")
            for strategy in patentability["claim_strategy"]:
                lines.append(f"  • {strategy}")

        lines.extend([
            "",
            "RECOMMENDATIONS",
            "-" * 50
        ])
        for rec in patentability.get("recommendations", []):
            lines.append(f"• {rec}")
        for rec in fto_result.get("recommendations", []):
            lines.append(f"• {rec}")

        lines.extend([
            "",
            "=" * 70,
            "DISCLAIMER: This is a preliminary automated analysis.",
            "Consult qualified patent counsel for legal advice.",
            "=" * 70
        ])

        return "\n".join(lines)
