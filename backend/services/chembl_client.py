"""
ChEMBL API Client

Fetches real bioactivity data (IC50, Ki, EC50) for targets and compounds.
This is critical for learning from actual SAR data.
"""
import asyncio
import logging
from typing import Dict, List, Optional, Any
from dataclasses import dataclass
import httpx

logger = logging.getLogger(__name__)

CHEMBL_API_BASE = "https://www.ebi.ac.uk/chembl/api/data"


@dataclass
class BioactivityData:
    """Bioactivity measurement from ChEMBL."""
    chembl_id: str
    compound_name: str
    smiles: str
    activity_type: str  # IC50, Ki, EC50, etc.
    value: float  # in nM
    units: str
    assay_description: str
    target_name: str


@dataclass
class TargetData:
    """Target information from ChEMBL."""
    chembl_id: str
    name: str
    target_type: str
    organism: str
    uniprot_id: Optional[str]


class ChEMBLClient:
    """Client for ChEMBL API - the largest public medicinal chemistry database."""

    def __init__(self):
        self.base_url = CHEMBL_API_BASE
        self.timeout = 30.0

    async def search_target(self, query: str) -> List[TargetData]:
        """Search for a target by name."""
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                # Search targets
                url = f"{self.base_url}/target/search.json"
                params = {"q": query, "limit": 10}

                response = await client.get(url, params=params)

                if response.status_code != 200:
                    logger.warning(f"ChEMBL target search failed: {response.status_code}")
                    return []

                data = response.json()
                targets = []

                for t in data.get("targets", []):
                    targets.append(TargetData(
                        chembl_id=t.get("target_chembl_id", ""),
                        name=t.get("pref_name", ""),
                        target_type=t.get("target_type", ""),
                        organism=t.get("organism", ""),
                        uniprot_id=t.get("target_components", [{}])[0].get("accession") if t.get("target_components") else None
                    ))

                logger.info(f"Found {len(targets)} targets for '{query}'")
                return targets

        except Exception as e:
            logger.error(f"ChEMBL target search error: {e}")
            return []

    async def get_bioactivities_for_target(
        self,
        target_chembl_id: str,
        activity_types: List[str] = None,
        max_results: int = 100
    ) -> List[BioactivityData]:
        """Get bioactivity data for a target (IC50, Ki, etc.)."""
        if activity_types is None:
            activity_types = ["IC50", "Ki", "EC50", "Kd"]

        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                url = f"{self.base_url}/activity.json"
                params = {
                    "target_chembl_id": target_chembl_id,
                    "standard_type__in": ",".join(activity_types),
                    "limit": max_results,
                    "order_by": "standard_value"  # Best binders first
                }

                response = await client.get(url, params=params)

                if response.status_code != 200:
                    logger.warning(f"ChEMBL activity fetch failed: {response.status_code}")
                    return []

                data = response.json()
                activities = []

                for a in data.get("activities", []):
                    # Only include entries with SMILES and valid values
                    smiles = a.get("canonical_smiles")
                    value = a.get("standard_value")

                    if smiles and value is not None:
                        try:
                            float_value = float(value)
                            # Filter out invalid values: 0, negative, or unrealistically small
                            # Real IC50/Ki values are typically > 0.001 nM
                            if float_value <= 0.001:
                                logger.debug(f"Skipping invalid bioactivity value: {float_value} nM")
                                continue

                            activities.append(BioactivityData(
                                chembl_id=a.get("molecule_chembl_id", ""),
                                compound_name=a.get("molecule_pref_name", "Unknown"),
                                smiles=smiles,
                                activity_type=a.get("standard_type", ""),
                                value=float_value,
                                units=a.get("standard_units", "nM"),
                                assay_description=a.get("assay_description", ""),
                                target_name=a.get("target_pref_name", "")
                            ))
                        except (ValueError, TypeError):
                            continue

                # Sort by potency (lower is better for IC50/Ki)
                activities.sort(key=lambda x: x.value)

                logger.info(f"Found {len(activities)} bioactivities for {target_chembl_id}")
                return activities

        except Exception as e:
            logger.error(f"ChEMBL activity fetch error: {e}")
            return []

    async def get_compound_info(self, smiles_or_name: str) -> Optional[Dict]:
        """Get compound information from ChEMBL."""
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                # Try compound search
                url = f"{self.base_url}/molecule/search.json"
                params = {"q": smiles_or_name, "limit": 5}

                response = await client.get(url, params=params)

                if response.status_code == 200:
                    data = response.json()
                    molecules = data.get("molecules", [])

                    if molecules:
                        mol = molecules[0]
                        return {
                            "chembl_id": mol.get("molecule_chembl_id"),
                            "name": mol.get("pref_name"),
                            "smiles": mol.get("molecule_structures", {}).get("canonical_smiles"),
                            "max_phase": mol.get("max_phase"),  # Clinical phase
                            "molecular_weight": mol.get("molecule_properties", {}).get("full_mwt"),
                            "alogp": mol.get("molecule_properties", {}).get("alogp"),
                        }

                return None

        except Exception as e:
            logger.error(f"ChEMBL compound lookup error: {e}")
            return None

    async def get_similar_compounds(
        self,
        smiles: str,
        similarity: float = 0.7,
        max_results: int = 20
    ) -> List[Dict]:
        """Find similar compounds in ChEMBL (for novelty checking)."""
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                # Similarity search endpoint
                url = f"{self.base_url}/similarity/{smiles}/{int(similarity * 100)}.json"
                params = {"limit": max_results}

                response = await client.get(url, params=params)

                if response.status_code == 200:
                    data = response.json()
                    similar = []

                    for mol in data.get("molecules", []):
                        similar.append({
                            "chembl_id": mol.get("molecule_chembl_id"),
                            "name": mol.get("pref_name"),
                            "smiles": mol.get("molecule_structures", {}).get("canonical_smiles"),
                            "similarity": mol.get("similarity"),
                            "max_phase": mol.get("max_phase")
                        })

                    return similar

                return []

        except Exception as e:
            logger.error(f"ChEMBL similarity search error: {e}")
            return []

    async def get_target_with_best_binders(self, target_name: str) -> Dict:
        """
        Full pipeline: search target, get bioactivities, return structured data.
        This is what the LLM should call to understand a target.
        """
        result = {
            "target_found": False,
            "target_info": None,
            "known_binders": [],
            "sar_insights": []
        }

        # Search for target
        targets = await self.search_target(target_name)

        if not targets:
            logger.warning(f"No ChEMBL target found for: {target_name}")
            return result

        # Use first matching target
        target = targets[0]
        result["target_found"] = True
        result["target_info"] = {
            "chembl_id": target.chembl_id,
            "name": target.name,
            "type": target.target_type,
            "organism": target.organism,
            "uniprot": target.uniprot_id
        }

        # Get bioactivities
        activities = await self.get_bioactivities_for_target(
            target.chembl_id,
            max_results=50
        )

        if activities:
            # Get top binders (lowest IC50/Ki)
            top_binders = activities[:20]
            result["known_binders"] = [
                {
                    "name": a.compound_name,
                    "smiles": a.smiles,
                    "activity_type": a.activity_type,
                    "value_nM": a.value,
                    "chembl_id": a.chembl_id
                }
                for a in top_binders
            ]

            # Generate SAR insights
            if len(top_binders) >= 5:
                # Calculate potency ranges
                values = [a.value for a in top_binders]
                result["sar_insights"].append(
                    f"Best known binder: {top_binders[0].compound_name} ({top_binders[0].activity_type} = {top_binders[0].value:.1f} nM)"
                )
                result["sar_insights"].append(
                    f"Top 20 compounds range: {min(values):.1f} - {max(values):.1f} nM"
                )

        return result


# Singleton instance
chembl_client = ChEMBLClient()
