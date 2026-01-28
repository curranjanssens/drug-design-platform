"""
PDB Structure Fetcher

Fetches protein structures from RCSB PDB and AlphaFold databases.
Provides PDB content for molecular docking.
"""

import asyncio
import json
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Dict, Any

import httpx

logger = logging.getLogger(__name__)


@dataclass
class PDBInfo:
    """Information about a PDB structure."""
    pdb_id: str
    title: str
    resolution: Optional[float] = None
    method: Optional[str] = None  # X-RAY, NMR, CRYO-EM
    uniprot_id: Optional[str] = None
    organism: Optional[str] = None
    source: str = "RCSB"  # RCSB or AlphaFold


@dataclass
class FetchResult:
    """Result of fetching a PDB structure."""
    success: bool
    pdb_content: Optional[str] = None
    pdb_info: Optional[PDBInfo] = None
    error_message: Optional[str] = None


class PDBFetcher:
    """
    Fetches protein structures from public databases.

    Priority order:
    1. RCSB PDB (experimental structures)
    2. AlphaFold DB (predicted structures)

    Usage:
        fetcher = PDBFetcher()
        result = await fetcher.fetch_by_target_name("EGFR")
        if result.success:
            pdb_content = result.pdb_content
    """

    RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
    RCSB_DATA_URL = "https://data.rcsb.org/rest/v1/core/entry"
    RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download"
    ALPHAFOLD_URL = "https://alphafold.ebi.ac.uk/files"
    UNIPROT_URL = "https://rest.uniprot.org/uniprotkb"

    def __init__(self, cache_dir: Optional[Path] = None):
        self.cache_dir = cache_dir
        if cache_dir:
            cache_dir.mkdir(parents=True, exist_ok=True)
        self._client: Optional[httpx.AsyncClient] = None

    async def _get_client(self) -> httpx.AsyncClient:
        """Get or create HTTP client."""
        if self._client is None:
            self._client = httpx.AsyncClient(timeout=60.0)
        return self._client

    async def close(self):
        """Close the HTTP client."""
        if self._client:
            await self._client.aclose()
            self._client = None

    async def fetch_by_target_name(
        self,
        target_name: str,
        organism: str = "Homo sapiens",
        prefer_resolution: bool = True,
    ) -> FetchResult:
        """
        Fetch PDB structure by target protein name.

        Args:
            target_name: Protein name (e.g., "EGFR", "FAAH", "CDK2")
            organism: Organism filter (default: human)
            prefer_resolution: If True, prefer higher resolution structures

        Returns:
            FetchResult with PDB content if successful
        """
        logger.info(f"Fetching PDB for target: {target_name}")

        # Check cache first
        if self.cache_dir:
            cache_file = self.cache_dir / f"{target_name.lower().replace(' ', '_')}.pdb"
            if cache_file.exists():
                logger.info(f"Using cached PDB: {cache_file}")
                return FetchResult(
                    success=True,
                    pdb_content=cache_file.read_text(),
                    pdb_info=PDBInfo(pdb_id="cached", title=target_name, source="cache"),
                )

        # Try RCSB PDB first
        pdb_ids = await self._search_rcsb(target_name, organism)

        if pdb_ids:
            # Get the best structure
            best_pdb = await self._select_best_structure(pdb_ids, prefer_resolution)
            if best_pdb:
                result = await self._download_pdb(best_pdb)
                if result.success and self.cache_dir:
                    # Cache the result
                    cache_file = self.cache_dir / f"{target_name.lower().replace(' ', '_')}.pdb"
                    cache_file.write_text(result.pdb_content)
                return result

        # Fall back to AlphaFold
        logger.info(f"No RCSB structure found, trying AlphaFold for {target_name}")
        uniprot_id = await self._get_uniprot_id(target_name, organism)

        if uniprot_id:
            result = await self._fetch_alphafold(uniprot_id, target_name)
            if result.success and self.cache_dir:
                cache_file = self.cache_dir / f"{target_name.lower().replace(' ', '_')}.pdb"
                cache_file.write_text(result.pdb_content)
            return result

        return FetchResult(
            success=False,
            error_message=f"Could not find structure for {target_name}",
        )

    async def fetch_by_pdb_id(self, pdb_id: str) -> FetchResult:
        """Fetch PDB structure by PDB ID."""
        return await self._download_pdb(pdb_id)

    async def fetch_by_uniprot(self, uniprot_id: str, target_name: str = "") -> FetchResult:
        """Fetch structure by UniProt ID (uses AlphaFold)."""
        return await self._fetch_alphafold(uniprot_id, target_name or uniprot_id)

    async def _search_rcsb(
        self,
        target_name: str,
        organism: str = "Homo sapiens",
        limit: int = 10,
    ) -> List[str]:
        """Search RCSB PDB for structures matching target name."""
        client = await self._get_client()

        # Build search query using RCSB v2 API
        # Search by gene name with organism filter
        query = {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
                            "operator": "exact_match",
                            "value": target_name.upper(),
                        },
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_entity_source_organism.ncbi_scientific_name",
                            "operator": "exact_match",
                            "value": organism,
                        },
                    },
                ],
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {"start": 0, "rows": limit},
            },
        }

        try:
            response = await client.post(
                self.RCSB_SEARCH_URL,
                json=query,
                headers={"Content-Type": "application/json"},
            )

            if response.status_code == 200:
                data = response.json()
                pdb_ids = [hit["identifier"] for hit in data.get("result_set", [])]
                if pdb_ids:
                    logger.info(f"Found {len(pdb_ids)} RCSB structures for {target_name} (gene name)")
                    return pdb_ids

            # If gene name search fails, try title search
            title_query = {
                "query": {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "struct.title",
                        "operator": "contains_phrase",
                        "value": target_name,
                    },
                },
                "return_type": "entry",
                "request_options": {
                    "paginate": {"start": 0, "rows": limit},
                },
            }

            response = await client.post(
                self.RCSB_SEARCH_URL,
                json=title_query,
                headers={"Content-Type": "application/json"},
            )

            if response.status_code == 200:
                data = response.json()
                pdb_ids = [hit["identifier"] for hit in data.get("result_set", [])]
                if pdb_ids:
                    logger.info(f"Found {len(pdb_ids)} RCSB structures for {target_name} (title)")
                    return pdb_ids

            logger.warning(f"RCSB search found no results for: {target_name}")

        except Exception as e:
            logger.warning(f"RCSB search error: {e}")

        return []

    async def _select_best_structure(
        self,
        pdb_ids: List[str],
        prefer_resolution: bool = True,
    ) -> Optional[str]:
        """Select the best structure from a list of PDB IDs."""
        if not pdb_ids:
            return None

        if not prefer_resolution:
            return pdb_ids[0]

        client = await self._get_client()
        best_pdb = pdb_ids[0]
        best_resolution = float("inf")

        for pdb_id in pdb_ids[:5]:  # Check top 5
            try:
                response = await client.get(f"{self.RCSB_DATA_URL}/{pdb_id}")
                if response.status_code == 200:
                    data = response.json()
                    resolution = data.get("rcsb_entry_info", {}).get("resolution_combined")
                    if resolution and resolution < best_resolution:
                        best_resolution = resolution
                        best_pdb = pdb_id
            except Exception:
                pass

        logger.info(f"Selected best PDB: {best_pdb} (resolution: {best_resolution})")
        return best_pdb

    async def _download_pdb(self, pdb_id: str) -> FetchResult:
        """Download PDB file by ID."""
        client = await self._get_client()

        try:
            # Try PDB format first
            url = f"{self.RCSB_DOWNLOAD_URL}/{pdb_id}.pdb"
            response = await client.get(url)

            if response.status_code == 200:
                pdb_content = response.text

                # Get metadata
                info = await self._get_pdb_info(pdb_id)

                return FetchResult(
                    success=True,
                    pdb_content=pdb_content,
                    pdb_info=info or PDBInfo(pdb_id=pdb_id, title=pdb_id, source="RCSB"),
                )
            else:
                logger.warning(f"Failed to download PDB {pdb_id}: {response.status_code}")

        except Exception as e:
            logger.error(f"Error downloading PDB {pdb_id}: {e}")

        return FetchResult(
            success=False,
            error_message=f"Failed to download PDB {pdb_id}",
        )

    async def _get_pdb_info(self, pdb_id: str) -> Optional[PDBInfo]:
        """Get metadata about a PDB structure."""
        client = await self._get_client()

        try:
            response = await client.get(f"{self.RCSB_DATA_URL}/{pdb_id}")
            if response.status_code == 200:
                data = response.json()
                entry_info = data.get("rcsb_entry_info", {})
                struct = data.get("struct", {})

                return PDBInfo(
                    pdb_id=pdb_id,
                    title=struct.get("title", pdb_id),
                    resolution=entry_info.get("resolution_combined"),
                    method=entry_info.get("experimental_method"),
                    source="RCSB",
                )
        except Exception:
            pass

        return None

    async def _get_uniprot_id(
        self,
        target_name: str,
        organism: str = "Homo sapiens",
    ) -> Optional[str]:
        """Get UniProt ID from target name."""
        client = await self._get_client()

        try:
            # Search UniProt
            query = f"(gene:{target_name} OR protein_name:{target_name}) AND organism_name:{organism}"
            url = f"{self.UNIPROT_URL}/search"
            params = {
                "query": query,
                "format": "json",
                "size": 1,
            }

            response = await client.get(url, params=params)

            if response.status_code == 200:
                data = response.json()
                results = data.get("results", [])
                if results:
                    uniprot_id = results[0].get("primaryAccession")
                    logger.info(f"Found UniProt ID: {uniprot_id} for {target_name}")
                    return uniprot_id

        except Exception as e:
            logger.warning(f"UniProt search error: {e}")

        return None

    async def _fetch_alphafold(
        self,
        uniprot_id: str,
        target_name: str,
    ) -> FetchResult:
        """Fetch AlphaFold predicted structure."""
        client = await self._get_client()

        try:
            # Try different AlphaFold versions
            for version in ["v4", "v3", "v2"]:
                url = f"{self.ALPHAFOLD_URL}/AF-{uniprot_id}-F1-model_{version}.pdb"
                response = await client.get(url)

                if response.status_code == 200:
                    logger.info(f"Found AlphaFold structure: {uniprot_id} ({version})")
                    return FetchResult(
                        success=True,
                        pdb_content=response.text,
                        pdb_info=PDBInfo(
                            pdb_id=f"AF-{uniprot_id}",
                            title=target_name,
                            uniprot_id=uniprot_id,
                            source="AlphaFold",
                        ),
                    )

            logger.warning(f"No AlphaFold structure found for {uniprot_id}")

        except Exception as e:
            logger.error(f"AlphaFold fetch error: {e}")

        return FetchResult(
            success=False,
            error_message=f"No AlphaFold structure for {uniprot_id}",
        )


# Global instance
pdb_fetcher = PDBFetcher()


def get_pdb_fetcher(cache_dir: Optional[Path] = None) -> PDBFetcher:
    """Get PDB fetcher with optional cache directory."""
    global pdb_fetcher
    if cache_dir and pdb_fetcher.cache_dir != cache_dir:
        pdb_fetcher = PDBFetcher(cache_dir=cache_dir)
    return pdb_fetcher
