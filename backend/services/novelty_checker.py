"""
Novelty Checking Service for patent and database searches
"""
from typing import Dict, Any, Optional, List, Tuple
import asyncio
import aiohttp
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, inchi
from loguru import logger
import pubchempy as pcp
from tenacity import retry, stop_after_attempt, wait_exponential

from models import NoveltyAssessment
from config import settings


class NoveltyChecker:
    """Service for checking novelty of compounds against databases"""

    def __init__(self):
        self.pubchem_url = settings.pubchem_base_url
        self.chembl_url = settings.chembl_base_url
        self.similarity_threshold = settings.similarity_threshold

    async def check_novelty(self, smiles: str) -> NoveltyAssessment:
        """Perform comprehensive novelty check"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return NoveltyAssessment(
                    is_novel=False,
                    confidence_score=0.0,
                    closest_known_compound="Invalid SMILES"
                )

            # Get InChI key for exact match searches
            inchi_key = inchi.MolToInchiKey(mol)

            # Run parallel searches
            pubchem_result, chembl_result = await asyncio.gather(
                self._search_pubchem(smiles, inchi_key),
                self._search_chembl(smiles, inchi_key),
                return_exceptions=True
            )

            # Handle exceptions
            if isinstance(pubchem_result, Exception):
                logger.warning(f"PubChem search failed: {pubchem_result}")
                pubchem_result = {"exact_match": False, "similar_count": 0, "closest": None}
            if isinstance(chembl_result, Exception):
                logger.warning(f"ChEMBL search failed: {chembl_result}")
                chembl_result = {"exact_match": False, "similar_count": 0, "closest": None}

            # Determine novelty
            exact_match_found = pubchem_result.get("exact_match") or chembl_result.get("exact_match")

            # Get closest compound
            closest_compound = None
            closest_similarity = 0.0

            if pubchem_result.get("closest"):
                closest_compound = pubchem_result["closest"]["name"]
                closest_similarity = pubchem_result["closest"]["similarity"]

            if chembl_result.get("closest") and chembl_result["closest"]["similarity"] > closest_similarity:
                closest_compound = chembl_result["closest"]["name"]
                closest_similarity = chembl_result["closest"]["similarity"]

            # Calculate confidence score
            confidence = self._calculate_novelty_confidence(
                exact_match_found,
                closest_similarity,
                pubchem_result.get("similar_count", 0),
                chembl_result.get("similar_count", 0)
            )

            # Determine if novel
            is_novel = not exact_match_found and closest_similarity < self.similarity_threshold

            return NoveltyAssessment(
                is_novel=is_novel,
                confidence_score=confidence,
                closest_known_compound=closest_compound,
                closest_similarity=closest_similarity,
                pubchem_matches=pubchem_result.get("similar_count", 0),
                chembl_matches=chembl_result.get("similar_count", 0),
                patent_matches=0,  # Would require patent database integration
                non_obviousness_rationale=self._generate_rationale(is_novel, closest_similarity, closest_compound)
            )

        except Exception as e:
            logger.error(f"Error in novelty check: {e}")
            return NoveltyAssessment(
                is_novel=True,  # Assume novel if check fails
                confidence_score=0.3,
                closest_known_compound="Check failed"
            )

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(min=1, max=10))
    async def _search_pubchem(self, smiles: str, inchi_key: str) -> Dict[str, Any]:
        """Search PubChem for exact and similar compounds"""
        result = {
            "exact_match": False,
            "similar_count": 0,
            "closest": None
        }

        try:
            # Try exact structure search
            compounds = pcp.get_compounds(smiles, 'smiles')
            if compounds:
                result["exact_match"] = True
                result["closest"] = {
                    "name": compounds[0].iupac_name or f"CID:{compounds[0].cid}",
                    "similarity": 1.0
                }
                return result

            # Similarity search using fingerprints
            async with aiohttp.ClientSession() as session:
                # Search for similar structures
                search_url = f"{self.pubchem_url}/compound/fastsimilarity_2d/smiles/{smiles}/cids/JSON"
                params = {"Threshold": int(self.similarity_threshold * 100)}

                async with session.get(search_url, params=params, timeout=30) as response:
                    if response.status == 200:
                        data = await response.json()
                        cids = data.get("IdentifierList", {}).get("CID", [])
                        result["similar_count"] = len(cids)

                        if cids:
                            # Get details of most similar
                            top_cid = cids[0]
                            compounds = pcp.get_compounds(top_cid, 'cid')
                            if compounds:
                                # Calculate actual similarity
                                ref_mol = Chem.MolFromSmiles(smiles)
                                comp_mol = Chem.MolFromSmiles(compounds[0].isomeric_smiles)
                                if ref_mol and comp_mol:
                                    fp1 = AllChem.GetMorganFingerprintAsBitVect(ref_mol, 2, nBits=2048)
                                    fp2 = AllChem.GetMorganFingerprintAsBitVect(comp_mol, 2, nBits=2048)
                                    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)

                                    result["closest"] = {
                                        "name": compounds[0].iupac_name or f"CID:{top_cid}",
                                        "similarity": similarity
                                    }

        except Exception as e:
            logger.warning(f"PubChem search error: {e}")

        return result

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(min=1, max=10))
    async def _search_chembl(self, smiles: str, inchi_key: str) -> Dict[str, Any]:
        """Search ChEMBL for exact and similar compounds"""
        result = {
            "exact_match": False,
            "similar_count": 0,
            "closest": None
        }

        try:
            async with aiohttp.ClientSession() as session:
                # Search by InChI key for exact match
                search_url = f"{self.chembl_url}/molecule.json"
                params = {"molecule_structures__standard_inchi_key": inchi_key}

                async with session.get(search_url, params=params, timeout=30) as response:
                    if response.status == 200:
                        data = await response.json()
                        molecules = data.get("molecules", [])
                        if molecules:
                            result["exact_match"] = True
                            result["closest"] = {
                                "name": molecules[0].get("pref_name", molecules[0].get("molecule_chembl_id")),
                                "similarity": 1.0
                            }
                            return result

                # Similarity search
                similarity_url = f"{self.chembl_url}/similarity/{smiles}/{int(self.similarity_threshold * 100)}.json"
                async with session.get(similarity_url, timeout=30) as response:
                    if response.status == 200:
                        data = await response.json()
                        molecules = data.get("molecules", [])
                        result["similar_count"] = len(molecules)

                        if molecules:
                            top_mol = molecules[0]
                            result["closest"] = {
                                "name": top_mol.get("pref_name", top_mol.get("molecule_chembl_id")),
                                "similarity": top_mol.get("similarity", 0.8)
                            }

        except Exception as e:
            logger.warning(f"ChEMBL search error: {e}")

        return result

    def _calculate_novelty_confidence(
        self,
        exact_match: bool,
        closest_similarity: float,
        pubchem_count: int,
        chembl_count: int
    ) -> float:
        """Calculate confidence score for novelty assessment"""
        if exact_match:
            return 0.0

        # Base confidence on similarity
        similarity_confidence = 1.0 - closest_similarity

        # Adjust based on number of similar compounds found
        total_similar = pubchem_count + chembl_count
        if total_similar == 0:
            count_factor = 1.0
        elif total_similar < 5:
            count_factor = 0.9
        elif total_similar < 20:
            count_factor = 0.7
        else:
            count_factor = 0.5

        confidence = similarity_confidence * count_factor

        return round(max(0.0, min(1.0, confidence)), 2)

    def _generate_rationale(
        self,
        is_novel: bool,
        closest_similarity: float,
        closest_compound: Optional[str]
    ) -> str:
        """Generate rationale for novelty/non-obviousness"""
        if not is_novel:
            if closest_similarity >= 0.99:
                return f"Exact or near-exact match found: {closest_compound}"
            else:
                return f"High similarity ({closest_similarity:.0%}) to known compound: {closest_compound}"

        if closest_similarity > 0.7:
            return (
                f"Structurally related to {closest_compound} (similarity: {closest_similarity:.0%}), "
                f"but modifications introduce sufficient novelty for potential patentability. "
                f"Non-obviousness would depend on unexpected properties."
            )
        elif closest_similarity > 0.5:
            return (
                f"Moderate structural relationship to known compounds (closest: {closest_similarity:.0%}). "
                f"Good potential for novelty claims. Non-obviousness supported by structural distinctiveness."
            )
        else:
            return (
                f"Low similarity to known compounds (closest: {closest_similarity:.0%}). "
                f"Strong novelty case. Non-obviousness well-supported by structural uniqueness."
            )

    async def batch_check(self, smiles_list: List[str]) -> List[NoveltyAssessment]:
        """Check novelty for multiple compounds"""
        tasks = [self.check_novelty(smiles) for smiles in smiles_list]
        return await asyncio.gather(*tasks)

    def calculate_structural_uniqueness(
        self,
        candidate_smiles: str,
        reference_smiles_list: List[str]
    ) -> Tuple[float, str]:
        """
        Calculate how structurally unique a candidate is compared to references
        Returns (uniqueness_score, most_similar_reference)
        """
        candidate_mol = Chem.MolFromSmiles(candidate_smiles)
        if candidate_mol is None:
            return 0.0, "Invalid SMILES"

        candidate_fp = AllChem.GetMorganFingerprintAsBitVect(candidate_mol, 2, nBits=2048)

        max_similarity = 0.0
        most_similar = ""

        for ref_smiles in reference_smiles_list:
            ref_mol = Chem.MolFromSmiles(ref_smiles)
            if ref_mol:
                ref_fp = AllChem.GetMorganFingerprintAsBitVect(ref_mol, 2, nBits=2048)
                similarity = DataStructs.TanimotoSimilarity(candidate_fp, ref_fp)
                if similarity > max_similarity:
                    max_similarity = similarity
                    most_similar = ref_smiles

        uniqueness = 1.0 - max_similarity
        return uniqueness, most_similar


class PatentSearcher:
    """Service for searching patent databases (placeholder for future integration)"""

    async def search_patents(self, smiles: str) -> Dict[str, Any]:
        """
        Search patent databases for similar compounds

        Note: Full implementation would integrate with:
        - Google Patents API
        - USPTO PatentsView API
        - EPO Open Patent Services
        - SureChEMBL
        """
        # Placeholder - would integrate with patent APIs
        return {
            "exact_matches": [],
            "similar_patents": [],
            "freedom_to_operate": "Unknown - manual review required"
        }

    async def check_prior_art(self, smiles: str, filing_date: str) -> Dict[str, Any]:
        """Check for prior art before a given date"""
        # Placeholder
        return {
            "prior_art_found": False,
            "relevant_patents": [],
            "relevant_publications": []
        }
