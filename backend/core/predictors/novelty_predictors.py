"""
Novelty and Patentability Predictors

Assess how novel a compound is compared to known chemical space.
"""

import asyncio
import logging
from typing import Dict, Any, List, Optional
from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import httpx

from ..prediction import Prediction, ConfidenceLevel, PredictionBasis
from ..predictor_registry import Predictor, PredictorCapability, registry

logger = logging.getLogger(__name__)


@registry.register
class FingerprintNoveltyPredictor(Predictor):
    """
    Assess novelty by comparing fingerprints to reference compounds.

    This is a LOCAL predictor that compares against provided reference
    compounds without external database queries.
    """

    capabilities = [
        PredictorCapability(
            predicts="reference_similarity",
            requires=["smiles", "reference_smiles"],
            confidence_tier=ConfidenceLevel.HIGH,
            description="Maximum Tanimoto similarity to reference compounds",
            tags=["novelty", "similarity"],
        ),
        PredictorCapability(
            predicts="novelty_score",
            requires=["smiles", "reference_smiles"],
            confidence_tier=ConfidenceLevel.HIGH,
            description="Novelty score (1 - max_similarity)",
            tags=["novelty", "score"],
        ),
    ]

    def _calculate_fingerprint(self, mol: Chem.Mol):
        """Calculate Morgan fingerprint."""
        return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)

    async def predict(
        self,
        inputs: Dict[str, Any],
        capability: PredictorCapability
    ) -> Prediction:
        smiles = inputs["smiles"]
        reference_smiles = inputs["reference_smiles"]

        if not isinstance(reference_smiles, list):
            reference_smiles = [reference_smiles]

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")

        query_fp = self._calculate_fingerprint(mol)

        # Calculate similarities to all references
        similarities = []
        closest_ref = None
        max_sim = 0.0

        for ref_smi in reference_smiles:
            ref_mol = Chem.MolFromSmiles(ref_smi)
            if ref_mol is None:
                continue

            ref_fp = self._calculate_fingerprint(ref_mol)
            sim = DataStructs.TanimotoSimilarity(query_fp, ref_fp)
            similarities.append(sim)

            if sim > max_sim:
                max_sim = sim
                closest_ref = ref_smi

        if not similarities:
            return Prediction.not_assessed(
                "No valid reference compounds provided",
                predictor_id=self.predictor_id
            )

        prop = capability.predicts

        if prop == "reference_similarity":
            # Determine confidence based on number of references
            if len(similarities) >= 5:
                conf = ConfidenceLevel.HIGH
                conf_score = 0.9
            elif len(similarities) >= 2:
                conf = ConfidenceLevel.MODERATE
                conf_score = 0.7
            else:
                conf = ConfidenceLevel.LOW
                conf_score = 0.5

            caveats = []
            if len(similarities) < 5:
                caveats.append(f"Only {len(similarities)} reference compounds")

            return Prediction(
                value=round(max_sim, 3),
                confidence=conf,
                confidence_score=conf_score,
                basis=PredictionBasis.SIMILARITY,
                caveats=caveats,
                metadata={
                    "closest_reference": closest_ref,
                    "num_references": len(similarities),
                    "mean_similarity": round(sum(similarities) / len(similarities), 3),
                },
            )

        elif prop == "novelty_score":
            novelty = 1.0 - max_sim

            # Determine if patentable (rough heuristic)
            patentable = max_sim < 0.85

            caveats = []
            if max_sim >= 0.85:
                caveats.append(f"High similarity ({max_sim:.0%}) to known compound - may not be patentable")
            elif max_sim >= 0.7:
                caveats.append(f"Moderate similarity ({max_sim:.0%}) to known compound")

            return Prediction(
                value=round(novelty, 3),
                confidence=ConfidenceLevel.HIGH,
                confidence_score=0.9,
                basis=PredictionBasis.SIMILARITY,
                caveats=caveats,
                metadata={
                    "closest_similarity": max_sim,
                    "closest_reference": closest_ref,
                    "patentable_estimate": patentable,
                },
            )

        raise ValueError(f"Unknown property: {prop}")


@registry.register
class DatabaseNoveltyPredictor(Predictor):
    """
    Check novelty against external databases (PubChem, ChEMBL).

    This predictor makes external API calls, so predictions may be slow
    and subject to API availability.
    """

    capabilities = [
        PredictorCapability(
            predicts="pubchem_novelty",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.MODERATE,
            description="Novelty check against PubChem database",
            tags=["novelty", "database", "pubchem"],
        ),
        PredictorCapability(
            predicts="chembl_novelty",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.MODERATE,
            description="Novelty check against ChEMBL database",
            tags=["novelty", "database", "chembl"],
        ),
        PredictorCapability(
            predicts="database_novelty",
            requires=["smiles"],
            confidence_tier=ConfidenceLevel.MODERATE,
            description="Combined database novelty check",
            tags=["novelty", "database"],
        ),
    ]

    def __init__(self):
        super().__init__()
        self._client: Optional[httpx.AsyncClient] = None

    async def initialize(self):
        """Initialize HTTP client."""
        self._client = httpx.AsyncClient(timeout=30.0)

    async def _search_pubchem(self, smiles: str) -> Dict[str, Any]:
        """Search PubChem for similar compounds."""
        try:
            # First convert SMILES to InChIKey for exact match
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return {"found": False, "error": "Invalid SMILES"}

            inchi = Chem.MolToInchi(mol)
            inchi_key = Chem.InchiToInchiKey(inchi) if inchi else None

            # Try exact match first
            if inchi_key:
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchi_key}/cids/JSON"
                response = await self._client.get(url)

                if response.status_code == 200:
                    data = response.json()
                    if "IdentifierList" in data and data["IdentifierList"].get("CID"):
                        return {
                            "found": True,
                            "exact_match": True,
                            "cid": data["IdentifierList"]["CID"][0],
                        }

            # No exact match - compound is novel to PubChem
            return {
                "found": False,
                "exact_match": False,
            }

        except Exception as e:
            logger.warning(f"PubChem search failed: {e}")
            return {"found": False, "error": str(e)}

    async def _search_chembl(self, smiles: str) -> Dict[str, Any]:
        """Search ChEMBL for similar compounds."""
        try:
            # Use similarity search
            url = "https://www.ebi.ac.uk/chembl/api/data/similarity.json"
            params = {
                "smiles": smiles,
                "similarity": 85,  # 85% threshold
            }

            response = await self._client.get(url, params=params)

            if response.status_code == 200:
                data = response.json()
                molecules = data.get("molecules", [])

                if molecules:
                    # Find most similar - handle similarity as string or number
                    def get_sim(m):
                        sim = m.get("similarity", 0)
                        try:
                            return float(sim)
                        except (ValueError, TypeError):
                            return 0

                    best = max(molecules, key=get_sim)
                    sim_value = get_sim(best)
                    return {
                        "found": True,
                        "most_similar_id": best.get("molecule_chembl_id"),
                        "similarity": sim_value / 100 if sim_value > 1 else sim_value,
                        "match_count": len(molecules),
                    }

            return {"found": False}

        except Exception as e:
            logger.warning(f"ChEMBL search failed: {e}")
            return {"found": False, "error": str(e)}

    async def predict(
        self,
        inputs: Dict[str, Any],
        capability: PredictorCapability
    ) -> Prediction:
        smiles = inputs["smiles"]
        prop = capability.predicts

        if prop == "pubchem_novelty":
            result = await self._search_pubchem(smiles)

            if "error" in result:
                return Prediction(
                    value=None,
                    confidence=ConfidenceLevel.LOW,
                    confidence_score=0.3,
                    basis=PredictionBasis.LITERATURE,
                    caveats=[f"PubChem search failed: {result['error']}"],
                    assessed=False,
                )

            is_novel = not result.get("exact_match", False)

            caveats = []
            if result.get("exact_match"):
                caveats.append(f"Exact match in PubChem (CID: {result.get('cid')})")

            return Prediction(
                value=is_novel,
                confidence=ConfidenceLevel.MODERATE,
                confidence_score=0.8,
                basis=PredictionBasis.LITERATURE,
                caveats=caveats,
                metadata=result,
            )

        elif prop == "chembl_novelty":
            result = await self._search_chembl(smiles)

            if "error" in result:
                return Prediction(
                    value=None,
                    confidence=ConfidenceLevel.LOW,
                    confidence_score=0.3,
                    basis=PredictionBasis.LITERATURE,
                    caveats=[f"ChEMBL search failed: {result['error']}"],
                    assessed=False,
                )

            is_novel = not result.get("found", False)
            similarity = result.get("similarity", 0)

            caveats = []
            if result.get("found"):
                caveats.append(
                    f"Similar compound in ChEMBL: {result.get('most_similar_id')} "
                    f"({similarity:.0%} similar)"
                )

            return Prediction(
                value=is_novel,
                confidence=ConfidenceLevel.MODERATE,
                confidence_score=0.75,
                basis=PredictionBasis.LITERATURE,
                caveats=caveats,
                metadata=result,
            )

        elif prop == "database_novelty":
            # Run both searches in parallel
            pubchem_result, chembl_result = await asyncio.gather(
                self._search_pubchem(smiles),
                self._search_chembl(smiles),
            )

            # Combine results
            pubchem_novel = not pubchem_result.get("exact_match", False)
            chembl_novel = not chembl_result.get("found", False)

            # Novel only if novel in both databases
            is_novel = pubchem_novel and chembl_novel

            caveats = []
            if not pubchem_novel:
                caveats.append(f"Found in PubChem (CID: {pubchem_result.get('cid')})")
            if not chembl_novel:
                sim = chembl_result.get("similarity", 0)
                caveats.append(
                    f"Similar to ChEMBL compound {chembl_result.get('most_similar_id')} ({sim:.0%})"
                )

            # Confidence based on whether both searches succeeded
            pubchem_ok = "error" not in pubchem_result
            chembl_ok = "error" not in chembl_result

            if pubchem_ok and chembl_ok:
                conf = ConfidenceLevel.MODERATE
                conf_score = 0.85
            elif pubchem_ok or chembl_ok:
                conf = ConfidenceLevel.LOW
                conf_score = 0.5
                caveats.append("One database search failed")
            else:
                conf = ConfidenceLevel.LOW
                conf_score = 0.2
                caveats.append("Both database searches failed")

            return Prediction(
                value=is_novel,
                confidence=conf,
                confidence_score=conf_score,
                basis=PredictionBasis.LITERATURE,
                caveats=caveats,
                metadata={
                    "pubchem": pubchem_result,
                    "chembl": chembl_result,
                },
            )

        raise ValueError(f"Unknown property: {prop}")
