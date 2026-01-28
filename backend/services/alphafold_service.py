"""
AlphaFold 3 Integration Service
"""
from typing import Dict, Any, Optional, List
from pathlib import Path
import json
import subprocess
import tempfile
import asyncio
from loguru import logger

from config import settings


class AlphaFoldService:
    """Service for AlphaFold 3 structure predictions"""

    def __init__(self):
        self.model_dir = Path(settings.alphafold_model_dir)
        self.db_dir = Path(settings.alphafold_db_dir)
        self.docker_image = settings.alphafold_docker_image
        self.output_dir = settings.output_dir / "alphafold"
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def is_available(self) -> bool:
        """Check if AlphaFold 3 is available"""
        try:
            # Check for Docker
            result = subprocess.run(
                ["docker", "images", "-q", self.docker_image],
                capture_output=True,
                text=True
            )
            return bool(result.stdout.strip())
        except Exception:
            return False

    async def predict_structure(
        self,
        sequences: List[Dict[str, Any]],
        name: str = "prediction",
        model_seeds: List[int] = [1]
    ) -> Dict[str, Any]:
        """
        Run AlphaFold 3 structure prediction

        Args:
            sequences: List of sequence definitions
            name: Name for the prediction
            model_seeds: Random seeds for predictions

        Returns:
            Prediction results including structure files
        """
        if not self.is_available():
            return await self._fallback_prediction(sequences, name)

        try:
            # Create input JSON
            input_json = {
                "name": name,
                "sequences": sequences,
                "modelSeeds": model_seeds,
                "dialect": "alphafold3",
                "version": 1
            }

            # Write input file
            input_path = self.output_dir / f"{name}_input.json"
            with open(input_path, 'w') as f:
                json.dump(input_json, f, indent=2)

            # Run AlphaFold 3 via Docker
            output_path = self.output_dir / name

            cmd = [
                "docker", "run", "-it", "--rm",
                "--volume", f"{input_path.parent}:/root/af_input",
                "--volume", f"{output_path}:/root/af_output",
                "--volume", f"{self.model_dir}:/root/models",
                "--volume", f"{self.db_dir}:/root/public_databases",
                "--gpus", "all",
                self.docker_image,
                "python", "run_alphafold.py",
                f"--json_path=/root/af_input/{input_path.name}",
                "--model_dir=/root/models",
                "--output_dir=/root/af_output"
            ]

            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE
            )

            stdout, stderr = await asyncio.wait_for(
                process.communicate(),
                timeout=3600  # 1 hour timeout
            )

            if process.returncode == 0:
                # Parse results
                results = self._parse_results(output_path)
                return {
                    "success": True,
                    "prediction_name": name,
                    "output_directory": str(output_path),
                    **results
                }
            else:
                logger.error(f"AlphaFold error: {stderr.decode()}")
                return {
                    "success": False,
                    "error": stderr.decode()
                }

        except asyncio.TimeoutError:
            logger.error("AlphaFold prediction timed out")
            return {"success": False, "error": "Prediction timed out"}
        except Exception as e:
            logger.error(f"AlphaFold error: {e}")
            return {"success": False, "error": str(e)}

    async def predict_peptide_receptor_complex(
        self,
        peptide_sequence: str,
        receptor_sequence: str,
        peptide_name: str = "peptide",
        receptor_name: str = "receptor"
    ) -> Dict[str, Any]:
        """
        Predict structure of peptide-receptor complex

        Useful for predicting binding mode of enkephalin analogues
        to delta opioid receptor
        """
        sequences = [
            {
                "protein": {
                    "id": [receptor_name],
                    "sequence": receptor_sequence
                }
            },
            {
                "protein": {
                    "id": [peptide_name],
                    "sequence": peptide_sequence
                }
            }
        ]

        return await self.predict_structure(
            sequences,
            name=f"{peptide_name}_{receptor_name}_complex"
        )

    async def _fallback_prediction(
        self,
        sequences: List[Dict[str, Any]],
        name: str
    ) -> Dict[str, Any]:
        """
        Fallback when AlphaFold 3 is not available
        Uses AlphaFold DB API for known sequences
        """
        logger.warning("AlphaFold 3 not available, using fallback")

        results = {
            "success": True,
            "method": "fallback",
            "predictions": [],
            "note": "AlphaFold 3 not available. Using AlphaFold DB lookup for known sequences."
        }

        for seq_def in sequences:
            if "protein" in seq_def:
                sequence = seq_def["protein"].get("sequence", "")
                seq_id = seq_def["protein"].get("id", ["unknown"])[0]

                # Try to find in AlphaFold DB
                db_result = await self._search_alphafold_db(sequence)
                if db_result:
                    results["predictions"].append({
                        "id": seq_id,
                        "source": "alphafold_db",
                        **db_result
                    })
                else:
                    results["predictions"].append({
                        "id": seq_id,
                        "source": "not_found",
                        "note": "Sequence not in AlphaFold DB"
                    })

        return results

    async def _search_alphafold_db(self, sequence: str) -> Optional[Dict[str, Any]]:
        """Search AlphaFold DB for a sequence"""
        import aiohttp

        # For short peptides like enkephalins, they won't be in the DB
        # This is mainly for receptor structures
        if len(sequence) < 50:
            return None

        try:
            # Search UniProt for the sequence first
            async with aiohttp.ClientSession() as session:
                # This is a simplified search - full implementation would use
                # BLAST or exact sequence matching
                uniprot_url = "https://rest.uniprot.org/uniprotkb/search"
                params = {
                    "query": f"sequence:{sequence[:50]}",  # First 50 residues
                    "format": "json",
                    "size": 1
                }

                async with session.get(uniprot_url, params=params) as response:
                    if response.status == 200:
                        data = await response.json()
                        results = data.get("results", [])
                        if results:
                            uniprot_id = results[0].get("primaryAccession")

                            # Get AlphaFold structure
                            af_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
                            async with session.get(af_url) as af_response:
                                if af_response.status == 200:
                                    af_data = await af_response.json()
                                    return {
                                        "uniprot_id": uniprot_id,
                                        "pdb_url": af_data[0].get("pdbUrl"),
                                        "confidence": af_data[0].get("confidenceAvgLocalScore"),
                                        "model_created": af_data[0].get("latestVersion")
                                    }

        except Exception as e:
            logger.warning(f"AlphaFold DB search error: {e}")

        return None

    def _parse_results(self, output_dir: Path) -> Dict[str, Any]:
        """Parse AlphaFold 3 output"""
        results = {
            "structures": [],
            "confidence_scores": [],
            "pae_scores": []
        }

        try:
            # Look for output files
            for pdb_file in output_dir.glob("*.pdb"):
                results["structures"].append(str(pdb_file))

            for json_file in output_dir.glob("*confidence*.json"):
                with open(json_file) as f:
                    confidence_data = json.load(f)
                    results["confidence_scores"].append(confidence_data)

            for json_file in output_dir.glob("*pae*.json"):
                with open(json_file) as f:
                    pae_data = json.load(f)
                    results["pae_scores"].append(pae_data)

        except Exception as e:
            logger.warning(f"Error parsing AlphaFold results: {e}")

        return results


class DORSequenceData:
    """Delta Opioid Receptor sequence and structure data"""

    # Human Delta Opioid Receptor sequence (UniProt P41143)
    HUMAN_DOR_SEQUENCE = """
MEPAPSAGAELQPPLFANASDAYPSACPSAGANASGPPGARSASSLALAIAITALYSAVC
AVGLLGNVLVMFGIVRYTKMKTATNIYIFNLALADALATSTLPFQSAKYLMETWPFGELL
CKAVLSIDYYNMFTSIFTLTMMSVDRYIAVCHPVKALDFRTPAKAKLINICIWVLASGVG
VPIMVMAVTRPRDGAVVCMLQFPSPSWYWDTVTKICVFLFAFVVPILIITVCYGLMLLRL
RSVRLLSGSKEKDRSLRRITRMVLVVVGAFVVCWAPIHIFVIVWTLVDIDRRDPLVVAAL
HLCIALGYANSSLNPVLYAFLDENFKRCFRQLCRKPCGRPDPSSFSRAREATARERVTAC
TPSDGPGGGAAA
""".replace('\n', '').replace(' ', '')

    # Key binding residues for DOR (from crystal structures)
    BINDING_RESIDUES = [
        ("D128", "Salt bridge with amine"),
        ("Y129", "Aromatic interaction"),
        ("M132", "Hydrophobic pocket"),
        ("W274", "Aromatic stacking"),
        ("H278", "H-bond acceptor"),
        ("Y308", "Phenol interaction"),
    ]

    @classmethod
    def get_receptor_sequence(cls) -> str:
        """Get the DOR sequence for structure prediction"""
        return cls.HUMAN_DOR_SEQUENCE

    @classmethod
    def get_binding_site_info(cls) -> Dict[str, Any]:
        """Get information about the binding site"""
        return {
            "residues": cls.BINDING_RESIDUES,
            "pocket_volume": "~600 Å³",
            "key_interactions": [
                "Salt bridge between ligand amine and D128",
                "π-π stacking with aromatic residues",
                "Hydrophobic contacts in subpockets"
            ]
        }
