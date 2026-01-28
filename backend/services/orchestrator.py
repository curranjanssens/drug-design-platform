"""
Claude-powered orchestrator for the drug design pipeline.
Coordinates all services and provides natural language interface.
"""
import asyncio
import os
from typing import Dict, List, Optional, Any
from datetime import datetime
import logging
import json
import re

from models import (
    DesignJob, DesignedMolecule, DesignType, JobStatus,
    MolecularProperties, ADMETProfile, BindingPrediction,
    NoveltyAssessment, SynthesisRoute
)
from .molecular_analyzer import MolecularAnalyzer
from .admet_predictor import ADMETPredictor
from .retrosynthesis import RetrosynthesisService
from .patent_analyzer import PatentAnalyzer
from .molecule_generator import MoleculeGenerator

logger = logging.getLogger(__name__)


class DrugDesignOrchestrator:
    """
    Main orchestrator that coordinates all drug design services.
    Uses Claude to interpret requests and guide the design process.
    """

    def __init__(self, anthropic_api_key: str = None):
        self.api_key = anthropic_api_key or os.getenv("ANTHROPIC_API_KEY")

        # Initialize services
        self.analyzer = MolecularAnalyzer()
        self.admet = ADMETPredictor()
        self.retrosynthesis = RetrosynthesisService()
        self.patent = PatentAnalyzer()
        self.generator = MoleculeGenerator()

        # Job storage (in production, use database)
        self.jobs: Dict[str, DesignJob] = {}

    async def process_request(self, prompt: str, user_id: str = None) -> DesignJob:
        """
        Process a natural language drug design request.

        Args:
            prompt: User's design request
            user_id: Optional user identifier

        Returns:
            DesignJob with results
        """
        # Create job
        job = DesignJob(
            prompt=prompt,
            user_id=user_id,
            status=JobStatus.RUNNING,
            started_at=datetime.utcnow()
        )
        self.jobs[job.id] = job

        try:
            # Parse request using Claude
            job.update_progress("Analyzing request", 0.05)
            parsed = await self._parse_request(prompt)

            job.design_type = parsed["design_type"]
            job.reference_smiles = parsed.get("reference_smiles")
            job.reference_name = parsed.get("reference_name")
            job.target_name = parsed.get("target_name")
            job.constraints = parsed.get("constraints", {})

            # Execute design pipeline
            molecules = await self._execute_pipeline(job, parsed)

            # Complete job
            job.complete(molecules)

        except Exception as e:
            logger.error(f"Job {job.id} failed: {e}")
            job.fail(str(e))

        return job

    async def _parse_request(self, prompt: str) -> Dict[str, Any]:
        """
        Use Claude to parse and understand the natural language request.
        Claude extracts design parameters and looks up compound structures.
        """
        import httpx
        import pubchempy as pcp

        result = {
            "design_type": DesignType.ANALOGUE,
            "reference_smiles": None,
            "reference_name": None,
            "target_name": None,
            "constraints": {}
        }

        # Use Claude to understand the request
        logger.info(f"API key available: {bool(self.api_key)}")
        if self.api_key:
            try:
                logger.info("Calling Claude API to parse request...")
                async with httpx.AsyncClient(timeout=30.0) as client:
                    response = await client.post(
                        "https://api.anthropic.com/v1/messages",
                        headers={
                            "x-api-key": self.api_key,
                            "anthropic-version": "2023-06-01",
                            "content-type": "application/json"
                        },
                        json={
                            "model": "claude-sonnet-4-20250514",
                            "max_tokens": 1024,
                            "messages": [{
                                "role": "user",
                                "content": f"""Analyze this drug design request and extract key information.

Request: "{prompt}"

Return a JSON object with:
- design_type: one of "analogue", "bivalent", "scaffold_hop", "de_novo", "optimization"
- compound_names: list of drug/compound names mentioned (e.g., ["PF-04457845", "ibuprofen"])
- target_names: list of biological targets mentioned (e.g., ["FAAH", "COX-2", "5HT1A"])
- constraints: object with any mentioned requirements like {{"bbb_penetration": true, "oral": true, "stability": true}}

Only return valid JSON, no explanation."""
                            }]
                        }
                    )

                    logger.info(f"Claude API response status: {response.status_code}")
                    if response.status_code == 200:
                        data = response.json()
                        content = data["content"][0]["text"]
                        logger.info(f"Claude response: {content[:200]}")

                        # Strip markdown code blocks if present
                        content = content.strip()
                        if content.startswith("```"):
                            # Remove ```json or ``` prefix and trailing ```
                            lines = content.split("\n")
                            content = "\n".join(lines[1:-1] if lines[-1].strip() == "```" else lines[1:])

                        # Extract JSON from response
                        try:
                            parsed = json.loads(content)

                            # Map design type
                            type_map = {
                                "analogue": DesignType.ANALOGUE,
                                "bivalent": DesignType.BIVALENT,
                                "scaffold_hop": DesignType.SCAFFOLD_HOP,
                                "de_novo": DesignType.DE_NOVO,
                                "optimization": DesignType.OPTIMIZATION
                            }
                            result["design_type"] = type_map.get(parsed.get("design_type", "analogue"), DesignType.ANALOGUE)
                            result["constraints"] = parsed.get("constraints", {})

                            if parsed.get("target_names"):
                                result["target_name"] = parsed["target_names"][0]

                            # Look up compound structures from PubChem
                            compound_names = parsed.get("compound_names", [])
                            for name in compound_names:
                                logger.info(f"Looking up compound: {name}")
                                smiles = await self._lookup_compound(name)
                                if smiles:
                                    result["reference_name"] = name
                                    result["reference_smiles"] = smiles
                                    logger.info(f"Found {name}: {smiles[:50]}...")
                                    break

                        except json.JSONDecodeError:
                            logger.warning("Failed to parse Claude response as JSON")
                    else:
                        logger.warning(f"Claude API returned {response.status_code}")

            except Exception as e:
                logger.warning(f"Claude API call failed: {e}, falling back to pattern matching")

        # Fallback: if no SMILES found yet, try direct PubChem lookup on any compound-like names in prompt
        if not result["reference_smiles"]:
            # Extract potential compound names (alphanumeric with hyphens like PF-04457845)
            potential_names = re.findall(r'\b[A-Z]{1,3}[-]?\d{2,}[-]?\d*\b', prompt, re.IGNORECASE)
            potential_names += re.findall(r'\b[a-zA-Z]+[-]?[a-zA-Z]*\b', prompt)

            # Filter common words
            stop_words = {'design', 'novel', 'new', 'version', 'analogue', 'analog', 'improved',
                         'patentable', 'for', 'with', 'the', 'treating', 'treatment', 'of', 'and',
                         'dual', 'agonist', 'antagonist', 'inhibitor', 'bivalent', 'ligand',
                         'cns', 'penetration', 'depression', 'anxiety', 'pain', 'chronic'}

            for name in potential_names:
                if name.lower() in stop_words or len(name) < 3:
                    continue
                smiles = await self._lookup_compound(name)
                if smiles:
                    result["reference_name"] = name
                    result["reference_smiles"] = smiles
                    logger.info(f"Fallback found {name}: {smiles[:50]}...")
                    break

        return result

    async def _lookup_compound(self, name: str) -> Optional[str]:
        """Look up a compound structure by name from PubChem."""
        import pubchempy as pcp
        from rdkit import Chem

        try:
            compounds = pcp.get_compounds(name, 'name')
            if not compounds:
                return None

            c = compounds[0]

            # Try SMILES first
            if c.isomeric_smiles:
                return c.isomeric_smiles
            if c.canonical_smiles:
                return c.canonical_smiles

            # If no SMILES, try to get full compound with InChI and convert
            full_compound = pcp.Compound.from_cid(c.cid)
            if full_compound.inchi:
                mol = Chem.MolFromInchi(full_compound.inchi)
                if mol:
                    return Chem.MolToSmiles(mol)

            # Last resort: try molecular formula lookup
            logger.warning(f"Could not get SMILES for {name} (CID: {c.cid})")
            return None

        except Exception as e:
            logger.warning(f"PubChem lookup failed for {name}: {e}")
            return None

    async def _execute_pipeline(
        self,
        job: DesignJob,
        parsed: Dict[str, Any]
    ) -> List[DesignedMolecule]:
        """Execute the full design pipeline."""
        molecules = []

        design_type = parsed["design_type"]

        # Step 1: Generate candidate molecules
        job.update_progress("Generating candidate molecules", 0.1)

        if design_type == DesignType.ANALOGUE:
            if parsed.get("reference_smiles"):
                candidates = await self.generator.generate_analogues(
                    parsed["reference_smiles"],
                    n_molecules=20,
                    constraints=parsed.get("constraints", {})
                )
            else:
                raise ValueError("Reference compound required for analogue design")

        elif design_type == DesignType.BIVALENT:
            # For bivalent, need two pharmacophores
            # This is simplified - real implementation would identify pharmacophores
            candidates = await self._generate_bivalent_candidates(parsed)

        elif design_type == DesignType.SCAFFOLD_HOP:
            if parsed.get("reference_smiles"):
                candidates = await self.generator.scaffold_hop(
                    parsed["reference_smiles"],
                    n_molecules=20
                )
            else:
                raise ValueError("Reference compound required for scaffold hopping")

        else:
            candidates = []

        if not candidates:
            logger.warning("No candidates generated")
            return []

        # Step 2: Analyze each candidate
        job.update_progress("Analyzing candidates", 0.2)

        for i, smiles in enumerate(candidates):
            progress = 0.2 + (0.6 * (i + 1) / len(candidates))
            job.update_progress(f"Analyzing molecule {i+1}/{len(candidates)}", progress)

            molecule = await self._analyze_molecule(smiles, parsed)
            if molecule:
                molecules.append(molecule)

        # Step 3: Rank molecules
        job.update_progress("Ranking molecules", 0.85)
        molecules = self._rank_molecules(molecules, parsed.get("constraints", {}))

        # Step 4: Get synthesis routes for top candidates
        job.update_progress("Planning synthesis routes", 0.9)
        for mol in molecules[:5]:  # Top 5
            routes = await self.retrosynthesis.plan_synthesis(mol.smiles)
            mol.synthesis_routes = routes
            if routes:
                mol.best_route = routes[0]

        return molecules

    async def _analyze_molecule(
        self,
        smiles: str,
        parsed: Dict[str, Any]
    ) -> Optional[DesignedMolecule]:
        """Perform full analysis of a candidate molecule."""
        try:
            # Validate
            validation = self.analyzer.validate_molecule(smiles)
            if not validation["valid"]:
                return None

            canonical_smiles = validation["canonical_smiles"]

            # Get InChI
            inchi, inchi_key = self.analyzer.smiles_to_inchi(canonical_smiles)

            # Calculate properties
            properties = self.analyzer.calculate_properties(canonical_smiles)

            # ADMET prediction
            admet = await self.admet.predict(canonical_smiles)

            # Novelty assessment
            known_compounds = []
            if parsed.get("reference_smiles"):
                known_compounds.append({
                    "smiles": parsed["reference_smiles"],
                    "name": parsed.get("reference_name", "Reference")
                })

            novelty = self.analyzer.assess_novelty(canonical_smiles, known_compounds)

            # Patent check
            patent_result = await self.patent.check_novelty(canonical_smiles)

            # Update novelty with patent info
            novelty.exact_patent_match = bool(patent_result.get("exact_matches"))
            novelty.patentable = novelty.novel and not novelty.exact_patent_match

            # Generate image
            image_svg = self.analyzer.molecule_to_svg(canonical_smiles)

            # Generate 3D structure
            mol_block = self.analyzer.generate_conformer(canonical_smiles)

            # Create molecule object
            molecule = DesignedMolecule(
                smiles=canonical_smiles,
                inchi=inchi,
                inchi_key=inchi_key,
                mol_block=mol_block,
                image_svg=image_svg,
                properties=properties,
                admet=admet,
                novelty=novelty,
            )

            # Calculate overall score
            molecule.overall_score = self._calculate_score(molecule, parsed)

            return molecule

        except Exception as e:
            logger.error(f"Error analyzing {smiles}: {e}")
            return None

    async def _generate_bivalent_candidates(
        self,
        parsed: Dict[str, Any]
    ) -> List[str]:
        """Generate bivalent ligand candidates."""
        # This is a simplified implementation
        # Real version would identify appropriate pharmacophores based on targets

        target = parsed.get("target_name", "").lower()

        # Example pharmacophores for 5HT1A-5HT2A
        if "5ht1a" in target or "5ht2a" in target:
            # 5HT1A agonist scaffold (simplified buspirone-like)
            pharm1 = "c1ccc2c(c1)CCCC2"  # Tetralin-like

            # 5HT2A agonist scaffold (simplified tryptamine)
            pharm2 = "NCCc1c[nH]c2ccccc12"  # Tryptamine

            return await self.generator.generate_bivalent_ligand(
                pharm1, pharm2,
                linker_length=6,
                linker_type="alkyl"
            )

        return []

    def _calculate_score(
        self,
        molecule: DesignedMolecule,
        parsed: Dict[str, Any]
    ) -> float:
        """Calculate overall score for a molecule."""
        score = 0.0
        weights = {
            "novelty": 0.25,
            "properties": 0.20,
            "admet": 0.25,
            "patentability": 0.30
        }

        # Novelty score
        if molecule.novelty:
            score += weights["novelty"] * molecule.novelty.novelty_score

        # Properties score (drug-likeness)
        if molecule.properties:
            props_score = molecule.properties.qed
            if molecule.properties.passes_lipinski:
                props_score += 0.2
            score += weights["properties"] * min(1.0, props_score)

        # ADMET score
        if molecule.admet:
            admet_score = molecule.admet.get_safety_score()
            score += weights["admet"] * admet_score

        # Patentability score
        if molecule.novelty:
            patent_score = 1.0 if molecule.novelty.patentable else 0.3
            score += weights["patentability"] * patent_score

        return round(score, 3)

    def _rank_molecules(
        self,
        molecules: List[DesignedMolecule],
        constraints: Dict[str, Any]
    ) -> List[DesignedMolecule]:
        """Rank molecules by overall score, applying constraint penalties."""
        # Apply constraint penalties instead of hard filtering
        for mol in molecules:
            penalty = 0.0

            if constraints.get("bbb_penetration"):
                if mol.admet and not mol.admet.bbb_penetration:
                    # Penalize but don't exclude - some bivalent ligands can still cross BBB
                    penalty += 0.1

            if constraints.get("max_mw") and mol.properties:
                if mol.properties.molecular_weight > constraints["max_mw"]:
                    # Penalize for exceeding MW
                    penalty += 0.15

            mol.overall_score = max(0, mol.overall_score - penalty)

        # Sort by score
        molecules.sort(key=lambda m: m.overall_score, reverse=True)

        # Assign ranks
        for i, mol in enumerate(molecules):
            mol.rank = i + 1

        return molecules

    async def get_job(self, job_id: str) -> Optional[DesignJob]:
        """Retrieve a job by ID."""
        return self.jobs.get(job_id)

    async def generate_report(self, job: DesignJob) -> str:
        """Generate comprehensive report for a completed job."""
        lines = [
            "=" * 80,
            "DRUG DESIGN REPORT",
            "=" * 80,
            f"Job ID: {job.id}",
            f"Request: {job.prompt}",
            f"Design Type: {job.design_type.value}",
            f"Status: {job.status.value}",
            f"Created: {job.created_at.isoformat()}",
            f"Completed: {job.completed_at.isoformat() if job.completed_at else 'N/A'}",
            "",
        ]

        if job.reference_name:
            lines.append(f"Reference Compound: {job.reference_name}")
        if job.reference_smiles:
            lines.append(f"Reference SMILES: {job.reference_smiles}")

        lines.extend([
            "",
            f"RESULTS: {len(job.molecules)} molecules designed",
            "-" * 60,
            ""
        ])

        for mol in job.molecules[:10]:  # Top 10
            lines.extend([
                f"Rank #{mol.rank}: Score = {mol.overall_score:.3f}",
                f"  SMILES: {mol.smiles}",
                f"  InChIKey: {mol.inchi_key or 'N/A'}",
            ])

            if mol.properties:
                lines.append(
                    f"  Properties: MW={mol.properties.molecular_weight:.1f}, "
                    f"LogP={mol.properties.logp:.2f}, QED={mol.properties.qed:.3f}"
                )

            if mol.novelty:
                lines.append(
                    f"  Novelty: {mol.novelty.tanimoto_similarity:.2f} Tanimoto to nearest, "
                    f"Patentable={'Yes' if mol.novelty.patentable else 'No'}"
                )

            if mol.admet:
                lines.append(
                    f"  ADMET Safety: {mol.admet.get_safety_score() * 100:.0f}%"
                )

            if mol.best_route:
                lines.append(
                    f"  Synthesis: {mol.best_route.total_steps} steps, "
                    f"{mol.best_route.overall_yield * 100:.0f}% yield"
                )

            lines.append("")

        # Best molecule details
        if job.best_molecule:
            lines.extend([
                "=" * 60,
                "BEST CANDIDATE DETAILS",
                "=" * 60,
                f"SMILES: {job.best_molecule.smiles}",
                "",
            ])

            if job.best_molecule.best_route:
                lines.append(
                    self.retrosynthesis.format_synthesis_report(job.best_molecule.best_route)
                )

            if job.best_molecule.admet:
                lines.append(
                    self.admet.format_admet_report(job.best_molecule.admet)
                )

            if job.best_molecule.novelty:
                patent_result = {
                    "novel": job.best_molecule.novelty.novel,
                    "exact_matches": [],
                    "confidence": job.best_molecule.novelty.confidence
                }
                fto_result = {
                    "risk_level": job.best_molecule.novelty.fto_risk,
                    "blocking_patents": job.best_molecule.novelty.blocking_patents,
                    "recommendations": []
                }
                patentability = self.patent.assess_patentability(patent_result)

                lines.append(
                    self.patent.generate_patent_report(
                        job.best_molecule.smiles,
                        patent_result,
                        fto_result,
                        patentability
                    )
                )

        lines.extend([
            "",
            "=" * 80,
            "END OF REPORT",
            "=" * 80
        ])

        return "\n".join(lines)
