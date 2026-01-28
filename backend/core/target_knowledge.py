"""
Target Knowledge Extraction and Management (Phase 2/3)

This module provides intelligent extraction and management of target-specific
knowledge, enabling truly target-agnostic drug design.

Design Philosophy:
- Extract knowledge dynamically from multiple sources
- Adapt design constraints to the specific target
- Handle uncertainty about target biology honestly
- Support any drug target without hardcoded classes
"""

import asyncio
import logging
import json
import os
from typing import Dict, List, Optional, Any, Set
from dataclasses import dataclass, field
from enum import Enum
from datetime import datetime

import httpx
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

logger = logging.getLogger(__name__)


class TargetType(Enum):
    """Classification of drug target types."""
    ENZYME = "enzyme"
    GPCR = "gpcr"
    ION_CHANNEL = "ion_channel"
    NUCLEAR_RECEPTOR = "nuclear_receptor"
    KINASE = "kinase"
    PROTEASE = "protease"
    TRANSPORTER = "transporter"
    PROTEIN_PROTEIN = "protein_protein_interaction"
    RNA = "rna"
    DNA = "dna"
    UNKNOWN = "unknown"


class BindingMechanism(Enum):
    """Types of binding mechanisms."""
    REVERSIBLE_COMPETITIVE = "reversible_competitive"
    REVERSIBLE_ALLOSTERIC = "reversible_allosteric"
    REVERSIBLE_UNCOMPETITIVE = "reversible_uncompetitive"
    COVALENT_IRREVERSIBLE = "covalent_irreversible"
    COVALENT_REVERSIBLE = "covalent_reversible"
    SUBSTRATE_COMPETITIVE = "substrate_competitive"
    UNKNOWN = "unknown"


@dataclass
class StructuralFeature:
    """A structural feature relevant to target binding."""
    name: str
    description: str
    smarts_pattern: Optional[str] = None
    importance: str = "recommended"  # "required", "recommended", "optional", "avoid"
    rationale: str = ""


@dataclass
class PropertyRange:
    """Optimal property range for a target."""
    property_name: str
    min_value: Optional[float] = None
    max_value: Optional[float] = None
    optimal_value: Optional[float] = None
    rationale: str = ""
    confidence: str = "moderate"  # "high", "moderate", "low"


@dataclass
class ReferenceCompound:
    """A known active compound for the target."""
    name: str
    smiles: str
    activity_type: str = ""  # "IC50", "Ki", "EC50", etc.
    activity_value: Optional[float] = None  # in nM
    activity_unit: str = "nM"
    source: str = ""  # "ChEMBL", "literature", "LLM"
    notes: str = ""

    @property
    def pIC50(self) -> Optional[float]:
        """Convert activity to pIC50 scale."""
        if self.activity_value and self.activity_value > 0:
            return -1 * (self.activity_value / 1e9).__log10__()
        return None


@dataclass
class TargetKnowledgeBase:
    """
    Comprehensive knowledge about a drug target.

    This aggregates information from multiple sources:
    - LLM reasoning about target biology
    - ChEMBL bioactivity data
    - PubChem compound information
    - Literature-derived insights
    """
    # Basic identification
    target_name: str = ""
    target_aliases: List[str] = field(default_factory=list)
    uniprot_id: Optional[str] = None
    chembl_id: Optional[str] = None

    # Classification
    target_type: TargetType = TargetType.UNKNOWN
    organism: str = "Homo sapiens"
    disease_relevance: List[str] = field(default_factory=list)

    # Binding mechanism
    binding_mechanism: BindingMechanism = BindingMechanism.UNKNOWN
    binding_site_residues: List[str] = field(default_factory=list)
    key_interactions: List[str] = field(default_factory=list)

    # Covalent inhibitor specifics
    is_covalent: bool = False
    warhead_types: List[str] = field(default_factory=list)
    nucleophile_residue: Optional[str] = None
    staying_portion_features: List[str] = field(default_factory=list)
    leaving_group_features: List[str] = field(default_factory=list)

    # Structural requirements
    required_features: List[StructuralFeature] = field(default_factory=list)
    recommended_features: List[StructuralFeature] = field(default_factory=list)
    features_to_avoid: List[StructuralFeature] = field(default_factory=list)

    # Property constraints
    property_ranges: List[PropertyRange] = field(default_factory=list)

    # Reference compounds
    reference_compounds: List[ReferenceCompound] = field(default_factory=list)
    scaffold_templates: List[str] = field(default_factory=list)  # SMILES

    # SAR insights
    sar_insights: List[str] = field(default_factory=list)

    # Design strategy
    design_strategy: str = ""
    novelty_strategies: List[str] = field(default_factory=list)

    # Metadata
    knowledge_sources: List[str] = field(default_factory=list)
    extraction_timestamp: datetime = field(default_factory=datetime.utcnow)
    confidence_assessment: str = ""

    def get_optimal_mw_range(self) -> tuple:
        """Get optimal molecular weight range."""
        for pr in self.property_ranges:
            if pr.property_name in ["molecular_weight", "mw"]:
                return (pr.min_value or 200, pr.max_value or 600)
        return (300, 550)  # Default

    def get_optimal_logp_range(self) -> tuple:
        """Get optimal LogP range."""
        for pr in self.property_ranges:
            if pr.property_name in ["logp", "clogp"]:
                return (pr.min_value or 0, pr.max_value or 5)
        return (1, 5)  # Default

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "target_name": self.target_name,
            "target_aliases": self.target_aliases,
            "uniprot_id": self.uniprot_id,
            "chembl_id": self.chembl_id,
            "target_type": self.target_type.value,
            "binding_mechanism": self.binding_mechanism.value,
            "is_covalent": self.is_covalent,
            "warhead_types": self.warhead_types,
            "nucleophile_residue": self.nucleophile_residue,
            "required_features": [
                {"name": f.name, "description": f.description, "smarts": f.smarts_pattern}
                for f in self.required_features
            ],
            "reference_compounds": [
                {"name": r.name, "smiles": r.smiles, "activity": r.activity_value}
                for r in self.reference_compounds
            ],
            "design_strategy": self.design_strategy,
            "sar_insights": self.sar_insights,
        }


class TargetKnowledgeExtractor:
    """
    Extracts comprehensive knowledge about a drug target.

    This is the core of the target-agnostic design system. It:
    1. Parses the user request to identify the target
    2. Queries multiple data sources (ChEMBL, PubChem, LLM)
    3. Synthesizes a comprehensive knowledge base
    4. Adapts design constraints to the specific target
    """

    def __init__(self, anthropic_api_key: str = None):
        self.api_key = anthropic_api_key or os.getenv("ANTHROPIC_API_KEY")
        self._chembl_client = None
        self._http_client = None

    async def _get_http_client(self) -> httpx.AsyncClient:
        """Get or create HTTP client."""
        if self._http_client is None:
            self._http_client = httpx.AsyncClient(timeout=60.0)
        return self._http_client

    async def extract_knowledge(
        self,
        request: str,
        existing_knowledge: Optional[TargetKnowledgeBase] = None
    ) -> TargetKnowledgeBase:
        """
        Extract comprehensive target knowledge from a design request.

        Args:
            request: User's drug design request
            existing_knowledge: Optional existing knowledge to augment

        Returns:
            TargetKnowledgeBase with extracted information
        """
        knowledge = existing_knowledge or TargetKnowledgeBase()

        # Step 1: LLM-based target identification and initial knowledge
        logger.info("Extracting target knowledge from request...")
        await self._extract_from_llm(request, knowledge)

        # Step 2: ChEMBL database lookup
        if knowledge.target_name:
            logger.info(f"Querying ChEMBL for {knowledge.target_name}...")
            await self._extract_from_chembl(knowledge)

        # Step 3: Mechanism-specific knowledge
        if knowledge.is_covalent:
            logger.info("Extracting covalent mechanism details...")
            await self._extract_covalent_mechanism(knowledge)

        # Step 4: Synthesize design constraints
        logger.info("Synthesizing design constraints...")
        await self._synthesize_design_constraints(knowledge)

        knowledge.extraction_timestamp = datetime.utcnow()
        return knowledge

    async def _extract_from_llm(
        self,
        request: str,
        knowledge: TargetKnowledgeBase
    ) -> None:
        """Extract target knowledge using LLM reasoning."""
        prompt = f"""You are an expert medicinal chemist. Analyze this drug design request and extract target knowledge.

REQUEST: "{request}"

Identify and provide detailed information about:

1. TARGET IDENTIFICATION
   - What is the biological target?
   - Any aliases or alternative names?
   - UniProt ID if known?

2. TARGET CLASSIFICATION
   - What type of target? (enzyme, GPCR, kinase, ion_channel, nuclear_receptor, protease, transporter, protein_protein_interaction, rna, dna, unknown)
   - What organism?
   - What diseases is it relevant to?

3. BINDING MECHANISM
   - How do inhibitors/ligands bind? (reversible_competitive, reversible_allosteric, covalent_irreversible, covalent_reversible, substrate_competitive, unknown)
   - What residues are in the binding site?
   - What molecular interactions are key? (H-bonds, pi-stacking, hydrophobic, covalent, etc.)

4. FOR COVALENT INHIBITORS (if applicable)
   - Which nucleophile residue attacks?
   - What warhead types work? (acrylamide, urea, carbamate, nitrile, etc.)
   - What stays bound after reaction?
   - What leaves as leaving group?

5. STRUCTURAL REQUIREMENTS
   - What features are REQUIRED for activity?
   - What features are RECOMMENDED?
   - What features should be AVOIDED?
   - Include SMARTS patterns if you know them

6. OPTIMAL PROPERTIES
   - What MW range is optimal and why?
   - What LogP range is optimal and why?
   - Other important properties?

7. KNOWN ACTIVE COMPOUNDS
   - Name clinical candidates or approved drugs
   - Include SMILES if you know them
   - Note their activity values if known

8. SAR INSIGHTS
   - What structural changes improve activity?
   - What are common SAR trends?

Return as JSON:
{{
    "target_name": "official name",
    "target_aliases": ["list of aliases"],
    "uniprot_id": "if known",
    "target_type": "enzyme|gpcr|kinase|...",
    "organism": "Homo sapiens",
    "disease_relevance": ["diseases"],
    "binding_mechanism": "reversible_competitive|covalent_irreversible|...",
    "binding_site_residues": ["Ser241", "His469", ...],
    "key_interactions": ["hydrophobic pocket", "H-bond to backbone", ...],
    "is_covalent": true/false,
    "warhead_types": ["urea", "carbamate", ...],
    "nucleophile_residue": "Ser241",
    "staying_portion_features": ["lipophilic tail", "aromatic for pi-stacking"],
    "leaving_group_features": ["low pKa amine", "aromatic"],
    "required_features": [
        {{"name": "feature", "description": "why needed", "smarts": "pattern or null"}}
    ],
    "recommended_features": [...],
    "features_to_avoid": [...],
    "property_ranges": [
        {{"property": "molecular_weight", "min": 300, "max": 550, "rationale": "why"}}
    ],
    "known_compounds": [
        {{"name": "PF-04457845", "smiles": "...", "activity_nM": 10, "activity_type": "IC50"}}
    ],
    "sar_insights": ["larger lipophilic groups improve potency", ...],
    "design_strategy": "overall approach for this target"
}}

Base your analysis on chemical first principles and specific target biology."""

        response = await self._call_claude(prompt)

        try:
            data = self._parse_json(response)

            # Populate knowledge base
            knowledge.target_name = data.get("target_name", "")
            knowledge.target_aliases = data.get("target_aliases", [])
            knowledge.uniprot_id = data.get("uniprot_id")

            # Target type
            target_type_str = data.get("target_type", "unknown").lower()
            for tt in TargetType:
                if tt.value == target_type_str:
                    knowledge.target_type = tt
                    break

            # Binding mechanism
            mechanism_str = data.get("binding_mechanism", "unknown").lower()
            for bm in BindingMechanism:
                if bm.value == mechanism_str:
                    knowledge.binding_mechanism = bm
                    break

            knowledge.organism = data.get("organism", "Homo sapiens")
            knowledge.disease_relevance = data.get("disease_relevance", [])
            knowledge.binding_site_residues = data.get("binding_site_residues", [])
            knowledge.key_interactions = data.get("key_interactions", [])

            # Covalent specifics
            knowledge.is_covalent = data.get("is_covalent", False)
            knowledge.warhead_types = data.get("warhead_types", [])
            knowledge.nucleophile_residue = data.get("nucleophile_residue")
            knowledge.staying_portion_features = data.get("staying_portion_features", [])
            knowledge.leaving_group_features = data.get("leaving_group_features", [])

            # Structural features
            for feat in data.get("required_features", []):
                knowledge.required_features.append(StructuralFeature(
                    name=feat.get("name", ""),
                    description=feat.get("description", ""),
                    smarts_pattern=feat.get("smarts"),
                    importance="required",
                ))

            for feat in data.get("recommended_features", []):
                knowledge.recommended_features.append(StructuralFeature(
                    name=feat.get("name", ""),
                    description=feat.get("description", ""),
                    smarts_pattern=feat.get("smarts"),
                    importance="recommended",
                ))

            for feat in data.get("features_to_avoid", []):
                knowledge.features_to_avoid.append(StructuralFeature(
                    name=feat.get("name", ""),
                    description=feat.get("description", ""),
                    smarts_pattern=feat.get("smarts"),
                    importance="avoid",
                ))

            # Property ranges
            for pr in data.get("property_ranges", []):
                knowledge.property_ranges.append(PropertyRange(
                    property_name=pr.get("property", ""),
                    min_value=pr.get("min"),
                    max_value=pr.get("max"),
                    rationale=pr.get("rationale", ""),
                ))

            # Reference compounds
            for comp in data.get("known_compounds", []):
                smiles = comp.get("smiles", "")
                if smiles and Chem.MolFromSmiles(smiles):
                    knowledge.reference_compounds.append(ReferenceCompound(
                        name=comp.get("name", ""),
                        smiles=Chem.MolToSmiles(Chem.MolFromSmiles(smiles)),
                        activity_type=comp.get("activity_type", "IC50"),
                        activity_value=comp.get("activity_nM"),
                        source="LLM",
                    ))

            knowledge.sar_insights = data.get("sar_insights", [])
            knowledge.design_strategy = data.get("design_strategy", "")

            knowledge.knowledge_sources.append("LLM_reasoning")

        except Exception as e:
            logger.warning(f"Failed to parse LLM response: {e}")

    async def _extract_from_chembl(self, knowledge: TargetKnowledgeBase) -> None:
        """Extract knowledge from ChEMBL database."""
        try:
            client = await self._get_http_client()

            # Search for target
            search_url = "https://www.ebi.ac.uk/chembl/api/data/target/search.json"
            params = {"q": knowledge.target_name, "limit": 5}

            response = await client.get(search_url, params=params)
            if response.status_code != 200:
                return

            data = response.json()
            targets = data.get("targets", [])

            if not targets:
                return

            # Find best match
            target = targets[0]
            knowledge.chembl_id = target.get("target_chembl_id")

            # Get bioactivity data
            if knowledge.chembl_id:
                activity_url = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
                activity_params = {
                    "target_chembl_id": knowledge.chembl_id,
                    "standard_type__in": "IC50,Ki,EC50",
                    "limit": 100,
                }

                activity_response = await client.get(activity_url, params=activity_params)
                if activity_response.status_code == 200:
                    activity_data = activity_response.json()
                    activities = activity_data.get("activities", [])

                    # Process activities
                    seen_smiles = set()
                    for act in activities[:20]:  # Top 20
                        smiles = act.get("canonical_smiles", "")
                        if not smiles or smiles in seen_smiles:
                            continue

                        mol = Chem.MolFromSmiles(smiles)
                        if mol is None:
                            continue

                        seen_smiles.add(smiles)

                        value = act.get("standard_value")
                        if value:
                            try:
                                value = float(value)
                            except:
                                value = None

                        knowledge.reference_compounds.append(ReferenceCompound(
                            name=act.get("molecule_pref_name", act.get("molecule_chembl_id", "")),
                            smiles=Chem.MolToSmiles(mol),
                            activity_type=act.get("standard_type", "IC50"),
                            activity_value=value,
                            activity_unit=act.get("standard_units", "nM"),
                            source="ChEMBL",
                        ))

            knowledge.knowledge_sources.append("ChEMBL")

        except Exception as e:
            logger.warning(f"ChEMBL extraction failed: {e}")

    async def _extract_covalent_mechanism(self, knowledge: TargetKnowledgeBase) -> None:
        """Extract detailed covalent mechanism information."""
        if not knowledge.is_covalent or not knowledge.reference_compounds:
            return

        # Analyze reference compounds for covalent patterns
        from services.mechanistic_analyzer import mechanistic_analyzer

        for ref in knowledge.reference_compounds[:5]:
            try:
                analysis = mechanistic_analyzer.analyze_mechanism(ref.smiles)
                if analysis.is_covalent:
                    # Update warhead types if found new one
                    warhead = analysis.warhead_type.value
                    if warhead and warhead not in knowledge.warhead_types:
                        knowledge.warhead_types.append(warhead)

                    ref.notes = f"Warhead: {warhead}"
            except Exception as e:
                logger.warning(f"Covalent analysis failed for {ref.name}: {e}")

    async def _synthesize_design_constraints(self, knowledge: TargetKnowledgeBase) -> None:
        """Synthesize all knowledge into actionable design constraints."""
        # Ensure we have property ranges
        if not any(pr.property_name == "molecular_weight" for pr in knowledge.property_ranges):
            # Add default based on target type
            if knowledge.target_type == TargetType.KINASE:
                knowledge.property_ranges.append(PropertyRange(
                    property_name="molecular_weight",
                    min_value=350,
                    max_value=550,
                    rationale="Kinase inhibitors typically 350-550 Da",
                ))
            elif knowledge.target_type == TargetType.GPCR:
                knowledge.property_ranges.append(PropertyRange(
                    property_name="molecular_weight",
                    min_value=300,
                    max_value=500,
                    rationale="GPCR ligands typically 300-500 Da",
                ))
            else:
                knowledge.property_ranges.append(PropertyRange(
                    property_name="molecular_weight",
                    min_value=300,
                    max_value=550,
                    rationale="Standard drug-like range",
                ))

        # Generate confidence assessment
        sources = len(knowledge.knowledge_sources)
        refs = len(knowledge.reference_compounds)

        if sources >= 2 and refs >= 5:
            knowledge.confidence_assessment = "HIGH: Multiple data sources with abundant reference compounds"
        elif sources >= 1 and refs >= 2:
            knowledge.confidence_assessment = "MODERATE: Some data available but could be more comprehensive"
        else:
            knowledge.confidence_assessment = "LOW: Limited data - predictions should be treated with caution"

    async def _call_claude(self, prompt: str) -> str:
        """Call Claude API."""
        if not self.api_key:
            return "{}"

        try:
            client = await self._get_http_client()
            response = await client.post(
                "https://api.anthropic.com/v1/messages",
                headers={
                    "x-api-key": self.api_key,
                    "anthropic-version": "2023-06-01",
                    "content-type": "application/json"
                },
                json={
                    "model": "claude-sonnet-4-20250514",
                    "max_tokens": 4096,
                    "messages": [{"role": "user", "content": prompt}]
                }
            )

            if response.status_code == 200:
                data = response.json()
                return data["content"][0]["text"]
            else:
                logger.warning(f"Claude API error: {response.status_code}")
                return "{}"
        except Exception as e:
            logger.error(f"Claude API call failed: {e}")
            return "{}"

    def _parse_json(self, response: str) -> Dict:
        """Parse JSON from response."""
        content = response.strip()

        # Remove markdown code blocks
        if "```json" in content:
            start = content.find("```json") + 7
            end = content.find("```", start)
            if end > start:
                content = content[start:end].strip()
        elif content.startswith("```"):
            lines = content.split("\n")[1:]
            if lines and lines[-1].strip() == "```":
                lines = lines[:-1]
            content = "\n".join(lines)

        # Find JSON start
        for i, char in enumerate(content):
            if char in "{[":
                content = content[i:]
                break

        try:
            decoder = json.JSONDecoder()
            result, _ = decoder.raw_decode(content)
            return result
        except:
            return {}


# Global instance
target_knowledge_extractor = TargetKnowledgeExtractor()
