"""
Claude AI Service for intelligent drug design decisions
"""
import anthropic
from typing import List, Dict, Any, Optional
from loguru import logger
import json

from config import settings


class ClaudeService:
    """Service for interacting with Claude AI for drug design intelligence"""

    def __init__(self):
        self.client = anthropic.Anthropic(api_key=settings.anthropic_api_key)
        self.model = settings.claude_model

    async def suggest_modifications(
        self,
        compound_info: Dict[str, Any],
        target_receptor: str,
        num_suggestions: int = 10
    ) -> List[Dict[str, Any]]:
        """
        Use Claude to suggest intelligent modifications for analogue generation
        """
        prompt = f"""You are an expert medicinal chemist specializing in peptide drug design and opioid receptor pharmacology.

Given the following compound information:
- Name: {compound_info.get('name', 'Unknown')}
- Structure/Sequence: {compound_info.get('structure', 'Unknown')}
- Target Receptor: {target_receptor}
- Known Properties: {json.dumps(compound_info.get('properties', {}), indent=2)}

Generate {num_suggestions} novel modifications that could:
1. Improve metabolic stability
2. Maintain or improve binding affinity to {target_receptor}
3. Result in a patentable (novel, non-obvious) compound
4. Be synthetically accessible

For each modification, provide:
1. modification_type: The category (n_terminal, c_terminal, amino_acid_substitution, backbone, cyclization)
2. description: Clear description of the modification
3. rationale: Scientific rationale for why this modification might work
4. expected_effect: Predicted effect on properties
5. synthetic_accessibility: Estimated difficulty (1-10)
6. novelty_assessment: Why this would be novel/patentable

Return your response as a JSON array of modification objects."""

        try:
            response = self.client.messages.create(
                model=self.model,
                max_tokens=settings.claude_max_tokens,
                messages=[{"role": "user", "content": prompt}]
            )

            # Parse the response
            content = response.content[0].text

            # Extract JSON from response
            json_start = content.find('[')
            json_end = content.rfind(']') + 1
            if json_start != -1 and json_end > json_start:
                modifications = json.loads(content[json_start:json_end])
                return modifications

            logger.warning("Could not parse Claude response as JSON")
            return []

        except Exception as e:
            logger.error(f"Error calling Claude API: {e}")
            return []

    async def assess_patentability(
        self,
        compound_smiles: str,
        modification_description: str,
        similar_compounds: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """
        Use Claude to assess patentability of a novel compound
        """
        prompt = f"""You are a patent attorney specializing in pharmaceutical patents.

Assess the patentability of this novel compound:
- SMILES: {compound_smiles}
- Modification from parent: {modification_description}

Similar known compounds found:
{json.dumps(similar_compounds, indent=2)}

Evaluate:
1. NOVELTY: Is this compound structurally distinct from known compounds?
2. NON-OBVIOUSNESS: Would this modification be obvious to a skilled medicinal chemist?
3. UTILITY: Does this have clear pharmaceutical utility?

Provide your assessment as JSON:
{{
    "novelty_score": 0.0-1.0,
    "non_obviousness_score": 0.0-1.0,
    "utility_score": 0.0-1.0,
    "overall_patentability": 0.0-1.0,
    "key_claims": ["claim1", "claim2"],
    "potential_issues": ["issue1"],
    "recommendation": "string"
}}"""

        try:
            response = self.client.messages.create(
                model=self.model,
                max_tokens=2048,
                messages=[{"role": "user", "content": prompt}]
            )

            content = response.content[0].text
            json_start = content.find('{')
            json_end = content.rfind('}') + 1
            if json_start != -1 and json_end > json_start:
                return json.loads(content[json_start:json_end])

            return {"overall_patentability": 0.5, "recommendation": "Manual review needed"}

        except Exception as e:
            logger.error(f"Error in patentability assessment: {e}")
            return {"overall_patentability": 0.5, "recommendation": "Assessment failed"}

    async def generate_synthesis_route(
        self,
        compound_smiles: str,
        compound_name: str,
        is_peptide: bool = True
    ) -> Dict[str, Any]:
        """
        Use Claude to suggest a synthesis route for a compound
        """
        prompt = f"""You are a synthetic organic chemist specializing in peptide synthesis.

Design a synthesis route for:
- Name: {compound_name}
- SMILES: {compound_smiles}
- Type: {"Peptide/Peptidomimetic" if is_peptide else "Small molecule"}

Provide:
1. Recommended synthesis strategy (solid-phase vs solution-phase for peptides)
2. Key building blocks needed
3. Step-by-step synthesis outline
4. Potential challenges and solutions
5. Estimated number of steps
6. Special equipment or conditions needed

Return as JSON:
{{
    "strategy": "string",
    "building_blocks": ["block1", "block2"],
    "steps": [
        {{"step_number": 1, "description": "string", "conditions": "string"}}
    ],
    "challenges": ["challenge1"],
    "estimated_steps": number,
    "special_requirements": ["req1"],
    "overall_difficulty": "Easy/Moderate/Difficult/Very Difficult"
}}"""

        try:
            response = self.client.messages.create(
                model=self.model,
                max_tokens=2048,
                messages=[{"role": "user", "content": prompt}]
            )

            content = response.content[0].text
            json_start = content.find('{')
            json_end = content.rfind('}') + 1
            if json_start != -1 and json_end > json_start:
                return json.loads(content[json_start:json_end])

            return {"overall_difficulty": "Unknown", "steps": []}

        except Exception as e:
            logger.error(f"Error generating synthesis route: {e}")
            return {"overall_difficulty": "Unknown", "steps": []}

    async def rank_candidates(
        self,
        candidates: List[Dict[str, Any]],
        optimization_goals: List[str]
    ) -> List[Dict[str, Any]]:
        """
        Use Claude to intelligently rank analogue candidates
        """
        prompt = f"""You are an expert in drug discovery and lead optimization.

Rank these analogue candidates based on the following optimization goals:
{json.dumps(optimization_goals, indent=2)}

Candidates:
{json.dumps(candidates, indent=2)}

For each candidate, provide:
1. Overall score (0-1)
2. Strengths
3. Weaknesses
4. Recommendation (Prioritize/Consider/Deprioritize)

Return as JSON array with added ranking fields for each candidate."""

        try:
            response = self.client.messages.create(
                model=self.model,
                max_tokens=4096,
                messages=[{"role": "user", "content": prompt}]
            )

            content = response.content[0].text
            json_start = content.find('[')
            json_end = content.rfind(']') + 1
            if json_start != -1 and json_end > json_start:
                return json.loads(content[json_start:json_end])

            return candidates

        except Exception as e:
            logger.error(f"Error ranking candidates: {e}")
            return candidates

    async def explain_modification(
        self,
        original_smiles: str,
        modified_smiles: str,
        property_changes: Dict[str, Any]
    ) -> str:
        """
        Generate a human-readable explanation of a modification
        """
        prompt = f"""Explain this molecular modification in clear, concise terms:

Original: {original_smiles}
Modified: {modified_smiles}
Property changes: {json.dumps(property_changes, indent=2)}

Provide a 2-3 sentence explanation suitable for a medicinal chemist, focusing on:
1. What structural change was made
2. Why this change might improve drug properties
3. Any trade-offs to consider"""

        try:
            response = self.client.messages.create(
                model=self.model,
                max_tokens=500,
                messages=[{"role": "user", "content": prompt}]
            )
            return response.content[0].text

        except Exception as e:
            logger.error(f"Error explaining modification: {e}")
            return "Modification explanation unavailable."
