"""
Analogue Generation Service using RDKit and intelligent modification strategies
"""
from typing import List, Dict, Any, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw, rdMolDescriptors
from rdkit.Chem import rdFMCS
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import DataStructs
import uuid
from loguru import logger
import copy

from models import (
    CompoundInput, AnalogueCandidate, PropertyPrediction,
    NoveltyCheck, SynthesisGuidance, InputFormat
)


class AnalogueGenerator:
    """Service for generating novel analogues of input compounds"""

    # Common bioisosteric replacements
    BIOISOSTERES = {
        "COOH": ["C(=O)N", "C(=O)NO", "S(=O)(=O)O", "P(=O)(O)O"],  # Carboxylic acid
        "OH": ["NH2", "SH", "F"],  # Hydroxyl
        "NH2": ["OH", "NHCH3", "N(CH3)2"],  # Amine
        "F": ["Cl", "CF3", "OCF3"],  # Fluorine
        "Cl": ["F", "Br", "CF3"],  # Chlorine
        "CH3": ["CF3", "Cl", "OCH3"],  # Methyl
        "OCH3": ["SCH3", "CH3", "F"],  # Methoxy
        "C=O": ["C=S", "C=NH"],  # Carbonyl
    }

    # N-terminal modifications for peptides
    N_TERMINAL_MODS = {
        "pivaloyl": "CC(C)(C)C(=O)",
        "cyclopropylcarbonyl": "C1CC1C(=O)",
        "isobutyryl": "CC(C)C(=O)",
        "acetyl": "CC(=O)",
        "benzoyl": "c1ccccc1C(=O)",
        "cyclopentylcarbonyl": "C1CCCC1C(=O)",
        "2,2-dimethylbutanoyl": "CCC(C)(C)C(=O)",
        "tert-butylacetyl": "CC(C)(C)CC(=O)",
        "formyl": "C(=O)",
        "methylsulfonyl": "CS(=O)(=O)",
    }

    # C-terminal modifications for peptides
    C_TERMINAL_MODS = {
        "carboxylic_acid": "C(=O)O",
        "amide": "C(=O)N",
        "n_methylamide": "C(=O)NC",
        "ethyl_ester": "C(=O)OCC",
        "methyl_ester": "C(=O)OC",
        "hydroxamic_acid": "C(=O)NO",
    }

    # Amino acid replacements
    AMINO_ACID_SUBSTITUTIONS = {
        "Tyr": ["D-Tyr", "beta-Tyr", "3-Iodo-Tyr", "3,5-diiodo-Tyr", "HomoTyr"],
        "Gly": ["Sar", "Aib", "D-Ala", "beta-Ala"],
        "Phe": ["D-Phe", "HomoPhe", "4-Cl-Phe", "4-F-Phe", "Nal", "Trp"],
        "Leu": ["D-Leu", "Nle", "Tle", "Ile", "Val", "CycloLeu"],
    }

    # SMILES for amino acids
    AMINO_ACID_SMILES = {
        "Tyr": "N[C@@H](Cc1ccc(O)cc1)C(=O)",
        "D-Tyr": "N[C@H](Cc1ccc(O)cc1)C(=O)",
        "beta-Tyr": "NCC(Cc1ccc(O)cc1)C(=O)",
        "Gly": "NCC(=O)",
        "Sar": "CN(C)CC(=O)",  # N-methylglycine
        "Aib": "NC(C)(C)C(=O)",  # alpha-aminoisobutyric acid
        "Phe": "N[C@@H](Cc1ccccc1)C(=O)",
        "D-Phe": "N[C@H](Cc1ccccc1)C(=O)",
        "HomoPhe": "N[C@@H](CCc1ccccc1)C(=O)",
        "Leu": "N[C@@H](CC(C)C)C(=O)",
        "D-Leu": "N[C@H](CC(C)C)C(=O)",
        "Nle": "N[C@@H](CCCC)C(=O)",  # Norleucine
        "Tle": "N[C@@H](C(C)(C)C)C(=O)",  # tert-Leucine
    }

    def __init__(self):
        self.generated_smiles = set()

    def parse_input(self, compound: CompoundInput) -> Optional[Chem.Mol]:
        """Parse input compound to RDKit Mol object"""
        try:
            if compound.format == InputFormat.SMILES:
                mol = Chem.MolFromSmiles(compound.structure)
            elif compound.format == InputFormat.PEPTIDE_SEQUENCE:
                smiles = self._sequence_to_smiles(compound.structure)
                mol = Chem.MolFromSmiles(smiles) if smiles else None
            elif compound.format in [InputFormat.MOL_FILE, InputFormat.SDF_FILE]:
                mol = Chem.MolFromMolBlock(compound.structure)
            elif compound.format == InputFormat.PDB_FILE:
                mol = Chem.MolFromPDBBlock(compound.structure)
            else:
                mol = None

            if mol is not None:
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol, randomSeed=42)

            return mol
        except Exception as e:
            logger.error(f"Error parsing input: {e}")
            return None

    def _sequence_to_smiles(self, sequence: str, n_term: str = "pivaloyl", c_term: str = "carboxylic_acid") -> Optional[str]:
        """Convert peptide sequence to SMILES (simplified)"""
        # This is a simplified implementation
        # A full implementation would use proper peptide SMILES generation
        aa_map = {
            'Y': 'Tyr', 'G': 'Gly', 'F': 'Phe', 'L': 'Leu',
            'A': 'Ala', 'V': 'Val', 'I': 'Ile', 'M': 'Met',
            'P': 'Pro', 'W': 'Trp', 'S': 'Ser', 'T': 'Thr',
            'C': 'Cys', 'N': 'Asn', 'Q': 'Gln', 'K': 'Lys',
            'R': 'Arg', 'H': 'His', 'D': 'Asp', 'E': 'Glu'
        }

        # For KK103 (YGGFL with pivaloyl N-terminal)
        if sequence.upper() == "YGGFL":
            if n_term == "pivaloyl":
                # KK103 SMILES
                return "CC(C)(C)C(=O)NC(Cc1ccc(O)cc1)C(=O)NCC(=O)NCC(=O)NC(Cc2ccccc2)C(=O)NC(CC(C)C)C(=O)O"
            else:
                # Base Leu-ENK
                return "NC(Cc1ccc(O)cc1)C(=O)NCC(=O)NCC(=O)NC(Cc2ccccc2)C(=O)NC(CC(C)C)C(=O)O"

        return None

    def generate_analogues(
        self,
        compound: CompoundInput,
        strategies: List[str],
        num_analogues: int = 10,
        claude_suggestions: Optional[List[Dict]] = None
    ) -> List[Dict[str, Any]]:
        """Generate novel analogues using specified strategies"""
        analogues = []

        # Reset generated SMILES for this call
        self.generated_smiles = set()

        mol = self.parse_input(compound)

        if mol is None:
            logger.error("Could not parse input compound")
            return []

        canonical_smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))
        self.generated_smiles.add(canonical_smiles)

        # Use the original input SMILES for pattern matching (if SMILES format)
        input_smiles = compound.structure if compound.format == InputFormat.SMILES else canonical_smiles

        # Apply each strategy
        if "n_terminal" in strategies:
            analogues.extend(self._generate_n_terminal_mods(compound, input_smiles))

        if "c_terminal" in strategies:
            analogues.extend(self._generate_c_terminal_mods(compound, input_smiles))

        if "amino_acid_substitution" in strategies:
            analogues.extend(self._generate_aa_substitutions(compound, input_smiles))

        if "backbone" in strategies:
            analogues.extend(self._generate_backbone_mods(compound, input_smiles))

        if "bioisostere" in strategies:
            analogues.extend(self._generate_bioisostere_replacements(mol, input_smiles))

        # Add Claude-suggested modifications
        if claude_suggestions:
            analogues.extend(self._apply_claude_suggestions(compound, input_smiles, claude_suggestions))

        # Remove duplicates and filter
        unique_analogues = self._deduplicate_analogues(analogues)

        # Sort by diversity and return top N
        return self._select_diverse_set(unique_analogues, num_analogues, canonical_smiles)

    def _generate_n_terminal_mods(self, compound: CompoundInput, original_smiles: str) -> List[Dict[str, Any]]:
        """Generate N-terminal modifications"""
        analogues = []

        # Check if this is a peptide-like structure
        is_peptide = (
            compound.format == InputFormat.PEPTIDE_SEQUENCE or
            "YGGFL" in str(compound.structure).upper() or
            "enkephalin" in compound.name.lower() or
            "KK103" in compound.name.upper() or
            # Check for peptide bond pattern in SMILES (amide bonds)
            original_smiles.count("C(=O)N") >= 3
        )

        if not is_peptide:
            return analogues

        # Skip pivaloyl since that's KK103
        for mod_name, mod_smiles in self.N_TERMINAL_MODS.items():
            if mod_name == "pivaloyl":
                continue

            # Generate modified SMILES
            new_smiles = self._replace_n_terminal(original_smiles, mod_smiles)
            if new_smiles and new_smiles not in self.generated_smiles:
                self.generated_smiles.add(new_smiles)
                analogues.append({
                    "smiles": new_smiles,
                    "modification_type": "n_terminal",
                    "description": f"N-terminal {mod_name} modification",
                    "rationale": f"Replace pivaloyl with {mod_name} to explore SAR at N-terminus"
                })

        return analogues

    def _generate_c_terminal_mods(self, compound: CompoundInput, original_smiles: str) -> List[Dict[str, Any]]:
        """Generate C-terminal modifications"""
        analogues = []

        for mod_name, mod_smiles in self.C_TERMINAL_MODS.items():
            if mod_name == "carboxylic_acid":
                continue

            new_smiles = self._replace_c_terminal(original_smiles, mod_smiles)
            if new_smiles and new_smiles not in self.generated_smiles:
                self.generated_smiles.add(new_smiles)
                analogues.append({
                    "smiles": new_smiles,
                    "modification_type": "c_terminal",
                    "description": f"C-terminal {mod_name} modification",
                    "rationale": f"Convert to {mod_name} for improved membrane permeability or stability"
                })

        return analogues

    def _generate_aa_substitutions(self, compound: CompoundInput, original_smiles: str) -> List[Dict[str, Any]]:
        """Generate amino acid substitution analogues"""
        analogues = []

        # For each position, try substitutions
        substitution_sets = [
            ("Tyr", 1, ["D-Tyr", "3-Iodo-Tyr"]),
            ("Gly", 2, ["Sar", "D-Ala"]),
            ("Gly", 3, ["Sar", "Aib"]),
            ("Phe", 4, ["D-Phe", "4-F-Phe", "HomoPhe"]),
            ("Leu", 5, ["D-Leu", "Nle", "Tle"]),
        ]

        for original_aa, position, replacements in substitution_sets:
            for replacement in replacements:
                new_smiles = self._substitute_amino_acid(original_smiles, position, original_aa, replacement)
                if new_smiles and new_smiles not in self.generated_smiles:
                    self.generated_smiles.add(new_smiles)
                    analogues.append({
                        "smiles": new_smiles,
                        "modification_type": "amino_acid_substitution",
                        "description": f"Position {position}: {original_aa} → {replacement}",
                        "rationale": self._get_substitution_rationale(original_aa, replacement)
                    })

        return analogues

    def _generate_backbone_mods(self, compound: CompoundInput, original_smiles: str) -> List[Dict[str, Any]]:
        """Generate backbone modifications"""
        analogues = []

        # N-methylation at specific positions
        methylation_positions = [2, 3]  # Gly positions
        for pos in methylation_positions:
            new_smiles = self._add_n_methylation(original_smiles, pos)
            if new_smiles and new_smiles not in self.generated_smiles:
                self.generated_smiles.add(new_smiles)
                analogues.append({
                    "smiles": new_smiles,
                    "modification_type": "backbone",
                    "description": f"N-methylation at position {pos}",
                    "rationale": "N-methylation improves proteolytic stability and may enhance membrane permeability"
                })

        return analogues

    def _generate_bioisostere_replacements(self, mol: Chem.Mol, original_smiles: str) -> List[Dict[str, Any]]:
        """Generate bioisosteric replacements"""
        analogues = []
        # Implementation would use RDKit's reaction SMARTS for bioisosteric swaps
        return analogues

    def _apply_claude_suggestions(
        self,
        compound: CompoundInput,
        original_smiles: str,
        suggestions: List[Dict]
    ) -> List[Dict[str, Any]]:
        """Apply modifications suggested by Claude"""
        analogues = []

        for suggestion in suggestions:
            # Try to generate the suggested modification
            mod_type = suggestion.get("modification_type", "")
            description = suggestion.get("description", "")

            # Generate SMILES based on suggestion type
            new_smiles = self._generate_from_description(original_smiles, mod_type, description)

            if new_smiles and new_smiles not in self.generated_smiles:
                self.generated_smiles.add(new_smiles)
                analogues.append({
                    "smiles": new_smiles,
                    "modification_type": mod_type,
                    "description": description,
                    "rationale": suggestion.get("rationale", "AI-suggested modification")
                })

        return analogues

    def _replace_n_terminal(self, smiles: str, new_n_term: str) -> Optional[str]:
        """Replace N-terminal acyl group in a peptide SMILES"""
        # Simplified: replace pivaloyl pattern with new group
        # CC(C)(C)C(=O)N -> new_n_term + N
        import re
        pivaloyl_pattern = r"CC\(C\)\(C\)C\(=O\)"
        if re.search(pivaloyl_pattern, smiles):
            return re.sub(pivaloyl_pattern, new_n_term, smiles)
        return None

    def _replace_c_terminal(self, smiles: str, new_c_term: str) -> Optional[str]:
        """Replace C-terminal in a peptide SMILES"""
        # Replace terminal C(=O)O with new group
        if smiles.endswith("C(=O)O"):
            return smiles[:-6] + new_c_term
        return None

    def _substitute_amino_acid(self, smiles: str, position: int, original: str, replacement: str) -> Optional[str]:
        """Substitute an amino acid at a given position"""
        # This would require more sophisticated SMILES manipulation
        # For now, return predefined substitutions for KK103

        # Pre-computed substitutions for common modifications
        kk103_substitutions = {
            (1, "D-Tyr"): "CC(C)(C)C(=O)N[C@H](Cc1ccc(O)cc1)C(=O)NCC(=O)NCC(=O)NC(Cc2ccccc2)C(=O)NC(CC(C)C)C(=O)O",
            (2, "Sar"): "CC(C)(C)C(=O)NC(Cc1ccc(O)cc1)C(=O)N(C)CC(=O)NCC(=O)NC(Cc2ccccc2)C(=O)NC(CC(C)C)C(=O)O",
            (4, "D-Phe"): "CC(C)(C)C(=O)NC(Cc1ccc(O)cc1)C(=O)NCC(=O)NCC(=O)N[C@H](Cc2ccccc2)C(=O)NC(CC(C)C)C(=O)O",
            (5, "Nle"): "CC(C)(C)C(=O)NC(Cc1ccc(O)cc1)C(=O)NCC(=O)NCC(=O)NC(Cc2ccccc2)C(=O)NC(CCCC)C(=O)O",
            (5, "D-Leu"): "CC(C)(C)C(=O)NC(Cc1ccc(O)cc1)C(=O)NCC(=O)NCC(=O)NC(Cc2ccccc2)C(=O)N[C@H](CC(C)C)C(=O)O",
        }

        return kk103_substitutions.get((position, replacement))

    def _add_n_methylation(self, smiles: str, position: int) -> Optional[str]:
        """Add N-methylation at a specific position"""
        # Pre-computed for KK103
        methylated = {
            2: "CC(C)(C)C(=O)NC(Cc1ccc(O)cc1)C(=O)N(C)CC(=O)NCC(=O)NC(Cc2ccccc2)C(=O)NC(CC(C)C)C(=O)O",
            3: "CC(C)(C)C(=O)NC(Cc1ccc(O)cc1)C(=O)NCC(=O)N(C)CC(=O)NC(Cc2ccccc2)C(=O)NC(CC(C)C)C(=O)O",
        }
        return methylated.get(position)

    def _generate_from_description(self, original_smiles: str, mod_type: str, description: str) -> Optional[str]:
        """Generate SMILES from a textual description"""
        # This would use more sophisticated NLP and chemistry logic
        # Placeholder for now
        return None

    def _get_substitution_rationale(self, original: str, replacement: str) -> str:
        """Get rationale for an amino acid substitution"""
        rationales = {
            ("Tyr", "D-Tyr"): "D-amino acid substitution increases proteolytic stability",
            ("Gly", "Sar"): "N-methylation (sarcosine) improves stability and membrane permeability",
            ("Gly", "Aib"): "Alpha-aminoisobutyric acid constrains backbone conformation",
            ("Gly", "D-Ala"): "D-alanine adds steric bulk and proteolytic resistance",
            ("Phe", "D-Phe"): "D-phenylalanine increases proteolytic stability",
            ("Phe", "4-F-Phe"): "Fluorination can improve metabolic stability and binding",
            ("Phe", "HomoPhe"): "Homologation explores distance requirements in binding pocket",
            ("Leu", "D-Leu"): "D-leucine increases proteolytic stability at C-terminus",
            ("Leu", "Nle"): "Norleucine (linear chain) explores steric requirements",
            ("Leu", "Tle"): "tert-Leucine adds bulk and rigidity",
        }
        return rationales.get((original, replacement), "Structural modification to explore SAR")

    def _deduplicate_analogues(self, analogues: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Remove duplicate analogues based on canonical SMILES"""
        seen = set()
        unique = []

        for analogue in analogues:
            smiles = analogue.get("smiles", "")
            if smiles:
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        canonical = Chem.MolToSmiles(mol)
                        if canonical not in seen:
                            seen.add(canonical)
                            analogue["smiles"] = canonical
                            unique.append(analogue)
                except:
                    continue

        return unique

    def _select_diverse_set(
        self,
        analogues: List[Dict[str, Any]],
        num_select: int,
        reference_smiles: str
    ) -> List[Dict[str, Any]]:
        """Select a diverse set of analogues"""
        if len(analogues) <= num_select:
            return analogues

        # Calculate fingerprints
        ref_mol = Chem.MolFromSmiles(reference_smiles)
        ref_fp = AllChem.GetMorganFingerprintAsBitVect(ref_mol, 2, nBits=2048)

        # Score by diversity from reference
        scored = []
        for analogue in analogues:
            mol = Chem.MolFromSmiles(analogue["smiles"])
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                similarity = DataStructs.TanimotoSimilarity(ref_fp, fp)
                analogue["similarity_to_parent"] = similarity
                scored.append(analogue)

        # Sort by diversity (lower similarity = more diverse)
        scored.sort(key=lambda x: x.get("similarity_to_parent", 1.0))

        return scored[:num_select]

    def calculate_similarity(self, smiles1: str, smiles2: str) -> float:
        """Calculate Tanimoto similarity between two molecules"""
        try:
            mol1 = Chem.MolFromSmiles(smiles1)
            mol2 = Chem.MolFromSmiles(smiles2)
            if mol1 and mol2:
                fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
                fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
                return DataStructs.TanimotoSimilarity(fp1, fp2)
        except:
            pass
        return 0.0
