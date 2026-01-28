"""
Mechanistic Analyzer for Covalent Inhibitor Design

This module solves the critical problem of inverted topology in covalent inhibitors.
It analyzes which portion of a molecule STAYS bound vs which portion LEAVES during
the covalent reaction, and validates that binding groups are correctly placed.

The key insight: For covalent inhibitors, what matters is what REMAINS bound after
the enzyme attacks. The lipophilic/binding groups must be on the STAYING portion,
not the LEAVING portion.

Example failure this prevents:
  Wrong:   hexylphenyl-NH-C(=O)-N-piperidine
           ^^^^^^^^^^^        ^^^^^^^^^^^
           This LEAVES!       This STAYS (empty pocket)

  Correct: pyridazinyl-N-C(=O)-N-piperidine-biaryl-CF3
           ^^^^^^^^^^          ^^^^^^^^^^^^^^^^^^^^^
           This LEAVES         This STAYS (fills pocket)
"""

import logging
import re
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field
from enum import Enum

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

logger = logging.getLogger(__name__)


class WarheadType(Enum):
    """Types of covalent warheads and their mechanisms."""
    CARBAMATE = "carbamate"           # O-C(=O)-N - serine attacks carbonyl
    UREA = "urea"                     # N-C(=O)-N - serine attacks carbonyl
    ACRYLAMIDE = "acrylamide"         # C=C-C(=O)-N - Michael addition
    NITRILE = "nitrile"               # C#N - reversible covalent
    ALDEHYDE = "aldehyde"             # C=O - reversible covalent
    BORONIC_ACID = "boronic_acid"     # B(OH)2 - reversible covalent
    EPOXIDE = "epoxide"               # C1OC1 - irreversible
    VINYL_SULFONE = "vinyl_sulfone"   # C=C-S(=O)2 - Michael addition
    BETA_LACTAM = "beta_lactam"       # 4-membered cyclic amide
    REVERSIBLE = "reversible"         # Non-covalent
    UNKNOWN = "unknown"


@dataclass
class CovalentMechanism:
    """Describes the covalent mechanism of an inhibitor."""
    warhead_type: WarheadType = WarheadType.UNKNOWN
    warhead_smarts: str = ""

    # The key distinction
    leaving_group_smiles: str = ""      # What departs after enzyme attack
    staying_portion_smiles: str = ""    # What remains bound (the ADDUCT)

    # Analysis
    leaving_group_pka: float = 0.0      # pKa determines leaving group ability
    staying_portion_mw: float = 0.0
    staying_portion_logp: float = 0.0

    # Validation flags
    is_covalent: bool = False
    topology_correct: bool = False      # True if binding groups on staying portion
    mechanism_valid: bool = False

    # Details
    analysis_notes: List[str] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)


# SMARTS patterns for identifying warhead types
WARHEAD_PATTERNS = {
    # Carbamate: O-C(=O)-N
    WarheadType.CARBAMATE: [
        "[O;!R]-C(=O)-[N;!R]",           # Acyclic carbamate
        "[O;!R]-C(=O)-[N;R]",            # O-aryl N-cyclic carbamate
    ],
    # Urea: N-C(=O)-N
    WarheadType.UREA: [
        "[N;!R]-C(=O)-[N;!R]",           # Acyclic urea
        "[N;!R]-C(=O)-[N;R]",            # Asymmetric urea
        "[NX3;!R]-C(=O)-[NX3;!R]",       # Explicit NH urea
    ],
    # Acrylamide: C=C-C(=O)-N (Michael acceptor)
    WarheadType.ACRYLAMIDE: [
        "C=CC(=O)N",                      # Simple acrylamide
        "C=CC(=O)[NH]",                   # Secondary amide
    ],
    # Nitrile: C#N
    WarheadType.NITRILE: [
        "[C]#N",
    ],
    # Aldehyde: C(=O)H
    WarheadType.ALDEHYDE: [
        "[CX3H1](=O)",
    ],
    # Boronic acid: B(O)(O)
    WarheadType.BORONIC_ACID: [
        "[B]([OH])([OH])",
        "[B](O)(O)",
    ],
    # Vinyl sulfone
    WarheadType.VINYL_SULFONE: [
        "C=C[SX4](=O)(=O)",
    ],
}

# pKa values for common leaving groups (lower = better leaving group)
# These are approximate conjugate acid pKa values
LEAVING_GROUP_PKA = {
    # Aromatic amines (good leaving groups for carbamates/ureas)
    "aniline": 4.6,
    "4-nitroaniline": 1.0,
    "4-methoxyaniline": 5.3,
    "4-fluoroaniline": 4.7,
    "pyridazine": 2.3,        # Excellent leaving group
    "pyrimidine": 1.3,        # Excellent leaving group

    # Aliphatic amines (poor leaving groups)
    "piperidine": 11.1,       # Terrible leaving group
    "morpholine": 8.4,
    "pyrrolidine": 11.3,      # Terrible leaving group
    "diethylamine": 10.9,     # Terrible leaving group

    # Alcohols/phenols (for carbamates)
    "phenol": 10.0,
    "4-nitrophenol": 7.1,     # Good leaving group
    "hexanol": 16.0,          # Poor leaving group
}


class MechanisticAnalyzer:
    """
    Analyzes covalent inhibitor mechanisms to ensure correct topology.

    The key function: validate_covalent_topology() ensures that:
    1. Binding groups (lipophilic chains, π-stacking aromatics) are on STAYING portion
    2. Leaving group has appropriate pKa (low enough to leave)
    3. The covalent adduct (what remains) has drug-like properties
    """

    def __init__(self):
        self._compiled_patterns = {}
        self._compile_patterns()

    def _compile_patterns(self):
        """Pre-compile SMARTS patterns for efficiency."""
        for warhead_type, patterns in WARHEAD_PATTERNS.items():
            compiled = []
            for smarts in patterns:
                try:
                    pattern = Chem.MolFromSmarts(smarts)
                    if pattern:
                        compiled.append((smarts, pattern))
                except:
                    logger.warning(f"Could not compile SMARTS: {smarts}")
            self._compiled_patterns[warhead_type] = compiled

    def analyze_mechanism(self, smiles: str, target_type: str = "") -> CovalentMechanism:
        """
        Full mechanistic analysis of a molecule.

        Args:
            smiles: SMILES string of the candidate compound
            target_type: Type of target (e.g., "serine hydrolase", "cysteine protease")

        Returns:
            CovalentMechanism object with full analysis
        """
        result = CovalentMechanism()

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            result.errors.append(f"Invalid SMILES: {smiles}")
            return result

        # Step 1: Identify warhead type
        warhead_type, warhead_match = self._identify_warhead(mol)
        result.warhead_type = warhead_type

        if warhead_type == WarheadType.UNKNOWN:
            result.is_covalent = False
            result.analysis_notes.append("No covalent warhead detected - treating as reversible inhibitor")
            return result

        if warhead_type == WarheadType.REVERSIBLE:
            result.is_covalent = False
            return result

        result.is_covalent = True
        result.analysis_notes.append(f"Identified warhead: {warhead_type.value}")

        # Step 2: Split molecule into leaving/staying portions
        leaving, staying = self._split_at_warhead(mol, warhead_type, warhead_match)

        if leaving:
            result.leaving_group_smiles = leaving
            result.analysis_notes.append(f"Leaving group: {leaving}")
        else:
            result.errors.append("Could not identify leaving group")

        if staying:
            result.staying_portion_smiles = staying
            result.analysis_notes.append(f"Staying portion (adduct): {staying}")

            # Calculate properties of staying portion
            staying_mol = Chem.MolFromSmiles(staying)
            if staying_mol:
                result.staying_portion_mw = Descriptors.MolWt(staying_mol)
                result.staying_portion_logp = Descriptors.MolLogP(staying_mol)
        else:
            result.errors.append("Could not identify staying portion")

        # Step 3: Estimate leaving group pKa
        result.leaving_group_pka = self._estimate_leaving_group_pka(leaving)
        result.analysis_notes.append(f"Estimated leaving group pKa: {result.leaving_group_pka:.1f}")

        # Step 4: Validate topology
        result.topology_correct = self._validate_topology(
            mol, leaving, staying, warhead_type
        )

        # Step 5: Overall mechanism validity
        result.mechanism_valid = (
            result.is_covalent and
            result.topology_correct and
            result.leaving_group_pka < 12  # Must be reasonable leaving group
        )

        return result

    def _identify_warhead(self, mol: Chem.Mol) -> Tuple[WarheadType, Optional[Tuple]]:
        """Identify the type of covalent warhead in the molecule."""
        for warhead_type, patterns in self._compiled_patterns.items():
            for smarts_str, pattern in patterns:
                matches = mol.GetSubstructMatches(pattern)
                if matches:
                    return warhead_type, matches[0]

        return WarheadType.UNKNOWN, None

    def _split_at_warhead(
        self,
        mol: Chem.Mol,
        warhead_type: WarheadType,
        warhead_match: Tuple
    ) -> Tuple[str, str]:
        """
        Split molecule into leaving group and staying portion at the warhead.

        For carbamate O-C(=O)-N:
          - Enzyme attacks carbonyl carbon
          - O-C(=O) stays bound to enzyme
          - N-[rest] leaves

        For urea N-C(=O)-N:
          - Enzyme attacks carbonyl carbon
          - One N-[rest] stays, other N-[rest] leaves
          - The BETTER leaving group (lower pKa N) leaves
        """
        if not warhead_match:
            return "", ""

        try:
            if warhead_type == WarheadType.CARBAMATE:
                return self._split_carbamate(mol, warhead_match)
            elif warhead_type == WarheadType.UREA:
                return self._split_urea(mol, warhead_match)
            elif warhead_type == WarheadType.ACRYLAMIDE:
                return self._split_acrylamide(mol, warhead_match)
            else:
                return "", ""
        except Exception as e:
            logger.warning(f"Error splitting molecule: {e}")
            return "", ""

    def _split_carbamate(self, mol: Chem.Mol, match: Tuple) -> Tuple[str, str]:
        """
        Split carbamate O-C(=O)-N

        The O side stays bound to enzyme (as carbamate adduct)
        The N side leaves
        """
        # match is (O_idx, C_idx, N_idx) from pattern [O]-C(=O)-[N]
        # We need to find what's attached to O (stays) vs N (leaves)

        smiles = Chem.MolToSmiles(mol)

        # Find carbamate pattern
        carbamate_pattern = Chem.MolFromSmarts("[OX2]-C(=O)-[NX3]")
        matches = mol.GetSubstructMatches(carbamate_pattern)

        if not matches:
            return "", ""

        # Use first match
        o_idx, c_idx, _, n_idx = matches[0][0], matches[0][1], None, matches[0][2]

        # The carbonyl oxygen is at index 2 (double-bonded to C)
        # We need to find all atoms on O side vs N side

        # Simple approach: fragment at C-N bond
        # The staying portion is O-C(=O)*
        # The leaving portion is *-N-[rest]

        # Get atoms connected to the nitrogen (excluding the carbonyl carbon)
        n_atom = mol.GetAtomWithIdx(n_idx)
        n_neighbors = [a.GetIdx() for a in n_atom.GetNeighbors() if a.GetIdx() != c_idx]

        # Get atoms connected to the oxygen (excluding carbonyl carbon)
        o_atom = mol.GetAtomWithIdx(o_idx)
        o_neighbors = [a.GetIdx() for a in o_atom.GetNeighbors() if a.GetIdx() != c_idx]

        # Build staying portion (O-attached + carbamate)
        staying_atoms = self._get_connected_atoms(mol, o_neighbors, exclude={n_idx} | set(n_neighbors))
        staying_atoms.add(o_idx)
        staying_atoms.add(c_idx)
        # Add carbonyl oxygen
        for neighbor in mol.GetAtomWithIdx(c_idx).GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetIdx() != o_idx:
                staying_atoms.add(neighbor.GetIdx())

        # Build leaving portion (N-attached)
        leaving_atoms = self._get_connected_atoms(mol, [n_idx] + n_neighbors, exclude=staying_atoms)

        # Convert to SMILES
        staying_smiles = self._atoms_to_smiles(mol, staying_atoms)
        leaving_smiles = self._atoms_to_smiles(mol, leaving_atoms)

        return leaving_smiles, staying_smiles

    def _split_urea(self, mol: Chem.Mol, match: Tuple) -> Tuple[str, str]:
        """
        Split urea N-C(=O)-N

        One nitrogen stays (as carbamoyl adduct), the other leaves.
        The nitrogen with LOWER pKa (better leaving group) leaves.

        Key insight: For FAAH inhibitors:
          - Aromatic amines (pKa ~4-5) are LEAVING groups
          - Aliphatic amines like piperidine (pKa ~11) STAY
          - Lipophilic groups should be attached to the STAYING nitrogen
        """
        # Find urea pattern
        urea_pattern = Chem.MolFromSmarts("[NX3]-C(=O)-[NX3]")
        matches = mol.GetSubstructMatches(urea_pattern)

        if not matches:
            return "", ""

        n1_idx, c_idx, n2_idx = matches[0][0], matches[0][1], matches[0][2]

        # Determine which nitrogen is the better leaving group
        n1_substituents = self._get_nitrogen_substituents(mol, n1_idx, c_idx)
        n2_substituents = self._get_nitrogen_substituents(mol, n2_idx, c_idx)

        n1_pka = self._estimate_nitrogen_pka(mol, n1_idx, n1_substituents)
        n2_pka = self._estimate_nitrogen_pka(mol, n2_idx, n2_substituents)

        # Lower pKa = better leaving group = LEAVES
        if n1_pka < n2_pka:
            leaving_n_idx = n1_idx
            staying_n_idx = n2_idx
        else:
            leaving_n_idx = n2_idx
            staying_n_idx = n1_idx

        # Get atoms for each portion
        leaving_atoms = self._get_connected_atoms(
            mol, [leaving_n_idx],
            exclude={c_idx, staying_n_idx}
        )

        staying_atoms = self._get_connected_atoms(
            mol, [staying_n_idx],
            exclude={c_idx, leaving_n_idx}
        )
        staying_atoms.add(c_idx)
        staying_atoms.add(staying_n_idx)
        # Add carbonyl oxygen
        for neighbor in mol.GetAtomWithIdx(c_idx).GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                staying_atoms.add(neighbor.GetIdx())

        leaving_smiles = self._atoms_to_smiles(mol, leaving_atoms)
        staying_smiles = self._atoms_to_smiles(mol, staying_atoms)

        return leaving_smiles, staying_smiles

    def _split_acrylamide(self, mol: Chem.Mol, match: Tuple) -> Tuple[str, str]:
        """
        Split acrylamide C=C-C(=O)-N

        Cysteine attacks the beta carbon (Michael addition)
        The whole amide portion stays, only H leaves
        """
        # For acrylamide, there's no "leaving group" per se
        # The whole molecule stays bound
        smiles = Chem.MolToSmiles(mol)
        return "", smiles

    def _get_connected_atoms(
        self,
        mol: Chem.Mol,
        start_atoms: List[int],
        exclude: set = None
    ) -> set:
        """BFS to get all atoms connected to start_atoms, excluding some atoms."""
        if exclude is None:
            exclude = set()

        visited = set()
        queue = list(start_atoms)

        while queue:
            atom_idx = queue.pop(0)
            if atom_idx in visited or atom_idx in exclude:
                continue
            visited.add(atom_idx)

            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                n_idx = neighbor.GetIdx()
                if n_idx not in visited and n_idx not in exclude:
                    queue.append(n_idx)

        return visited

    def _atoms_to_smiles(self, mol: Chem.Mol, atom_indices: set) -> str:
        """Convert a set of atom indices to a SMILES string."""
        if not atom_indices:
            return ""

        try:
            # Create editable molecule
            edit_mol = Chem.RWMol(mol)

            # Get atoms to remove (not in our set)
            all_atoms = set(range(mol.GetNumAtoms()))
            atoms_to_remove = all_atoms - atom_indices

            # Remove atoms in reverse order (so indices don't shift)
            for idx in sorted(atoms_to_remove, reverse=True):
                edit_mol.RemoveAtom(idx)

            # Get SMILES
            frag_mol = edit_mol.GetMol()

            # Handle fragments (take largest)
            frags = Chem.GetMolFrags(frag_mol, asMols=True)
            if frags:
                largest = max(frags, key=lambda x: x.GetNumAtoms())
                return Chem.MolToSmiles(largest)

            return Chem.MolToSmiles(frag_mol)
        except Exception as e:
            logger.warning(f"Error converting atoms to SMILES: {e}")
            return ""

    def _get_nitrogen_substituents(
        self,
        mol: Chem.Mol,
        n_idx: int,
        exclude_idx: int
    ) -> List[Dict]:
        """Get information about substituents on a nitrogen."""
        n_atom = mol.GetAtomWithIdx(n_idx)
        substituents = []

        for neighbor in n_atom.GetNeighbors():
            if neighbor.GetIdx() == exclude_idx:
                continue

            sub_info = {
                "idx": neighbor.GetIdx(),
                "symbol": neighbor.GetSymbol(),
                "is_aromatic": neighbor.GetIsAromatic(),
                "is_in_ring": neighbor.IsInRing(),
                "ring_size": 0,
            }

            if neighbor.IsInRing():
                ring_info = mol.GetRingInfo()
                for ring in ring_info.AtomRings():
                    if neighbor.GetIdx() in ring:
                        sub_info["ring_size"] = len(ring)
                        break

            substituents.append(sub_info)

        return substituents

    def _estimate_nitrogen_pka(
        self,
        mol: Chem.Mol,
        n_idx: int,
        substituents: List[Dict]
    ) -> float:
        """
        Estimate the pKa of a nitrogen (as conjugate acid).

        Lower pKa = more stable as neutral = better leaving group
        - Aromatic amines: pKa ~4-5 (good leaving groups)
        - Aliphatic amines: pKa ~10-11 (poor leaving groups)
        """
        n_atom = mol.GetAtomWithIdx(n_idx)

        # Check if nitrogen is part of aromatic ring
        if n_atom.GetIsAromatic():
            return 3.0  # Like pyridine nitrogen - excellent leaving group

        # Check substituents
        has_aromatic = any(s["is_aromatic"] for s in substituents)
        all_aromatic = all(s["is_aromatic"] for s in substituents) if substituents else False

        has_ring = any(s["is_in_ring"] for s in substituents)
        ring_sizes = [s["ring_size"] for s in substituents if s["ring_size"] > 0]

        if all_aromatic:
            # Diarylamine or arylamine - good leaving group
            return 4.5

        if has_aromatic:
            # Arylamine - moderate leaving group
            return 5.5

        # Check for saturated ring (piperidine, morpholine, etc.)
        if has_ring and ring_sizes:
            ring_size = ring_sizes[0]
            if ring_size == 6:
                return 11.0  # Piperidine-like - poor leaving group
            elif ring_size == 5:
                return 11.3  # Pyrrolidine-like - poor leaving group

        # Generic aliphatic amine
        return 10.5

    def _estimate_leaving_group_pka(self, leaving_smiles: str) -> float:
        """Estimate pKa of the leaving group."""
        if not leaving_smiles:
            return 15.0  # No leaving group = not a good inhibitor

        mol = Chem.MolFromSmiles(leaving_smiles)
        if mol is None:
            return 15.0

        # Check for known patterns
        smiles_lower = leaving_smiles.lower()

        # Aromatic nitrogen heterocycles (excellent leaving groups)
        pyridazine = Chem.MolFromSmarts("c1cnncc1")
        pyrimidine = Chem.MolFromSmarts("c1ncncn1")
        pyridine = Chem.MolFromSmarts("c1ccncc1")

        if pyridazine and mol.HasSubstructMatch(pyridazine):
            return 2.3
        if pyrimidine and mol.HasSubstructMatch(pyrimidine):
            return 1.3
        if pyridine and mol.HasSubstructMatch(pyridine):
            return 5.2

        # Anilines
        aniline = Chem.MolFromSmarts("Nc1ccccc1")
        if aniline and mol.HasSubstructMatch(aniline):
            # Check for electron-withdrawing groups
            nitro = Chem.MolFromSmarts("[N+](=O)[O-]")
            cf3 = Chem.MolFromSmarts("C(F)(F)F")
            fluoro = Chem.MolFromSmarts("F")

            if nitro and mol.HasSubstructMatch(nitro):
                return 1.0  # Nitroaniline - excellent
            if cf3 and mol.HasSubstructMatch(cf3):
                return 3.5  # CF3-aniline - good
            if fluoro and mol.HasSubstructMatch(fluoro):
                return 4.7  # Fluoroaniline - good
            return 4.6  # Plain aniline - acceptable

        # Saturated amines (poor leaving groups)
        piperidine = Chem.MolFromSmarts("C1CCNCC1")
        morpholine = Chem.MolFromSmarts("C1COCCN1")
        pyrrolidine = Chem.MolFromSmarts("C1CCNC1")

        if piperidine and mol.HasSubstructMatch(piperidine):
            return 11.1
        if morpholine and mol.HasSubstructMatch(morpholine):
            return 8.4
        if pyrrolidine and mol.HasSubstructMatch(pyrrolidine):
            return 11.3

        # Default based on molecular features
        # Count aromatic atoms
        aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())

        if aromatic_atoms > 4:
            return 5.0  # Aromatic amine
        else:
            return 10.0  # Aliphatic amine

    def _validate_topology(
        self,
        mol: Chem.Mol,
        leaving_smiles: str,
        staying_smiles: str,
        warhead_type: WarheadType
    ) -> bool:
        """
        Validate that binding groups are on the STAYING portion.

        For covalent inhibitors, the pharmacophoric features (lipophilic chains,
        aromatic groups for π-stacking) MUST be on the portion that remains
        bound to the enzyme after covalent modification.

        Returns True if topology is correct, False if inverted.
        """
        if not staying_smiles or not leaving_smiles:
            return True  # Can't validate without both portions

        staying_mol = Chem.MolFromSmiles(staying_smiles)
        leaving_mol = Chem.MolFromSmiles(leaving_smiles)

        if staying_mol is None or leaving_mol is None:
            return True  # Can't validate

        # Calculate properties
        staying_mw = Descriptors.MolWt(staying_mol)
        leaving_mw = Descriptors.MolWt(leaving_mol)

        staying_logp = Descriptors.MolLogP(staying_mol)
        leaving_logp = Descriptors.MolLogP(leaving_mol)

        staying_carbons = sum(1 for atom in staying_mol.GetAtoms() if atom.GetAtomicNum() == 6)
        leaving_carbons = sum(1 for atom in leaving_mol.GetAtoms() if atom.GetAtomicNum() == 6)

        staying_aromatics = sum(1 for atom in staying_mol.GetAtoms() if atom.GetIsAromatic())
        leaving_aromatics = sum(1 for atom in leaving_mol.GetAtoms() if atom.GetIsAromatic())

        # Check for long alkyl chains (lipophilic tails)
        long_chain_pattern = Chem.MolFromSmarts("[CH2][CH2][CH2][CH2][CH2]")  # 5+ CH2
        staying_has_long_chain = bool(staying_mol.HasSubstructMatch(long_chain_pattern)) if long_chain_pattern else False
        leaving_has_long_chain = bool(leaving_mol.HasSubstructMatch(long_chain_pattern)) if long_chain_pattern else False

        # CRITICAL CHECK: Lipophilic tail should be on STAYING portion
        # This is the exact error we're trying to catch
        if leaving_has_long_chain and not staying_has_long_chain:
            logger.warning(f"TOPOLOGY ERROR: Long alkyl chain on leaving group!")
            logger.warning(f"Leaving: {leaving_smiles} (LogP: {leaving_logp:.1f})")
            logger.warning(f"Staying: {staying_smiles} (LogP: {staying_logp:.1f})")
            return False

        # Check if leaving portion is too large/lipophilic
        if leaving_logp > staying_logp + 2.0:
            logger.warning(f"TOPOLOGY WARNING: Leaving group more lipophilic than staying portion")
            logger.warning(f"Leaving LogP: {leaving_logp:.1f}, Staying LogP: {staying_logp:.1f}")
            return False

        # Check if leaving portion has most of the mass
        if leaving_mw > staying_mw * 1.5:
            logger.warning(f"TOPOLOGY WARNING: Most of molecule leaves!")
            logger.warning(f"Leaving MW: {leaving_mw:.0f}, Staying MW: {staying_mw:.0f}")
            return False

        return True

    def validate_for_target(
        self,
        smiles: str,
        target_type: str,
        mechanism: str
    ) -> Dict[str, Any]:
        """
        Validate a candidate for a specific target.

        Args:
            smiles: Candidate SMILES
            target_type: Type of target (e.g., "serine hydrolase", "GPCR")
            mechanism: Mechanism of action (e.g., "covalent inhibitor", "reversible")

        Returns:
            Dictionary with validation results
        """
        result = {
            "smiles": smiles,
            "valid": True,
            "mechanism_analysis": None,
            "errors": [],
            "warnings": [],
            "recommendations": []
        }

        # Analyze mechanism
        mechanism_info = self.analyze_mechanism(smiles, target_type)
        result["mechanism_analysis"] = {
            "warhead_type": mechanism_info.warhead_type.value,
            "is_covalent": mechanism_info.is_covalent,
            "leaving_group": mechanism_info.leaving_group_smiles,
            "staying_portion": mechanism_info.staying_portion_smiles,
            "leaving_group_pka": mechanism_info.leaving_group_pka,
            "topology_correct": mechanism_info.topology_correct,
            "notes": mechanism_info.analysis_notes
        }

        # Check for covalent mechanism targets
        covalent_keywords = ["covalent", "irreversible", "serine hydrolase", "faah", "magl"]
        is_covalent_target = any(kw in mechanism.lower() or kw in target_type.lower() for kw in covalent_keywords)

        if is_covalent_target:
            if not mechanism_info.is_covalent:
                result["warnings"].append("Target requires covalent mechanism but no warhead detected")

            if mechanism_info.is_covalent and not mechanism_info.topology_correct:
                result["valid"] = False
                result["errors"].append(
                    f"INVERTED TOPOLOGY: Binding groups on leaving portion! "
                    f"Leaving group ({mechanism_info.leaving_group_smiles}) has lipophilic features "
                    f"that should be on the staying portion ({mechanism_info.staying_portion_smiles})"
                )
                result["recommendations"].append(
                    "Redesign: Move lipophilic tail to the nitrogen that STAYS bound"
                )

            if mechanism_info.leaving_group_pka > 10:
                result["warnings"].append(
                    f"Poor leaving group (pKa ~{mechanism_info.leaving_group_pka:.0f}). "
                    f"Consider aromatic amine or better leaving group (pKa < 8)"
                )

        result["errors"].extend(mechanism_info.errors)

        return result


# Singleton instance
mechanistic_analyzer = MechanisticAnalyzer()
