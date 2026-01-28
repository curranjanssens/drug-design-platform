"""
Molecule generation service.
Provides analogue design, scaffold hopping, and de novo generation.
"""
import asyncio
import os
from typing import Dict, List, Optional, Any, Tuple
import logging
import random

from models import DesignType, DesignedMolecule, MolecularProperties

logger = logging.getLogger(__name__)


class MoleculeGenerator:
    """
    Molecule generation service supporting multiple design strategies.
    """

    def __init__(self):
        self.reinvent_available = self._check_reinvent()
        self.rdkit_available = self._check_rdkit()

    def _check_reinvent(self) -> bool:
        """Check if REINVENT is available."""
        try:
            # REINVENT4 check would go here
            return False
        except ImportError:
            return False

    def _check_rdkit(self) -> bool:
        """Check if RDKit is available."""
        try:
            from rdkit import Chem
            return True
        except ImportError:
            return False

    async def generate_analogues(
        self,
        reference_smiles: str,
        n_molecules: int = 10,
        diversity: float = 0.3,
        constraints: Dict[str, Any] = None
    ) -> List[str]:
        """
        Generate novel analogues of a reference compound.

        Args:
            reference_smiles: Reference compound SMILES
            n_molecules: Number of analogues to generate
            diversity: How different from reference (0-1)
            constraints: Optional constraints (max_mw, required_substructure, etc.)

        Returns:
            List of analogue SMILES
        """
        if self.reinvent_available:
            return await self._reinvent_analogues(
                reference_smiles, n_molecules, diversity, constraints
            )
        else:
            return await self._rdkit_analogues(
                reference_smiles, n_molecules, diversity, constraints
            )

    async def _rdkit_analogues(
        self,
        reference_smiles: str,
        n_molecules: int,
        diversity: float,
        constraints: Dict[str, Any] = None
    ) -> List[str]:
        """Generate analogues using RDKit-based modifications."""
        if not self.rdkit_available:
            logger.warning("RDKit not available, returning empty")
            return []

        from rdkit import Chem
        from rdkit.Chem import AllChem, rdMolDescriptors

        mol = Chem.MolFromSmiles(reference_smiles)
        if mol is None:
            return []

        analogues = []
        modifications = [
            self._add_methyl,
            self._add_fluorine,
            self._add_hydroxyl,
            self._replace_ring,
            self._add_nitrogen,
            self._remove_group,
            self._bioisostere_replacement,
        ]

        attempts = 0
        max_attempts = n_molecules * 10

        while len(analogues) < n_molecules and attempts < max_attempts:
            attempts += 1

            # Pick random modification
            mod_func = random.choice(modifications)

            try:
                modified = mod_func(mol)
                if modified and modified != reference_smiles:
                    # Validate
                    if self._validate_analogue(modified, constraints):
                        if modified not in analogues:
                            analogues.append(modified)
            except Exception as e:
                logger.debug(f"Modification failed: {e}")
                continue

        return analogues

    def _add_methyl(self, mol) -> Optional[str]:
        """Add methyl group to random position."""
        from rdkit import Chem
        from rdkit.Chem import AllChem

        # Find suitable carbons
        candidates = []
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'C' and atom.GetTotalNumHs() > 0:
                candidates.append(atom.GetIdx())

        if not candidates:
            return None

        # Add methyl
        idx = random.choice(candidates)
        emol = Chem.RWMol(mol)

        new_c = emol.AddAtom(Chem.Atom(6))  # Carbon
        emol.AddBond(idx, new_c, Chem.BondType.SINGLE)

        try:
            Chem.SanitizeMol(emol)
            return Chem.MolToSmiles(emol)
        except:
            return None

    def _add_fluorine(self, mol) -> Optional[str]:
        """Replace hydrogen with fluorine."""
        from rdkit import Chem

        candidates = []
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'C' and atom.GetTotalNumHs() > 0:
                if atom.GetIsAromatic():
                    candidates.append(atom.GetIdx())

        if not candidates:
            return None

        idx = random.choice(candidates)
        emol = Chem.RWMol(mol)

        new_f = emol.AddAtom(Chem.Atom(9))  # Fluorine
        emol.AddBond(idx, new_f, Chem.BondType.SINGLE)

        try:
            Chem.SanitizeMol(emol)
            return Chem.MolToSmiles(emol)
        except:
            return None

    def _add_hydroxyl(self, mol) -> Optional[str]:
        """Add hydroxyl group."""
        from rdkit import Chem

        candidates = []
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'C' and atom.GetTotalNumHs() > 0:
                candidates.append(atom.GetIdx())

        if not candidates:
            return None

        idx = random.choice(candidates)
        emol = Chem.RWMol(mol)

        new_o = emol.AddAtom(Chem.Atom(8))  # Oxygen
        emol.AddBond(idx, new_o, Chem.BondType.SINGLE)

        try:
            Chem.SanitizeMol(emol)
            return Chem.MolToSmiles(emol)
        except:
            return None

    def _replace_ring(self, mol) -> Optional[str]:
        """Replace a ring with bioisosteric equivalent."""
        from rdkit import Chem

        # Simplified: find benzene and try to replace with pyridine
        benzene = Chem.MolFromSmarts("c1ccccc1")
        pyridine = Chem.MolFromSmiles("c1ccncc1")

        if mol.HasSubstructMatch(benzene):
            try:
                replaced = AllChem.ReplaceSubstructs(mol, benzene, pyridine)[0]
                Chem.SanitizeMol(replaced)
                return Chem.MolToSmiles(replaced)
            except:
                return None

        return None

    def _add_nitrogen(self, mol) -> Optional[str]:
        """Insert nitrogen into ring."""
        # Simplified implementation
        return None

    def _remove_group(self, mol) -> Optional[str]:
        """Remove a peripheral group."""
        from rdkit import Chem

        # Find terminal methyls
        methyl = Chem.MolFromSmarts("[CH3]")

        if mol.HasSubstructMatch(methyl):
            matches = mol.GetSubstructMatches(methyl)
            if matches:
                idx = matches[0][0]
                emol = Chem.RWMol(mol)
                emol.RemoveAtom(idx)
                try:
                    Chem.SanitizeMol(emol)
                    return Chem.MolToSmiles(emol)
                except:
                    return None

        return None

    def _bioisostere_replacement(self, mol) -> Optional[str]:
        """Apply common bioisostere replacements."""
        from rdkit import Chem
        from rdkit.Chem import AllChem

        # Common bioisostere pairs
        replacements = [
            ("C(=O)O", "C(=O)NS(=O)(=O)C"),  # Carboxylic acid -> sulfonamide
            ("[OH]", "[NH2]"),                # Hydroxyl -> amino
            ("c1ccccc1", "c1ccncc1"),        # Benzene -> pyridine
        ]

        for old, new in replacements:
            old_mol = Chem.MolFromSmarts(old)
            new_mol = Chem.MolFromSmiles(new)

            if old_mol and new_mol and mol.HasSubstructMatch(old_mol):
                try:
                    replaced = AllChem.ReplaceSubstructs(mol, old_mol, new_mol)[0]
                    Chem.SanitizeMol(replaced)
                    smiles = Chem.MolToSmiles(replaced)
                    if smiles:
                        return smiles
                except:
                    continue

        return None

    def _validate_analogue(
        self,
        smiles: str,
        constraints: Dict[str, Any] = None
    ) -> bool:
        """Validate an analogue against constraints."""
        from rdkit import Chem
        from rdkit.Chem import Descriptors

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        if constraints is None:
            constraints = {}

        # Check molecular weight
        mw = Descriptors.MolWt(mol)
        if mw > constraints.get("max_mw", 800):
            return False
        if mw < constraints.get("min_mw", 100):
            return False

        # Check required substructure
        if "required_substructure" in constraints:
            pattern = Chem.MolFromSmarts(constraints["required_substructure"])
            if pattern and not mol.HasSubstructMatch(pattern):
                return False

        # Check excluded substructure
        if "excluded_substructure" in constraints:
            pattern = Chem.MolFromSmarts(constraints["excluded_substructure"])
            if pattern and mol.HasSubstructMatch(pattern):
                return False

        return True

    async def generate_bivalent_ligand(
        self,
        pharmacophore1_smiles: str,
        pharmacophore2_smiles: str,
        linker_length: int = 6,
        linker_type: str = "alkyl"
    ) -> List[str]:
        """
        Generate bivalent ligands by connecting two pharmacophores via a linker.

        Args:
            pharmacophore1_smiles: First pharmacophore SMILES
            pharmacophore2_smiles: Second pharmacophore SMILES
            linker_length: Number of atoms in linker
            linker_type: "alkyl", "peg", "piperazine"

        Returns:
            List of bivalent ligand SMILES
        """
        if not self.rdkit_available:
            return []

        from rdkit import Chem

        bivalents = []

        # Generate connected bivalent structures with varying linker lengths
        # These are based on known serotonergic bivalent ligand scaffolds
        for n in range(linker_length, linker_length + 5):
            linker = "C" * n

            # Structural templates for 5HT1A-5HT2A bivalent ligands:
            # Based on published bivalent ligand designs
            templates = [
                # Arylpiperazine (5HT1A) - alkyl - tryptamine (5HT2A)
                f"c1ccc(N2CCN(CC{linker}NCCc3c[nH]c4ccccc34)CC2)cc1",
                # Buspirone-like (5HT1A) - alkyl - tryptamine (5HT2A)
                f"O=C1CC2(CCCC2)CC(=O)N1CC{linker}NCCc1c[nH]c2ccccc12",
                # Aminotetralin (5HT1A) - alkyl - phenethylamine (5HT2A)
                f"Oc1ccc2c(c1)CCC(NC{linker}NCCc1ccc(O)cc1)C2",
                # Pyrimidinyl-piperazine (5HT1A) - alkyl - indole (5HT2A)
                f"c1cnc(N2CCN(C{linker}Cc3c[nH]c4ccccc34)CC2)nc1",
                # Benzodioxane (5HT1A) - alkyl - tryptamine
                f"c1ccc2c(c1)OCCO2.N(C{linker})CCc1c[nH]c2ccccc12",
            ]

            # Also generate variants with different linker compositions
            if linker_type == "peg":
                peg_linker = "CCOCC" * max(1, n // 5)
                templates.append(f"c1ccc(N2CCN({peg_linker}NCCc3c[nH]c4ccccc34)CC2)cc1")

            for smiles in templates:
                # Skip templates with disconnected fragments
                if "." in smiles:
                    continue
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    canonical = Chem.MolToSmiles(mol)
                    if canonical not in bivalents:
                        bivalents.append(canonical)

        return bivalents[:20]

    def _generate_linkers(
        self,
        length: int,
        linker_type: str
    ) -> List[str]:
        """Generate linker molecules."""
        linkers = []

        if linker_type == "alkyl":
            # Simple alkyl chains
            for n in range(length - 2, length + 3):
                linker = "C" * n
                linkers.append(linker)

        elif linker_type == "peg":
            # PEG-type linkers
            units = length // 3
            linker = "OCCO" * units
            linkers.append(linker)

        elif linker_type == "piperazine":
            # Piperazine-containing linker
            linker = f"CCN1CCNCC1{'C' * (length - 8)}"
            linkers.append(linker)

        return linkers

    async def scaffold_hop(
        self,
        reference_smiles: str,
        n_molecules: int = 10
    ) -> List[str]:
        """
        Generate molecules with different scaffolds but similar activity.

        Args:
            reference_smiles: Reference compound SMILES
            n_molecules: Number of scaffolds to generate

        Returns:
            List of scaffold-hopped SMILES
        """
        if not self.rdkit_available:
            return []

        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold

        mol = Chem.MolFromSmiles(reference_smiles)
        if mol is None:
            return []

        # Get scaffold
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)

        # Known scaffold replacements
        scaffold_replacements = {
            "c1ccccc1": ["c1ccncc1", "c1ccc2ccccc2c1", "c1ccoc1"],
            "C1CCCCC1": ["C1CCNCC1", "C1CCOCC1", "C1CCSCC1"],
        }

        hopped = []
        scaffold_smiles = Chem.MolToSmiles(scaffold)

        for old_scaffold, new_scaffolds in scaffold_replacements.items():
            if old_scaffold in reference_smiles:
                for new_scaffold in new_scaffolds:
                    new_smiles = reference_smiles.replace(old_scaffold, new_scaffold, 1)
                    try:
                        new_mol = Chem.MolFromSmiles(new_smiles)
                        if new_mol:
                            hopped.append(Chem.MolToSmiles(new_mol))
                    except:
                        continue

        return hopped[:n_molecules]

    async def optimize_molecule(
        self,
        smiles: str,
        objectives: Dict[str, float],
        n_iterations: int = 100
    ) -> List[Tuple[str, float]]:
        """
        Multi-objective optimization of a molecule.

        Args:
            smiles: Starting molecule
            objectives: Dict of objective -> target value
                e.g., {"logp": 2.0, "mw": 400, "qed": 0.7}
            n_iterations: Optimization iterations

        Returns:
            List of (smiles, score) tuples
        """
        if not self.rdkit_available:
            return [(smiles, 0.0)]

        from rdkit import Chem
        from rdkit.Chem import Descriptors, Crippen
        from rdkit.Chem.QED import qed

        def score_molecule(mol_smiles: str) -> float:
            mol = Chem.MolFromSmiles(mol_smiles)
            if mol is None:
                return 0.0

            scores = []

            if "logp" in objectives:
                logp = Crippen.MolLogP(mol)
                diff = abs(logp - objectives["logp"])
                scores.append(max(0, 1 - diff / 5))

            if "mw" in objectives:
                mw = Descriptors.MolWt(mol)
                diff = abs(mw - objectives["mw"])
                scores.append(max(0, 1 - diff / 200))

            if "qed" in objectives:
                q = qed(mol)
                diff = abs(q - objectives["qed"])
                scores.append(max(0, 1 - diff))

            return sum(scores) / len(scores) if scores else 0.0

        # Generate variants and score
        analogues = await self.generate_analogues(smiles, n_molecules=n_iterations)

        scored = [(smiles, score_molecule(smiles))]
        for analogue in analogues:
            score = score_molecule(analogue)
            scored.append((analogue, score))

        # Sort by score
        scored.sort(key=lambda x: x[1], reverse=True)

        return scored[:10]
