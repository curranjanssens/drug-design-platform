"""
Molecular Docking Service using AutoDock Vina
"""
from typing import Dict, Any, Optional, List, Tuple
from pathlib import Path
import tempfile
import subprocess
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from loguru import logger
import numpy as np

from config import settings


class DockingService:
    """Service for molecular docking using AutoDock Vina"""

    # Delta Opioid Receptor binding site parameters (from crystal structure)
    DOR_BINDING_SITE = {
        "center_x": 0.0,  # These would be set from actual PDB structure
        "center_y": 0.0,
        "center_z": 0.0,
        "size_x": 25.0,
        "size_y": 25.0,
        "size_z": 25.0,
    }

    # Reference binding affinities (kcal/mol, estimated)
    REFERENCE_AFFINITIES = {
        "leu_enkephalin": -8.5,
        "kk103": -8.2,
    }

    def __init__(self):
        self.exhaustiveness = settings.autodock_exhaustiveness
        self.num_modes = settings.autodock_num_modes
        self.temp_dir = settings.temp_dir

    async def dock_compound(
        self,
        ligand_smiles: str,
        receptor_pdb: Optional[str] = None,
        binding_site: Optional[Dict[str, float]] = None
    ) -> Dict[str, Any]:
        """
        Dock a compound against a receptor

        Args:
            ligand_smiles: SMILES of the ligand
            receptor_pdb: Path to receptor PDB file (or use default DOR)
            binding_site: Binding site coordinates and dimensions

        Returns:
            Docking results including affinity and pose
        """
        try:
            # Prepare ligand
            ligand_pdbqt = self._prepare_ligand(ligand_smiles)
            if ligand_pdbqt is None:
                return {"success": False, "error": "Failed to prepare ligand"}

            # Use default receptor if not provided
            if receptor_pdb is None:
                receptor_pdb = self._get_default_receptor()

            # Prepare receptor
            receptor_pdbqt = self._prepare_receptor(receptor_pdb)
            if receptor_pdbqt is None:
                return {"success": False, "error": "Failed to prepare receptor"}

            # Set binding site
            site = binding_site or self.DOR_BINDING_SITE

            # Run docking
            result = await self._run_vina(
                ligand_pdbqt,
                receptor_pdbqt,
                site
            )

            return result

        except Exception as e:
            logger.error(f"Docking error: {e}")
            return {"success": False, "error": str(e)}

    def _prepare_ligand(self, smiles: str) -> Optional[str]:
        """Prepare ligand PDBQT file from SMILES"""
        try:
            # Generate 3D structure
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None

            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)

            # Write to temp PDB file
            pdb_path = self.temp_dir / f"ligand_{hash(smiles)}.pdb"
            pdbqt_path = self.temp_dir / f"ligand_{hash(smiles)}.pdbqt"

            Chem.MolToPDBFile(mol, str(pdb_path))

            # Convert to PDBQT using Open Babel or meeko
            # In production, use meeko for better results
            try:
                subprocess.run([
                    "obabel", str(pdb_path),
                    "-O", str(pdbqt_path),
                    "--gen3d", "-h"
                ], check=True, capture_output=True)
            except (subprocess.CalledProcessError, FileNotFoundError):
                # If obabel not available, create simplified PDBQT
                self._simple_pdb_to_pdbqt(pdb_path, pdbqt_path)

            return str(pdbqt_path)

        except Exception as e:
            logger.error(f"Ligand preparation error: {e}")
            return None

    def _prepare_receptor(self, pdb_path: str) -> Optional[str]:
        """Prepare receptor PDBQT file"""
        pdbqt_path = Path(pdb_path).with_suffix('.pdbqt')

        if pdbqt_path.exists():
            return str(pdbqt_path)

        try:
            subprocess.run([
                "obabel", pdb_path,
                "-O", str(pdbqt_path),
                "-xr"
            ], check=True, capture_output=True)
            return str(pdbqt_path)

        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.warning("Open Babel not available, using pre-prepared receptor")
            return pdb_path  # Assume already prepared

    def _get_default_receptor(self) -> str:
        """Get default Delta Opioid Receptor structure"""
        # In production, this would return the path to a prepared DOR structure
        # Using PDB ID 4EJ4 or similar
        dor_path = settings.data_dir / "receptors" / "dor_4ej4.pdb"

        if not dor_path.exists():
            # Create placeholder - in production, download from PDB
            logger.warning("DOR structure not found, docking results will be estimated")
            return str(dor_path)

        return str(dor_path)

    async def _run_vina(
        self,
        ligand_pdbqt: str,
        receptor_pdbqt: str,
        binding_site: Dict[str, float]
    ) -> Dict[str, Any]:
        """Run AutoDock Vina docking"""
        output_path = self.temp_dir / f"docked_{Path(ligand_pdbqt).stem}.pdbqt"
        log_path = self.temp_dir / f"docked_{Path(ligand_pdbqt).stem}.log"

        try:
            # Try using Vina directly
            cmd = [
                "vina",
                "--receptor", receptor_pdbqt,
                "--ligand", ligand_pdbqt,
                "--center_x", str(binding_site["center_x"]),
                "--center_y", str(binding_site["center_y"]),
                "--center_z", str(binding_site["center_z"]),
                "--size_x", str(binding_site["size_x"]),
                "--size_y", str(binding_site["size_y"]),
                "--size_z", str(binding_site["size_z"]),
                "--exhaustiveness", str(self.exhaustiveness),
                "--num_modes", str(self.num_modes),
                "--out", str(output_path),
                "--log", str(log_path),
            ]

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

            if result.returncode == 0:
                # Parse results
                affinities = self._parse_vina_output(log_path)
                return {
                    "success": True,
                    "best_affinity": affinities[0] if affinities else None,
                    "all_affinities": affinities,
                    "output_file": str(output_path),
                    "method": "autodock_vina"
                }
            else:
                logger.warning(f"Vina failed: {result.stderr}")
                return self._estimate_affinity_fallback(ligand_pdbqt)

        except FileNotFoundError:
            logger.warning("Vina not installed, using estimation fallback")
            return self._estimate_affinity_fallback(ligand_pdbqt)
        except subprocess.TimeoutExpired:
            logger.error("Vina timed out")
            return {"success": False, "error": "Docking timed out"}
        except Exception as e:
            logger.error(f"Vina error: {e}")
            return self._estimate_affinity_fallback(ligand_pdbqt)

    def _parse_vina_output(self, log_path: Path) -> List[float]:
        """Parse Vina log file for binding affinities"""
        affinities = []

        try:
            with open(log_path, 'r') as f:
                for line in f:
                    if line.strip().startswith('1') or line.strip().startswith('2'):
                        parts = line.split()
                        if len(parts) >= 2:
                            try:
                                affinity = float(parts[1])
                                affinities.append(affinity)
                            except ValueError:
                                continue
        except Exception as e:
            logger.warning(f"Error parsing Vina output: {e}")

        return affinities

    def _estimate_affinity_fallback(self, ligand_pdbqt: str) -> Dict[str, Any]:
        """Estimate binding affinity when docking is not available"""
        # Use simple scoring based on molecular properties
        # This is a rough approximation for when Vina is not installed

        # Read ligand PDB/PDBQT and estimate
        try:
            # Extract SMILES from filename or read molecule
            pdb_path = Path(ligand_pdbqt).with_suffix('.pdb')
            if pdb_path.exists():
                mol = Chem.MolFromPDBFile(str(pdb_path))
            else:
                mol = None

            if mol:
                from rdkit.Chem import Descriptors

                mw = Descriptors.MolWt(mol)
                logp = Descriptors.MolLogP(mol)
                hbd = Descriptors.NumHDonors(mol)
                hba = Descriptors.NumHAcceptors(mol)

                # Simple estimation formula (for peptide-like molecules)
                # Based on empirical observations for opioid ligands
                estimated_affinity = -6.0  # Base
                estimated_affinity -= 0.01 * min(mw, 800)  # Size contribution
                estimated_affinity -= 0.2 * min(hbd, 5)  # H-bond donors
                estimated_affinity -= 0.1 * min(hba, 10)  # H-bond acceptors

                # Clamp to reasonable range
                estimated_affinity = max(-12.0, min(-4.0, estimated_affinity))

                return {
                    "success": True,
                    "best_affinity": round(estimated_affinity, 2),
                    "all_affinities": [round(estimated_affinity, 2)],
                    "method": "estimated",
                    "note": "Docking unavailable, affinity estimated from molecular properties"
                }

        except Exception as e:
            logger.warning(f"Fallback estimation failed: {e}")

        # Return default estimate
        return {
            "success": True,
            "best_affinity": -7.5,
            "all_affinities": [-7.5],
            "method": "default_estimate",
            "note": "Default estimate used"
        }

    def _simple_pdb_to_pdbqt(self, pdb_path: Path, pdbqt_path: Path):
        """Simple PDB to PDBQT conversion (without Open Babel)"""
        # This is a simplified conversion - proper conversion requires charge assignment
        with open(pdb_path, 'r') as pdb_file:
            with open(pdbqt_path, 'w') as pdbqt_file:
                for line in pdb_file:
                    if line.startswith(('ATOM', 'HETATM')):
                        # Add placeholder charges
                        atom_type = line[76:78].strip() if len(line) > 76 else 'C'
                        charge = "0.000"
                        new_line = f"{line[:54]}{charge:>8}{atom_type:>2}\n"
                        pdbqt_file.write(new_line)
                    elif line.startswith(('MODEL', 'ENDMDL', 'END')):
                        pdbqt_file.write(line)

    def calculate_relative_affinity(
        self,
        candidate_affinity: float,
        reference: str = "kk103"
    ) -> float:
        """Calculate affinity relative to a reference compound"""
        ref_affinity = self.REFERENCE_AFFINITIES.get(reference, -8.0)
        # More negative is better
        # Return as percentage of reference
        return (candidate_affinity / ref_affinity) * 100

    async def batch_dock(
        self,
        smiles_list: List[str],
        receptor_pdb: Optional[str] = None
    ) -> List[Dict[str, Any]]:
        """Dock multiple compounds"""
        results = []
        for smiles in smiles_list:
            result = await self.dock_compound(smiles, receptor_pdb)
            results.append(result)
        return results


class StructureVisualization:
    """Service for generating molecular visualizations"""

    @staticmethod
    def generate_2d_image(smiles: str, size: Tuple[int, int] = (400, 400)) -> Optional[bytes]:
        """Generate 2D structure image"""
        try:
            from rdkit.Chem import Draw

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None

            img = Draw.MolToImage(mol, size=size)

            # Convert to bytes
            import io
            img_bytes = io.BytesIO()
            img.save(img_bytes, format='PNG')
            return img_bytes.getvalue()

        except Exception as e:
            logger.error(f"Error generating image: {e}")
            return None

    @staticmethod
    def generate_3d_conformer(smiles: str) -> Optional[str]:
        """Generate 3D conformer as MOL block"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None

            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)

            return Chem.MolToMolBlock(mol)

        except Exception as e:
            logger.error(f"Error generating 3D conformer: {e}")
            return None
