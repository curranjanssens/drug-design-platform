"""
Molecular analysis service using RDKit.
Provides structure handling, property calculation, and similarity analysis.
"""
from typing import Dict, List, Optional, Tuple, Any
import io
import base64

try:
    from rdkit import Chem
    from rdkit.Chem import (
        AllChem, Descriptors, Lipinski, Crippen,
        rdMolDescriptors, Draw, FilterCatalog,
        inchi as rdkit_inchi
    )
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit.DataStructs import TanimotoSimilarity
    from rdkit.Chem.QED import qed
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("Warning: RDKit not installed. Molecular analysis will be limited.")

from models import MolecularProperties, NoveltyAssessment


class MolecularAnalyzer:
    """
    Service for molecular structure analysis and property calculation.
    """

    def __init__(self):
        if not RDKIT_AVAILABLE:
            raise RuntimeError("RDKit is required for molecular analysis")

        # Initialize PAINS filter
        self.pains_catalog = self._init_pains_filter()

    def _init_pains_filter(self):
        """Initialize PAINS (Pan Assay Interference Compounds) filter."""
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        return FilterCatalog.FilterCatalog(params)

    def parse_smiles(self, smiles: str) -> Optional[Chem.Mol]:
        """Parse SMILES string to RDKit molecule."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            # Try to sanitize
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol is not None:
                try:
                    Chem.SanitizeMol(mol)
                except:
                    return None
        return mol

    def canonicalize_smiles(self, smiles: str) -> Optional[str]:
        """Convert SMILES to canonical form."""
        mol = self.parse_smiles(smiles)
        if mol:
            return Chem.MolToSmiles(mol, canonical=True)
        return None

    def smiles_to_inchi(self, smiles: str) -> Tuple[Optional[str], Optional[str]]:
        """Convert SMILES to InChI and InChIKey."""
        mol = self.parse_smiles(smiles)
        if mol:
            inchi = rdkit_inchi.MolToInchi(mol)
            inchi_key = rdkit_inchi.MolToInchiKey(mol)
            return inchi, inchi_key
        return None, None

    def calculate_properties(self, smiles: str) -> Optional[MolecularProperties]:
        """Calculate molecular properties from SMILES."""
        mol = self.parse_smiles(smiles)
        if mol is None:
            return None

        # Basic properties
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        tpsa = rdMolDescriptors.CalcTPSA(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        num_rings = rdMolDescriptors.CalcNumRings(mol)
        num_arom_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        formula = rdMolDescriptors.CalcMolFormula(mol)

        # Drug-likeness
        qed_score = qed(mol)

        # Lipinski violations
        lipinski_violations = 0
        if mw > 500:
            lipinski_violations += 1
        if logp > 5:
            lipinski_violations += 1
        if hbd > 5:
            lipinski_violations += 1
        if hba > 10:
            lipinski_violations += 1

        # Synthetic accessibility score
        sa_score = self._calculate_sa_score(mol)

        return MolecularProperties(
            molecular_weight=round(mw, 2),
            logp=round(logp, 2),
            tpsa=round(tpsa, 2),
            hbd=hbd,
            hba=hba,
            rotatable_bonds=rot_bonds,
            num_rings=num_rings,
            num_aromatic_rings=num_arom_rings,
            formula=formula,
            qed=round(qed_score, 3),
            lipinski_violations=lipinski_violations,
            sa_score=round(sa_score, 2)
        )

    def _calculate_sa_score(self, mol: Chem.Mol) -> float:
        """
        Calculate synthetic accessibility score (1-10).
        Lower is easier to synthesize.
        Simplified implementation.
        """
        # This is a simplified SA score
        # Full implementation would use Ertl's SA_Score
        num_rings = rdMolDescriptors.CalcNumRings(mol)
        num_stereo = len(Chem.FindMolChiralCenters(mol))
        num_spiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
        num_bridgehead = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
        mw = Descriptors.MolWt(mol)

        # Simple heuristic
        score = 1.0
        score += num_stereo * 0.5
        score += num_spiro * 1.0
        score += num_bridgehead * 1.0
        score += max(0, num_rings - 3) * 0.3
        score += max(0, (mw - 300) / 100) * 0.2

        return min(10.0, score)

    def calculate_fingerprint(
        self,
        smiles: str,
        fp_type: str = "morgan",
        radius: int = 2,
        n_bits: int = 2048
    ):
        """Calculate molecular fingerprint."""
        mol = self.parse_smiles(smiles)
        if mol is None:
            return None

        if fp_type == "morgan":
            return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
        elif fp_type == "rdkit":
            return Chem.RDKFingerprint(mol)
        elif fp_type == "maccs":
            return AllChem.GetMACCSKeysFingerprint(mol)
        else:
            raise ValueError(f"Unknown fingerprint type: {fp_type}")

    def calculate_similarity(
        self,
        smiles1: str,
        smiles2: str,
        fp_type: str = "morgan"
    ) -> Optional[float]:
        """Calculate Tanimoto similarity between two molecules."""
        fp1 = self.calculate_fingerprint(smiles1, fp_type)
        fp2 = self.calculate_fingerprint(smiles2, fp_type)

        if fp1 is None or fp2 is None:
            return None

        return TanimotoSimilarity(fp1, fp2)

    def find_nearest_neighbor(
        self,
        query_smiles: str,
        database: List[str],
        top_k: int = 5
    ) -> List[Tuple[str, float]]:
        """Find most similar molecules in a database."""
        query_fp = self.calculate_fingerprint(query_smiles)
        if query_fp is None:
            return []

        similarities = []
        for db_smiles in database:
            db_fp = self.calculate_fingerprint(db_smiles)
            if db_fp:
                sim = TanimotoSimilarity(query_fp, db_fp)
                similarities.append((db_smiles, sim))

        return sorted(similarities, key=lambda x: x[1], reverse=True)[:top_k]

    def get_scaffold(self, smiles: str, generic: bool = False) -> Optional[str]:
        """Extract Murcko scaffold from molecule."""
        mol = self.parse_smiles(smiles)
        if mol is None:
            return None

        try:
            if generic:
                scaffold = MurckoScaffold.MakeScaffoldGeneric(
                    MurckoScaffold.GetScaffoldForMol(mol)
                )
            else:
                scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            return Chem.MolToSmiles(scaffold)
        except:
            return None

    def has_pains(self, smiles: str) -> Tuple[bool, List[str]]:
        """Check if molecule contains PAINS patterns."""
        mol = self.parse_smiles(smiles)
        if mol is None:
            return False, []

        matches = self.pains_catalog.GetMatches(mol)
        if matches:
            patterns = [match.GetDescription() for match in matches]
            return True, patterns
        return False, []

    def has_substructure(self, smiles: str, substructure_smarts: str) -> bool:
        """Check if molecule contains a substructure."""
        mol = self.parse_smiles(smiles)
        pattern = Chem.MolFromSmarts(substructure_smarts)

        if mol is None or pattern is None:
            return False

        return mol.HasSubstructMatch(pattern)

    def generate_conformer(self, smiles: str, num_confs: int = 1) -> Optional[str]:
        """Generate 3D conformer and return as MOL block."""
        mol = self.parse_smiles(smiles)
        if mol is None:
            return None

        mol = Chem.AddHs(mol)

        # Generate conformers
        result = AllChem.EmbedMolecule(mol, randomSeed=42)
        if result != 0:
            # Fallback to random coords
            AllChem.EmbedMolecule(mol, useRandomCoords=True)

        # Optimize
        try:
            AllChem.MMFFOptimizeMolecule(mol)
        except:
            try:
                AllChem.UFFOptimizeMolecule(mol)
            except:
                pass

        return Chem.MolToMolBlock(mol)

    def molecule_to_image(
        self,
        smiles: str,
        size: Tuple[int, int] = (300, 300),
        highlight_atoms: List[int] = None
    ) -> Optional[str]:
        """Generate molecule image as base64 encoded PNG."""
        mol = self.parse_smiles(smiles)
        if mol is None:
            return None

        # Generate 2D coordinates
        AllChem.Compute2DCoords(mol)

        # Draw
        drawer = Draw.MolDraw2DCairo(size[0], size[1])
        if highlight_atoms:
            drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms)
        else:
            drawer.DrawMolecule(mol)
        drawer.FinishDrawing()

        # Convert to base64
        png_data = drawer.GetDrawingText()
        return base64.b64encode(png_data).decode('utf-8')

    def molecule_to_svg(self, smiles: str, size: Tuple[int, int] = (300, 300)) -> Optional[str]:
        """Generate molecule image as SVG."""
        mol = self.parse_smiles(smiles)
        if mol is None:
            return None

        AllChem.Compute2DCoords(mol)

        drawer = Draw.MolDraw2DSVG(size[0], size[1])
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()

        return drawer.GetDrawingText()

    def assess_novelty(
        self,
        smiles: str,
        known_compounds: List[Dict[str, str]],
        patent_hits: List[Dict[str, Any]] = None
    ) -> NoveltyAssessment:
        """
        Assess novelty of a compound against known compounds and patents.

        Args:
            smiles: SMILES of query compound
            known_compounds: List of {"smiles": "...", "name": "..."} dicts
            patent_hits: Optional list of patent search results
        """
        # Find nearest known compound
        similarities = []
        for known in known_compounds:
            sim = self.calculate_similarity(smiles, known["smiles"])
            if sim is not None:
                similarities.append((known, sim))

        similarities.sort(key=lambda x: x[1], reverse=True)

        if similarities:
            nearest = similarities[0]
            nearest_smiles = nearest[0]["smiles"]
            nearest_name = nearest[0].get("name")
            nearest_sim = nearest[1]
        else:
            nearest_smiles = ""
            nearest_name = None
            nearest_sim = 0.0

        # Check scaffold novelty
        query_scaffold = self.get_scaffold(smiles, generic=True)
        scaffold_novel = True
        for known in known_compounds:
            known_scaffold = self.get_scaffold(known["smiles"], generic=True)
            if query_scaffold == known_scaffold:
                scaffold_novel = False
                break

        # Analyze patent hits
        exact_match = False
        markush_coverage = []
        blocking_patents = []
        similar_patents = []
        expiry_dates = {}

        if patent_hits:
            for hit in patent_hits:
                if hit.get("exact_match"):
                    exact_match = True
                    blocking_patents.append(hit["patent_id"])
                elif hit.get("markush_match"):
                    markush_coverage.append(hit["patent_id"])
                elif hit.get("similarity", 0) > 0.7:
                    similar_patents.append(hit)

                if hit.get("expiry_date"):
                    expiry_dates[hit["patent_id"]] = hit["expiry_date"]

        # Determine FTO risk
        if blocking_patents:
            fto_risk = "high"
        elif markush_coverage:
            fto_risk = "medium"
        elif similar_patents:
            fto_risk = "low"
        else:
            fto_risk = "low"

        # Conclusions
        novel = not exact_match and nearest_sim < 0.9
        patentable = novel and not markush_coverage and scaffold_novel

        # Confidence
        confidence = 0.7  # Base confidence
        if not patent_hits:
            confidence -= 0.2  # Lower confidence without patent data
        if len(known_compounds) < 10:
            confidence -= 0.1  # Lower confidence with small comparison set

        return NoveltyAssessment(
            nearest_known_smiles=nearest_smiles,
            nearest_known_name=nearest_name,
            tanimoto_similarity=round(nearest_sim, 3),
            scaffold_novel=scaffold_novel,
            exact_patent_match=exact_match,
            markush_coverage=markush_coverage,
            similar_patents=similar_patents,
            fto_risk=fto_risk,
            blocking_patents=blocking_patents,
            expiry_dates=expiry_dates,
            novel=novel,
            patentable=patentable,
            confidence=confidence,
            notes=self._generate_novelty_notes(novel, patentable, nearest_sim, scaffold_novel)
        )

    def _generate_novelty_notes(
        self,
        novel: bool,
        patentable: bool,
        similarity: float,
        scaffold_novel: bool
    ) -> str:
        """Generate human-readable novelty notes."""
        notes = []

        if similarity > 0.9:
            notes.append(f"Very similar to known compound (Tanimoto={similarity:.2f})")
        elif similarity > 0.7:
            notes.append(f"Moderately similar to known compounds (Tanimoto={similarity:.2f})")
        else:
            notes.append(f"Structurally distinct from known compounds (Tanimoto={similarity:.2f})")

        if scaffold_novel:
            notes.append("Novel scaffold identified")
        else:
            notes.append("Scaffold present in prior art")

        if novel and patentable:
            notes.append("Compound appears novel and potentially patentable")
        elif novel:
            notes.append("Compound appears novel but patentability uncertain")
        else:
            notes.append("Compound may not meet novelty requirements")

        return ". ".join(notes) + "."

    def validate_molecule(self, smiles: str) -> Dict[str, Any]:
        """Comprehensive validation of a molecule."""
        result = {
            "valid": False,
            "canonical_smiles": None,
            "errors": [],
            "warnings": []
        }

        # Try to parse
        mol = self.parse_smiles(smiles)
        if mol is None:
            result["errors"].append("Invalid SMILES string")
            return result

        result["valid"] = True
        result["canonical_smiles"] = Chem.MolToSmiles(mol, canonical=True)

        # Check for PAINS
        has_pains, pains_patterns = self.has_pains(smiles)
        if has_pains:
            result["warnings"].append(f"Contains PAINS patterns: {', '.join(pains_patterns)}")

        # Check properties
        props = self.calculate_properties(smiles)
        if props:
            if props.lipinski_violations > 1:
                result["warnings"].append(f"Lipinski violations: {props.lipinski_violations}")
            if props.sa_score > 6:
                result["warnings"].append(f"Difficult to synthesize (SA score: {props.sa_score})")
            if props.molecular_weight > 800:
                result["warnings"].append(f"Very high molecular weight: {props.molecular_weight}")

        return result
