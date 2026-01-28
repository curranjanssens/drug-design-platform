"""
Pydantic models for the Drug Design Platform API
"""
from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any
from enum import Enum
from datetime import datetime


class InputFormat(str, Enum):
    """Supported input formats"""
    SMILES = "smiles"
    PEPTIDE_SEQUENCE = "peptide_sequence"
    MOL_FILE = "mol_file"
    SDF_FILE = "sdf_file"
    PDB_FILE = "pdb_file"


class CompoundInput(BaseModel):
    """Input compound specification"""
    name: str = Field(..., description="Name or identifier for the compound")
    structure: str = Field(..., description="Molecular structure (SMILES, sequence, or file content)")
    format: InputFormat = Field(default=InputFormat.SMILES, description="Format of the structure")
    target_receptor: Optional[str] = Field(None, description="Target receptor (e.g., 'DOR' for delta opioid receptor)")
    target_pdb_id: Optional[str] = Field(None, description="PDB ID of target receptor structure")
    constraints: Optional[Dict[str, Any]] = Field(None, description="Property constraints for analogues")


class PropertyPrediction(BaseModel):
    """Predicted properties for a compound"""
    molecular_weight: float
    logp: float
    hbd: int  # Hydrogen bond donors
    hba: int  # Hydrogen bond acceptors
    tpsa: float  # Topological polar surface area
    rotatable_bonds: int
    qed: float  # Quantitative Estimate of Drug-likeness
    synthetic_accessibility: float
    lipinski_violations: int
    predicted_solubility: Optional[str] = None
    predicted_half_life: Optional[str] = None
    predicted_binding_affinity: Optional[float] = None


class NoveltyAssessment(BaseModel):
    """Novelty and patentability assessment"""
    is_novel: bool
    confidence_score: float = Field(..., ge=0, le=1)
    closest_known_compound: Optional[str] = None
    closest_similarity: Optional[float] = None
    pubchem_matches: int = 0
    chembl_matches: int = 0
    patent_matches: int = 0
    non_obviousness_rationale: Optional[str] = None


class SynthesisGuidance(BaseModel):
    """Synthesis feasibility and guidance"""
    synthetic_accessibility_score: float = Field(..., ge=1, le=10)
    feasibility_rating: str  # "Easy", "Moderate", "Difficult", "Very Difficult"
    key_building_blocks: List[str] = []
    suggested_route: Optional[str] = None
    estimated_steps: Optional[int] = None
    potential_challenges: List[str] = []


class AnalogueCandidate(BaseModel):
    """A generated analogue candidate"""
    id: str
    name: str
    smiles: str
    inchi: Optional[str] = None
    inchi_key: Optional[str] = None
    modification_description: str
    modification_type: str  # e.g., "N-terminal", "amino_acid_substitution", "backbone"
    properties: PropertyPrediction
    novelty: NoveltyAssessment
    synthesis: SynthesisGuidance
    rank_score: float = Field(..., ge=0, le=1)
    mol_file_content: Optional[str] = None
    binding_pose_file: Optional[str] = None
    rationale: str


class AnalogueGenerationRequest(BaseModel):
    """Request for analogue generation"""
    compound: CompoundInput
    num_analogues: int = Field(default=10, ge=1, le=50)
    strategies: List[str] = Field(
        default=["n_terminal", "amino_acid_substitution", "c_terminal", "backbone"],
        description="Modification strategies to apply"
    )
    include_docking: bool = Field(default=True, description="Perform molecular docking")
    include_alphafold: bool = Field(default=False, description="Use AlphaFold for structure prediction")


class AnalogueGenerationResponse(BaseModel):
    """Response from analogue generation"""
    request_id: str
    input_compound: CompoundInput
    generated_at: datetime
    total_candidates: int
    passed_filters: int
    candidates: List[AnalogueCandidate]
    generation_summary: str
    warnings: List[str] = []


class PipelineStatus(BaseModel):
    """Status of the generation pipeline"""
    request_id: str
    status: str  # "pending", "running", "completed", "failed"
    current_step: str
    progress: float = Field(..., ge=0, le=1)
    steps_completed: List[str] = []
    errors: List[str] = []
    estimated_time_remaining: Optional[int] = None  # seconds


# KK103 Reference Data
KK103_REFERENCE = {
    "name": "KK103",
    "full_name": "N-pivaloyl-Leu-Enkephalin",
    "sequence": "YGGFL",
    "n_terminal": "pivaloyl",
    "c_terminal": "carboxylic_acid",
    "estimated_smiles": "CC(C)(C)C(=O)NC(Cc1ccc(O)cc1)C(=O)NCC(=O)NCC(=O)NC(Cc2ccccc2)C(=O)NC(CC(C)C)C(=O)O",
    "molecular_formula": "C33H45N5O8",
    "target_receptor": "DOR",
    "dor_binding_affinity_percent": 68,  # relative to Leu-ENK
    "plasma_half_life_hours": 37,
    "therapeutic_use": "Antinociception"
}

# Leu-Enkephalin Reference Data
LEU_ENKEPHALIN_REFERENCE = {
    "name": "Leu-Enkephalin",
    "sequence": "YGGFL",
    "smiles": "CC(C)CC(NC(=O)C(Cc1ccccc1)NC(=O)CNC(=O)CNC(=O)C(Cc2ccc(O)cc2)N)C(=O)O",
    "molecular_formula": "C28H37N5O7",
    "molecular_weight": 555.62,
    "target_receptor": "DOR",
    "plasma_half_life_minutes": 2
}
