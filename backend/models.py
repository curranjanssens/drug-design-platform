"""
Pydantic models for the Drug Design Platform API
"""
from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any
from enum import Enum
from datetime import datetime


class MolecularProperties(BaseModel):
    """Molecular properties calculated by RDKit"""
    molecular_weight: float
    logp: float
    tpsa: float
    hbd: int
    hba: int
    rotatable_bonds: int
    num_rings: int
    num_aromatic_rings: int
    formula: str
    qed: float
    lipinski_violations: int
    sa_score: float


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
    """Novelty and patentability assessment from molecular analyzer"""
    nearest_known_smiles: str = ""
    nearest_known_name: Optional[str] = None
    tanimoto_similarity: float = 0.0
    scaffold_novel: bool = True
    exact_patent_match: bool = False
    markush_coverage: List[str] = []
    similar_patents: List[Dict[str, Any]] = []
    fto_risk: str = "low"
    blocking_patents: List[str] = []
    expiry_dates: Dict[str, str] = {}
    novel: bool = True
    patentable: bool = True
    confidence: float = 0.7
    notes: str = ""


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


# ============== Additional Models for Services ==============

class JobStatus(str, Enum):
    """Status of a design job"""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


class DesignType(str, Enum):
    """Types of molecular design"""
    ANALOGUE = "analogue"
    BIVALENT = "bivalent"
    SCAFFOLD_HOP = "scaffold_hop"
    DE_NOVO = "de_novo"


class ADMETProfile(BaseModel):
    """ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity) profile"""
    # Absorption
    caco2_permeability: float = -5.5
    human_intestinal_absorption: float = 80.0
    pgp_substrate: bool = False
    pgp_inhibitor: bool = False

    # Distribution
    vdss: float = 0.5
    bbb_penetration: bool = False
    cns_penetration: bool = False
    plasma_protein_binding: float = 90.0

    # Metabolism
    cyp1a2_inhibitor: bool = False
    cyp2c9_inhibitor: bool = False
    cyp2c19_inhibitor: bool = False
    cyp2d6_inhibitor: bool = False
    cyp3a4_inhibitor: bool = False
    cyp2d6_substrate: bool = False
    cyp3a4_substrate: bool = True

    # Excretion
    clearance: float = 10.0
    half_life: float = 6.0

    # Toxicity
    herg_inhibitor: bool = False
    herg_pic50: Optional[float] = None
    ames_toxicity: bool = False
    hepatotoxicity: bool = False
    skin_sensitization: bool = False
    ld50: Optional[float] = None

    warnings: List[str] = []

    def get_safety_score(self) -> float:
        """Calculate overall safety score (0-1, higher is safer)"""
        score = 1.0
        if self.herg_inhibitor:
            score -= 0.3
        if self.ames_toxicity:
            score -= 0.3
        if self.hepatotoxicity:
            score -= 0.2
        if len(self.warnings) > 3:
            score -= 0.1
        return max(0.0, score)


class ReactionStep(BaseModel):
    """A single step in a synthesis route"""
    step_number: int
    reaction_smiles: str
    reaction_name: str
    reactants: List[str] = []
    reagents: List[str] = []
    solvent: Optional[str] = None
    temperature: Optional[str] = None
    time: Optional[str] = None
    atmosphere: Optional[str] = None
    expected_yield: float = 0.8
    confidence: float = 0.8
    reference: Optional[str] = None


class SynthesisRoute(BaseModel):
    """A complete synthesis route"""
    steps: List[ReactionStep] = []
    starting_materials: List[Dict[str, Any]] = []
    source: str = "unknown"
    total_steps: int = 0
    overall_yield: float = 0.0
    estimated_cost: Optional[float] = None
    confidence: float = 0.5

    def calculate_metrics(self):
        """Calculate route metrics"""
        self.total_steps = len(self.steps)
        if self.steps:
            self.overall_yield = 1.0
            for step in self.steps:
                self.overall_yield *= step.expected_yield
            self.confidence = sum(s.confidence for s in self.steps) / len(self.steps)


class DesignedMolecule(BaseModel):
    """A designed molecule with all associated data"""
    rank: int = 0
    smiles: str
    inchi_key: Optional[str] = None
    overall_score: float = 0.0
    properties: Optional[MolecularProperties] = None
    novelty: Optional[NoveltyAssessment] = None
    admet: Optional[ADMETProfile] = None
    best_route: Optional[SynthesisRoute] = None
    rationale: str = ""


class BindingPrediction(BaseModel):
    """Binding affinity prediction"""
    value: float = 0.0
    unit: str = "pKi"
    confidence: float = 0.5
    method: str = "unknown"


class DesignJob(BaseModel):
    """A drug design job tracking state and results"""
    id: str = ""
    prompt: str = ""
    user_id: Optional[str] = None
    status: JobStatus = JobStatus.PENDING
    design_type: DesignType = DesignType.ANALOGUE
    current_step: str = ""
    progress: float = 0.0

    # Reference compound info
    reference_smiles: Optional[str] = None
    reference_name: Optional[str] = None
    target_name: Optional[str] = None
    constraints: Dict[str, Any] = {}

    # Timestamps
    created_at: datetime = None
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None

    # Results
    molecules: List[DesignedMolecule] = []
    error_message: Optional[str] = None

    class Config:
        arbitrary_types_allowed = True

    def __init__(self, **data):
        if "id" not in data or not data["id"]:
            import uuid
            data["id"] = str(uuid.uuid4())
        if "created_at" not in data or data["created_at"] is None:
            data["created_at"] = datetime.utcnow()
        super().__init__(**data)

    def update_progress(self, step: str, progress: float):
        """Update job progress"""
        self.current_step = step
        self.progress = progress

    def complete(self, molecules: List[DesignedMolecule]):
        """Mark job as completed"""
        self.status = JobStatus.COMPLETED
        self.molecules = molecules
        self.completed_at = datetime.utcnow()
        self.progress = 1.0

    def fail(self, error: str):
        """Mark job as failed"""
        self.status = JobStatus.FAILED
        self.error_message = error
        self.completed_at = datetime.utcnow()
