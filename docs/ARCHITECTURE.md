# Automated Drug Design Platform - System Architecture

## Overview

A fully automated platform that takes a target compound/therapeutic goal as input and produces chemist-ready novel, patentable molecule designs with complete synthesis routes.

## System Flow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              USER INPUT                                       │
│  "Design a novel analogue of KK103 with improved stability"                  │
│  "Create a dual 5HT1A-5HT2A agonist bivalent ligand"                        │
└─────────────────────────────────┬───────────────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         CLAUDE ORCHESTRATOR                                   │
│  • Parses natural language request                                           │
│  • Determines target type (analogue, de novo, bivalent, etc.)               │
│  • Plans multi-step workflow                                                 │
│  • Coordinates all downstream services                                       │
└─────────────────────────────────┬───────────────────────────────────────────┘
                                  │
          ┌───────────────────────┼───────────────────────┐
          ▼                       ▼                       ▼
┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│  TARGET ANALYSIS │    │  LITERATURE     │    │  PATENT         │
│  SERVICE         │    │  SEARCH         │    │  LANDSCAPE      │
│                 │    │                 │    │                 │
│  • PDB lookup   │    │  • PubMed API   │    │  • SureChEMBL   │
│  • UniProt      │    │  • SAR data     │    │  • EPO/USPTO    │
│  • ChEMBL       │    │  • Known drugs  │    │  • FTO analysis │
└────────┬────────┘    └────────┬────────┘    └────────┬────────┘
         │                      │                      │
         └──────────────────────┼──────────────────────┘
                                ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                      MOLECULE GENERATION ENGINE                              │
│                                                                             │
│  ┌─────────────────┐  ┌─────────────────┐  ┌─────────────────┐             │
│  │   REINVENT4     │  │   Scaffold      │  │   De Novo       │             │
│  │   Optimization  │  │   Hopping       │  │   Generation    │             │
│  └─────────────────┘  └─────────────────┘  └─────────────────┘             │
│                                                                             │
│  Multi-objective optimization:                                              │
│  • Binding affinity (predicted)                                             │
│  • Novelty (Tanimoto distance from known)                                  │
│  • Drug-likeness (Lipinski, QED)                                           │
│  • Synthetic accessibility (SAScore)                                       │
│  • Selectivity (off-target predictions)                                    │
└─────────────────────────────────┬───────────────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                      STRUCTURE PREDICTION                                    │
│                                                                             │
│  ┌─────────────────────────────────────────────────────────────┐           │
│  │  AlphaFold 3 / Boltz-1 / Chai-1                             │           │
│  │  • Protein-ligand complex prediction                        │           │
│  │  • Binding pose generation                                  │           │
│  │  • Interface confidence (ipTM) scoring                      │           │
│  └─────────────────────────────────────────────────────────────┘           │
│                                                                             │
│  ┌─────────────────────────────────────────────────────────────┐           │
│  │  Molecular Docking (AutoDock Vina / Gnina)                  │           │
│  │  • Binding affinity estimation                              │           │
│  │  • Multiple pose generation                                 │           │
│  │  • Cross-validation with AF3                                │           │
│  └─────────────────────────────────────────────────────────────┘           │
└─────────────────────────────────┬───────────────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                      ADMET & PROPERTY PREDICTION                             │
│                                                                             │
│  • Absorption: Caco-2 permeability, HIA, P-gp substrate                    │
│  • Distribution: VDss, BBB penetration, PPB                                │
│  • Metabolism: CYP450 inhibition/substrate                                 │
│  • Excretion: Clearance, half-life                                         │
│  • Toxicity: hERG, AMES, hepatotoxicity, LD50                              │
│                                                                             │
│  Tools: ADMETlab 2.0, pkCSM, DeepChem, SwissADME                           │
└─────────────────────────────────┬───────────────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                      NOVELTY & PATENTABILITY CHECK                           │
│                                                                             │
│  ┌─────────────────────────────────────────────────────────────┐           │
│  │  Structural Novelty                                         │           │
│  │  • Tanimoto similarity to known compounds                   │           │
│  │  • Substructure search in ChEMBL/PubChem                    │           │
│  │  • Scaffold novelty assessment                              │           │
│  └─────────────────────────────────────────────────────────────┘           │
│                                                                             │
│  ┌─────────────────────────────────────────────────────────────┐           │
│  │  Patent Analysis                                            │           │
│  │  • SureChEMBL structure search                              │           │
│  │  • Markush structure analysis                               │           │
│  │  • Freedom-to-operate preliminary check                     │           │
│  └─────────────────────────────────────────────────────────────┘           │
└─────────────────────────────────┬───────────────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                      RETROSYNTHESIS PLANNING                                 │
│                                                                             │
│  Primary: IBM RXN for Chemistry (API)                                       │
│  Secondary: AiZynthFinder (self-hosted)                                     │
│  Validation: ASKCOS (cross-reference)                                       │
│                                                                             │
│  Output:                                                                    │
│  • Multiple synthesis routes (ranked)                                       │
│  • Starting materials with availability/cost                                │
│  • Reaction conditions for each step                                        │
│  • Confidence scores                                                        │
└─────────────────────────────────┬───────────────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                      OUTPUT GENERATOR                                        │
│                                                                             │
│  Chemist-Ready Deliverables:                                                │
│  ┌─────────────────────────────────────────────────────────────┐           │
│  │  1. MOLECULE SPECIFICATION                                  │           │
│  │     • SMILES, InChI, InChIKey                               │           │
│  │     • 2D structure image (SVG/PNG)                          │           │
│  │     • 3D structure file (SDF/MOL2)                          │           │
│  │     • Molecular formula, MW, key properties                 │           │
│  └─────────────────────────────────────────────────────────────┘           │
│  ┌─────────────────────────────────────────────────────────────┐           │
│  │  2. BINDING ANALYSIS                                        │           │
│  │     • Predicted binding pose (PDB file)                     │           │
│  │     • Key interactions (H-bonds, hydrophobic, etc.)         │           │
│  │     • Predicted binding affinity                            │           │
│  │     • Comparison to reference compound                      │           │
│  └─────────────────────────────────────────────────────────────┘           │
│  ┌─────────────────────────────────────────────────────────────┐           │
│  │  3. SYNTHESIS PLAN                                          │           │
│  │     • Step-by-step synthesis route                          │           │
│  │     • Reaction schemes (image)                              │           │
│  │     • Reagents, conditions, expected yields                 │           │
│  │     • Starting materials with vendors/pricing               │           │
│  │     • Alternative routes ranked                             │           │
│  └─────────────────────────────────────────────────────────────┘           │
│  ┌─────────────────────────────────────────────────────────────┐           │
│  │  4. ADMET PROFILE                                           │           │
│  │     • Full property predictions table                       │           │
│  │     • Flags/warnings for problematic properties             │           │
│  │     • Comparison to drug-likeness criteria                  │           │
│  └─────────────────────────────────────────────────────────────┘           │
│  ┌─────────────────────────────────────────────────────────────┐           │
│  │  5. NOVELTY REPORT                                          │           │
│  │     • Structural novelty score                              │           │
│  │     • Nearest known compounds                               │           │
│  │     • Patent landscape summary                              │           │
│  │     • Preliminary FTO assessment                            │           │
│  └─────────────────────────────────────────────────────────────┘           │
└─────────────────────────────────────────────────────────────────────────────┘
```

## Technology Stack

### Backend
- **Framework**: FastAPI (Python)
- **AI Orchestration**: Claude API (claude-sonnet-4-20250514 or claude-opus-4-5-20251101)
- **Task Queue**: Celery + Redis (for long-running computations)
- **Database**: PostgreSQL (results, projects) + Redis (caching)

### Cheminformatics Core
- **RDKit**: Molecular manipulation, fingerprints, similarity
- **Open Babel**: Format conversion, 3D generation
- **Meeko**: AutoDock preparation

### ML/AI Services
- **Structure Prediction**: Boltz-1 (MIT license) or Chai-1
- **Molecule Generation**: REINVENT4
- **Property Prediction**: DeepChem, ChemBERTa

### External APIs
- **Retrosynthesis**: IBM RXN for Chemistry
- **ADMET**: ADMETlab 2.0 API, pkCSM
- **Patent Search**: SureChEMBL, EPO OPS
- **Chemical Data**: PubChem, ChEMBL, UniProt

### Frontend
- **Framework**: React + TypeScript
- **Visualization**: 3Dmol.js (molecular), Ketcher (drawing)
- **UI**: Tailwind CSS

### Infrastructure
- **Containerization**: Docker + Docker Compose
- **GPU**: NVIDIA GPU for ML models (optional cloud GPU)

## API Endpoints

```
POST /api/v1/design
  - Input: { prompt: string, constraints?: object }
  - Output: { job_id: string }

GET /api/v1/design/{job_id}
  - Output: { status, progress, results? }

GET /api/v1/design/{job_id}/molecules
  - Output: List of designed molecules with properties

GET /api/v1/design/{job_id}/synthesis/{molecule_id}
  - Output: Synthesis routes for specific molecule

POST /api/v1/analyze
  - Input: { smiles: string }
  - Output: { properties, admet, novelty }

POST /api/v1/retrosynthesis
  - Input: { smiles: string, options? }
  - Output: { routes: Route[] }
```

## Data Models

```python
class DesignJob:
    id: str
    prompt: str
    status: Enum['pending', 'running', 'completed', 'failed']
    created_at: datetime
    molecules: List[DesignedMolecule]

class DesignedMolecule:
    id: str
    smiles: str
    inchi: str
    properties: MolecularProperties
    binding_prediction: BindingPrediction
    admet: ADMETProfile
    novelty: NoveltyAssessment
    synthesis_routes: List[SynthesisRoute]

class SynthesisRoute:
    id: str
    steps: List[ReactionStep]
    total_steps: int
    confidence: float
    estimated_cost: float
    starting_materials: List[Chemical]

class ReactionStep:
    reaction_smiles: str
    conditions: ReactionConditions
    expected_yield: float
    reference: str
```

## Deployment Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                      Load Balancer                          │
└─────────────────────────┬───────────────────────────────────┘
                          │
        ┌─────────────────┼─────────────────┐
        ▼                 ▼                 ▼
┌───────────────┐ ┌───────────────┐ ┌───────────────┐
│   Web Server  │ │   Web Server  │ │   Web Server  │
│   (FastAPI)   │ │   (FastAPI)   │ │   (FastAPI)   │
└───────┬───────┘ └───────┬───────┘ └───────┬───────┘
        │                 │                 │
        └─────────────────┼─────────────────┘
                          │
┌─────────────────────────▼───────────────────────────────────┐
│                      Redis (Queue/Cache)                     │
└─────────────────────────┬───────────────────────────────────┘
                          │
        ┌─────────────────┼─────────────────┐
        ▼                 ▼                 ▼
┌───────────────┐ ┌───────────────┐ ┌───────────────┐
│ Celery Worker │ │ Celery Worker │ │ GPU Worker    │
│ (CPU tasks)   │ │ (CPU tasks)   │ │ (ML models)   │
└───────────────┘ └───────────────┘ └───────────────┘
                          │
                          ▼
              ┌───────────────────────┐
              │    PostgreSQL DB      │
              └───────────────────────┘
```

## Environment Variables

```bash
# Claude API
ANTHROPIC_API_KEY=

# IBM RXN
IBM_RXN_API_KEY=

# Database
DATABASE_URL=postgresql://...
REDIS_URL=redis://...

# External Services
PUBCHEM_API_KEY=
CHEMBL_API_KEY=

# GPU Config
CUDA_VISIBLE_DEVICES=0
```
