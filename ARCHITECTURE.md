# Drug Design Platform Architecture

## System Overview

This platform transforms natural language prompts into validated, novel drug candidate SMILES strings through an AI-driven iterative design process.

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         DRUG DESIGN PLATFORM                                 │
│                     Prompt → SMILES Pipeline                                 │
└─────────────────────────────────────────────────────────────────────────────┘

┌──────────────┐     ┌──────────────┐     ┌──────────────┐     ┌──────────────┐
│   Frontend   │────▶│   FastAPI    │────▶│   Agentic    │────▶│    SMILES    │
│   (Browser)  │     │   Backend    │     │   Designer   │     │   Output     │
└──────────────┘     └──────────────┘     └──────────────┘     └──────────────┘
                            │                    │
                            ▼                    ▼
                     ┌──────────────┐     ┌──────────────┐
                     │   Claude     │     │   External   │
                     │   Opus 4.5   │     │   APIs       │
                     └──────────────┘     └──────────────┘
                                                │
                            ┌───────────────────┼───────────────────┐
                            ▼                   ▼                   ▼
                     ┌──────────────┐   ┌──────────────┐   ┌──────────────┐
                     │   ChEMBL     │   │   PubChem    │   │   RDKit      │
                     │   Database   │   │   Database   │   │   Library    │
                     └──────────────┘   └──────────────┘   └──────────────┘
```

---

## Complete Data Flow: Prompt → SMILES

```
USER PROMPT
    │
    │  "Design a novel NOP/MOP dual agonist like cebranopadol"
    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ PHASE 1: REQUEST PARSING                                            [5%]   │
├─────────────────────────────────────────────────────────────────────────────┤
│ • Claude extracts: design_type, compound_names, target_names               │
│ • PubChem lookup: Get reference SMILES for mentioned compounds             │
│ • ChEMBL lookup: Find target IDs for mentioned receptors                   │
└─────────────────────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ PHASE 2: DEEP TARGET RESEARCH                                       [15%]  │
├─────────────────────────────────────────────────────────────────────────────┤
│ Claude Opus 4.5 researches:                                                 │
│ • Target biology & disease relevance                                        │
│ • Binding mechanism (reversible vs covalent)                                │
│ • Key binding site residues                                                 │
│ • Essential structural features                                             │
│ • Optimal property ranges (MW, LogP)                                        │
│                                                                             │
│ OUTPUT: TargetKnowledge object populated                                    │
└─────────────────────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ PHASE 3: ChEMBL BIOACTIVITY DATA                                    [25%]  │
├─────────────────────────────────────────────────────────────────────────────┤
│ Query: GET https://www.ebi.ac.uk/chembl/api/data/activity.json             │
│                                                                             │
│ Returns real bioactivity data:                                              │
│ • Known active compounds with SMILES                                        │
│ • IC50, Ki, EC50, Kd values (in nM)                                        │
│ • Assay descriptions                                                        │
│                                                                             │
│ Example: "PF-04457845 Ki=0.7nM", "URB597 IC50=4.6nM"                       │
└─────────────────────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ PHASE 4: MECHANISTIC ANALYSIS                                       [45%]  │
├─────────────────────────────────────────────────────────────────────────────┤
│ Claude determines mechanism from first principles:                          │
│                                                                             │
│ IF COVALENT INHIBITOR:                                                      │
│ ┌─────────────────────────────────────────────────────────────────────────┐│
│ │ • Identify warhead type (urea, carbamate, acrylamide, etc.)            ││
│ │ • Determine nucleophile (Ser241, Cys797, etc.)                         ││
│ │ • Define STAYING portion requirements (must have binding groups)       ││
│ │ • Define LEAVING group requirements (good leaving group ability)       ││
│ │ • Extract mechanistic constraints                                       ││
│ └─────────────────────────────────────────────────────────────────────────┘│
│                                                                             │
│ IF REVERSIBLE: Standard binding principles apply                            │
└─────────────────────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ PHASE 5: DESIGN STRATEGY GENERATION                                 [50%]  │
├─────────────────────────────────────────────────────────────────────────────┤
│ Claude synthesizes research into actionable strategy:                       │
│ • Core scaffold classes appropriate for target                              │
│ • Required functional groups (with SMARTS patterns)                         │
│ • Modification sites for novelty                                            │
│ • Property optimization targets                                             │
│ • Design rules specific to this target                                      │
└─────────────────────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ PHASE 6: SCORING CRITERIA GENERATION                                [55%]  │
├─────────────────────────────────────────────────────────────────────────────┤
│ Claude creates multi-objective scoring function:                            │
│                                                                             │
│ ┌─────────────────────────────────────────────────────────────────────────┐│
│ │  Binding Score    │████████████████████│  35%  (most important)        ││
│ │  Selectivity      │██████████████      │  25%                          ││
│ │  ADMET Safety     │████████████        │  20%                          ││
│ │  Novelty          │████████████        │  20%                          ││
│ └─────────────────────────────────────────────────────────────────────────┘│
│                                                                             │
│ Plus target-specific bonus/penalty patterns                                 │
└─────────────────────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ PHASE 7: ITERATIVE DESIGN LOOP                                   [55-90%]  │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   ┌──────────────────────────────────────────────────────────────────┐     │
│   │                    ITERATION N (up to 10)                         │     │
│   └──────────────────────────────────────────────────────────────────┘     │
│                              │                                              │
│   ┌──────────────────────────▼──────────────────────────┐                  │
│   │ A. GENERATE CANDIDATES                               │                  │
│   │    Claude generates 5-10 SMILES per iteration        │                  │
│   │    Following strategy + mechanistic constraints      │                  │
│   └──────────────────────────┬──────────────────────────┘                  │
│                              │                                              │
│   ┌──────────────────────────▼──────────────────────────┐                  │
│   │ B. VALIDATE STRUCTURES (RDKit)                       │                  │
│   │    • Parse SMILES → Mol object                       │                  │
│   │    • Check valid chemistry                           │                  │
│   │    • Apply PAINS filter                              │                  │
│   │    • For covalent: Check topology (MechanisticAnalyzer) │              │
│   └──────────────────────────┬──────────────────────────┘                  │
│                              │                                              │
│   ┌──────────────────────────▼──────────────────────────┐                  │
│   │ C. PREDICT PROPERTIES                                │                  │
│   │    MW, LogP, TPSA, HBD, HBA, QED, SA Score          │                  │
│   │    Lipinski violations, solubility                   │                  │
│   └──────────────────────────┬──────────────────────────┘                  │
│                              │                                              │
│   ┌──────────────────────────▼──────────────────────────┐                  │
│   │ D. CHECK NOVELTY (Parallel)                          │                  │
│   │    • Morgan fingerprint similarity                   │                  │
│   │    • PubChem similarity search                       │                  │
│   │    • ChEMBL similarity search                        │                  │
│   │    • Tanimoto < 0.85 = novel                        │                  │
│   └──────────────────────────┬──────────────────────────┘                  │
│                              │                                              │
│   ┌──────────────────────────▼──────────────────────────┐                  │
│   │ E. PREDICT ADMET                                     │                  │
│   │    Absorption, Distribution, Metabolism,             │                  │
│   │    Excretion, Toxicity → Safety Score               │                  │
│   └──────────────────────────┬──────────────────────────┘                  │
│                              │                                              │
│   ┌──────────────────────────▼──────────────────────────┐                  │
│   │ F. SCORE CANDIDATES                                  │                  │
│   │    overall = 0.35*binding + 0.25*selectivity        │                  │
│   │            + 0.20*admet + 0.20*novelty              │                  │
│   └──────────────────────────┬──────────────────────────┘                  │
│                              │                                              │
│   ┌──────────────────────────▼──────────────────────────┐                  │
│   │ G. CONVERGENCE CHECK                                 │                  │
│   │    If best_score > 0.85 OR no improvement for 3 iter │                  │
│   │    → BREAK                                           │                  │
│   │    Else → Next iteration with refined strategy       │                  │
│   └──────────────────────────────────────────────────────┘                  │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ PHASE 8: FINAL SELECTION                                            [92%]  │
├─────────────────────────────────────────────────────────────────────────────┤
│ Claude reviews all candidates:                                              │
│ • Sort by overall_score                                                     │
│ • Ensure diversity (not 10 similar molecules)                               │
│ • Rank top 10-20 candidates                                                 │
│ • Generate rationale for each                                               │
└─────────────────────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│ OUTPUT: VALIDATED SMILES WITH ANALYSIS                             [100%]  │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  #1  CC1CCN(C(=O)Nc2ccc(C(F)(F)F)cc2)CC1c1ccccc1                          │
│      ├── Score: 87%                                                         │
│      ├── MW: 378.4  LogP: 3.2  QED: 0.72                                   │
│      ├── Novelty: 0.73 (patentable)                                        │
│      ├── ADMET: Safe (no hERG, no hepatotox)                               │
│      └── Rationale: "Piperidine-biaryl scaffold with CF3..."              │
│                                                                             │
│  #2  ...                                                                    │
│  #3  ...                                                                    │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Key Components

### 1. Entry Points

| Endpoint | Purpose |
|----------|---------|
| `POST /api/v1/design/advanced/stream` | Main streaming endpoint (SSE) |
| `POST /api/v1/design/advanced` | Non-streaming advanced design |
| `POST /api/v1/design` | Quick design mode |

### 2. Core Services

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                            SERVICE LAYER                                     │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌─────────────────────┐    ┌─────────────────────┐                        │
│  │  AgenticDrugDesigner│───▶│    LiveNarrator     │                        │
│  │  (Main Engine)      │    │  (Real-time Updates)│                        │
│  └──────────┬──────────┘    └─────────────────────┘                        │
│             │                                                               │
│  ┌──────────┴──────────────────────────────────────────────────────┐       │
│  │                                                                  │       │
│  ▼                    ▼                    ▼                    ▼   │       │
│  ┌──────────────┐ ┌──────────────┐ ┌──────────────┐ ┌──────────────┐│       │
│  │ ChEMBL       │ │ Mechanistic  │ │ Property     │ │ Novelty      ││       │
│  │ Client       │ │ Analyzer     │ │ Predictor    │ │ Checker      ││       │
│  └──────────────┘ └──────────────┘ └──────────────┘ └──────────────┘│       │
│         │                │                │                │        │       │
│         ▼                ▼                ▼                ▼        │       │
│  ┌──────────────┐ ┌──────────────┐ ┌──────────────┐ ┌──────────────┐│       │
│  │ ChEMBL API   │ │ RDKit        │ │ RDKit        │ │ PubChem API  ││       │
│  │ (Bioactivity)│ │ (Validation) │ │ (Properties) │ │ (Similarity) ││       │
│  └──────────────┘ └──────────────┘ └──────────────┘ └──────────────┘│       │
│                                                                      │       │
└──────────────────────────────────────────────────────────────────────┘       │
                                                                               │
└─────────────────────────────────────────────────────────────────────────────┘
```

### 3. Data Structures

```python
DesignSession
├── request: str                    # User's prompt
├── target_knowledge: TargetKnowledge
│   ├── target_name: str            # "NOP/MOP Receptors"
│   ├── target_type: str            # "GPCR"
│   ├── mechanism: str              # "Dual agonist"
│   ├── is_covalent: bool           # False
│   ├── essential_features: []      # Required structural elements
│   ├── avoid_features: []          # What NOT to include
│   ├── reference_drugs: []         # Known actives from ChEMBL
│   └── scoring_criteria: {}        # Multi-objective weights
├── candidates: List[DesignCandidate]
│   └── DesignCandidate
│       ├── smiles: str             # "CC1CCN(C(=O)Nc2ccc..."
│       ├── name: str               # "Spiropiperidine-CF3-1"
│       ├── rationale: str          # Why this was designed
│       ├── molecular_weight: float # 378.4
│       ├── logp: float             # 3.2
│       ├── qed: float              # 0.72
│       ├── binding_score: float    # 0.85
│       ├── novelty_score: float    # 0.73
│       ├── overall_score: float    # 0.87
│       └── patentable: bool        # True
├── iterations: int                 # 5
├── strategy: str                   # Design approach used
└── design_log: List[str]           # Phase-by-phase log
```

### 4. External APIs Used

| API | Purpose | Data Retrieved |
|-----|---------|----------------|
| **Claude Opus 4.5** | LLM reasoning | Target analysis, SMILES generation, strategy |
| **Claude Haiku** | Fast narration | Real-time progress summaries |
| **ChEMBL** | Bioactivity data | IC50, Ki values, known active SMILES |
| **PubChem** | Compound lookup | Reference SMILES, similarity search |
| **RDKit** | Cheminformatics | Property calculation, validation, fingerprints |

---

## Covalent Inhibitor Handling

For covalent inhibitors (FAAH, EGFR, etc.), special topology validation ensures correct design:

```
WRONG TOPOLOGY (rejected):                CORRECT TOPOLOGY (accepted):

  hexyl─phenyl─NH─C(=O)─N─piperidine       pyridazinyl─N─C(=O)─N─piperidine─biaryl─CF3
  ^^^^^^^^^^^^                ^^^^^        ^^^^^^^^^^              ^^^^^^^^^^^^^^^
  LIPOPHILIC    WARHEAD    BASIC          LEAVING      WARHEAD     LIPOPHILIC
  (LEAVES!)                (STAYS)        GROUP                    (STAYS!)

  After reaction:                         After reaction:
  Only bare piperidine bound              Piperidine-biaryl-CF3 bound
  = NO BINDING = INACTIVE                 = FILLS POCKET = ACTIVE
```

The `MechanisticAnalyzer` validates:
- Warhead type detection (urea, carbamate, acrylamide, etc.)
- pKa-based leaving group identification
- Binding group placement verification

---

## Real-Time Streaming

```
Frontend ◄──────── SSE ──────── Backend
                    │
    data: {"phase": "research", "message": "Identified NOP/MOP receptors..."}
    data: {"phase": "chembl", "message": "Found 47 active compounds..."}
    data: {"phase": "design", "message": "Generated 10 candidates, 7 passed..."}
    data: {"phase": "completed", "data": {candidates: [...]}}
```

The `LiveNarrator` service uses Claude Haiku to generate human-readable summaries every ~8 seconds, transforming raw work content into insights like:
- "Found cebranopadol (Ki=0.1nM) as lead reference"
- "Generated spirocyclic scaffold with CF3 substituent"
- "Validation: 7/10 passed topology check"

---

## File Structure

```
drug-design-platform/
├── backend/
│   ├── api/
│   │   └── main.py              # FastAPI endpoints (streaming, advanced)
│   ├── services/
│   │   ├── agentic_designer.py  # Core design engine (8 phases)
│   │   ├── live_narrator.py     # Real-time Haiku summaries
│   │   ├── chembl_client.py     # ChEMBL bioactivity queries
│   │   ├── mechanistic_analyzer.py  # Covalent inhibitor topology
│   │   ├── property_predictor.py    # RDKit property calculation
│   │   ├── novelty_checker.py   # PubChem/ChEMBL similarity
│   │   └── admet_predictor.py   # ADMET safety predictions
│   └── models/
│       └── __init__.py          # Data structures
└── frontend/
    └── app.html                 # UI with live feed
```

---

## Summary

**Input**: Natural language prompt describing desired drug properties

**Process**:
1. Parse request → identify targets, compounds, constraints
2. Research target biology via Claude
3. Fetch real bioactivity data from ChEMBL
4. Analyze mechanism (covalent vs reversible)
5. Generate target-specific design strategy
6. Iteratively generate & validate SMILES (5-10 iterations)
7. Score on binding, selectivity, ADMET, novelty
8. Select top candidates with diversity

**Output**: 10-20 validated SMILES with:
- Full property profiles (MW, LogP, QED, etc.)
- Novelty assessment (Tanimoto similarity)
- ADMET predictions (safety score)
- Design rationale
- Real-time progress throughout
