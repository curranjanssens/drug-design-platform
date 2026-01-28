# Drug Design Platform: First Principles Philosophy

## THE GOAL
**>50% chance of each designed drug hitting its desired endpoint**

Current drug discovery: ~10% success rate. We aim for 5x improvement through intelligent LLM orchestration.

---

## WHY DRUGS FAIL (in order of frequency)

1. **Wrong target hypothesis** - The biology was wrong
2. **Poor binding affinity** - Doesn't bind strongly enough
3. **Off-target toxicity** - Binds things it shouldn't
4. **Poor pharmacokinetics** - Doesn't get to target
5. **Can't be synthesized** - Too complex to make

## WHAT WE MUST DO FOR EACH

| Failure Mode | Solution | Tool Required |
|--------------|----------|---------------|
| Wrong target | Deep literature/database research | PubMed, ChEMBL, UniProt |
| Poor binding | Structure-based design + docking | AlphaFold, AutoDock/GNINA |
| Off-target toxicity | Multi-target docking screen | Dock to off-target panel |
| Poor PK | ADMET prediction + optimization | ADMETlab, pkCSM, DeepChem |
| Synthesis issues | Retrosynthesis planning | ASKCOS, IBM RXN |

---

## THE IDEAL ARCHITECTURE

```
USER PROMPT
     ↓
┌─────────────────────────────────────────────────────────────┐
│  PHASE 1: DEEP TARGET RESEARCH                              │
│  ├── PubChem: Compound lookup (SMILES, properties)          │
│  ├── UniProt: Protein sequence, function, diseases          │
│  ├── PDB/AlphaFold: 3D protein structure                    │
│  ├── ChEMBL: Known binders with IC50/Ki values              │
│  ├── BindingDB: Additional binding data                     │
│  └── LLM: Synthesize all data, understand mechanism         │
└─────────────────────────────────────────────────────────────┘
     ↓
┌─────────────────────────────────────────────────────────────┐
│  PHASE 2: BINDING SITE ANALYSIS                             │
│  ├── Identify binding pocket from structure                 │
│  ├── Extract key interacting residues                       │
│  ├── Generate pharmacophore model                           │
│  ├── Analyze known ligand binding poses                     │
│  └── LLM: Define structural requirements                    │
└─────────────────────────────────────────────────────────────┘
     ↓
┌─────────────────────────────────────────────────────────────┐
│  PHASE 3: CANDIDATE GENERATION                              │
│  ├── LLM generates based on knowledge                       │
│  ├── Generative model (REINVENT/DrugGPT) for exploration    │
│  ├── Fragment-based growth from known binders               │
│  ├── Scaffold hopping for novelty                           │
│  └── Each candidate validated for essential features        │
└─────────────────────────────────────────────────────────────┘
     ↓
┌─────────────────────────────────────────────────────────────┐
│  PHASE 4: MULTI-OBJECTIVE SCORING                           │
│  ├── DOCKING: Binding affinity score (kcal/mol)             │
│  ├── SELECTIVITY: Dock to off-target panel                  │
│  ├── ADMET: Absorption, metabolism, toxicity                │
│  ├── SYNTHESIS: Retrosynthesis feasibility                  │
│  ├── NOVELTY: Tanimoto distance from known drugs            │
│  └── DRUGLIKENESS: QED, Lipinski, Veber                     │
└─────────────────────────────────────────────────────────────┘
     ↓
┌─────────────────────────────────────────────────────────────┐
│  PHASE 5: VALIDATION & ITERATION                            │
│  ├── Verify structures are chemically valid                 │
│  ├── Verify essential features are present                  │
│  ├── LLM critiques candidates against requirements          │
│  ├── Reject failures, learn why they failed                 │
│  └── Iterate with improved strategy                         │
└─────────────────────────────────────────────────────────────┘
     ↓
┌─────────────────────────────────────────────────────────────┐
│  PHASE 6: FINAL DELIVERABLES                                │
│  ├── Top 10-20 candidates ranked by overall score           │
│  ├── Full synthesis routes for each                         │
│  ├── Risk assessment and mitigation strategies              │
│  ├── Recommended assays for validation                      │
│  └── Patent landscape analysis                              │
└─────────────────────────────────────────────────────────────┘
```

---

## KEY TOOLS TO INTEGRATE

### Already Have
- [x] PubChem lookup
- [x] RDKit properties
- [x] Basic ADMET (fallback)
- [x] LLM orchestration
- [x] Validation loop

### Need to Add (Priority Order)

1. **ChEMBL Integration** (HIGH)
   - Get known binders with actual Ki/IC50 values
   - Learn from real SAR data
   - API: https://www.ebi.ac.uk/chembl/api/data/

2. **Protein Structure** (HIGH)
   - AlphaFold DB for predicted structures
   - PDB for experimental structures
   - Needed for docking

3. **Docking** (HIGH)
   - GNINA (ML-enhanced, free)
   - Or AutoDock Vina
   - Actual binding affinity prediction

4. **Better ADMET** (MEDIUM)
   - ADMETlab 2.0 API
   - Or SwissADME
   - More accurate toxicity prediction

5. **Off-Target Panel** (MEDIUM)
   - Dock to common off-targets (hERG, CYPs, etc.)
   - Predict selectivity issues early

6. **Synthesis Planning** (MEDIUM)
   - ASKCOS from MIT
   - Score synthetic accessibility
   - Full retrosynthesis routes

---

## LLM ORCHESTRATION PRINCIPLES

1. **LLM decides what tools to use** based on available data
2. **LLM has access to all information** gathered by tools
3. **LLM validates its own outputs** before accepting
4. **LLM learns from failures** and adjusts strategy
5. **LLM explains reasoning** for every decision
6. **No hardcoded drug classes** - all knowledge from research
7. **Iterate until quality threshold met** - no arbitrary limits

---

## SUCCESS METRICS

For each candidate, we need:

| Metric | Target | How to Measure |
|--------|--------|----------------|
| Binding affinity | < -8 kcal/mol | Docking score |
| Selectivity | >100x vs off-targets | Multi-target docking |
| Drug-likeness | QED > 0.5 | RDKit |
| ADMET | No red flags | Prediction tools |
| Synthesizable | < 6 steps | Retrosynthesis |
| Novel | Tanimoto < 0.7 | Fingerprint comparison |
| Essential features | All present | Validation loop |

**Overall target: At least 50% of final candidates should meet ALL criteria**

---

## REFERENCE THIS DOCUMENT

When making design decisions, always ask:
1. Does this increase our chance of hitting >50% success?
2. Are we using the best available tool for this task?
3. Is the LLM being given enough information to make good decisions?
4. Are we validating outputs before accepting them?
5. Are we learning from failures?

---

*Last updated: 2026-01-10*
