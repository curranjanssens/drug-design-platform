# Drug Design Platform - Complete Remodel Plan

## THE PROBLEM

Current system generates molecules via pattern matching:
- Identifies components (urea, piperidine, lipophilic tail) ✓
- Assembles them WITHOUT understanding mechanism ✗
- Result: Lipophilic tail on LEAVING GROUP = inactive compound

**Example failure (FAAH inhibitor):**
```
AI generated: hexyloxy-phenyl-NH-C(=O)-N-piperidine
                ^^^^^^^^^^^^          ^^^^^^^^^
                This LEAVES!          This STAYS (empty pocket)

Should be:   pyridazinyl-N=C(=O)-N-piperidine-biaryl-CF₃
             ^^^^^^^^^           ^^^^^^^^^^^^^^^^^^^^^^^^
             This LEAVES         This STAYS (fills pocket)
```

## ROOT CAUSE

LLMs pattern-match from reference structures but don't reason about:
1. **Covalent mechanism**: Which bond breaks, what stays vs leaves
2. **pKa/leaving group ability**: Which nitrogen is better leaving group
3. **Binding topology**: Where must lipophilic groups be AFTER reaction
4. **3D binding pocket**: What shape must the REMAINING fragment have

## THE SOLUTION: MECHANISM-FIRST DESIGN

### Architecture Overview

```
USER PROMPT
     ↓
┌─────────────────────────────────────────────────────────────┐
│  PHASE 1: MECHANISTIC ANALYSIS                              │
│  ├── Identify inhibitor type (reversible vs covalent)       │
│  ├── For covalent: identify warhead, leaving group          │
│  ├── Determine what STAYS bound to enzyme                   │
│  ├── Determine what LEAVES                                  │
│  └── Binding groups MUST be on the STAYING portion          │
└─────────────────────────────────────────────────────────────┘
     ↓
┌─────────────────────────────────────────────────────────────┐
│  PHASE 2: TEMPLATE EXTRACTION                               │
│  ├── Get protein structure (AlphaFold/PDB)                  │
│  ├── Dock reference compounds                               │
│  ├── Identify the COVALENT ADDUCT structure                 │
│  ├── Extract what remains bound after reaction              │
│  └── This becomes the TEMPLATE for new designs              │
└─────────────────────────────────────────────────────────────┘
     ↓
┌─────────────────────────────────────────────────────────────┐
│  PHASE 3: CONSTRAINED GENERATION                            │
│  ├── Use validated scaffold (not free-form SMILES)          │
│  ├── Enumerate substitutions at ALLOWED positions           │
│  ├── CONSTRAINT: Binding groups on staying portion          │
│  ├── CONSTRAINT: Leaving group is appropriate pKa           │
│  └── Generate 3D conformers                                 │
└─────────────────────────────────────────────────────────────┘
     ↓
┌─────────────────────────────────────────────────────────────┐
│  PHASE 4: MECHANISTIC VALIDATION                            │
│  ├── Simulate the covalent reaction                         │
│  ├── Check: Does the ADDUCT occupy the binding pocket?      │
│  ├── Check: Is leaving group pKa appropriate?               │
│  ├── Dock the ADDUCT (not the pro-drug)                     │
│  └── Reject if mechanism is wrong                           │
└─────────────────────────────────────────────────────────────┘
     ↓
┌─────────────────────────────────────────────────────────────┐
│  PHASE 5: MULTI-OBJECTIVE SCORING                           │
│  ├── Docking score of COVALENT ADDUCT                       │
│  ├── Kinact/KI prediction (covalent efficiency)             │
│  ├── ADMET properties                                       │
│  ├── Selectivity (dock to off-targets)                      │
│  └── Synthetic accessibility                                │
└─────────────────────────────────────────────────────────────┘
```

## KEY INSIGHT: DOCK THE ADDUCT, NOT THE PRO-DRUG

For covalent inhibitors, what matters is what REMAINS bound after the reaction:

```
Pro-drug:     LEAVING_GROUP-C(=O)-STAYING_PORTION
                    ↓ [enzyme attacks]
Adduct:       ENZYME-O-C(=O)-STAYING_PORTION  +  LEAVING_GROUP (gone)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
              THIS is what you dock and optimize!
```

## IMPLEMENTATION PLAN

### Step 1: Mechanistic Analyzer (NEW)
```python
class MechanisticAnalyzer:
    def analyze_inhibitor_type(self, target, mechanism):
        """Determine if reversible or covalent, warhead type, etc."""

    def identify_leaving_group(self, smiles, warhead_type):
        """For covalent: which portion leaves after enzyme attack?"""

    def extract_covalent_adduct(self, smiles, warhead_type):
        """Return SMILES of what remains bound to enzyme"""

    def validate_topology(self, smiles, target):
        """Check: Are binding groups on staying portion?"""
```

### Step 2: Template-Based Generator (REPLACE free-form)
```python
class TemplateGenerator:
    def get_validated_scaffold(self, target, mechanism):
        """Return scaffold with marked modification sites"""
        # For FAAH: piperidine-CF3-biaryl core with leaving group attachment point

    def enumerate_modifications(self, scaffold, positions):
        """Generate variants at allowed positions only"""

    def ensure_correct_topology(self, candidate, mechanism):
        """Verify binding groups on staying portion"""
```

### Step 3: Covalent Docking (NEW)
```python
class CovalentDocker:
    def dock_adduct(self, adduct_smiles, protein_structure):
        """Dock the COVALENT ADDUCT, not the pro-drug"""

    def score_binding(self, pose):
        """Score how well adduct occupies binding pocket"""
```

### Step 4: Mechanism Validator (NEW)
```python
class MechanismValidator:
    def validate_leaving_group_pka(self, leaving_group):
        """Check pKa is appropriate (8-12 for most warheads)"""

    def validate_adduct_binding(self, adduct, pocket):
        """Check adduct fills the binding pocket properly"""

    def reject_inverted_topology(self, candidate):
        """Catch the FAAH error: lipophilic on leaving group"""
```

## SPECIFIC FIX FOR FAAH

For FAAH urea inhibitors:
1. Template: `[LEAVING]-N=C(=O)-N(piperidine)-[BIARYL]-[LIPOPHILIC]`
2. Leaving group: Must be good leaving group (pyridazinyl, pKa ~8)
3. Staying portion: Piperidine-biaryl-lipophilic tail
4. Dock the carbamate adduct: `Enzyme-O-C(=O)-N(piperidine)-biaryl-CF3`

**NOT**: `[LIPOPHILIC]-phenyl-NH-C(=O)-N-piperidine`
(This puts lipophilic on the leaving group!)

## SUCCESS CRITERIA

1. **Mechanistic correctness**: 100% of covalent inhibitors have correct topology
2. **Binding validation**: All candidates docked as covalent adducts
3. **Activity prediction**: Predicted Ki within 10x of actual
4. **Overall hit rate**: >50% of designs show measurable activity

## TOOLS NEEDED

1. **Structure tools**: RDKit, OpenBabel
2. **Docking**: GNINA or AutoDock for covalent docking
3. **pKa prediction**: ChemAxon or RDKit
4. **3D conformers**: RDKit ETKDG
5. **Protein structures**: AlphaFold DB, PDB

## TIMELINE

This is a fundamental rewrite, not a patch. The current "generate SMILES from features" approach cannot work for covalent inhibitors.

Priority order:
1. Mechanistic analyzer (catch inverted topology)
2. Template extraction from reference compounds
3. Constrained generation (not free-form)
4. Covalent adduct docking
5. Multi-objective optimization

---

*The core insight: LLMs can identify components but cannot reason about reaction mechanisms. We must encode mechanistic constraints explicitly.*
