# Drug Design Platform - Acceptance Criteria

## Overview

This document defines the acceptance criteria for the automated drug design platform. The platform must be capable of accepting natural language requests and producing chemist-ready molecule designs with complete synthesis routes, novelty assessments, and patentability analysis.

---

## Primary Acceptance Tests

### Test 1: Novel KK103 Analogue Design

**Objective:** Design a novel, patentable analogue of KK103 (a precursor of leucine-enkephalin).

**Input Request:**
```
"Design a novel, patentable analogue of KK103 with improved metabolic stability"
```

**Reference Compound:**
- Name: KK103
- Structure: Leucine-enkephalin precursor
- Target: Opioid receptors (mu, delta)
- SMILES: `CC(C)C[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)CNC(=O)CNC(=O)[C@@H](N)Cc2ccc(O)cc2)C(=O)O`

**Acceptance Criteria:**

| Criterion | Requirement | Metric |
|-----------|-------------|--------|
| **Novelty** | Tanimoto similarity < 0.85 to reference | Pass/Fail |
| **Patentability** | No exact patent matches | Boolean |
| **Drug-likeness** | QED > 0.4 | Numeric |
| **Lipinski Compliance** | ≤ 1 violation | Count |
| **ADMET Safety** | Safety score > 0.6 | Numeric |
| **Metabolic Stability** | No obvious metabolic hotspots | Assessment |
| **Synthesis Feasibility** | ≤ 8 steps, ≥ 5% overall yield | Route provided |
| **Starting Materials** | All purchasable | Boolean |

**Deliverables:**
1. Top 5 ranked analogue structures (SMILES, InChI, MOL files)
2. 2D and 3D structure visualizations
3. Complete property profiles (MW, LogP, TPSA, HBD/HBA)
4. Full ADMET predictions
5. Novelty assessment with nearest known compounds
6. Patent landscape analysis
7. Step-by-step synthesis routes with conditions
8. Starting material list with vendors and pricing
9. Comprehensive report suitable for medicinal chemistry team

---

### Test 2: 5HT1A-5HT2A Dual Agonist Bivalent Ligand

**Objective:** Design a novel dual agonist bivalent ligand for 5HT1A-5HT2A receptors.

**Input Request:**
```
"Design a novel dual 5HT1A-5HT2A agonist bivalent ligand with CNS penetration"
```

**Design Parameters:**
- Target 1: 5HT1A receptor (Gi-coupled)
- Target 2: 5HT2A receptor (Gq-coupled)
- Design type: Bivalent ligand (two pharmacophores + linker)
- Constraint: Blood-brain barrier penetration required

**Acceptance Criteria:**

| Criterion | Requirement | Metric |
|-----------|-------------|--------|
| **Structure** | Two distinct pharmacophores + linker | Structural analysis |
| **5HT1A Pharmacophore** | Contains validated agonist scaffold | Substructure match |
| **5HT2A Pharmacophore** | Contains validated agonist scaffold | Substructure match |
| **Linker Length** | 12-24 atoms | Count |
| **Novelty** | No exact prior art | Patent search |
| **CNS Penetration** | BBB+ prediction | ADMET model |
| **MW** | < 800 Da | Numeric |
| **LogP** | 2-5 (CNS optimal range) | Numeric |
| **Safety** | hERG negative, Ames negative | ADMET predictions |
| **Synthesis** | Route identified | Pass/Fail |

**Deliverables:**
1. Top 3 bivalent ligand designs
2. Pharmacophore identification and rationale
3. Linker chemistry explanation
4. Full ADMET profile with CNS-specific metrics
5. Synthesis route for the bivalent molecule
6. Novelty and patentability assessment
7. Design rationale document

---

## Platform Functional Requirements

### Input Processing
- Accept natural language prompts (Claude-powered parsing)
- Accept SMILES input (validated with RDKit)
- Accept compound names (lookup against known database)
- Accept target specifications (PDB ID, UniProt, or name)
- Accept design constraints (MW limits, required substructures)

### Molecule Generation
- Generate novel analogues (RDKit-based modifications)
- Scaffold hopping (bioisostere replacements)
- Bivalent ligand design (pharmacophore + linker assembly)
- Multi-objective optimization (balance affinity, ADMET, novelty)
- Generate ≥ 10 candidates per request

### Property Prediction
- Molecular weight, LogP, TPSA, H-bond donors/acceptors
- Rotatable bonds, ring count
- QED (drug-likeness)
- Synthetic accessibility score

### ADMET Prediction
- Absorption (Caco-2, HIA)
- Distribution (BBB, VDss, PPB)
- Metabolism (CYP inhibition/substrate)
- Excretion (clearance, half-life)
- Toxicity (hERG, Ames, hepatotoxicity)

### Novelty Assessment
- Similarity to known compounds (Tanimoto)
- Scaffold novelty check
- Patent search integration

### Retrosynthesis Planning
- Multi-step route generation
- Starting material identification
- Reaction conditions
- Yield estimates and cost estimation

---

## Chemist-Ready Deliverable Checklist

For each designed molecule:

### Molecule Specification
- Canonical SMILES
- InChI and InChIKey
- Molecular formula
- 2D structure image (SVG/PNG)
- 3D structure file (SDF/MOL2)
- Key descriptors

### ADMET Profile
- Full predictions table
- Safety score
- Flagged warnings

### Novelty Report
- Tanimoto similarity to nearest known
- Scaffold novelty assessment
- Patent landscape summary
- FTO risk level

### Synthesis Plan
- Complete step-by-step route
- Reagents and conditions
- Expected yields
- Starting materials with vendors/pricing
- Route confidence score

---

## API Endpoints

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/api/v1/design` | POST | Submit design job |
| `/api/v1/design/{job_id}` | GET | Get job status/results |
| `/api/v1/design/{job_id}/report` | GET | Get comprehensive report |
| `/api/v1/analyze` | POST | Analyze single molecule |
| `/api/v1/analyze/admet` | POST | Get ADMET predictions |
| `/api/v1/retrosynthesis` | POST | Plan synthesis |
| `/api/v1/patent/novelty` | POST | Check novelty |
| `/api/v1/patent/fto` | POST | Assess FTO |
| `/health` | GET | Health check |

---

## Performance Requirements

| Metric | Target |
|--------|--------|
| Design job completion | < 5 minutes |
| Single molecule analysis | < 10 seconds |
| Retrosynthesis planning | < 60 seconds |
| Patent search | < 30 seconds |
| API response time | < 500ms |

---

## Sign-off Criteria

The platform is production-ready when:

1. All functional requirements implemented
2. Both primary acceptance tests pass
3. All API endpoints operational
4. Performance targets met
5. Documentation complete
6. Error handling robust
