# Drug Design Platform

An automated platform for generating novel, patentable analogues of drug compounds. This platform combines computational chemistry tools (RDKit), AI-powered design suggestions (Claude), molecular docking (AutoDock Vina), and structure prediction (AlphaFold 3) to deliver chemist-ready molecule designs.

## Features

- **Automated Analogue Generation**: Generate structurally diverse analogues using multiple modification strategies
- **AI-Powered Design**: Claude AI suggests intelligent modifications based on medicinal chemistry principles
- **Property Prediction**: Predict ADMET, drug-likeness, and pharmacokinetic properties
- **Novelty Checking**: Search PubChem and ChEMBL databases to verify compound novelty
- **Molecular Docking**: Predict binding affinity using AutoDock Vina
- **Structure Prediction**: Integration with AlphaFold 3 for protein-ligand complex prediction
- **Synthesis Guidance**: AI-generated synthesis routes and feasibility assessment
- **Web Interface**: Interactive frontend for compound input and results visualization

## Quick Start

### Prerequisites

- Python 3.9+
- pip
- (Optional) AutoDock Vina for docking
- (Optional) Docker for AlphaFold 3

### Installation

```bash
# Clone or navigate to the project
cd drug-design-platform

# Run setup script
chmod +x setup.sh
./setup.sh

# Or manual setup:
python3 -m venv venv
source venv/bin/activate
cd backend
pip install -r requirements.txt
```

### Configuration

Create a `.env` file in the `backend` directory:

```env
ANTHROPIC_API_KEY=your_api_key_here
```

### Running the Platform

```bash
# Start backend server
cd backend
python main.py

# In another terminal, serve the frontend
python -m http.server 3000 --directory frontend
```

- Backend API: http://localhost:8000
- API Docs: http://localhost:8000/docs
- Frontend: http://localhost:3000

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/generate` | POST | Generate novel analogues |
| `/api/demo/kk103` | GET | Get KK103 demo data |
| `/api/demo/kk103/full` | POST | Run full pipeline on KK103 |
| `/api/predict/properties` | POST | Predict compound properties |
| `/api/check/novelty` | POST | Check compound novelty |
| `/api/dock` | POST | Dock compound to receptor |
| `/api/visualize/2d/{smiles}` | GET | Get 2D structure image |
| `/api/visualize/3d/{smiles}` | GET | Get 3D conformer |

## Acceptance Criteria (KK103 Test Case)

The platform is validated against KK103 (N-pivaloyl-Leu-Enkephalin):

| Criterion | Threshold |
|-----------|-----------|
| Analogues generated | ≥ 5 |
| Tanimoto similarity to KK103 | < 0.85 |
| Predicted DOR affinity | > 50% of KK103 |
| Drug-likeness (QED) | > 0.3 |
| Synthetic accessibility | ≤ 6 (1-10 scale) |
| Novelty confidence | > 70% |

## Architecture

```
drug-design-platform/
├── backend/
│   ├── main.py              # FastAPI application
│   ├── pipeline.py          # Main orchestrator
│   ├── models.py            # Pydantic models
│   ├── config.py            # Configuration
│   ├── services/
│   │   ├── analogue_generator.py   # Structural modifications
│   │   ├── property_predictor.py   # ADMET & properties
│   │   ├── novelty_checker.py      # Database searches
│   │   ├── docking_service.py      # AutoDock Vina
│   │   ├── alphafold_service.py    # AlphaFold 3
│   │   └── claude_service.py       # Claude AI integration
│   └── tests/
│       └── test_acceptance.py      # Acceptance tests
├── frontend/
│   └── index.html           # Web interface
├── data/                    # Input data directory
├── output/                  # Generated outputs
└── docs/
    └── ACCEPTANCE_CRITERIA.md
```

## Modification Strategies

### N-Terminal Modifications
- Cyclopropylcarbonyl, isobutyryl, benzoyl, etc.
- Replaces pivaloyl group to explore SAR

### Amino Acid Substitutions
- D-amino acids for proteolytic stability
- Non-natural amino acids (Sar, Aib, Nle, etc.)
- Halogenated variants (4-F-Phe, etc.)

### C-Terminal Modifications
- Amide instead of carboxylic acid
- Ester prodrugs
- Hydroxamic acids

### Backbone Modifications
- N-methylation at specific positions
- Reduced amide bonds
- Beta-amino acid insertions

## Pre-built KK103 Analogues

| Name | Modification | Rationale |
|------|--------------|-----------|
| KK103-A1 | N-cyclopropylcarbonyl | Reduced MW, maintained steric bulk |
| KK103-A2 | D-Tyr1 | Proteolytic resistance |
| KK103-A3 | Sar2 (N-methylglycine) | Improved stability |
| KK103-A4 | C-terminal amide | Enhanced permeability |
| KK103-A5 | Nle5 (Norleucine) | Linear chain variant |
| KK103-A6 | 4-Fluoro-Phe4 | Metabolic stability |
| KK103-A7 | D-Tyr1 + Sar2 | Dual protection |
| KK103-A8 | N-isobutyryl | Smaller N-cap |

## Running Tests

```bash
cd backend
python -m pytest tests/test_acceptance.py -v
```

## Technology Stack

- **Backend**: FastAPI, Python 3.9+
- **Chemistry**: RDKit, Open Babel
- **AI**: Claude (Anthropic API)
- **Docking**: AutoDock Vina
- **Structure Prediction**: AlphaFold 3
- **Database Search**: PubChem, ChEMBL
- **Frontend**: HTML5, TailwindCSS, 3Dmol.js

## References

1. KK103 Research: "An Effective and Safe Enkephalin Analog for Antinociception" (PMC8308721)
2. RDKit: https://www.rdkit.org/
3. AlphaFold 3: https://github.com/google-deepmind/alphafold3
4. AutoDock Vina: https://vina.scripps.edu/

## License

For research and educational purposes only.

## Contributing

Contributions welcome! Please read the acceptance criteria and ensure all tests pass.
