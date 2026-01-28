#!/bin/bash

# Drug Design Platform Setup Script
# This script sets up the complete environment for the automated drug design platform

set -e

echo "========================================"
echo "Drug Design Platform Setup"
echo "========================================"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check Python version
echo -e "\n${YELLOW}Checking Python version...${NC}"
python3 --version || { echo -e "${RED}Python 3 is required${NC}"; exit 1; }

# Create virtual environment
echo -e "\n${YELLOW}Creating virtual environment...${NC}"
if [ ! -d "venv" ]; then
    python3 -m venv venv
    echo -e "${GREEN}Virtual environment created${NC}"
else
    echo -e "${GREEN}Virtual environment already exists${NC}"
fi

# Activate virtual environment
source venv/bin/activate

# Upgrade pip
echo -e "\n${YELLOW}Upgrading pip...${NC}"
pip install --upgrade pip

# Install backend dependencies
echo -e "\n${YELLOW}Installing backend dependencies...${NC}"
cd backend
pip install -r requirements.txt

# Check for optional dependencies
echo -e "\n${YELLOW}Checking optional dependencies...${NC}"

# Check for Open Babel (for molecular file conversion)
if command -v obabel &> /dev/null; then
    echo -e "${GREEN}Open Babel found${NC}"
else
    echo -e "${YELLOW}Open Babel not found - some features may be limited${NC}"
    echo "Install with: brew install open-babel (macOS) or apt-get install openbabel (Linux)"
fi

# Check for AutoDock Vina
if command -v vina &> /dev/null; then
    echo -e "${GREEN}AutoDock Vina found${NC}"
else
    echo -e "${YELLOW}AutoDock Vina not found - docking will use estimation fallback${NC}"
    echo "Install from: https://vina.scripps.edu/downloads/"
fi

# Check for Docker (for AlphaFold)
if command -v docker &> /dev/null; then
    echo -e "${GREEN}Docker found${NC}"
else
    echo -e "${YELLOW}Docker not found - AlphaFold 3 integration unavailable${NC}"
fi

# Create necessary directories
echo -e "\n${YELLOW}Creating directories...${NC}"
mkdir -p ../data/receptors
mkdir -p ../output
mkdir -p ../temp

# Create .env file if not exists
if [ ! -f ".env" ]; then
    echo -e "\n${YELLOW}Creating .env file...${NC}"
    cat > .env << EOF
# Drug Design Platform Configuration

# Anthropic API Key (required for AI-powered suggestions)
ANTHROPIC_API_KEY=your_api_key_here

# AlphaFold settings (optional)
ALPHAFOLD_MODEL_DIR=/models/alphafold3
ALPHAFOLD_DB_DIR=/databases/alphafold

# AutoDock Vina settings
AUTODOCK_EXHAUSTIVENESS=32
EOF
    echo -e "${GREEN}.env file created - please add your API keys${NC}"
fi

cd ..

# Run tests
echo -e "\n${YELLOW}Running acceptance tests...${NC}"
cd backend
python -m pytest tests/test_acceptance.py -v --tb=short || {
    echo -e "${YELLOW}Some tests may fail without API keys or optional dependencies${NC}"
}

cd ..

echo -e "\n========================================"
echo -e "${GREEN}Setup Complete!${NC}"
echo "========================================"
echo ""
echo "To start the platform:"
echo "  1. Activate the virtual environment:"
echo "     source venv/bin/activate"
echo ""
echo "  2. Set your Anthropic API key in backend/.env"
echo ""
echo "  3. Start the backend server:"
echo "     cd backend && python main.py"
echo ""
echo "  4. Open frontend/index.html in your browser"
echo "     Or serve it: python -m http.server 3000 --directory frontend"
echo ""
echo "API Documentation: http://localhost:8000/docs"
echo "========================================"
