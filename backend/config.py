"""
Configuration settings for the Drug Design Platform
"""
import os
from pathlib import Path
from typing import Optional

from dotenv import load_dotenv

# Load .env file from the backend directory FIRST before any imports
_backend_dir = Path(__file__).parent
_env_file = _backend_dir / ".env"
if _env_file.exists():
    load_dotenv(_env_file, override=True)

# Now get API keys after dotenv has loaded
_api_key = os.getenv("ANTHROPIC_API_KEY", "")


class Settings:
    """Application settings loaded from environment variables"""

    def __init__(self):
        # API Keys
        self.anthropic_api_key: str = _api_key

        # AlphaFold settings
        self.alphafold_model_dir: str = "/models/alphafold3"
        self.alphafold_db_dir: str = "/databases/alphafold"
        self.alphafold_docker_image: str = "alphafold3:latest"

        # Docking settings
        self.autodock_exhaustiveness: int = 32
        self.autodock_num_modes: int = 9

        # Paths
        self.base_dir: Path = Path(__file__).parent.parent
        self.data_dir: Path = self.base_dir / "data"
        self.output_dir: Path = self.base_dir / "output"
        self.temp_dir: Path = self.base_dir / "temp"

        # Database URLs
        self.pubchem_base_url: str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        self.chembl_base_url: str = "https://www.ebi.ac.uk/chembl/api/data"

        # Generation settings
        self.max_analogues: int = 20
        self.similarity_threshold: float = 0.85
        self.min_qed_score: float = 0.3
        self.max_synthetic_accessibility: float = 6.0

        # Claude settings
        self.claude_model: str = "claude-sonnet-4-20250514"
        self.claude_max_tokens: int = 4096


settings = Settings()

# Ensure directories exist
for dir_path in [settings.data_dir, settings.output_dir, settings.temp_dir]:
    dir_path.mkdir(parents=True, exist_ok=True)
