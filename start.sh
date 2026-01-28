#!/bin/bash
# Activate conda environment and start the application
source /opt/conda/etc/profile.d/conda.sh
conda activate drugdesign
exec uvicorn backend.api.main:app --host 0.0.0.0 --port ${PORT:-8000}
