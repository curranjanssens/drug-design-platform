# Drug Design Platform - Backend Dockerfile
# Using conda for RDKit which requires specialized builds
FROM continuumio/miniconda3:latest

# Set working directory
WORKDIR /app

# Install RDKit via conda (this is the reliable way)
RUN conda install -c conda-forge rdkit=2023.09.1 -y && \
    conda clean -afy

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements (simplified)
COPY requirements-railway.txt .

# Install Python dependencies via pip (non-chemistry ones)
RUN pip install --no-cache-dir -r requirements-railway.txt

# Copy application code
COPY . .

# Create data directories
RUN mkdir -p /app/data /app/output /app/temp

# Expose port
EXPOSE 8000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8000/health || exit 1

# Run the application
CMD ["uvicorn", "backend.api.main:app", "--host", "0.0.0.0", "--port", "8000"]
