# Drug Design Platform - Backend Dockerfile
# Using conda for RDKit which requires specialized builds
FROM continuumio/miniconda3:latest

# Set working directory
WORKDIR /app

# Create conda environment with Python 3.11 and RDKit
RUN conda create -n drugdesign python=3.11 rdkit -c conda-forge -y && \
    conda clean -afy

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "drugdesign", "/bin/bash", "-c"]

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements (simplified)
COPY requirements-railway.txt .

# Install Python dependencies via pip
RUN pip install --no-cache-dir -r requirements-railway.txt

# Copy application code
COPY . .

# Create data directories
RUN mkdir -p /app/data /app/output /app/temp

# Make startup script executable
RUN chmod +x /app/start.sh

# Default port (Railway will override via PORT env var)
ENV PORT=8000

# Expose port (informational)
EXPOSE 8000

# Use bash to run the startup script
CMD ["/bin/bash", "/app/start.sh"]
