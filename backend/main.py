"""
FastAPI Main Application for Drug Design Platform
"""
from fastapi import FastAPI, HTTPException, BackgroundTasks, WebSocket
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, FileResponse
from fastapi.staticfiles import StaticFiles
from contextlib import asynccontextmanager
import asyncio
from typing import Optional
from datetime import datetime
import json
from pathlib import Path
from loguru import logger

from models import (
    CompoundInput, AnalogueGenerationRequest, AnalogueGenerationResponse,
    PipelineStatus, InputFormat, KK103_REFERENCE
)
from pipeline import DrugDesignPipeline, PREBUILT_KK103_ANALOGUES
from services.property_predictor import PropertyPredictor
from services.novelty_checker import NoveltyChecker
from services.docking_service import DockingService, StructureVisualization
from config import settings


# Initialize services
pipeline = DrugDesignPipeline()
property_predictor = PropertyPredictor()
novelty_checker = NoveltyChecker()
docking_service = DockingService()


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application lifespan handler"""
    logger.info("Starting Drug Design Platform")
    yield
    logger.info("Shutting down Drug Design Platform")


app = FastAPI(
    title="Drug Design Platform",
    description="Automated platform for generating novel, patentable drug analogues",
    version="1.0.0",
    lifespan=lifespan
)

# CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# ============= API Endpoints =============

@app.get("/")
async def root():
    """Root endpoint"""
    return {
        "name": "Drug Design Platform",
        "version": "1.0.0",
        "status": "running",
        "endpoints": {
            "generate": "/api/generate",
            "status": "/api/status/{request_id}",
            "predict_properties": "/api/predict/properties",
            "check_novelty": "/api/check/novelty",
            "kk103_demo": "/api/demo/kk103"
        }
    }


@app.post("/api/generate", response_model=AnalogueGenerationResponse)
async def generate_analogues(request: AnalogueGenerationRequest):
    """
    Generate novel analogues for an input compound

    This endpoint runs the full drug design pipeline:
    1. Parse and validate input
    2. Generate structural modifications
    3. Predict properties
    4. Check novelty
    5. Rank and filter candidates
    """
    try:
        result = await pipeline.run_pipeline(request)
        return result
    except Exception as e:
        logger.error(f"Generation error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/api/status/{request_id}", response_model=PipelineStatus)
async def get_pipeline_status(request_id: str):
    """Get the status of a running pipeline"""
    status = pipeline.get_status(request_id)
    if status is None:
        raise HTTPException(status_code=404, detail="Pipeline not found")
    return status


@app.post("/api/predict/properties")
async def predict_properties(smiles: str):
    """Predict properties for a single compound"""
    try:
        properties = property_predictor.predict_properties(smiles)
        if properties is None:
            raise HTTPException(status_code=400, detail="Invalid SMILES")

        admet = property_predictor.predict_admet(smiles)
        druglikeness = property_predictor.calculate_druglikeness_profile(smiles)

        return {
            "smiles": smiles,
            "basic_properties": properties.dict(),
            "admet": admet,
            "druglikeness": druglikeness
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/check/novelty")
async def check_novelty(smiles: str):
    """Check novelty of a compound against databases"""
    try:
        result = await novelty_checker.check_novelty(smiles)
        return result.dict()
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/dock")
async def dock_compound(
    smiles: str,
    receptor_pdb_id: Optional[str] = None
):
    """Dock a compound against a receptor"""
    try:
        result = await docking_service.dock_compound(smiles)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/api/visualize/2d/{smiles}")
async def get_2d_structure(smiles: str):
    """Get 2D structure image for a SMILES"""
    try:
        image_data = StructureVisualization.generate_2d_image(smiles)
        if image_data is None:
            raise HTTPException(status_code=400, detail="Invalid SMILES")

        return JSONResponse(
            content={"image": image_data.hex()},
            media_type="application/json"
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/api/visualize/3d/{smiles}")
async def get_3d_structure(smiles: str):
    """Get 3D conformer for a SMILES"""
    try:
        mol_block = StructureVisualization.generate_3d_conformer(smiles)
        if mol_block is None:
            raise HTTPException(status_code=400, detail="Invalid SMILES")

        return {"mol_block": mol_block}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# ============= Demo Endpoints =============

@app.get("/api/demo/kk103")
async def kk103_demo():
    """
    Get information about KK103 and pre-built analogues for demonstration

    This endpoint provides immediate access to KK103 data without running
    the full pipeline.
    """
    return {
        "reference_compound": KK103_REFERENCE,
        "prebuilt_analogues": PREBUILT_KK103_ANALOGUES,
        "note": "These are pre-computed analogues for demonstration. Use /api/generate for full pipeline."
    }


@app.post("/api/demo/kk103/full")
async def run_kk103_full_pipeline():
    """
    Run the full pipeline for KK103

    This is the acceptance test endpoint that demonstrates the complete
    capability of the platform.
    """
    request = AnalogueGenerationRequest(
        compound=CompoundInput(
            name="KK103",
            structure="CC(C)(C)C(=O)NC(Cc1ccc(O)cc1)C(=O)NCC(=O)NCC(=O)NC(Cc2ccccc2)C(=O)NC(CC(C)C)C(=O)O",
            format=InputFormat.SMILES,
            target_receptor="DOR",
            constraints={
                "min_qed": 0.3,
                "max_sa": 6.0,
                "maintain_binding": True
            }
        ),
        num_analogues=10,
        strategies=["n_terminal", "amino_acid_substitution", "c_terminal", "backbone"],
        include_docking=True,
        include_alphafold=False
    )

    try:
        result = await pipeline.run_pipeline(request)
        return result
    except Exception as e:
        logger.error(f"KK103 demo error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# ============= Health Check =============

@app.get("/health")
async def health_check():
    """Health check endpoint"""
    return {
        "status": "healthy",
        "timestamp": datetime.now().isoformat(),
        "services": {
            "pipeline": "ready",
            "property_predictor": "ready",
            "novelty_checker": "ready",
            "docking": "ready" if docking_service else "unavailable"
        }
    }


# ============= WebSocket for Progress Updates =============

@app.websocket("/ws/progress/{request_id}")
async def websocket_progress(websocket: WebSocket, request_id: str):
    """WebSocket endpoint for real-time progress updates"""
    await websocket.accept()

    try:
        while True:
            status = pipeline.get_status(request_id)
            if status:
                await websocket.send_json(status.dict())
                if status.status in ["completed", "failed"]:
                    break
            await asyncio.sleep(1)
    except Exception as e:
        logger.error(f"WebSocket error: {e}")
    finally:
        await websocket.close()


# ============= Main Entry Point =============

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(
        "main:app",
        host="0.0.0.0",
        port=8000,
        reload=True,
        log_level="info"
    )
