"""
FastAPI application for the Drug Design Platform.
Provides REST API for molecular design, analysis, and synthesis planning.
"""
from fastapi import FastAPI, BackgroundTasks, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, StreamingResponse
from pydantic import BaseModel, Field
from typing import Dict, List, Optional, Any
from datetime import datetime
import asyncio
import logging
import io
import json
import os

from services.orchestrator import DrugDesignOrchestrator
from services.molecular_analyzer import MolecularAnalyzer
from services.admet_predictor import ADMETPredictor
from services.retrosynthesis import RetrosynthesisService
from services.patent_analyzer import PatentAnalyzer
from services.agentic_designer import AgenticDrugDesigner
from models import JobStatus, DesignType
from config import settings

# Import new prediction infrastructure
from core import (
    prediction_service,
    honest_formatter,
    applicability_checker,
    report_generator,
    DomainStatus,
    TargetKnowledgeBase,
    TargetType,
    BindingMechanism,
)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Initialize FastAPI app
app = FastAPI(
    title="Drug Design Platform API",
    description="Automated platform for novel, patentable molecule design",
    version="1.0.0",
    docs_url="/docs",
    redoc_url="/redoc"
)

# CORS middleware - configurable via ALLOWED_ORIGINS env var
allowed_origins = os.getenv("ALLOWED_ORIGINS", "*").split(",")
app.add_middleware(
    CORSMiddleware,
    allow_origins=allowed_origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Initialize services
orchestrator = DrugDesignOrchestrator()

try:
    analyzer = MolecularAnalyzer()
except Exception as e:
    logger.warning(f"MolecularAnalyzer initialization failed: {e}")
    analyzer = None

admet_predictor = ADMETPredictor()
retrosynthesis_service = RetrosynthesisService()
patent_analyzer = PatentAnalyzer()
agentic_designer = AgenticDrugDesigner(anthropic_api_key=settings.anthropic_api_key)


# Request/Response Models
class DesignRequest(BaseModel):
    """Request for molecule design."""
    prompt: str = Field(..., description="Natural language design request")
    reference_smiles: Optional[str] = Field(None, description="Reference compound SMILES")
    reference_name: Optional[str] = Field(None, description="Reference compound name")
    target_name: Optional[str] = Field(None, description="Target protein/receptor name")
    design_type: Optional[str] = Field("analogue", description="Design type: analogue, bivalent, scaffold_hop, de_novo")
    constraints: Optional[Dict[str, Any]] = Field(default_factory=dict, description="Design constraints")

    class Config:
        json_schema_extra = {
            "example": {
                "prompt": "Design a novel, patentable analogue of KK103 with improved metabolic stability",
                "reference_name": "KK103",
                "design_type": "analogue",
                "constraints": {"max_mw": 600, "bbb_penetration": False}
            }
        }


class AnalyzeRequest(BaseModel):
    """Request for molecule analysis."""
    smiles: str = Field(..., description="SMILES string to analyze")

    class Config:
        json_schema_extra = {
            "example": {
                "smiles": "CC(=O)Oc1ccccc1C(=O)O"
            }
        }


class RetrosynthesisRequest(BaseModel):
    """Request for retrosynthesis planning."""
    smiles: str = Field(..., description="Target molecule SMILES")
    max_steps: int = Field(6, description="Maximum synthesis steps")
    use_ibm_rxn: bool = Field(True, description="Use IBM RXN backend")
    use_aizynthfinder: bool = Field(True, description="Use AiZynthFinder backend")

    class Config:
        json_schema_extra = {
            "example": {
                "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                "max_steps": 5
            }
        }


class PatentCheckRequest(BaseModel):
    """Request for patent/novelty check."""
    smiles: str = Field(..., description="Compound SMILES to check")
    thorough: bool = Field(True, description="Perform thorough search")
    jurisdictions: List[str] = Field(["US", "EP"], description="Target jurisdictions for FTO")


class JobResponse(BaseModel):
    """Response for job submission."""
    job_id: str
    status: str
    message: str


class MoleculeResponse(BaseModel):
    """Response for molecule data."""
    smiles: str
    canonical_smiles: Optional[str]
    inchi: Optional[str]
    inchi_key: Optional[str]
    properties: Optional[Dict[str, Any]]
    valid: bool
    warnings: List[str]


# API Endpoints

@app.get("/")
async def root():
    """API root endpoint."""
    return {
        "name": "Drug Design Platform API",
        "version": "1.0.0",
        "status": "operational",
        "endpoints": {
            "design": "/api/v1/design",
            "analyze": "/api/v1/analyze",
            "retrosynthesis": "/api/v1/retrosynthesis",
            "patent": "/api/v1/patent",
            "docs": "/docs"
        }
    }


@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {
        "status": "healthy",
        "timestamp": datetime.utcnow().isoformat(),
        "services": {
            "molecular_analyzer": analyzer is not None,
            "admet_predictor": True,
            "retrosynthesis": True,
            "patent_analyzer": True,
            "prediction_service": True,
        },
        "prediction_capabilities": len(prediction_service.list_available_predictions()),
    }


# Design endpoints

@app.post("/api/v1/design", response_model=JobResponse)
async def submit_design_job(request: DesignRequest, background_tasks: BackgroundTasks):
    """
    Submit a molecular design job.

    This endpoint accepts a natural language description of the desired molecule
    and returns a job ID for tracking progress.
    """
    try:
        # Start job in background
        job = await orchestrator.process_request(request.prompt)

        return JobResponse(
            job_id=job.id,
            status=job.status.value,
            message="Design job submitted successfully"
        )

    except Exception as e:
        logger.error(f"Design job submission failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/api/v1/design/{job_id}")
async def get_design_job(job_id: str):
    """
    Get the status and results of a design job.
    """
    job = await orchestrator.get_job(job_id)

    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    result = {
        "job_id": job.id,
        "status": job.status.value,
        "progress": job.progress,
        "current_step": job.current_step,
        "design_type": job.design_type.value,
        "created_at": job.created_at.isoformat(),
        "started_at": job.started_at.isoformat() if job.started_at else None,
        "completed_at": job.completed_at.isoformat() if job.completed_at else None,
    }

    if job.status == JobStatus.COMPLETED:
        result["molecules"] = [
            {
                "rank": mol.rank,
                "smiles": mol.smiles,
                "inchi_key": mol.inchi_key,
                "score": mol.overall_score,
                "properties": {
                    "molecular_weight": mol.properties.molecular_weight if mol.properties else None,
                    "logp": mol.properties.logp if mol.properties else None,
                    "qed": mol.properties.qed if mol.properties else None,
                },
                "novelty": {
                    "tanimoto_similarity": mol.novelty.tanimoto_similarity if mol.novelty else None,
                    "patentable": mol.novelty.patentable if mol.novelty else None,
                    "fto_risk": mol.novelty.fto_risk if mol.novelty else None,
                },
                "admet_safety_score": mol.admet.get_safety_score() if mol.admet else None,
                "synthesis_steps": mol.best_route.total_steps if mol.best_route else None,
            }
            for mol in job.molecules
        ]

    elif job.status == JobStatus.FAILED:
        result["error"] = job.error_message

    return result


class AdvancedDesignRequest(BaseModel):
    """Request for advanced agentic design."""
    prompt: str = Field(..., description="Detailed design request")

    class Config:
        json_schema_extra = {
            "example": {
                "prompt": "Design a selective kinase inhibitor targeting EGFR T790M mutation for non-small cell lung cancer"
            }
        }


@app.post("/api/v1/design/advanced/stream")
async def stream_advanced_design(request: AdvancedDesignRequest):
    """
    Stream advanced design with real-time progress updates via Server-Sent Events.

    Each phase sends progress updates so the UI can show what's happening.
    """
    import asyncio
    from services.agentic_designer import AgenticDrugDesigner

    # Create a queue to communicate between the designer and the stream
    progress_queue = asyncio.Queue()

    async def run_design():
        """Run the design process and push progress to the queue."""
        designer = AgenticDrugDesigner(
            anthropic_api_key=settings.anthropic_api_key,
        )

        async def progress_callback(phase: str, message: str, data: dict = None):
            await progress_queue.put({
                "phase": phase,
                "message": message,
                "data": data or {}
            })

        try:
            session = await designer.design_with_progress(
                request=request.prompt,
                progress_callback=progress_callback
            )

            # Push final results
            tk = session.target_knowledge
            await progress_queue.put({
                "phase": "completed",
                "message": "Design complete!",
                "data": {
                    "target": tk.target_name,
                    "target_type": tk.target_type,
                    "mechanism": tk.mechanism,
                    "is_covalent": tk.is_covalent,
                    "candidates": [
                        {
                            "rank": i + 1,
                            "name": c.name,
                            "smiles": c.smiles,
                            "score": c.overall_score,
                            "rationale": c.rationale,
                        }
                        for i, c in enumerate(session.candidates[:10])
                    ]
                }
            })
        except Exception as e:
            await progress_queue.put({
                "phase": "error",
                "message": str(e),
                "data": {}
            })

        # Signal end of stream
        await progress_queue.put(None)

    async def generate_events():
        # Start design task
        design_task = asyncio.create_task(run_design())

        # Send initial event
        yield f"data: {json.dumps({'phase': 'starting', 'message': 'Initializing drug design...', 'data': {}})}\n\n"

        # Stream progress events
        while True:
            try:
                event = await asyncio.wait_for(progress_queue.get(), timeout=300)  # 5 min timeout
                if event is None:
                    break
                yield f"data: {json.dumps(event)}\n\n"
                if event.get("phase") in ["completed", "error"]:
                    break
            except asyncio.TimeoutError:
                yield f"data: {json.dumps({'phase': 'error', 'message': 'Design timed out', 'data': {}})}\n\n"
                break

        # Ensure design task is done
        if not design_task.done():
            design_task.cancel()

    return StreamingResponse(
        generate_events(),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "Connection": "keep-alive",
            "X-Accel-Buffering": "no"
        }
    )


@app.post("/api/v1/design/advanced")
async def submit_advanced_design(request: AdvancedDesignRequest):
    """
    Submit an advanced agentic design job.

    This endpoint uses Claude as an intelligent medicinal chemist to:
    1. Research and fully understand the target
    2. Generate design guidelines and scoring criteria dynamically
    3. Generate and iteratively improve candidates
    4. Return the best molecules with full analysis

    This is slower but produces higher quality results suitable for animal testing.
    """
    try:
        logger.info(f"Starting advanced design: {request.prompt[:100]}...")

        # Create a new designer
        designer = AgenticDrugDesigner(
            anthropic_api_key=settings.anthropic_api_key,
        )

        # Run the agentic design process
        session = await designer.design(
            request=request.prompt
        )

        # Get target knowledge
        tk = session.target_knowledge

        # Format response
        return {
            "status": "completed",
            "target": tk.target_name,
            "target_type": tk.target_type,
            "mechanism": tk.mechanism,
            "essential_features": tk.essential_features,
            "avoid_features": tk.avoid_features,

            # MECHANISTIC ANALYSIS (new)
            "mechanistic_analysis": {
                "is_covalent": tk.is_covalent,
                "warhead_type": tk.warhead_type if tk.is_covalent else None,
                "nucleophile": tk.nucleophile if tk.is_covalent else None,
                "staying_portion_requirements": tk.staying_portion_requirements if tk.is_covalent else [],
                "leaving_group_requirements": tk.leaving_group_requirements if tk.is_covalent else [],
                "mechanistic_constraints": tk.mechanistic_constraints if tk.is_covalent else [],
            } if tk.is_covalent else None,

            "strategy": session.strategy,
            "iterations": session.iterations,
            "converged": session.converged,
            "design_log": session.design_log,
            "reference_compounds": [
                {"name": ref.get("name"), "smiles": ref.get("smiles")}
                for ref in tk.reference_drugs
            ],
            "candidates": [
                {
                    "rank": i + 1,
                    "name": c.name,
                    "smiles": c.smiles,
                    "rationale": c.rationale,
                    "scores": {
                        "overall": c.overall_score,
                        "binding": c.binding_score,
                        "selectivity": c.selectivity_score,
                        "admet": c.admet_score,
                        "novelty": c.novelty_score,
                        "synthesis": c.synthesis_score,
                    },
                    "properties": {
                        "molecular_weight": round(c.molecular_weight, 1),
                        "logp": round(c.logp, 2),
                        "tpsa": round(c.tpsa, 1),
                        "hbd": c.hbd,
                        "hba": c.hba,
                        "qed": round(c.qed, 3),
                    },
                    "confidence": c.overall_confidence,
                    "critical_issues": c.critical_issues,
                    "warnings": c.warnings,
                    "patentable": c.patentable,
                    "iteration": c.iteration,
                    # Include the comprehensive design report if available
                    "design_report": c.design_report,
                }
                for i, c in enumerate(session.candidates)
            ],
            "runtime_seconds": (session.completed_at - session.started_at).total_seconds() if session.completed_at else None,
        }

    except Exception as e:
        logger.error(f"Advanced design failed: {e}")
        import traceback
        logger.error(traceback.format_exc())
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/api/v1/design/{job_id}/molecules/{molecule_index}")
async def get_molecule_details(job_id: str, molecule_index: int):
    """
    Get detailed information about a specific designed molecule.
    """
    job = await orchestrator.get_job(job_id)

    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    if molecule_index < 0 or molecule_index >= len(job.molecules):
        raise HTTPException(status_code=404, detail="Molecule not found")

    mol = job.molecules[molecule_index]

    return {
        "smiles": mol.smiles,
        "inchi": mol.inchi,
        "inchi_key": mol.inchi_key,
        "image_svg": mol.image_svg,
        "mol_block": mol.mol_block,
        "properties": mol.properties.__dict__ if mol.properties else None,
        "admet": mol.admet.__dict__ if mol.admet else None,
        "novelty": mol.novelty.__dict__ if mol.novelty else None,
        "synthesis_routes": [
            {
                "route_id": route.route_id,
                "total_steps": route.total_steps,
                "overall_yield": route.overall_yield,
                "confidence": route.confidence,
                "estimated_cost": route.estimated_cost,
                "source": route.source,
                "steps": [
                    {
                        "step_number": step.step_number,
                        "reaction_name": step.reaction_name,
                        "reaction_smiles": step.reaction_smiles,
                        "conditions": {
                            "solvent": step.solvent,
                            "temperature": step.temperature,
                            "time": step.time,
                        },
                        "expected_yield": step.expected_yield,
                    }
                    for step in route.steps
                ],
                "starting_materials": route.starting_materials,
            }
            for route in mol.synthesis_routes
        ],
        "overall_score": mol.overall_score,
        "rank": mol.rank,
    }


@app.get("/api/v1/design/{job_id}/report")
async def get_design_report(job_id: str):
    """
    Get a comprehensive report for a completed design job.
    """
    job = await orchestrator.get_job(job_id)

    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    if job.status != JobStatus.COMPLETED:
        raise HTTPException(status_code=400, detail="Job not yet completed")

    report = await orchestrator.generate_report(job)

    return {
        "job_id": job_id,
        "report": report
    }


# Analysis endpoints

@app.post("/api/v1/analyze")
async def analyze_molecule(request: AnalyzeRequest):
    """
    Analyze a molecule and return its properties.
    """
    if analyzer is None:
        raise HTTPException(status_code=503, detail="Molecular analyzer not available")

    try:
        # Validate and canonicalize
        validation = analyzer.validate_molecule(request.smiles)

        if not validation["valid"]:
            return MoleculeResponse(
                smiles=request.smiles,
                canonical_smiles=None,
                inchi=None,
                inchi_key=None,
                properties=None,
                valid=False,
                warnings=validation["errors"]
            )

        canonical = validation["canonical_smiles"]
        inchi, inchi_key = analyzer.smiles_to_inchi(canonical)
        properties = analyzer.calculate_properties(canonical)

        return MoleculeResponse(
            smiles=request.smiles,
            canonical_smiles=canonical,
            inchi=inchi,
            inchi_key=inchi_key,
            properties=properties.__dict__ if properties else None,
            valid=True,
            warnings=validation.get("warnings", [])
        )

    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/v1/analyze/admet")
async def predict_admet(request: AnalyzeRequest):
    """
    Predict ADMET properties for a molecule.
    """
    try:
        profile = await admet_predictor.predict(request.smiles)

        return {
            "smiles": request.smiles,
            "admet": profile.__dict__,
            "safety_score": profile.get_safety_score(),
            "report": admet_predictor.format_admet_report(profile)
        }

    except Exception as e:
        logger.error(f"ADMET prediction failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/v1/analyze/similarity")
async def calculate_similarity(
    smiles1: str = Query(..., description="First molecule SMILES"),
    smiles2: str = Query(..., description="Second molecule SMILES")
):
    """
    Calculate Tanimoto similarity between two molecules.
    """
    if analyzer is None:
        raise HTTPException(status_code=503, detail="Molecular analyzer not available")

    try:
        similarity = analyzer.calculate_similarity(smiles1, smiles2)

        if similarity is None:
            raise HTTPException(status_code=400, detail="Invalid SMILES provided")

        return {
            "smiles1": smiles1,
            "smiles2": smiles2,
            "tanimoto_similarity": round(similarity, 4)
        }

    except Exception as e:
        logger.error(f"Similarity calculation failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/api/v1/analyze/image/{smiles}")
async def get_molecule_image(smiles: str, size: int = 300):
    """
    Get a 2D image of a molecule.
    """
    if analyzer is None:
        raise HTTPException(status_code=503, detail="Molecular analyzer not available")

    try:
        svg = analyzer.molecule_to_svg(smiles, size=(size, size))

        if svg is None:
            raise HTTPException(status_code=400, detail="Invalid SMILES")

        return StreamingResponse(
            io.BytesIO(svg.encode()),
            media_type="image/svg+xml"
        )

    except Exception as e:
        logger.error(f"Image generation failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# Retrosynthesis endpoints

@app.post("/api/v1/retrosynthesis")
async def plan_retrosynthesis(request: RetrosynthesisRequest):
    """
    Plan synthesis routes for a target molecule.
    """
    try:
        routes = await retrosynthesis_service.plan_synthesis(
            smiles=request.smiles,
            max_steps=request.max_steps,
            use_ibm_rxn=request.use_ibm_rxn,
            use_aizynthfinder=request.use_aizynthfinder
        )

        return {
            "smiles": request.smiles,
            "routes": [
                {
                    "route_id": route.route_id,
                    "total_steps": route.total_steps,
                    "overall_yield": route.overall_yield,
                    "confidence": route.confidence,
                    "estimated_cost": route.estimated_cost,
                    "source": route.source,
                    "steps": [
                        {
                            "step_number": step.step_number,
                            "reaction_name": step.reaction_name,
                            "reaction_smiles": step.reaction_smiles,
                            "reactants": step.reactants,
                            "reagents": step.reagents,
                            "solvent": step.solvent,
                            "temperature": step.temperature,
                            "time": step.time,
                            "expected_yield": step.expected_yield,
                            "confidence": step.confidence,
                        }
                        for step in route.steps
                    ],
                    "starting_materials": route.starting_materials,
                    "report": retrosynthesis_service.format_synthesis_report(route)
                }
                for route in routes
            ],
            "total_routes": len(routes)
        }

    except Exception as e:
        logger.error(f"Retrosynthesis planning failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# Patent endpoints

@app.post("/api/v1/patent/novelty")
async def check_novelty(request: PatentCheckRequest):
    """
    Check novelty of a compound against patent literature.
    """
    try:
        result = await patent_analyzer.check_novelty(
            smiles=request.smiles,
            thorough=request.thorough
        )

        return {
            "smiles": request.smiles,
            "novelty": result
        }

    except Exception as e:
        logger.error(f"Novelty check failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/v1/patent/fto")
async def assess_fto(request: PatentCheckRequest):
    """
    Assess freedom-to-operate for a compound.
    """
    try:
        result = await patent_analyzer.assess_fto(
            smiles=request.smiles,
            jurisdictions=request.jurisdictions
        )

        return {
            "smiles": request.smiles,
            "fto": result
        }

    except Exception as e:
        logger.error(f"FTO assessment failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/v1/patent/patentability")
async def assess_patentability(request: PatentCheckRequest):
    """
    Assess patentability of a compound.
    """
    try:
        # Get novelty first
        novelty_result = await patent_analyzer.check_novelty(
            smiles=request.smiles,
            thorough=request.thorough
        )

        # Get FTO
        fto_result = await patent_analyzer.assess_fto(
            smiles=request.smiles,
            jurisdictions=request.jurisdictions
        )

        # Assess patentability
        patentability = patent_analyzer.assess_patentability(novelty_result)

        # Generate report
        report = patent_analyzer.generate_patent_report(
            request.smiles,
            novelty_result,
            fto_result,
            patentability
        )

        return {
            "smiles": request.smiles,
            "novelty": novelty_result,
            "fto": fto_result,
            "patentability": patentability,
            "report": report
        }

    except Exception as e:
        logger.error(f"Patentability assessment failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# ========================================
# NEW PREDICTION INFRASTRUCTURE ENDPOINTS
# ========================================

class ScoreRequest(BaseModel):
    """Request for comprehensive molecule scoring."""
    smiles: str = Field(..., description="SMILES string to score")
    target: Optional[str] = Field(None, description="Target protein for binding predictions")
    reference_smiles: Optional[List[str]] = Field(None, description="Reference compounds for novelty")

    class Config:
        json_schema_extra = {
            "example": {
                "smiles": "O=C(Nc1ccc(OC(F)(F)F)cc1)N1CCC(c2ccnc3ccccc23)CC1",
                "target": "FAAH",
            }
        }


class BatchScoreRequest(BaseModel):
    """Request for batch scoring of multiple molecules."""
    molecules: List[str] = Field(..., description="List of SMILES strings")
    target: Optional[str] = Field(None, description="Target protein")
    reference_smiles: Optional[List[str]] = Field(None, description="Reference compounds")


@app.post("/api/v2/score")
async def score_molecule(request: ScoreRequest):
    """
    Comprehensive molecule scoring with honest uncertainty.

    This endpoint uses the new prediction infrastructure to provide:
    - Multi-dimensional scoring (drug-likeness, safety, synthesis, novelty, binding)
    - Confidence levels for each prediction
    - Applicability domain checking
    - Actionable recommendations

    Returns a score card formatted for clear, honest communication of results.
    """
    try:
        # Check applicability domain first
        domain_check = applicability_checker.check_domain(
            request.smiles,
            reference_smiles=request.reference_smiles
        )

        # Score the molecule
        score = await prediction_service.score_molecule(
            smiles=request.smiles,
            target=request.target,
            reference_compounds=request.reference_smiles,
        )

        # Format for honest output
        score_card = honest_formatter.format_score_card(score)

        # Add domain information
        result = score_card.to_dict()
        result["applicability_domain"] = {
            "status": domain_check.status.value,
            "confidence_modifier": domain_check.confidence_modifier,
            "reasons": domain_check.reasons,
        }

        return result

    except Exception as e:
        logger.error(f"Scoring failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/v2/score/batch")
async def batch_score_molecules(request: BatchScoreRequest):
    """
    Score multiple molecules in parallel with honest uncertainty.

    Returns a summary of results with top candidates and common issues.
    """
    try:
        # Score all molecules
        scores = await prediction_service.batch_score(
            molecules=request.molecules,
            target=request.target,
            reference_compounds=request.reference_smiles,
        )

        # Format batch summary
        summary = honest_formatter.format_batch_summary(scores)

        return summary

    except Exception as e:
        logger.error(f"Batch scoring failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/v2/predict")
async def predict_properties(
    smiles: str = Query(..., description="SMILES string"),
    properties: Optional[List[str]] = Query(None, description="Properties to predict")
):
    """
    Predict specific properties for a molecule with confidence levels.

    If no properties specified, returns all available predictions.
    """
    try:
        predictions = await prediction_service.predict(
            smiles=smiles,
            properties=properties,
        )

        # Format each prediction
        formatted = {}
        for prop, pred in predictions.items():
            formatted_pred = honest_formatter.format_prediction(prop, pred)
            formatted[prop] = formatted_pred.to_dict()

        return {
            "smiles": smiles,
            "predictions": formatted,
        }

    except Exception as e:
        logger.error(f"Prediction failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/api/v2/capabilities")
async def list_prediction_capabilities():
    """
    List all available prediction capabilities.

    Returns the properties that can be predicted along with their
    confidence tiers and requirements.
    """
    try:
        predictors = prediction_service.list_predictors()
        properties = prediction_service.list_available_predictions()

        return {
            "available_properties": properties,
            "predictors": predictors,
            "total_predictors": len(predictors),
            "total_properties": len(properties),
        }

    except Exception as e:
        logger.error(f"Failed to list capabilities: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/v2/domain-check")
async def check_applicability_domain(
    smiles: str = Query(..., description="SMILES string to check"),
    reference_smiles: Optional[List[str]] = Query(None, description="Reference compounds")
):
    """
    Check if a molecule is within the applicability domain of predictors.

    This helps determine how reliable predictions will be for this molecule.
    """
    try:
        result = applicability_checker.check_domain(
            smiles=smiles,
            reference_smiles=reference_smiles
        )

        return {
            "smiles": smiles,
            "status": result.status.value,
            "is_reliable": result.is_reliable,
            "needs_caution": result.needs_caution,
            "confidence_modifier": result.confidence_modifier,
            "reasons": result.reasons,
            "details": result.details,
        }

    except Exception as e:
        logger.error(f"Domain check failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


class DesignReportRequest(BaseModel):
    """Request for generating an actionable design report."""
    smiles: str = Field(..., description="SMILES string to evaluate")
    candidate_name: Optional[str] = Field("Candidate", description="Name for the candidate")
    target: Optional[str] = Field(None, description="Target protein name")
    target_type: Optional[str] = Field(None, description="Type: enzyme, kinase, gpcr, etc.")
    mechanism: Optional[str] = Field(None, description="Binding mechanism: covalent, competitive, etc.")
    reference_smiles: Optional[List[str]] = Field(None, description="Reference compounds for novelty")

    class Config:
        json_schema_extra = {
            "example": {
                "smiles": "O=C(Nc1ccc(OC(F)(F)F)cc1)N1CCC(c2ccnc3ccccc23)CC1",
                "candidate_name": "FAAH-001",
                "target": "FAAH",
                "target_type": "enzyme",
                "mechanism": "covalent",
            }
        }


@app.post("/api/v2/design-report")
async def generate_design_report(request: DesignReportRequest):
    """
    Generate a comprehensive, actionable design report for a drug candidate.

    This endpoint provides decision-support for medicinal chemistry:
    - Clear recommendation: ADVANCE, OPTIMIZE, DEPRIORITIZE, or REJECT
    - Executive summary with key strengths and concerns
    - Target-contextualized analysis
    - Prioritized next steps (immediate, optimization, experimental)
    - Go/no-go criteria
    - Honest documentation of limitations and unknowns

    This is what a GREAT drug design output looks like.
    """
    try:
        # Score the molecule first
        score = await prediction_service.score_molecule(
            smiles=request.smiles,
            target=request.target,
            reference_compounds=request.reference_smiles,
        )

        # Build target knowledge if provided
        target_kb = None
        if request.target:
            # Map target type string to enum
            target_type = TargetType.UNKNOWN
            if request.target_type:
                type_map = {
                    "enzyme": TargetType.ENZYME,
                    "kinase": TargetType.KINASE,
                    "gpcr": TargetType.GPCR,
                    "ion_channel": TargetType.ION_CHANNEL,
                    "nuclear_receptor": TargetType.NUCLEAR_RECEPTOR,
                    "transporter": TargetType.TRANSPORTER,
                    "ppi": TargetType.PROTEIN_PROTEIN,
                    "protein-protein": TargetType.PROTEIN_PROTEIN,
                }
                target_type = type_map.get(request.target_type.lower(), TargetType.UNKNOWN)

            # Map mechanism to enum
            binding_mechanism = BindingMechanism.UNKNOWN
            if request.mechanism:
                mech_map = {
                    "covalent": BindingMechanism.COVALENT_IRREVERSIBLE,
                    "irreversible": BindingMechanism.COVALENT_IRREVERSIBLE,
                    "competitive": BindingMechanism.COMPETITIVE,
                    "allosteric": BindingMechanism.ALLOSTERIC,
                    "uncompetitive": BindingMechanism.UNCOMPETITIVE,
                }
                binding_mechanism = mech_map.get(request.mechanism.lower(), BindingMechanism.UNKNOWN)

            target_kb = TargetKnowledgeBase(
                target_name=request.target,
                target_type=target_type,
                binding_mechanism=binding_mechanism,
                is_covalent=binding_mechanism == BindingMechanism.COVALENT_IRREVERSIBLE,
            )

        # Generate the design report
        report = report_generator.generate_report(
            score=score,
            target_knowledge=target_kb,
            candidate_name=request.candidate_name or "Candidate",
        )

        return report.to_dict()

    except Exception as e:
        logger.error(f"Design report generation failed: {e}")
        import traceback
        logger.error(traceback.format_exc())
        raise HTTPException(status_code=500, detail=str(e))


# ========================================
# PDB STRUCTURE ENDPOINTS
# ========================================

@app.get("/api/v2/pdb/{target}")
async def get_pdb_structure(target: str):
    """
    Fetch PDB structure for a target protein.

    Searches RCSB PDB first, falls back to AlphaFold if no experimental structure found.
    """
    try:
        from services.pdb_fetcher import get_pdb_fetcher

        pdb_fetcher = get_pdb_fetcher(cache_dir=settings.data_dir / "pdb_cache")
        result = await pdb_fetcher.fetch_by_target_name(target)

        if not result.success:
            raise HTTPException(
                status_code=404,
                detail=f"Could not find structure for: {target}. {result.error_message}"
            )

        return {
            "target": target,
            "pdb_id": result.pdb_info.pdb_id,
            "source": result.pdb_info.source,
            "title": result.pdb_info.title,
            "resolution": result.pdb_info.resolution,
            "pdb_content_length": len(result.pdb_content),
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"PDB fetch failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# Run with: uvicorn backend.api.main:app --reload
if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
