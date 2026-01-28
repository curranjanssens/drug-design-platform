"""
Honest Output Formatting (Phase 4)

This module ensures all predictions and results are communicated with
appropriate epistemic humility. It transforms raw predictions into
user-facing outputs that honestly convey uncertainty.

Design Philosophy:
- Never present uncertain predictions as definitive
- Distinguish between computed facts and predictions
- Communicate confidence levels clearly
- Highlight what we don't know
- Provide actionable recommendations based on uncertainty
"""

from typing import Dict, List, Any, Optional
from dataclasses import dataclass, field
from enum import Enum

from .prediction import Prediction, ConfidenceLevel, PredictionBasis
from .prediction_service import MoleculeScore


class OutputTone(Enum):
    """Tone for output messages based on confidence."""
    DEFINITIVE = "definitive"  # For computed facts (MW, LogP)
    CONFIDENT = "confident"     # For high-confidence predictions
    TENTATIVE = "tentative"     # For moderate-confidence predictions
    CAUTIOUS = "cautious"       # For low-confidence predictions
    UNKNOWN = "unknown"         # For unassessed properties


@dataclass
class FormattedPrediction:
    """A prediction formatted for human consumption."""
    property_name: str
    display_name: str
    value: Any
    value_formatted: str
    unit: str = ""

    # Confidence communication
    confidence_label: str = ""  # "High", "Moderate", "Low", "Unknown"
    confidence_explanation: str = ""

    # Tone and caveats
    tone: OutputTone = OutputTone.TENTATIVE
    caveats: List[str] = field(default_factory=list)
    recommendation: str = ""

    # Source information
    basis: str = ""  # How the prediction was made
    predictor: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "property": self.property_name,
            "display_name": self.display_name,
            "value": self.value,
            "value_formatted": self.value_formatted,
            "unit": self.unit,
            "confidence": {
                "label": self.confidence_label,
                "explanation": self.confidence_explanation,
            },
            "tone": self.tone.value,
            "caveats": self.caveats,
            "recommendation": self.recommendation,
            "basis": self.basis,
            "predictor": self.predictor,
        }


@dataclass
class FormattedScoreCard:
    """A complete molecule score card with honest uncertainty."""
    smiles: str
    molecule_name: str = ""

    # Overall assessment
    overall_score: float = 0.0
    overall_grade: str = ""  # A, B, C, D, F
    overall_confidence: str = ""
    overall_summary: str = ""

    # Component scores with uncertainty
    component_scores: Dict[str, Dict[str, Any]] = field(default_factory=dict)

    # Individual predictions
    predictions: List[FormattedPrediction] = field(default_factory=list)

    # Issues and recommendations
    critical_issues: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    recommendations: List[str] = field(default_factory=list)

    # What we don't know
    unassessed_properties: List[str] = field(default_factory=list)
    limitations: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "smiles": self.smiles,
            "molecule_name": self.molecule_name,
            "overall": {
                "score": round(self.overall_score, 3),
                "grade": self.overall_grade,
                "confidence": self.overall_confidence,
                "summary": self.overall_summary,
            },
            "component_scores": self.component_scores,
            "predictions": [p.to_dict() for p in self.predictions],
            "issues": {
                "critical": self.critical_issues,
                "warnings": self.warnings,
            },
            "recommendations": self.recommendations,
            "limitations": {
                "unassessed": self.unassessed_properties,
                "caveats": self.limitations,
            },
        }


class HonestOutputFormatter:
    """
    Formats predictions and scores for honest, clear communication.

    This is the last stage before results go to the user. It ensures:
    1. Confidence levels are communicated clearly
    2. Caveats are not buried
    3. Recommendations are actionable
    4. Uncertainty is acknowledged, not hidden
    """

    # Display names for properties
    PROPERTY_DISPLAY_NAMES = {
        "molecular_weight": "Molecular Weight",
        "logp": "LogP (Lipophilicity)",
        "tpsa": "Polar Surface Area",
        "hbd": "H-Bond Donors",
        "hba": "H-Bond Acceptors",
        "rotatable_bonds": "Rotatable Bonds",
        "qed": "Drug-likeness (QED)",
        "lipinski_violations": "Lipinski Violations",
        "synthetic_accessibility": "Synthetic Accessibility",
        "synthesis_feasibility_score": "Synthesis Feasibility",
        "oral_absorption_estimate": "Oral Absorption",
        "bbb_penetration_estimate": "BBB Penetration",
        "herg_liability_estimate": "hERG Liability",
        "ames_estimate": "Mutagenicity (Ames)",
        "pains_alerts": "PAINS Alerts",
        "brenk_alerts": "Structural Alerts",
        "binding_affinity_estimate": "Binding Affinity",
        "similarity_to_known_actives": "Similarity to Actives",
        "novelty_score": "Novelty Score",
        "admet_safety_score": "ADMET Safety",
    }

    PROPERTY_UNITS = {
        "molecular_weight": "Da",
        "logp": "",
        "tpsa": "Å²",
        "synthetic_accessibility": "(1-10 scale)",
        "binding_affinity_estimate": "pIC50",
    }

    # Confidence label mapping
    CONFIDENCE_LABELS = {
        ConfidenceLevel.HIGH: ("High", "This is a reliable calculation or well-validated prediction."),
        ConfidenceLevel.MODERATE: ("Moderate", "This is a reasonable estimate but should be verified experimentally."),
        ConfidenceLevel.LOW: ("Low", "This is a rough estimate based on heuristics. Treat with caution."),
        ConfidenceLevel.OUTSIDE_DOMAIN: ("Outside Domain", "The model is extrapolating. This prediction may be unreliable."),
        ConfidenceLevel.NOT_ASSESSED: ("Unknown", "This property could not be assessed."),
    }

    # Basis explanations
    BASIS_EXPLANATIONS = {
        PredictionBasis.HEURISTIC: "Calculated using rule-based methods",
        PredictionBasis.VALIDATED_MODEL: "Predicted by validated ML model",
        PredictionBasis.ML_MODEL: "Predicted by ML model (not prospectively validated)",
        PredictionBasis.SIMILARITY: "Estimated from similar known compounds",
        PredictionBasis.DOCKING: "Estimated from molecular docking",
        PredictionBasis.LITERATURE: "Derived from database/literature",
        PredictionBasis.LLM_REASONING: "Estimated by AI reasoning",
        PredictionBasis.AGGREGATED: "Combined from multiple predictions",
    }

    def format_prediction(self, prop_name: str, prediction: Prediction) -> FormattedPrediction:
        """Format a single prediction for display."""
        display_name = self.PROPERTY_DISPLAY_NAMES.get(prop_name, prop_name.replace("_", " ").title())
        unit = self.PROPERTY_UNITS.get(prop_name, "")

        # Format value
        if not prediction.assessed:
            value_formatted = "Not assessed"
        elif isinstance(prediction.value, bool):
            value_formatted = "Yes" if prediction.value else "No"
        elif isinstance(prediction.value, float):
            value_formatted = f"{prediction.value:.2f}"
        elif isinstance(prediction.value, int):
            value_formatted = str(prediction.value)
        else:
            value_formatted = str(prediction.value)

        if unit:
            value_formatted = f"{value_formatted} {unit}"

        # Confidence label
        conf_label, conf_explanation = self.CONFIDENCE_LABELS.get(
            prediction.confidence,
            ("Unknown", "Confidence level not specified")
        )

        # Determine tone
        if prediction.confidence == ConfidenceLevel.HIGH:
            tone = OutputTone.DEFINITIVE if prediction.basis == PredictionBasis.HEURISTIC else OutputTone.CONFIDENT
        elif prediction.confidence == ConfidenceLevel.MODERATE:
            tone = OutputTone.TENTATIVE
        elif prediction.confidence in (ConfidenceLevel.LOW, ConfidenceLevel.OUTSIDE_DOMAIN):
            tone = OutputTone.CAUTIOUS
        else:
            tone = OutputTone.UNKNOWN

        # Basis explanation
        basis = self.BASIS_EXPLANATIONS.get(prediction.basis, "Method unspecified")

        # Recommendation based on confidence
        if prediction.confidence == ConfidenceLevel.LOW:
            recommendation = "Experimental verification strongly recommended"
        elif prediction.confidence == ConfidenceLevel.OUTSIDE_DOMAIN:
            recommendation = "Consider experimental measurement - prediction may be unreliable"
        elif prediction.confidence == ConfidenceLevel.MODERATE:
            recommendation = "Consider experimental verification for critical decisions"
        else:
            recommendation = ""

        return FormattedPrediction(
            property_name=prop_name,
            display_name=display_name,
            value=prediction.value,
            value_formatted=value_formatted,
            unit=unit,
            confidence_label=conf_label,
            confidence_explanation=conf_explanation,
            tone=tone,
            caveats=prediction.caveats.copy(),
            recommendation=recommendation,
            basis=basis,
            predictor=prediction.predictor_id or "unknown",
        )

    def format_score_card(
        self,
        score: MoleculeScore,
        molecule_name: str = ""
    ) -> FormattedScoreCard:
        """Format a complete molecule score card."""
        # Overall grade
        if score.overall_score >= 0.8:
            grade = "A"
            summary = "Excellent drug-like profile"
        elif score.overall_score >= 0.7:
            grade = "B"
            summary = "Good drug-like profile with minor concerns"
        elif score.overall_score >= 0.6:
            grade = "C"
            summary = "Moderate drug-like profile - optimization recommended"
        elif score.overall_score >= 0.5:
            grade = "D"
            summary = "Poor drug-like profile - significant issues"
        else:
            grade = "F"
            summary = "Unsuitable - major issues detected"

        # Format confidence
        conf_label, _ = self.CONFIDENCE_LABELS.get(
            score.overall_confidence,
            ("Unknown", "")
        )

        # Component scores with uncertainty
        component_scores = {
            "drug_likeness": {
                "score": round(score.drug_likeness_score, 3),
                "label": self._score_to_label(score.drug_likeness_score),
                "description": "Overall drug-like properties",
            },
            "safety": {
                "score": round(score.safety_score, 3),
                "label": self._score_to_label(score.safety_score),
                "description": "ADMET and toxicity assessment",
            },
            "synthesis": {
                "score": round(score.synthesis_score, 3),
                "label": self._score_to_label(score.synthesis_score),
                "description": "Synthetic accessibility",
            },
            "novelty": {
                "score": round(score.novelty_score, 3),
                "label": self._score_to_label(score.novelty_score),
                "description": "Structural novelty vs known compounds",
            },
        }

        if score.binding_score is not None:
            component_scores["binding"] = {
                "score": round(score.binding_score, 3),
                "label": self._score_to_label(score.binding_score),
                "description": "Predicted target binding",
            }

        # Format individual predictions
        formatted_predictions = []
        for prop_name, pred in score.predictions.items():
            # Reconstruct Prediction from dict if needed
            if isinstance(pred, dict):
                pred_obj = Prediction(
                    value=pred.get("value"),
                    confidence=ConfidenceLevel(pred.get("confidence", "moderate")),
                    confidence_score=pred.get("confidence_score"),
                    basis=PredictionBasis(pred.get("basis", "unknown")),
                    caveats=pred.get("caveats", []),
                    assessed=pred.get("assessed", True),
                )
            else:
                pred_obj = pred

            formatted_predictions.append(self.format_prediction(prop_name, pred_obj))

        # Generate recommendations
        recommendations = self._generate_recommendations(score)

        # Limitations
        limitations = [
            "All ADMET predictions are estimates - experimental validation required",
            "Binding affinity predictions are similarity-based, not structure-based",
        ]

        if score.overall_confidence in (ConfidenceLevel.LOW, ConfidenceLevel.OUTSIDE_DOMAIN):
            limitations.append("Overall confidence is low - treat all predictions with caution")

        return FormattedScoreCard(
            smiles=score.smiles,
            molecule_name=molecule_name,
            overall_score=score.overall_score,
            overall_grade=grade,
            overall_confidence=conf_label,
            overall_summary=summary,
            component_scores=component_scores,
            predictions=formatted_predictions,
            critical_issues=score.critical_issues.copy(),
            warnings=score.warnings.copy(),
            recommendations=recommendations,
            unassessed_properties=score.unassessed_properties.copy(),
            limitations=limitations,
        )

    def _score_to_label(self, score: float) -> str:
        """Convert numeric score to label."""
        if score >= 0.8:
            return "Excellent"
        elif score >= 0.7:
            return "Good"
        elif score >= 0.6:
            return "Moderate"
        elif score >= 0.5:
            return "Poor"
        else:
            return "Very Poor"

    def _generate_recommendations(self, score: MoleculeScore) -> List[str]:
        """Generate actionable recommendations based on score."""
        recommendations = []

        # Drug-likeness recommendations
        if score.drug_likeness_score < 0.6:
            recommendations.append(
                "Consider reducing molecular weight or LogP to improve drug-likeness"
            )

        # Safety recommendations
        if score.safety_score < 0.7:
            recommendations.append(
                "Run experimental ADMET assays early - predicted safety concerns"
            )

        if score.critical_issues:
            recommendations.append(
                "Address critical structural alerts before proceeding"
            )

        # Synthesis recommendations
        if score.synthesis_score < 0.5:
            recommendations.append(
                "Consult with synthetic chemists - this structure may be challenging to make"
            )

        # Novelty recommendations
        if score.novelty_score < 0.3:
            recommendations.append(
                "Structure is very similar to known compounds - consider patentability"
            )
        elif score.novelty_score > 0.8:
            recommendations.append(
                "Highly novel structure - validate binding hypothesis carefully"
            )

        # General
        if not recommendations:
            recommendations.append(
                "Proceed with standard lead optimization and experimental validation"
            )

        return recommendations

    def format_batch_summary(
        self,
        scores: List[MoleculeScore],
        top_n: int = 5
    ) -> Dict[str, Any]:
        """Format a summary of batch scoring results."""
        if not scores:
            return {"error": "No molecules scored"}

        # Sort by overall score
        sorted_scores = sorted(scores, key=lambda s: s.overall_score, reverse=True)

        # Top candidates
        top_candidates = []
        for score in sorted_scores[:top_n]:
            card = self.format_score_card(score)
            top_candidates.append({
                "smiles": score.smiles,
                "score": round(score.overall_score, 3),
                "grade": card.overall_grade,
                "summary": card.overall_summary,
                "critical_issues": len(score.critical_issues),
            })

        # Statistics
        all_scores = [s.overall_score for s in scores]
        avg_score = sum(all_scores) / len(all_scores)
        max_score = max(all_scores)
        min_score = min(all_scores)

        # Grade distribution
        grades = {}
        for score in scores:
            if score.overall_score >= 0.8:
                g = "A"
            elif score.overall_score >= 0.7:
                g = "B"
            elif score.overall_score >= 0.6:
                g = "C"
            elif score.overall_score >= 0.5:
                g = "D"
            else:
                g = "F"
            grades[g] = grades.get(g, 0) + 1

        # Common issues
        all_issues = []
        for score in scores:
            all_issues.extend(score.critical_issues)

        issue_counts = {}
        for issue in all_issues:
            issue_counts[issue] = issue_counts.get(issue, 0) + 1

        common_issues = sorted(issue_counts.items(), key=lambda x: x[1], reverse=True)[:5]

        return {
            "total_molecules": len(scores),
            "top_candidates": top_candidates,
            "statistics": {
                "average_score": round(avg_score, 3),
                "max_score": round(max_score, 3),
                "min_score": round(min_score, 3),
            },
            "grade_distribution": grades,
            "common_issues": [{"issue": i, "count": c} for i, c in common_issues],
            "limitations": [
                "All predictions are computational estimates",
                "Experimental validation is required for all candidates",
                "Scores are relative within this batch, not absolute",
            ],
        }


# Global instance
honest_formatter = HonestOutputFormatter()
