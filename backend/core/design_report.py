"""
Design Report Generator

Creates actionable, contextual reports for drug design decisions.
This is what a GREAT drug design output should look like.

Design Philosophy:
- Lead with the decision: Is this worth pursuing?
- Provide context: Why do these scores matter for THIS target?
- Be actionable: What should a medicinal chemist do next?
- Be honest: What don't we know?
"""

import logging
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field
from enum import Enum

from .prediction import Prediction, ConfidenceLevel
from .prediction_service import MoleculeScore
from .target_knowledge import TargetKnowledgeBase, TargetType, BindingMechanism

logger = logging.getLogger(__name__)


class Recommendation(Enum):
    """Decision recommendations for a candidate."""
    ADVANCE = "advance"           # Advance to next stage
    OPTIMIZE = "optimize"         # Good potential, needs optimization
    DEPRIORITIZE = "deprioritize" # Significant issues, lower priority
    REJECT = "reject"             # Critical issues, do not pursue


@dataclass
class DesignInsight:
    """A specific insight about a design candidate."""
    category: str           # "strength", "concern", "opportunity", "risk"
    title: str
    description: str
    severity: str = "info"  # "info", "warning", "critical"
    action: str = ""        # What to do about it


@dataclass
class DesignReport:
    """
    Comprehensive design report that answers:
    - Should we advance this candidate?
    - What are the key strengths and risks?
    - What should we do next?
    """
    smiles: str
    candidate_name: str = ""

    # The Decision
    recommendation: Recommendation = Recommendation.OPTIMIZE
    recommendation_rationale: str = ""
    confidence_in_recommendation: str = "moderate"

    # Executive Summary
    one_liner: str = ""
    key_strengths: List[str] = field(default_factory=list)
    key_concerns: List[str] = field(default_factory=list)

    # Scores (simplified for decision-making)
    overall_score: float = 0.0
    overall_grade: str = ""

    # Target-Contextualized Analysis
    target_context: str = ""
    mechanism_implications: List[str] = field(default_factory=list)

    # Detailed Insights
    insights: List[DesignInsight] = field(default_factory=list)

    # Actionable Next Steps (prioritized)
    immediate_actions: List[str] = field(default_factory=list)
    optimization_opportunities: List[str] = field(default_factory=list)
    experimental_priorities: List[str] = field(default_factory=list)

    # Risk Assessment
    go_no_go_factors: Dict[str, Any] = field(default_factory=dict)

    # Honest Limitations
    what_we_dont_know: List[str] = field(default_factory=list)
    assumptions_made: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "smiles": self.smiles,
            "candidate_name": self.candidate_name,
            "decision": {
                "recommendation": self.recommendation.value,
                "rationale": self.recommendation_rationale,
                "confidence": self.confidence_in_recommendation,
            },
            "executive_summary": {
                "one_liner": self.one_liner,
                "key_strengths": self.key_strengths,
                "key_concerns": self.key_concerns,
            },
            "scores": {
                "overall": round(self.overall_score, 3),
                "grade": self.overall_grade,
            },
            "target_analysis": {
                "context": self.target_context,
                "mechanism_implications": self.mechanism_implications,
            },
            "insights": [
                {
                    "category": i.category,
                    "title": i.title,
                    "description": i.description,
                    "severity": i.severity,
                    "action": i.action,
                }
                for i in self.insights
            ],
            "next_steps": {
                "immediate": self.immediate_actions,
                "optimization": self.optimization_opportunities,
                "experiments": self.experimental_priorities,
            },
            "go_no_go": self.go_no_go_factors,
            "limitations": {
                "unknowns": self.what_we_dont_know,
                "assumptions": self.assumptions_made,
            },
        }


class DesignReportGenerator:
    """
    Generates actionable design reports.

    This is what transforms raw predictions into decision-support.
    """

    def generate_report(
        self,
        score: MoleculeScore,
        target_knowledge: Optional[TargetKnowledgeBase] = None,
        candidate_name: str = "",
        comparison_compounds: Optional[List[str]] = None,
    ) -> DesignReport:
        """
        Generate a comprehensive design report.

        Args:
            score: MoleculeScore from prediction service
            target_knowledge: Knowledge about the target (optional)
            candidate_name: Name for this candidate
            comparison_compounds: Compounds to compare against (for novelty)
        """
        report = DesignReport(
            smiles=score.smiles,
            candidate_name=candidate_name,
            overall_score=score.overall_score,
        )

        # 1. Analyze and grade
        report.overall_grade = self._calculate_grade(score)

        # 2. Generate target-contextualized analysis
        if target_knowledge:
            self._add_target_context(report, score, target_knowledge)

        # 3. Extract insights
        self._extract_insights(report, score, target_knowledge)

        # 4. Make recommendation
        self._determine_recommendation(report, score, target_knowledge)

        # 5. Generate actionable next steps
        self._generate_next_steps(report, score, target_knowledge)

        # 6. Document limitations
        self._document_limitations(report, score)

        return report

    def _calculate_grade(self, score: MoleculeScore) -> str:
        """
        Calculate grade with more nuanced criteria.

        A: Excellent - ready for advancement
        B: Good - minor optimization needed
        C: Moderate - significant optimization needed
        D: Poor - major issues, unlikely to advance
        F: Fail - critical issues, reject
        """
        base_score = score.overall_score

        # Adjust for critical issues
        critical_penalty = len(score.critical_issues) * 0.05

        # Adjust for confidence
        conf_modifier = {
            ConfidenceLevel.HIGH: 0.0,
            ConfidenceLevel.MODERATE: 0.0,
            ConfidenceLevel.LOW: -0.05,
            ConfidenceLevel.OUTSIDE_DOMAIN: -0.1,
        }.get(score.overall_confidence, 0)

        adjusted = base_score - critical_penalty + conf_modifier

        # More nuanced grading
        if adjusted >= 0.75 and score.safety_score >= 0.7:
            return "A"
        elif adjusted >= 0.65 and score.safety_score >= 0.6:
            return "B"
        elif adjusted >= 0.55:
            return "C"
        elif adjusted >= 0.45:
            return "D"
        else:
            return "F"

    def _add_target_context(
        self,
        report: DesignReport,
        score: MoleculeScore,
        tk: TargetKnowledgeBase
    ):
        """Add target-specific context to the report."""
        # Build context string
        context_parts = []

        if tk.target_name:
            context_parts.append(f"Target: {tk.target_name}")

        if tk.target_type != TargetType.UNKNOWN:
            context_parts.append(f"Type: {tk.target_type.value}")

        if tk.binding_mechanism != BindingMechanism.UNKNOWN:
            context_parts.append(f"Mechanism: {tk.binding_mechanism.value}")

        report.target_context = " | ".join(context_parts) if context_parts else ""

        # Mechanism implications
        if tk.is_covalent:
            report.mechanism_implications.append(
                f"Covalent inhibitor: Must have appropriate {tk.warhead_types[0] if tk.warhead_types else 'warhead'}"
            )
            if tk.nucleophile_residue:
                report.mechanism_implications.append(
                    f"Targets {tk.nucleophile_residue} for covalent bond formation"
                )
            if tk.staying_portion_features:
                report.mechanism_implications.append(
                    f"Binding portion must include: {', '.join(tk.staying_portion_features[:2])}"
                )

        # Target-specific property requirements
        mw_min, mw_max = tk.get_optimal_mw_range()
        logp_min, logp_max = tk.get_optimal_logp_range()

        # Check if within optimal ranges
        mw_pred = score.predictions.get("molecular_weight")
        logp_pred = score.predictions.get("logp")

        if mw_pred and mw_pred.assessed:
            mw = mw_pred.value
            if mw_min <= mw <= mw_max:
                report.mechanism_implications.append(
                    f"MW ({mw:.0f} Da) within optimal range for {tk.target_name or 'target'}"
                )
            else:
                report.mechanism_implications.append(
                    f"MW ({mw:.0f} Da) outside optimal range ({mw_min}-{mw_max} Da)"
                )

    def _extract_insights(
        self,
        report: DesignReport,
        score: MoleculeScore,
        tk: Optional[TargetKnowledgeBase]
    ):
        """Extract key insights from the scoring."""
        insights = []

        # Strengths
        if score.drug_likeness_score >= 0.6:
            insights.append(DesignInsight(
                category="strength",
                title="Good Drug-likeness",
                description=f"Drug-likeness score of {score.drug_likeness_score:.2f} suggests favorable oral drug properties",
                severity="info",
            ))

        if score.synthesis_score >= 0.8:
            insights.append(DesignInsight(
                category="strength",
                title="Easily Synthesizable",
                description="Low synthetic complexity means rapid analog preparation is feasible",
                severity="info",
                action="Prioritize for synthesis if other properties are acceptable",
            ))

        if score.safety_score >= 0.8:
            insights.append(DesignInsight(
                category="strength",
                title="Clean Safety Profile",
                description="No major safety flags detected in computational assessment",
                severity="info",
            ))

        # Concerns
        if score.safety_score < 0.6:
            insights.append(DesignInsight(
                category="concern",
                title="Safety Concerns",
                description=f"Safety score of {score.safety_score:.2f} indicates potential issues",
                severity="warning",
                action="Run experimental safety assays before advancing",
            ))

        # Check specific predictions
        herg = score.predictions.get("herg_liability_estimate")
        if herg and herg.assessed and herg.value == "high":
            insights.append(DesignInsight(
                category="risk",
                title="hERG Liability Risk",
                description="Predicted high hERG liability (rule-based, LOW confidence)",
                severity="warning",
                action="Consider hERG patch clamp assay. Reduce LogP if possible.",
            ))

        pains = score.predictions.get("pains_alerts")
        if pains and pains.assessed and pains.value and pains.value > 0:
            insights.append(DesignInsight(
                category="concern",
                title="PAINS Alert",
                description=f"{pains.value} PAINS pattern(s) detected - may be promiscuous binder",
                severity="critical",
                action="Investigate specificity. Consider removing problematic substructure.",
            ))

        # Novelty
        if score.novelty_score < 0.2:
            insights.append(DesignInsight(
                category="concern",
                title="Low Novelty",
                description="Very similar to known compounds - patentability may be challenging",
                severity="warning",
                action="Explore scaffold modifications for novelty",
            ))
        elif score.novelty_score > 0.7:
            insights.append(DesignInsight(
                category="opportunity",
                title="Novel Scaffold",
                description="Significantly different from known compounds",
                severity="info",
                action="Validate mechanism of action carefully for novel scaffold",
            ))

        report.insights = insights

        # Populate key strengths/concerns
        report.key_strengths = [
            i.title for i in insights if i.category == "strength"
        ][:3]
        report.key_concerns = [
            i.title for i in insights if i.category in ("concern", "risk")
        ][:3]

    def _determine_recommendation(
        self,
        report: DesignReport,
        score: MoleculeScore,
        tk: Optional[TargetKnowledgeBase]
    ):
        """Determine go/no-go recommendation."""
        # Critical failures
        pains = score.predictions.get("pains_alerts")
        if pains and pains.assessed and pains.value and pains.value >= 2:
            report.recommendation = Recommendation.REJECT
            report.recommendation_rationale = "Multiple PAINS alerts indicate likely promiscuous binding"
            report.confidence_in_recommendation = "high"
            return

        if score.safety_score < 0.4:
            report.recommendation = Recommendation.REJECT
            report.recommendation_rationale = "Serious predicted safety issues"
            report.confidence_in_recommendation = "moderate"
            return

        # Strong candidates
        if score.overall_score >= 0.75 and score.safety_score >= 0.7:
            report.recommendation = Recommendation.ADVANCE
            report.recommendation_rationale = "Strong profile across all dimensions"
            report.confidence_in_recommendation = "moderate"

        elif score.overall_score >= 0.6 and score.safety_score >= 0.6:
            report.recommendation = Recommendation.OPTIMIZE
            report.recommendation_rationale = "Good potential with optimization opportunities"
            report.confidence_in_recommendation = "moderate"

        elif score.overall_score >= 0.5:
            report.recommendation = Recommendation.DEPRIORITIZE
            report.recommendation_rationale = "Significant issues that may be difficult to address"
            report.confidence_in_recommendation = "low"

        else:
            report.recommendation = Recommendation.REJECT
            report.recommendation_rationale = "Overall profile does not support advancement"
            report.confidence_in_recommendation = "low"

        # Generate one-liner
        if report.recommendation == Recommendation.ADVANCE:
            report.one_liner = f"Strong candidate for advancement - {report.key_strengths[0] if report.key_strengths else 'good overall profile'}"
        elif report.recommendation == Recommendation.OPTIMIZE:
            report.one_liner = f"Promising candidate needing optimization - address {report.key_concerns[0] if report.key_concerns else 'identified issues'}"
        elif report.recommendation == Recommendation.DEPRIORITIZE:
            report.one_liner = f"Lower priority - multiple concerns including {', '.join(report.key_concerns[:2])}"
        else:
            report.one_liner = f"Not recommended - {report.recommendation_rationale}"

    def _generate_next_steps(
        self,
        report: DesignReport,
        score: MoleculeScore,
        tk: Optional[TargetKnowledgeBase]
    ):
        """Generate prioritized next steps."""
        immediate = []
        optimization = []
        experiments = []

        # Based on recommendation
        if report.recommendation == Recommendation.ADVANCE:
            immediate.append("Initiate synthesis")
            experiments.append("Primary activity assay")
            experiments.append("Basic ADMET panel (solubility, permeability, metabolic stability)")

        elif report.recommendation == Recommendation.OPTIMIZE:
            # What needs fixing?
            if score.safety_score < 0.7:
                optimization.append("Address safety liabilities before synthesis")

            herg = score.predictions.get("herg_liability_estimate")
            if herg and herg.assessed and herg.value == "high":
                optimization.append("Reduce LogP to mitigate hERG risk (target LogP < 4)")
                experiments.append("hERG patch clamp if advancing")

            if score.novelty_score < 0.3:
                optimization.append("Explore scaffold modifications for improved patentability")

            if score.drug_likeness_score < 0.5:
                optimization.append("Optimize molecular properties (MW, LogP, TPSA)")

            # Default experiments
            experiments.append("Biochemical assay to confirm target engagement")
            experiments.append("Microsomal stability assay")

        elif report.recommendation == Recommendation.DEPRIORITIZE:
            immediate.append("Review if compound fills critical portfolio gap")
            immediate.append("Consider as backup if higher-priority compounds fail")
            optimization.append("Major optimization required across multiple dimensions")

        else:  # REJECT
            immediate.append("Do not synthesize")
            immediate.append("Document learnings for SAR")

        report.immediate_actions = immediate
        report.optimization_opportunities = optimization
        report.experimental_priorities = experiments

        # Go/no-go factors
        report.go_no_go_factors = {
            "advance_if": [
                "Primary assay confirms target activity",
                "ADMET panel shows no major liabilities",
            ],
            "stop_if": [
                "Primary assay IC50 > 1 µM",
                "Major metabolic instability (t1/2 < 15 min)",
                "hERG IC50 < 10 µM",
            ],
        }

    def _document_limitations(self, report: DesignReport, score: MoleculeScore):
        """Document what we don't know."""
        unknowns = []
        assumptions = []

        # General unknowns
        unknowns.append("Actual binding affinity (predictions are similarity-based)")
        unknowns.append("Off-target activity profile")
        unknowns.append("In vivo pharmacokinetics")

        # Confidence-based unknowns
        low_conf_predictions = [
            prop for prop, pred in score.predictions.items()
            if pred.assessed and pred.confidence in (ConfidenceLevel.LOW, ConfidenceLevel.OUTSIDE_DOMAIN)
        ]
        if low_conf_predictions:
            unknowns.append(f"Low-confidence predictions: {', '.join(low_conf_predictions[:3])}")

        # Unassessed properties
        if score.unassessed_properties:
            unknowns.append(f"Could not assess: {', '.join(score.unassessed_properties[:3])}")

        # Assumptions
        assumptions.append("ADMET predictions assume standard oral administration")
        assumptions.append("Binding predictions based on known compound similarity")
        assumptions.append("Synthetic feasibility assumes standard medicinal chemistry resources")

        report.what_we_dont_know = unknowns
        report.assumptions_made = assumptions


# Global instance
report_generator = DesignReportGenerator()
