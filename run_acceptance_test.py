#!/usr/bin/env python3
"""
Standalone Acceptance Test Runner for KK103

This script validates the drug design platform against acceptance criteria
without requiring external API keys or optional dependencies.

Usage:
    python run_acceptance_test.py
"""

import sys
import os
from pathlib import Path
from datetime import datetime

# Add backend to path
backend_path = Path(__file__).parent / "backend"
sys.path.insert(0, str(backend_path))


def print_header(text):
    print("\n" + "=" * 60)
    print(text)
    print("=" * 60)


def print_result(test_name, passed, details=""):
    status = "PASS" if passed else "FAIL"
    symbol = "[✓]" if passed else "[✗]"
    print(f"  {symbol} {test_name}: {status}")
    if details:
        print(f"      {details}")


def run_acceptance_test():
    """Run the complete acceptance test suite"""
    print_header("Drug Design Platform - Acceptance Test")
    print(f"Test Date: {datetime.now().isoformat()}")
    print(f"Target Compound: KK103 (N-pivaloyl-Leu-Enkephalin)")

    results = {
        "passed": 0,
        "failed": 0,
        "skipped": 0
    }

    # Test 1: Import core modules
    print_header("1. Module Import Tests")

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors
        print_result("RDKit import", True)
        results["passed"] += 1
    except ImportError as e:
        print_result("RDKit import", False, str(e))
        print("\nCRITICAL: RDKit is required. Install with: pip install rdkit")
        results["failed"] += 1
        return results

    try:
        from models import CompoundInput, KK103_REFERENCE, InputFormat
        from services.analogue_generator import AnalogueGenerator
        from services.property_predictor import PropertyPredictor
        from pipeline import PREBUILT_KK103_ANALOGUES
        print_result("Platform modules import", True)
        results["passed"] += 1
    except ImportError as e:
        print_result("Platform modules import", False, str(e))
        results["failed"] += 1
        return results

    # Test 2: KK103 Structure Validation
    print_header("2. KK103 Structure Validation")

    KK103_SMILES = "CC(C)(C)C(=O)NC(Cc1ccc(O)cc1)C(=O)NCC(=O)NCC(=O)NC(Cc2ccccc2)C(=O)NC(CC(C)C)C(=O)O"

    mol = Chem.MolFromSmiles(KK103_SMILES)
    if mol is not None:
        print_result("KK103 SMILES valid", True)
        results["passed"] += 1

        mw = Descriptors.MolWt(mol)
        expected_mw_range = (630, 650)
        mw_valid = expected_mw_range[0] < mw < expected_mw_range[1]
        print_result(f"KK103 MW in range {expected_mw_range}", mw_valid, f"Actual MW: {mw:.2f}")
        results["passed" if mw_valid else "failed"] += 1

        # Check for pivaloyl group
        pivaloyl_pattern = Chem.MolFromSmarts("CC(C)(C)C(=O)N")
        has_pivaloyl = mol.HasSubstructMatch(pivaloyl_pattern)
        print_result("Contains pivaloyl N-terminal", has_pivaloyl)
        results["passed" if has_pivaloyl else "failed"] += 1
    else:
        print_result("KK103 SMILES valid", False, "Could not parse SMILES")
        results["failed"] += 1

    # Test 3: Analogue Generation
    print_header("3. Analogue Generation Tests")

    generator = AnalogueGenerator()
    predictor = PropertyPredictor()

    compound = CompoundInput(
        name="KK103",
        structure=KK103_SMILES,
        format=InputFormat.SMILES,
        target_receptor="DOR"
    )

    # Test N-terminal modifications
    n_term_analogues = generator.generate_analogues(
        compound=compound,
        strategies=["n_terminal"],
        num_analogues=5
    )
    n_term_count = len(n_term_analogues)
    print_result("N-terminal analogues generated", n_term_count > 0, f"Count: {n_term_count}")
    results["passed" if n_term_count > 0 else "failed"] += 1

    # Test amino acid substitutions
    aa_analogues = generator.generate_analogues(
        compound=compound,
        strategies=["amino_acid_substitution"],
        num_analogues=5
    )
    aa_count = len(aa_analogues)
    print_result("AA substitution analogues generated", aa_count > 0, f"Count: {aa_count}")
    results["passed" if aa_count > 0 else "failed"] += 1

    # Test C-terminal modifications
    c_term_analogues = generator.generate_analogues(
        compound=compound,
        strategies=["c_terminal"],
        num_analogues=5
    )
    c_term_count = len(c_term_analogues)
    print_result("C-terminal analogues generated", c_term_count > 0, f"Count: {c_term_count}")
    results["passed" if c_term_count > 0 else "failed"] += 1

    # Combined generation
    all_analogues = generator.generate_analogues(
        compound=compound,
        strategies=["n_terminal", "amino_acid_substitution", "c_terminal"],
        num_analogues=10
    )
    total_count = len(all_analogues)
    min_required = 5
    meets_minimum = total_count >= min_required
    print_result(f"Total analogues ≥ {min_required}", meets_minimum, f"Count: {total_count}")
    results["passed" if meets_minimum else "failed"] += 1

    # Test 4: Structural Diversity
    print_header("4. Structural Diversity Tests")

    from rdkit.Chem import DataStructs

    ref_mol = Chem.MolFromSmiles(KK103_SMILES)
    ref_fp = AllChem.GetMorganFingerprintAsBitVect(ref_mol, 2, nBits=2048)

    diverse_count = 0
    max_similarity_threshold = 0.85

    for analogue in all_analogues:
        try:
            anal_mol = Chem.MolFromSmiles(analogue["smiles"])
            if anal_mol:
                anal_fp = AllChem.GetMorganFingerprintAsBitVect(anal_mol, 2, nBits=2048)
                similarity = DataStructs.TanimotoSimilarity(ref_fp, anal_fp)
                if similarity < max_similarity_threshold:
                    diverse_count += 1
        except:
            pass

    diversity_pass = diverse_count >= min_required
    print_result(f"Analogues with Tanimoto < {max_similarity_threshold}", diversity_pass,
                 f"{diverse_count}/{total_count} pass threshold")
    results["passed" if diversity_pass else "failed"] += 1

    # Test 5: Property Prediction
    print_header("5. Property Prediction Tests")

    # Test KK103 properties
    kk103_props = predictor.predict_properties(KK103_SMILES)
    if kk103_props:
        print_result("KK103 properties predicted", True)
        results["passed"] += 1

        print(f"      MW: {kk103_props.molecular_weight:.2f}")
        print(f"      LogP: {kk103_props.logp:.2f}")
        print(f"      QED: {kk103_props.qed:.3f}")
        print(f"      SA Score: {kk103_props.synthetic_accessibility:.2f}")
    else:
        print_result("KK103 properties predicted", False)
        results["failed"] += 1

    # Test analogue properties
    druglike_count = 0
    synthesizable_count = 0
    min_qed = 0.3
    max_sa = 6.0

    for analogue in all_analogues:
        props = predictor.predict_properties(analogue["smiles"])
        if props:
            if props.qed >= min_qed:
                druglike_count += 1
            if props.synthetic_accessibility <= max_sa:
                synthesizable_count += 1

    druglike_pass = druglike_count >= min_required
    print_result(f"Drug-like analogues (QED ≥ {min_qed})", druglike_pass,
                 f"{druglike_count}/{total_count}")
    results["passed" if druglike_pass else "failed"] += 1

    synthesizable_pass = synthesizable_count >= min_required
    print_result(f"Synthesizable analogues (SA ≤ {max_sa})", synthesizable_pass,
                 f"{synthesizable_count}/{total_count}")
    results["passed" if synthesizable_pass else "failed"] += 1

    # Test 6: Pre-built Analogues Validation
    print_header("6. Pre-built Analogues Validation")

    valid_prebuilt = 0
    for analogue in PREBUILT_KK103_ANALOGUES:
        mol = Chem.MolFromSmiles(analogue["smiles"])
        if mol:
            valid_prebuilt += 1

    prebuilt_pass = valid_prebuilt == len(PREBUILT_KK103_ANALOGUES)
    print_result("All pre-built analogues valid", prebuilt_pass,
                 f"{valid_prebuilt}/{len(PREBUILT_KK103_ANALOGUES)}")
    results["passed" if prebuilt_pass else "failed"] += 1

    # Print pre-built analogues summary
    print("\n  Pre-built KK103 Analogues:")
    for analogue in PREBUILT_KK103_ANALOGUES[:5]:
        print(f"    - {analogue['name']}: {analogue['modification']}")

    # Test 7: Reference Data Completeness
    print_header("7. Reference Data Completeness")

    required_fields = ["name", "sequence", "target_receptor", "plasma_half_life_hours"]
    missing = [f for f in required_fields if f not in KK103_REFERENCE]

    ref_complete = len(missing) == 0
    print_result("KK103 reference data complete", ref_complete,
                 f"Missing: {missing}" if missing else "All fields present")
    results["passed" if ref_complete else "failed"] += 1

    # Final Summary
    print_header("ACCEPTANCE TEST SUMMARY")

    total_tests = results["passed"] + results["failed"]
    pass_rate = (results["passed"] / total_tests * 100) if total_tests > 0 else 0

    print(f"\nTotal Tests: {total_tests}")
    print(f"Passed: {results['passed']}")
    print(f"Failed: {results['failed']}")
    print(f"Pass Rate: {pass_rate:.1f}%")

    # Acceptance decision
    print("\n" + "-" * 60)
    if results["failed"] == 0:
        print("ACCEPTANCE STATUS: PASSED")
        print("\nThe platform meets all acceptance criteria for KK103 analogue generation.")
    else:
        print("ACCEPTANCE STATUS: FAILED")
        print(f"\n{results['failed']} test(s) failed. Review and fix before deployment.")

    print("-" * 60)

    # Generate sample output
    if all_analogues:
        print("\n" + "=" * 60)
        print("SAMPLE OUTPUT: Top 5 Generated Analogues")
        print("=" * 60)

        for i, analogue in enumerate(all_analogues[:5], 1):
            props = predictor.predict_properties(analogue["smiles"])
            print(f"\n{i}. {analogue.get('description', 'Unnamed')}")
            print(f"   Type: {analogue.get('modification_type', 'unknown')}")
            print(f"   SMILES: {analogue['smiles'][:60]}...")
            if props:
                print(f"   MW: {props.molecular_weight:.1f}, QED: {props.qed:.3f}, SA: {props.synthetic_accessibility:.1f}")
            print(f"   Rationale: {analogue.get('rationale', 'N/A')[:80]}")

    return results


if __name__ == "__main__":
    try:
        results = run_acceptance_test()
        sys.exit(0 if results["failed"] == 0 else 1)
    except Exception as e:
        print(f"\nFATAL ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
