#!/usr/bin/env python3
"""
Acceptance Test for 5HT1A-5HT2A Bivalent Ligand Design

This script validates the drug design platform against the acceptance criteria
for designing novel dual 5HT1A-5HT2A agonist bivalent ligands with CNS penetration.

Usage:
    python run_5ht_bivalent_test.py
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
    """Run the complete 5HT1A-5HT2A bivalent ligand acceptance test suite"""
    print_header("Drug Design Platform - 5HT1A-5HT2A Bivalent Ligand Test")
    print(f"Test Date: {datetime.now().isoformat()}")
    print(f"Target: Novel dual 5HT1A-5HT2A agonist bivalent ligand with CNS penetration")

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
        from services.bivalent_designer import BivalentLigandDesigner
        print_result("Bivalent Designer import", True)
        results["passed"] += 1
    except ImportError as e:
        print_result("Bivalent Designer import", False, str(e))
        results["failed"] += 1
        return results

    # Test 2: Pharmacophore Validation
    print_header("2. Pharmacophore Validation")

    designer = BivalentLigandDesigner()

    # Validate 5-HT1A pharmacophores
    ht1a_pharmacophores = [k for k, v in designer.PHARMACOPHORES.items() if v["target"] == "5-HT1A"]
    ht1a_valid = len(ht1a_pharmacophores) >= 1
    print_result("5-HT1A pharmacophores available", ht1a_valid, f"Found: {ht1a_pharmacophores}")
    results["passed" if ht1a_valid else "failed"] += 1

    # Validate 5-HT2A pharmacophores
    ht2a_pharmacophores = [k for k, v in designer.PHARMACOPHORES.items() if v["target"] == "5-HT2A"]
    ht2a_valid = len(ht2a_pharmacophores) >= 1
    print_result("5-HT2A pharmacophores available", ht2a_valid, f"Found: {ht2a_pharmacophores}")
    results["passed" if ht2a_valid else "failed"] += 1

    # Validate pharmacophore SMILES
    for name, data in designer.PHARMACOPHORES.items():
        mol = Chem.MolFromSmiles(data["smiles"])
        if mol:
            print_result(f"Pharmacophore '{name}' valid", True)
            results["passed"] += 1
        else:
            print_result(f"Pharmacophore '{name}' valid", False, "Invalid SMILES")
            results["failed"] += 1

    # Test 3: Linker Validation
    print_header("3. Linker Validation")

    for linker_name, linker_data in designer.LINKERS.items():
        mol = Chem.MolFromSmiles(linker_data["smiles"])
        linker_valid = mol is not None
        linker_atoms = linker_data.get("atoms", 0)
        print_result(
            f"Linker '{linker_name}' valid",
            linker_valid,
            f"Atoms: {linker_atoms}, Length: {linker_data.get('length_A', 'N/A')}Å"
        )
        results["passed" if linker_valid else "failed"] += 1

    # Test 4: Pre-designed Bivalent Ligands
    print_header("4. Pre-designed 5HT1A-5HT2A Bivalent Ligands")

    designs = designer.get_5ht1a_5ht2a_designs()

    valid_designs = [d for d in designs if d.get("valid", False)]
    print_result(
        "Pre-designed compounds available",
        len(valid_designs) >= 3,
        f"{len(valid_designs)}/{len(designs)} valid"
    )
    results["passed" if len(valid_designs) >= 3 else "failed"] += 1

    # Test 5: Structure Requirements (Two pharmacophores + linker)
    print_header("5. Bivalent Structure Requirements")

    for design in designs[:5]:  # Check first 5
        if not design.get("valid"):
            continue

        mol = Chem.MolFromSmiles(design["smiles"])
        if not mol:
            continue

        # Check for two distinct pharmacophores
        has_p1 = design.get("pharmacophore_1") is not None
        has_p2 = design.get("pharmacophore_2") is not None
        has_linker = design.get("linker") is not None

        structure_valid = has_p1 and has_p2 and has_linker
        print_result(
            f"{design['id']}: Two pharmacophores + linker",
            structure_valid,
            f"P1: {design.get('pharmacophore_1')}, P2: {design.get('pharmacophore_2')}, Linker: {design.get('linker')}"
        )
        results["passed" if structure_valid else "failed"] += 1

    # Test 6: Linker Length (12-24 atoms)
    print_header("6. Linker Length Requirements (12-24 atoms)")

    linker_length_pass = 0
    linker_length_total = 0

    for design in designs:
        if not design.get("valid"):
            continue

        linker_atoms = design.get("linker_atoms", 0)
        linker_length_total += 1

        # Acceptance criteria: 12-24 atoms
        in_range = 12 <= linker_atoms <= 24
        if in_range:
            linker_length_pass += 1

        print_result(
            f"{design['id']}: Linker length",
            in_range,
            f"{linker_atoms} atoms (target: 12-24)"
        )
        results["passed" if in_range else "failed"] += 1

    # Test 7: Molecular Weight Requirements (< 800 Da)
    print_header("7. Molecular Weight Requirements (< 800 Da)")

    mw_pass = 0
    for design in valid_designs:
        props = design.get("properties", {})
        mw = props.get("molecular_weight", 0)
        mw_valid = mw < 800
        if mw_valid:
            mw_pass += 1
        print_result(
            f"{design['id']}: MW < 800",
            mw_valid,
            f"MW: {mw:.1f} Da"
        )
        results["passed" if mw_valid else "failed"] += 1

    # Test 8: LogP Requirements (2-5 for CNS)
    print_header("8. LogP Requirements (2-5 for CNS optimal)")

    logp_pass = 0
    for design in valid_designs:
        props = design.get("properties", {})
        logp = props.get("logP", 0)
        logp_valid = 2 <= logp <= 5
        if logp_valid:
            logp_pass += 1
        print_result(
            f"{design['id']}: LogP 2-5",
            logp_valid,
            f"LogP: {logp:.2f}"
        )
        results["passed" if logp_valid else "failed"] += 1

    # Test 9: CNS/BBB Penetration
    print_header("9. CNS/BBB Penetration Prediction")

    bbb_pass = 0
    for design in valid_designs:
        props = design.get("properties", {})
        bbb_permeable = props.get("bbb_permeable", False)
        if bbb_permeable:
            bbb_pass += 1
        print_result(
            f"{design['id']}: BBB penetration",
            bbb_permeable,
            "BBB+" if bbb_permeable else "BBB-"
        )
        results["passed" if bbb_permeable else "failed"] += 1

    # Test 10: CNS MPO Score
    print_header("10. CNS Multi-Parameter Optimization (MPO) Score")

    for design in valid_designs:
        props = design.get("properties", {})
        cns_mpo = props.get("cns_mpo", 0)
        mpo_valid = cns_mpo >= 3.0  # Good CNS drug candidates have MPO >= 3
        print_result(
            f"{design['id']}: CNS MPO score",
            mpo_valid,
            f"MPO: {cns_mpo:.2f} (target: >= 3.0)"
        )
        results["passed" if mpo_valid else "failed"] += 1

    # Test 11: Drug-likeness
    print_header("11. Drug-likeness Properties")

    for design in valid_designs[:3]:  # Check first 3
        props = design.get("properties", {})
        hbd = props.get("hbd", 0)
        hba = props.get("hba", 0)
        rotatable = props.get("rotatable_bonds", 0)

        # Extended Lipinski for CNS drugs
        # HBD <= 3, HBA <= 7, Rotatable bonds <= 10
        hbd_ok = hbd <= 3
        hba_ok = hba <= 7
        rot_ok = rotatable <= 15  # More lenient for bivalent

        print(f"\n  {design['id']}:")
        print_result("    H-bond donors <= 3", hbd_ok, f"HBD: {hbd}")
        print_result("    H-bond acceptors <= 7", hba_ok, f"HBA: {hba}")
        print_result("    Rotatable bonds <= 15", rot_ok, f"Rotatable: {rotatable}")

        results["passed"] += sum([hbd_ok, hba_ok, rot_ok])
        results["failed"] += sum([not hbd_ok, not hba_ok, not rot_ok])

    # Test 12: Structural Diversity
    print_header("12. Structural Diversity of Designs")

    if len(valid_designs) >= 2:
        from rdkit.DataStructs import TanimotoSimilarity

        ref_mol = Chem.MolFromSmiles(valid_designs[0]["smiles"])
        ref_fp = AllChem.GetMorganFingerprintAsBitVect(ref_mol, 2, nBits=2048)

        diverse_count = 0
        for design in valid_designs[1:]:
            mol = Chem.MolFromSmiles(design["smiles"])
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                similarity = TanimotoSimilarity(ref_fp, fp)
                is_diverse = similarity < 0.9
                if is_diverse:
                    diverse_count += 1
                print_result(
                    f"{design['id']}: Structurally diverse",
                    is_diverse,
                    f"Similarity to lead: {similarity:.2%}"
                )
                results["passed" if is_diverse else "failed"] += 1

    # Test 13: Synthesis Information
    print_header("13. Synthesis Route Information")

    for design in valid_designs[:3]:
        synthesis = design.get("synthesis", {})
        has_synthesis = bool(synthesis)
        print_result(
            f"{design['id']}: Synthesis info available",
            has_synthesis,
            f"Route provided: {synthesis.get('route_summary', 'N/A')[:50] if synthesis else 'None'}..."
        )
        results["passed" if has_synthesis else "failed"] += 1

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
    if pass_rate >= 70:
        print("ACCEPTANCE STATUS: PASSED")
        print(f"\nThe platform meets acceptance criteria for 5HT1A-5HT2A bivalent ligand design.")
    else:
        print("ACCEPTANCE STATUS: FAILED")
        print(f"\n{results['failed']} test(s) failed. Review and fix before deployment.")

    print("-" * 60)

    # Show sample designs
    if valid_designs:
        print("\n" + "=" * 60)
        print("SAMPLE OUTPUT: Top Bivalent Ligand Designs")
        print("=" * 60)

        for i, design in enumerate(valid_designs[:3], 1):
            props = design.get("properties", {})
            print(f"\n{i}. {design['name']} ({design['id']})")
            print(f"   Pharmacophore 1: {design.get('pharmacophore_1')} (5-HT1A)")
            print(f"   Pharmacophore 2: {design.get('pharmacophore_2')} (5-HT2A)")
            print(f"   Linker: {design.get('linker')} ({design.get('linker_atoms', '?')} atoms)")
            print(f"   MW: {props.get('molecular_weight', 'N/A')}, LogP: {props.get('logP', 'N/A')}")
            print(f"   BBB: {'Permeable' if props.get('bbb_permeable') else 'Not permeable'}")
            print(f"   CNS MPO: {props.get('cns_mpo', 'N/A')}")
            print(f"   SMILES: {design['smiles'][:70]}...")
            print(f"   Rationale: {design.get('rationale', 'N/A')[:80]}")

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
