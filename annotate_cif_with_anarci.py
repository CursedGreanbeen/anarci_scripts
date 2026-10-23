#!/usr/bin/env python
"""
Annotate antibody sequences using ANARCI to identify chain types.

This script analyzes single-seed base model CIF files and produces a TSV
report mapping chains to antibody (AB) or antigen (AG) status for use with
DockQ analysis.

Usage:
    /opt/conda/bin/python3.10 annotate_cif_with_anarci.py

Output:
    - TSV report: /mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/results/reports/anarci_report.tsv

Report format (TSV with columns):
    - pdb_id: PDB identifier extracted from filename
    - chain: Chain ID (e.g., A, B, C)
    - status: ANTIBODY (H/K/L), NOT_RECOGNIZED (antigen), or NO_ANTIBODY_DOMAIN
"""

import argparse
import glob as globmod
import os
import sys
import traceback
from pathlib import Path

# Add ANARCI library to path for direct Python API access
ANARCI_LIB = Path("/mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/ANARCI/lib/python")
sys.path.insert(0, str(ANARCI_LIB))

# Add HMMER to PATH (required by ANARCI)
HMMER_PATH = "/data/user/mullagaliamova/.conda/envs/hmmer_env/bin"
os.environ["PATH"] = f"{HMMER_PATH}:{os.environ.get('PATH', '')}"

from anarci import anarci
from model_cif_to_fasta import extract_chain_sequences

# Configuration
INPUT_DIR = Path("/mnt/f0771a75-815c-4553-ba74-9aa70f57f995/single-seed")
REPORT_FILE = Path("/mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/results/reports/anarci_report.tsv")
RESIDUE_FILE = REPORT_FILE.with_name("anarci_residues.tsv")

ap = argparse.ArgumentParser()
ap.add_argument("--numbering", action="store_true")
ap.add_argument("--input", default=str(INPUT_DIR / "*/*_model.cif"),
                help="glob of CIF files (default: single-seed model files)")
ap.add_argument("--report", default=str(REPORT_FILE),
                help="chain-level TSV report output path")
ap.add_argument("--residues-out", default=str(RESIDUE_FILE),
                help="per-residue IMGT numbering TSV output path (--numbering)")
args = ap.parse_args()

REPORT_FILE = Path(args.report)
RESIDUE_FILE = Path(args.residues_out)


def extract_pdb_id(filename: str) -> str:
    """Extract PDB ID from filename (e.g., '8iuz.fasta' -> '8iuz')."""
    return Path(filename).stem


def run_anarci_annotation(sequence: str) -> tuple:
    """
    Run ANARCI Python API on a single sequence.

    Args:
        sequence: Amino acid sequence

    Returns:
        (chain_type, num_domains, residues) where:
        - chain_type: 'H', 'K', 'L' for antibody chains, or None if not recognized
        - num_domains: Number of antibody domains found (0 if none)
        - residues: list of (seq_index, imgt_num, ins_code) for the first domain
    """
    sequences = [("seq", sequence)]

    try:
        results = anarci(sequences, scheme="imgt", output=False)
        numbering, alignment_details, hit_tables = results

        if numbering[0] is None:
            return (None, 0, [])
        
        chain_type = alignment_details[0][0].get('chain_type')
        domain_numbering, start_index, end_index = numbering[0][0]

        residues = []                      # (seq_index, imgt_num, ins_code)
        idx = start_index
        for (num, ins), aa in domain_numbering:
            if aa == '-':                  # numbering gap, no residue
                continue
            assert sequence[idx] == aa, f"mapping mismatch at {idx}"
            residues.append((idx, num, ins.strip()))
            idx += 1

        return (chain_type, len(numbering[0]), residues)

    except Exception as e:
        print(f"    ERROR: {e}")
        traceback.print_exc()
        return (None, 0, [])


def imgt_region(n):
    if n <= 26:  return "FR1"
    if n <= 38:  return "CDR1"
    if n <= 55:  return "FR2"
    if n <= 65:  return "CDR2"
    if n <= 104: return "FR3"
    if n <= 117: return "CDR3"
    return "FR4"

H3_ANCHORS = {103, 104, 118, 119}   # flag, not a region


def is_antibody_chain_type(chain_type: str) -> bool:
    """Check if ANARCI chain type is an antibody chain (H, K, L)."""
    return chain_type in ('H', 'K', 'L')


def determine_status(chain_type: str, num_domains: int) -> str:
    """
    Determine the status for a chain based on ANARCI results.

    Returns:
        - 'ANTIBODY' for antibody chains (H, K, L)
        - 'NOT_RECOGNIZED' for non-antibody chains (antigen)
        - 'NO_ANTIBODY_DOMAIN' if ANARCI found no domains
    """
    if chain_type is None:
        return "NO_ANTIBODY_DOMAIN"

    if is_antibody_chain_type(chain_type):
        return "ANTIBODY"

    return "NOT_RECOGNIZED"


def main():
    print("ANARCI Annotation Script")
    print("=" * 60)
    print(f"Input:  {INPUT_DIR}")
    print(f"Report: {REPORT_FILE}")
    print()

    # Create reports directory if needed
    REPORT_FILE.parent.mkdir(parents=True, exist_ok=True)

    # Report data: list of dicts with pdb_id, chain, status
    report_data = []
    residue_rows = []

    # Counters
    total_files = 0
    total_chains = 0
    antibody_chains = 0
    antigen_chains = 0
    unrecognized_chains = 0
    errors = 0

    cif_files = sorted(Path(p) for p in globmod.glob(args.input))
    total_files = len(cif_files)

    for cif_file in cif_files:
        # e.g. '8jow' — subdir name for <pdb>/<model>.cif layout, file stem for flat *.cif
        pdb_id = cif_file.parent.name if cif_file.stem.endswith("_model") else cif_file.stem
        try:
            records = extract_chain_sequences(str(cif_file))   # [(chain_id, seq), ...]
        except Exception as e:
            print(f"  ERROR reading file: {e}")
            errors += 1
            report_data.append({
                "pdb_id": pdb_id,
                "chain": "-",
                "status": "FILE_ERROR"
            })
            continue

        file_has_antibody = False

        for chain_id, sequence in records:                  # was: records.items() / record.sequence
            total_chains += 1
            chain_type, num_domains, residues = run_anarci_annotation(sequence)

            # Determine status
            status = determine_status(chain_type, num_domains)

            if args.numbering:
                numbered = {i: (n, ins) for i, n, ins in residues}   # empty for non-antibody
                for i, aa in enumerate(sequence):
                    if status == "ANTIBODY" and i in numbered:
                        n, ins = numbered[i]
                        region, anchor = imgt_region(n), int(n in H3_ANCHORS)
                        imgt = n
                    elif status == "ANTIBODY":
                        region, anchor, imgt, ins = "OUT", 0, "", ""   # outside the V domain
                    else:
                        region, anchor, imgt, ins = "AG", 0, "", ""
                    residue_rows.append((pdb_id, chain_id, chain_type or "-", i, aa, imgt, ins, region, anchor))

            if status == "ANTIBODY":
                antibody_chains += 1
                file_has_antibody = True
                print(f"  Chain {chain_id}: ANTIBODY (type {chain_type})")
            elif status == "NOT_RECOGNIZED":
                antigen_chains += 1
                print(f"  Chain {chain_id}: NOT_RECOGNIZED (antigen)")
            else:
                unrecognized_chains += 1
                print(f"  Chain {chain_id}: {status}")

            report_data.append({
                "pdb_id": pdb_id,
                "chain": chain_id,
                "status": status
            })

        # If file has no antibody chains, add a summary entry
        if not file_has_antibody and records:
            report_data.append({
                "pdb_id": pdb_id,
                "chain": "-",
                "status": "NO_ANTIBODY_CHAINS"
            })
            print(f"  → No antibody chains in {pdb_id}")

        print()

    # Write report in format expected by mapping.py
    print(f"Writing report: {REPORT_FILE}")
    with open(REPORT_FILE, "w") as f:
        f.write("pdb_id\tchain\tstatus\n")
        for row in report_data:
            f.write(f"{row['pdb_id']}\t{row['chain']}\t{row['status']}\n")

    # Write per-residue IMGT numbering if requested
    if args.numbering:
        print(f"Writing residues: {RESIDUE_FILE}")
        with open(RESIDUE_FILE, "w") as f:
            f.write("pdb_id\tchain\tchain_type\tres_index\taa\timgt\tins\tregion\th3_anchor\n")
            for row in residue_rows:
                f.write("\t".join(str(v) for v in row) + "\n")

    # Summary
    print()
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Files processed:     {total_files}")
    print(f"Total chains:        {total_chains}")
    print(f"  - Antibody (AB):   {antibody_chains}")
    print(f"  - Antigen (AG):    {antigen_chains}")
    print(f"  - Unrecognized:    {unrecognized_chains}")
    print(f"  - Errors:          {errors}")
    print()
    print(f"Report file: {REPORT_FILE}")


if __name__ == "__main__":
    main()
