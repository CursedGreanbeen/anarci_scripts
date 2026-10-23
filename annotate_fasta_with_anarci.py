#!/usr/bin/env python
"""
Annotate antibody sequences using ANARCI to identify chain types.

This script analyzes FASTA files and produces a TSV report mapping chains
to antibody (AB) or antigen (AG) status for use with DockQ analysis.

Usage:
    python annotate_fasta_with_anarci.py

Output:
    - TSV report: /mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/results/reports/anarci_report.tsv

Report format (TSV with columns):
    - pdb_id: PDB identifier extracted from filename
    - chain: Chain ID (e.g., A, B, C)
    - status: ANTIBODY (H/K/L), NOT_RECOGNIZED (antigen), or NO_ANTIBODY_DOMAIN
"""

import os
import sys
from pathlib import Path

# Add ANARCI library to path for direct Python API access
ANARCI_LIB = Path("/mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/ANARCI/lib/python")
sys.path.insert(0, str(ANARCI_LIB))

# Add HMMER to PATH (required by ANARCI)
HMMER_PATH = "/data/user/mullagaliamova/.conda/envs/hmmer_env/bin"
os.environ["PATH"] = f"{HMMER_PATH}:{os.environ.get('PATH', '')}"

from anarci import anarci
from fasta_manager import read_fasta_file

# Configuration
INPUT_DIR = Path("/mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/blasted-fasta")
REPORT_FILE = Path("/mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/results/reports/anarci_report.tsv")


def extract_pdb_id(filename: str) -> str:
    """Extract PDB ID from filename (e.g., '8iuz.fasta' -> '8iuz')."""
    return Path(filename).stem


def run_anarci_annotation(sequence: str) -> tuple:
    """
    Run ANARCI Python API on a single sequence.

    Args:
        sequence: Amino acid sequence

    Returns:
        (chain_type, num_domains) where:
        - chain_type: 'H', 'K', 'L' for antibody chains, or None if not recognized
        - num_domains: Number of antibody domains found (0 if none)
    """
    sequences = [("seq", sequence)]

    try:
        results = anarci(sequences, scheme="imgt", output=False)
        numbering, alignment_details, hit_tables = results

        if numbering[0] is None:
            return (None, 0)

        # Get the first (most significant) domain
        domain_numbering, start_index, end_index = numbering[0][0]

        # Get chain type from alignment details
        if alignment_details[0]:
            chain_type = alignment_details[0][0].get('chain_type')
            return (chain_type, len(numbering[0]))

        return (None, 0)

    except Exception as e:
        print(f"    ERROR: {e}")
        return (None, 0)


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

    # Counters
    total_files = 0
    total_chains = 0
    antibody_chains = 0
    antigen_chains = 0
    unrecognized_chains = 0
    errors = 0

    # Process each FASTA file
    fasta_files = sorted(INPUT_DIR.glob("*.fasta"))
    total_files = len(fasta_files)

    print(f"Found {total_files} FASTA files\n")

    for fasta_file in fasta_files:
        print(f"Processing: {fasta_file.name}")
        pdb_id = extract_pdb_id(fasta_file.name)

        try:
            records = read_fasta_file(fasta_file)
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

        for chain_id, record in records.items():
            total_chains += 1

            # Run ANARCI on ALL chains - it's the authoritative source
            chain_type, num_domains = run_anarci_annotation(record.sequence)

            # Determine status
            status = determine_status(chain_type, num_domains)

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
