#!/usr/bin/env python3
"""Extract PDB IDs that are 2024+ only (not in pre-2024 antibody hits)."""

import csv
import os
import sys
import re

def extract_pdb_id(hit_id):
    """Extract PDB ID from BLAST hit_id."""
    match = re.search(r'pdb\|([A-Z0-9]{4})\|', hit_id)
    return match.group(1).upper() if match else None

def main():
    if len(sys.argv) < 2:
        print("Usage: python extract_2024_only_pdbs.py <results_dir>")
        print("  results_dir: Directory containing combined_filtered_hits.csv and pre_2024_antibody_hits.csv")
        sys.exit(1)

    results_dir = sys.argv[1]

    all_hits_file = os.path.join(results_dir, "combined_filtered_hits.csv")
    pre_2024_file = os.path.join(results_dir, "pre_2024_antibody_hits.csv")
    output_file = os.path.join(results_dir, "2024_only_pdbs.txt")
    detailed_output = os.path.join(results_dir, "2024_only_hits.csv")

    print(f"Processing: {results_dir}")

    # Load pre-2024 PDB IDs
    pre_2024_pdbs = set()
    with open(pre_2024_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            pdb_id = row.get('pdb_id', '').upper()
            if pdb_id:
                pre_2024_pdbs.add(pdb_id)

    print(f"Pre-2024 antibody PDBs: {len(pre_2024_pdbs)}")

    # Load all hits and find 2024+ only
    all_pdbs = set()
    pdb_to_details = {}  # pdb_id -> list of (hit_id, description)

    with open(all_hits_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            pdb_id = extract_pdb_id(row['hit_id'])
            if pdb_id:
                all_pdbs.add(pdb_id)
                if pdb_id not in pdb_to_details:
                    pdb_to_details[pdb_id] = []
                pdb_to_details[pdb_id].append({
                    'hit_id': row['hit_id'],
                    'description': row['hit_description']
                })

    print(f"Total unique PDBs: {len(all_pdbs)}")

    # Subtract pre-2024 from all
    only_pdbs = sorted(all_pdbs - pre_2024_pdbs)

    print(f"2024+ only (or no date match) PDBs: {len(only_pdbs)}")

    # Write simple list
    with open(output_file, 'w') as f:
        for pdb_id in only_pdbs:
            f.write(f"{pdb_id}\n")

    # Write detailed CSV with hit info
    with open(detailed_output, 'w', newline='') as f:
        # Get fieldnames from original file
        with open(all_hits_file, 'r') as orig:
            reader = csv.DictReader(orig)
            fieldnames = reader.fieldnames

        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        with open(all_hits_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                pdb_id = extract_pdb_id(row['hit_id'])
                if pdb_id and pdb_id in only_pdbs:
                    writer.writerow(row)

    print(f"\nOutputs:")
    print(f"  PDB list: {output_file}")
    print(f"  Detailed hits: {detailed_output}")

    return only_pdbs

if __name__ == '__main__':
    main()
