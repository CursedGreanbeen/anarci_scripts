#!/usr/bin/env python3
"""
Extract unique queries (max h3_pct_identity < 90%) from pre_2024_antibody_hits.csv
in each result subdirectory and output to pre_2024_unique.csv in each folder

A query is "unique" only if ALL its hits have h3_pct_identity < 90%.
For each unique query, keep the hit with the HIGHEST identity (best match).
"""

import csv
from collections import defaultdict
from pathlib import Path

RESULTS_DIR = Path("/mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/results")

SUBDIRS = [
    "fasta-fv-only-0-100",
    "fasta-fv-only-cdr-80-cutoff",
    "fasta-fv-only-levenshtein"
]

def process_subdirectory(subdir):
    """Extract unique queries from a single subdirectory."""

    csv_path = RESULTS_DIR / subdir / "pre_2024_antibody_hits.csv"
    output_path = RESULTS_DIR / subdir / "pre_2024_unique.csv"

    if not csv_path.exists():
        print(f"Warning: {csv_path} not found, skipping...")
        return

    print(f"Processing {csv_path}...")

    # Group all hits by query (source_file + chain_id)
    query_hits = defaultdict(list)
    fieldnames = None

    with open(csv_path, 'r', newline='') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames

        for row in reader:
            try:
                h3_pct_identity = float(row['h3_pct_identity'])
            except (ValueError, KeyError) as e:
                continue

            key = (row['source_file'], row['chain_id'])
            query_hits[key].append((h3_pct_identity, row))

    # For each query, check if ALL hits have h3_pct < 90%
    # If so, keep the hit with the HIGHEST identity (best match)
    unique_rows = {}

    for key, hits in query_hits.items():
        max_h3_pct = max(h[0] for h in hits)

        # Only keep query if ALL hits are < 90% (i.e., max < 90%)
        if max_h3_pct < 90:
            # Keep the hit with highest identity (best match among low-identity hits)
            best_hit = max(hits, key=lambda x: x[0])
            unique_rows[key] = best_hit[1]

    # Add pdb_id column to output (extracted from source_file, not hit)
    # release_date already exists in the data - preserve it
    output_fieldnames = list(fieldnames) + ['pdb_id'] if 'pdb_id' not in fieldnames else fieldnames

    print(f"  Found {len(unique_rows)} unique queries (all hits < 90%)")
    print(f"  Writing to {output_path}...")

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=output_fieldnames)
        writer.writeheader()
        for row in unique_rows.values():
            # Extract PDB ID from source_file (e.g., "7sn1.fasta" -> "7SN1")
            source_file = row.get('source_file', '')
            row['pdb_id'] = source_file.replace('.fasta', '').upper() if source_file else ''
            writer.writerow(row)

    print(f"  Done!")

def main():
    for subdir in SUBDIRS:
        process_subdirectory(subdir)
    print("\nAll done!")

if __name__ == "__main__":
    main()
