#!/usr/bin/env python3
"""Find input queries that have ONLY 2024+ hits (no pre-2024 antibody matches)."""

import csv
import os
import sys
import re
from collections import defaultdict

def extract_pdb_id(hit_id):
    """Extract PDB ID from BLAST hit_id."""
    match = re.search(r'pdb\|([A-Z0-9]{4})\|', hit_id)
    return match.group(1).upper() if match else None

def main():
    if len(sys.argv) < 2:
        print("Usage: python extract_2024_only_queries.py <results_dir>")
        sys.exit(1)

    results_dir = sys.argv[1]

    all_hits_file = os.path.join(results_dir, "combined_filtered_hits.csv")
    pre_2024_file = os.path.join(results_dir, "pre_2024_antibody_hits.csv")

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

    # Group all hits by query (source_file + chain_id)
    query_all_hits = defaultdict(set)  # query -> set of all hit PDBs
    query_details = defaultdict(list)  # query -> list of hit details

    with open(all_hits_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            query = f"{row['source_file']}:{row['chain_id']}"
            pdb_id = extract_pdb_id(row['hit_id'])
            if pdb_id:
                query_all_hits[query].add(pdb_id)
                query_details[query].append(row)

    print(f"Total queries with hits: {len(query_all_hits)}")

    # Find queries where ALL hits are 2024+ (none are in pre-2024 set)
    queries_2024_only = []
    for query, hit_pdbs in query_all_hits.items():
        # Check if any hit is in pre-2024 set
        has_pre_2024 = bool(hit_pdbs & pre_2024_pdbs)
        if not has_pre_2024:
            queries_2024_only.append(query)

    print(f"\nQueries with ONLY 2024+ hits (no pre-2024 matches): {len(queries_2024_only)}")

    # Write results
    output_list = os.path.join(results_dir, "queries_2024_only.txt")
    output_csv = os.path.join(results_dir, "queries_2024_only_details.csv")

    with open(output_list, 'w') as f:
        for query in sorted(queries_2024_only):
            f.write(f"{query}\n")

    # Write detailed CSV - get fieldnames from input file
    with open(all_hits_file, 'r') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames

    with open(output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for query in sorted(queries_2024_only):
            for row in query_details[query]:
                writer.writerow(row)

    print(f"\nOutputs:")
    print(f"  Query list: {output_list}")
    print(f"  Detailed hits: {output_csv}")

    # Show breakdown by hit count
    hit_counts = defaultdict(int)
    for query in queries_2024_only:
        count = len(query_all_hits[query])
        hit_counts[count] += 1

    print(f"\nBreakdown by number of 2024+ hits per query:")
    for count in sorted(hit_counts.keys()):
        print(f"  {count} hits: {hit_counts[count]} queries")

if __name__ == '__main__':
    main()
