#!/usr/bin/env python3
"""
Detect water, ions, and small molecules in CIF files.
Generates a report with counts per PDB ID.
"""

import os
import re
from pathlib import Path


# Common ions (monatomic)
ION_CODES = {
    'NA', 'K', 'CL', 'CA', 'MG', 'ZN', 'FE', 'MN', 'CO', 'NI', 'CU', 'HG',
    'PB', 'CD', 'CR', 'AL', 'LI', 'BE', 'BA', 'SR', 'CS', 'RB', 'GA', 'IN',
    'TL', 'GE', 'AS', 'SE', 'BR', 'I', 'AU', 'AG', 'PT', 'PD', 'IR', 'OS',
    'RH', 'RU', 'RE', 'MO', 'W', 'V', 'TI', 'ZR', 'NB', 'TA', 'Hf', 'LA',
    'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM',
    'YB', 'LU', 'Y', 'SC', 'AC', 'TH', 'PA', 'U', 'NP', 'PU', 'AM', 'CM'
}

# Water codes
WATER_CODES = {'HOH', 'WAT', 'WATER', 'DOD', 'H2O'}


def parse_chem_comp_section(content):
    """
    Parse the _chem_comp section to extract compound information.
    Returns dict: {comp_id: {'type': ..., 'name': ..., 'formula': ...}}
    """
    compounds = {}

    # Find _chem_comp loop
    chem_comp_match = re.search(
        r'loop_\s*\n(_chem_comp\.id\s*\n(?:_chem_comp\.\S+\s*\n)*)',
        content
    )

    if not chem_comp_match:
        return compounds

    # Get column headers
    headers_block = chem_comp_match.group(1)
    headers = re.findall(r'_chem_comp\.(\S+)', headers_block)

    # Find data lines after the headers
    after_headers = content[chem_comp_match.end():]

    # Split into lines and collect data until next section
    lines = after_headers.split('\n')
    data_lines = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith('_') or line.startswith('loop_') or line.startswith('#'):
            break
        data_lines.append(line)

    if not data_lines or not headers:
        return compounds

    # Parse each data line
    for data_line in data_lines:
        parts = data_line.split()
        if len(parts) >= len(headers):
            row = dict(zip(headers, parts))
            comp_id = row.get('id', '')
            if comp_id:
                compounds[comp_id] = {
                    'type': row.get('type', ''),
                    'name': row.get('name', ''),
                    'formula': row.get('formula', '')
                }

    return compounds


def count_entities(compounds):
    """
    Count water, ions, and small molecules from parsed compounds.
    """
    water_count = 0
    ion_count = 0
    small_mol_count = 0

    for comp_id, info in compounds.items():
        comp_type = info.get('type', '').lower()
        comp_name = info.get('name', '').upper()

        # Check for water
        if comp_id in WATER_CODES or 'water' in comp_name.lower():
            water_count += 1
        # Check for ions (non-polymer with simple formula or known ion)
        elif comp_id in ION_CODES or 'ion' in comp_name.lower():
            ion_count += 1
        # Check for small molecules (non-polymer that's not water or ion)
        elif 'non-polymer' in comp_type or 'nonpoly' in comp_type:
            small_mol_count += 1
        # Also check polymer types for ligands
        elif 'non-peptide' in comp_type or 'linking' in comp_type:
            small_mol_count += 1

    return water_count, ion_count, small_mol_count


def process_cif_file(file_path):
    """
    Process a single CIF file and return counts.
    """
    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
    except Exception as e:
        return None, None, None

    compounds = parse_chem_comp_section(content)
    water, ions, small_mols = count_entities(compounds)

    return water, ions, small_mols


def main():
    # Paths
    cif_folder = Path('/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/data/CIFs-filtered-new')
    results_folder = Path('/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/results/reports')

    # Create results folder if it doesn't exist
    results_folder.mkdir(parents=True, exist_ok=True)

    # Get all CIF files
    cif_files = sorted(cif_folder.glob('*.cif'))

    # Process each file and collect results
    results = []

    for cif_file in cif_files:
        pdb_id = cif_file.stem  # filename without extension
        water, ions, small_mols = process_cif_file(cif_file)

        if water is not None:
            results.append({
                'pdb_id': pdb_id,
                'water': water,
                'ions': ions,
                'small_molecules': small_mols
            })
            print(f"Processed {pdb_id}: water={water}, ions={ions}, small_mols={small_mols}")
        else:
            print(f"Failed to process {pdb_id}")

    # Write report
    report_path = results_folder / 'solvent_ligands_report.tsv'
    with open(report_path, 'w') as f:
        f.write('pdb_id\twater\tions\tsmall_molecules\n')
        for r in results:
            f.write(f"{r['pdb_id']}\t{r['water']}\t{r['ions']}\t{r['small_molecules']}\n")

    print(f"\nReport written to: {report_path}")
    print(f"Total files processed: {len(results)}")


if __name__ == '__main__':
    main()
