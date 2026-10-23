#!/usr/bin/env python3
"""Extract FASTA sequences from AlphaFold3 model CIF files without chain splitting.

Usage:
    /opt/conda/bin/python3.10 model_cif_to_fasta.py <cif_dir> <output_dir>

Note: Requires gemmi, installed for Python 3.10 at /opt/conda/bin/python3.10
"""

import gemmi
from pathlib import Path


def residue_to_one_letter(residue: gemmi.Residue) -> str:
    """Convert residue to one-letter code. Unknown residues become 'X'."""
    info = gemmi.find_tabulated_residue(residue.name)
    if info is None:
        return 'X'
    code = info.one_letter_code
    if not code or code.strip() == '':
        return 'X'
    return code.upper()


def extract_chain_sequences(cif_path: str) -> list[tuple[str, str]]:
    """Extract protein chain sequences from a CIF file.

    Returns:
        List of (chain_id, sequence) tuples.
    """
    structure = gemmi.read_structure(cif_path)
    model = structure[0]

    chains = []
    for chain in model:
        seq_parts = []
        for residue in chain:
            if residue.het_flag == 'A':  # Only standard residues
                seq_parts.append(residue_to_one_letter(residue))

        if seq_parts:
            chains.append((chain.name, ''.join(seq_parts)))

    return chains


def format_fasta_entry(pdb_id: str, chain_id: str, sequence: str) -> str:
    """Format a single FASTA entry."""
    header = f">{pdb_id}|Chain {chain_id}"
    return f"{header}\n{sequence}"


def process_folder(cif_dir: str, output_dir: str) -> None:
    """Process all CIF files in a directory and write FASTA files."""
    cif_dir_p = Path(cif_dir)
    output_dir_p = Path(output_dir)
    output_dir_p.mkdir(parents=True, exist_ok=True)

    cif_files = sorted(cif_dir_p.glob("*.cif"))
    print(f"Found {len(cif_files)} CIF files")

    for cif_path in cif_files:
        pdb_id = cif_path.stem
        try:
            chains = extract_chain_sequences(str(cif_path))
        except Exception as e:
            print(f"[error] Failed to process {cif_path.name}: {e}")
            continue

        if not chains:
            print(f"[warn] {cif_path.name}: no protein chains found")
            continue

        entries = [
            format_fasta_entry(pdb_id, chain_id, seq)
            for chain_id, seq in chains
        ]

        out_path = output_dir_p / f"{pdb_id}.fasta"
        out_path.write_text("\n".join(entries) + "\n")
        print(f"[ok] {cif_path.name} -> {out_path.name} ({len(chains)} chains)")


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <cif_dir> <output_dir>")
        sys.exit(1)

    process_folder(sys.argv[1], sys.argv[2])
