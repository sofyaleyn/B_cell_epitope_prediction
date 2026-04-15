"""Extract target chains from downloaded mmCIF files using the chains.csv table.

Reads chains.csv produced by fetch_pdbs.py, then for each PDB entry extracts
the target chain from the corresponding *_updated.cif and writes a cleaned
single-chain .pdb file into the output directory.

Output format is always PDB (.pdb) — required by GraphBepi and DiscoTope.

Usage:
    uv run python scripts/split_pdb_target_chain.py --uniprot-id P07992
    uv run python scripts/split_pdb_target_chain.py --uniprot-id P07992 \\
        --raw-pdb-dir data/my_target/raw_pdb --out-dir data/my_target/pdb
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from Bio import PDB as biopdb

REPO_ROOT = Path(__file__).parents[1]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract target chains from mmCIF files using chains.csv"
    )
    parser.add_argument("--uniprot-id", required=True, metavar="ACC",
                        help="UniProt accession (e.g. P07992) — used to derive default paths")
    parser.add_argument("--raw-pdb-dir", default=None, type=Path, metavar="DIR",
                        help="Directory with *_updated.cif files and chains.csv "
                             "(default: data/<uniprot-id>/raw_pdb)")
    parser.add_argument("--out-dir", default=None, type=Path, metavar="DIR",
                        help="Directory for extracted chain files "
                             "(default: data/<uniprot-id>/pdb)")
    args = parser.parse_args()

    raw_dir = args.raw_pdb_dir or REPO_ROOT / "data" / args.uniprot_id / "raw_pdb"
    out_dir = args.out_dir or REPO_ROOT / "data" / args.uniprot_id / "pdb"

    chains_csv = raw_dir / "chains.csv"
    if not chains_csv.exists():
        raise SystemExit(f"chains.csv not found: {chains_csv}")

    out_dir.mkdir(parents=True, exist_ok=True)

    with chains_csv.open(newline="") as f:
        rows = list(csv.DictReader(f))

    mmcif_parser = biopdb.MMCIFParser(QUIET=True)

    for row in rows:
        pdb_id = row["pdb_entry"]
        chain_field = row["chain"].strip()  # e.g. "A" or "A/B"
        pid = pdb_id.lower()
        cif_path = raw_dir / f"{pid}_updated.cif"

        if not cif_path.exists():
            print(f"  MISSING  {cif_path.name} — skipping")
            continue

        # "A/B" → extract each chain letter individually
        chains = [c.strip() for c in chain_field.split("/") if c.strip()]
        if not chains:
            print(f"  SKIP  {pdb_id} — no chain info")
            continue

        print(f"{pdb_id}  chains {'+'.join(chains)}")
        structure = mmcif_parser.get_structure(pid, str(cif_path))

        for chain in chains:
            dest = out_dir / f"{pid}_{chain}.pdb"
            if dest.exists():
                print(f"  skip  {dest.name}")
                continue

            class _Select(biopdb.Select):
                def accept_chain(self, c):
                    return c.id == chain
                def accept_atom(self, atom):
                    return not atom.is_disordered() or atom.get_altloc() == "A"

            io = biopdb.PDBIO()
            io.set_structure(structure)
            try:
                io.save(str(dest), select=_Select())
                print(f"  ok    {dest.name}")
            except Exception as e:
                print(f"  FAIL  {pdb_id} chain {chain}: {e}")
                dest.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
