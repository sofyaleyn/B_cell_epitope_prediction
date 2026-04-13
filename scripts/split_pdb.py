"""Split a multi-chain PDB or mmCIF into one file per chain, with optional grouping.

Without --groups, each chain goes into its own file.
With --groups, comma-separated chains within a group are written together;
any chains not mentioned in any group are each written individually.

Usage:
    # split all chains individually (PDB)
    uv run python scripts/split_pdb.py input.pdb

    # split all chains individually (mmCIF)
    uv run python scripts/split_pdb.py input.cif

    # keep A+B together, split C individually
    uv run python scripts/split_pdb.py input.pdb --groups A,B

    # explicit groups for all chains
    uv run python scripts/split_pdb.py input.pdb --groups A,B C

    # write to a custom directory
    uv run python scripts/split_pdb.py input.pdb --groups A,B --out-dir split/
"""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio import PDB as biopdb


def split_pdb(
    pdb_path: Path,
    groups: list[list[str]] | None,
    out_dir: Path,
) -> list[Path]:
    """Write one file per group of chains from pdb_path into out_dir.

    Supports PDB (.pdb, .ent) and mmCIF (.cif, .mmcif) input; output files
    use the same extension as the input.

    Args:
        pdb_path: input PDB or mmCIF file.
        groups:   list of chain-ID groups; each group is written to one file.
                  Chains not mentioned are each treated as their own group.
                  None means split every chain individually.
        out_dir:  destination directory (created if needed).

    Returns:
        List of paths to the written files.
    """
    suffix = pdb_path.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        parser = biopdb.MMCIFParser(QUIET=True)
        io_cls = biopdb.MMCIFIO
        out_ext = suffix
    else:
        parser = biopdb.PDBParser(QUIET=True)
        io_cls = biopdb.PDBIO
        out_ext = ".pdb"

    structure = parser.get_structure("target", str(pdb_path))

    present = {chain.id for model in structure for chain in model}

    # Build final list of groups, adding ungrouped chains individually
    final_groups: list[list[str]] = []
    if groups:
        grouped = {c for g in groups for c in g}
        missing = grouped - present
        if missing:
            raise ValueError(
                f"Chains {sorted(missing)} not found in {pdb_path.name}. "
                f"Present: {sorted(present)}"
            )
        for g in groups:
            final_groups.append(g)
        for chain_id in sorted(present - grouped):
            final_groups.append([chain_id])
    else:
        for chain_id in sorted(present):
            final_groups.append([chain_id])

    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    for group in final_groups:
        group_set = set(group)
        label = "".join(sorted(group))

        class _Select(biopdb.Select):
            def accept_chain(self, chain):
                return chain.id in group_set
            def accept_atom(self, atom):
                return not atom.is_disordered() or atom.get_altloc() == "A"

        out_path = out_dir / f"{pdb_path.stem}_{label}{out_ext}"
        io = io_cls()
        io.set_structure(structure)
        io.save(str(out_path), select=_Select())
        written.append(out_path)
        print(f"  chains {'+'.join(sorted(group))} → {out_path}")

    return written


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Split a PDB file into one file per chain or group of chains",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  split all chains individually:\n"
            "    split_pdb.py input.pdb\n\n"
            "  keep chains A+B together, split C individually:\n"
            "    split_pdb.py input.pdb --groups A,B\n\n"
            "  explicit groups for all chains:\n"
            "    split_pdb.py input.pdb --groups A,B C\n"
        ),
    )
    parser.add_argument("pdb", type=Path, help="Input PDB (.pdb, .ent) or mmCIF (.cif, .mmcif) file")
    parser.add_argument(
        "--groups", nargs="+", metavar="CHAINS",
        help=(
            "Comma-separated chain IDs to keep together in one file "
            "(e.g. A,B keeps heavy+light chain; remaining chains split individually)"
        ),
    )
    parser.add_argument(
        "--out-dir", type=Path, default=None, metavar="DIR",
        help="Output directory (default: same directory as input PDB)",
    )
    args = parser.parse_args()

    pdb_path = args.pdb.resolve()
    out_dir = args.out_dir.resolve() if args.out_dir else pdb_path.parent
    groups = [g.split(",") for g in args.groups] if args.groups else None

    print(f"Splitting {pdb_path.name} ...")
    split_pdb(pdb_path, groups, out_dir)


if __name__ == "__main__":
    main()
