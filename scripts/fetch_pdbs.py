"""Download all PDB structures associated with a UniProt entry.

Usage:
    uv run python scripts/fetch_pdbs.py --uniprot-id P07992
    uv run python scripts/fetch_pdbs.py --uniprot-id P07992 --out-dir data/my_target/pdb
"""
from __future__ import annotations

import argparse
import csv
import json
import urllib.request
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlretrieve

REPO_ROOT = Path(__file__).parents[1]

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{uid}.json"
FASTA_URL   = "https://www.ebi.ac.uk/pdbe/entry-files/download/{id}.fasta"
MMCIF_URL   = "https://www.ebi.ac.uk/pdbe/entry-files/download/{id}_updated.cif"


def fetch_pdb_entries(uniprot_id: str) -> list[tuple[str, str]]:
    """Return list of (pdb_id, chain) for all PDB cross-references in UniProt.

    Chain is parsed from the 'Chains' property, e.g. 'A=1-100' → 'A',
    'A/B=50-200' → 'A/B'. Empty string when the property is absent.
    """
    url = UNIPROT_URL.format(uid=uniprot_id)
    with urllib.request.urlopen(url) as r:
        data = json.load(r)
    entries = []
    for ref in data.get("uniProtKBCrossReferences", []):
        if ref.get("database") != "PDB":
            continue
        pdb_id = ref["id"]
        chains_val = next(
            (p["value"] for p in ref.get("properties", []) if p["key"] == "Chains"),
            "",
        )
        # "A=1-100" or "A/B=50-200" → keep only the part before "="
        chain = chains_val.split("=")[0] if chains_val else ""
        entries.append((pdb_id, chain))
    return entries


def download(url: str, dest: Path) -> bool:
    """Download url → dest. Return True if downloaded, False if skipped."""
    if dest.exists():
        print(f"  skip  {dest.name}")
        return False
    try:
        urlretrieve(url, dest)
        print(f"  ok    {dest.name}")
        return True
    except (HTTPError, URLError) as e:
        print(f"  FAIL  {dest.name}  ({e})")
        dest.unlink(missing_ok=True)
        return False


def main() -> None:
    parser = argparse.ArgumentParser(description="Download PDB structures for a UniProt entry")
    parser.add_argument("--uniprot-id", required=True, metavar="ACC",
                        help="UniProt accession (e.g. P07992)")
    parser.add_argument("--out-dir", default=None, type=Path, metavar="DIR",
                        help="Directory to save downloaded files (default: data/<uniprot-id>/raw_pdb)")
    args = parser.parse_args()

    out_dir = args.out_dir or REPO_ROOT / "data" / args.uniprot_id / "raw_pdb"
    out_dir.mkdir(parents=True, exist_ok=True)

    entries = fetch_pdb_entries(args.uniprot_id)
    print(f"{len(entries)} PDB entries for {args.uniprot_id}: {', '.join(e[0] for e in entries)}")

    for pdb_id, chain in entries:
        pid = pdb_id.lower()
        download(FASTA_URL.format(id=pid), out_dir / f"{pid}.fasta")
        download(MMCIF_URL.format(id=pid), out_dir / f"{pid}_updated.cif")

    chains_csv = out_dir / "chains.csv"
    with chains_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["pdb_entry", "chain"])
        writer.writerows(entries)
    print(f"chains table → {chains_csv}")


if __name__ == "__main__":
    main()
