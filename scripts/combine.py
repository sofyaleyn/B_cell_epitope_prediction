"""Stage 2 — load available predictor results and merge them.

Reads whatever raw CSVs exist under outputs/<target>/ and data/<target>/,
merges them into a single combined_scores.csv.

Usage:
    uv run python scripts/combine.py --target ERCC1
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from epitope_pipeline.integration.combine import load_and_combine


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge predictor results for a target")
    parser.add_argument("--target", required=True, metavar="NAME",
                        help="Target name matching a folder under data/")
    args = parser.parse_args()

    target_dir = REPO_ROOT / "data" / args.target
    if not target_dir.exists():
        raise SystemExit(f"Data directory not found: {target_dir}")

    out_dir = REPO_ROOT / "outputs" / args.target
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{args.target}] Combining results...")
    df = load_and_combine(target_dir, out_dir=out_dir)

    out_path = out_dir / "combined_scores.csv"
    df.to_csv(out_path, index=False)
    print(f"[{args.target}] Saved: {out_path}")

    for col in ["is_epitope_bepipred", "is_epitope_discotope", "is_epitope_graphbepi", "is_epitope_AND"]:
        if col in df.columns:
            n = df[col].sum()
            print(f"  {col}: {n} / {len(df)} ({100*n/len(df):.1f}%)")


if __name__ == "__main__":
    main()
