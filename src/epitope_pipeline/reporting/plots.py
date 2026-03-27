"""Epitope score visualisation.

Produces a two-panel per-residue score plot (BepiPred + DiscoTope)
saved to the specified output path.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def plot_epitope_scores(
    df: pd.DataFrame,
    target_name: str,
    output_path: Path,
) -> None:
    """Save a three-panel BepiPred + DiscoTope + GraphBepi score plot to output_path."""
    has_graphbepi = "graphbepi_score" in df.columns and df["graphbepi_score"].any()
    n_panels = 3 if has_graphbepi else 2

    sns.set_style("whitegrid")
    fig = plt.figure(figsize=(18, 4 * n_panels))
    gs  = fig.add_gridspec(n_panels, 1, hspace=0.5)

    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)

    bp_epitope = df.get("is_epitope_bepipred", False)
    dt_epitope = df.get("is_epitope_discotope", False)

    # BepiPred panel
    ax1.plot(df["res_id"], df["bepipred_score"], color="mediumpurple", linewidth=2, label="Score")
    ax1.axhline(y=0.5, color="gray", linestyle="--", linewidth=1, label="Threshold (0.50)")
    ax1.fill_between(df["res_id"], 0, df["bepipred_score"],
                     where=bp_epitope, alpha=0.8, color="lavender", label="Epitope")
    ax1.set_ylabel("BepiPred Score", fontsize=11, fontweight="bold")
    ax1.set_title(f"Linear B-cell Epitopes (BepiPred 3.0) — {target_name}", fontsize=13, fontweight="bold")
    ax1.set_ylim(-0.05, max(0.65, df["bepipred_score"].max() * 1.1))
    ax1.legend()

    # DiscoTope panel
    ax2.plot(df["res_id"], df["discotope_score"], color="green", linewidth=2, label="Score")
    ax2.axhline(y=0.9, color="gray", linestyle="--", linewidth=1, label="Threshold (0.90)")
    ax2.fill_between(df["res_id"], 0, df["discotope_score"],
                     where=dt_epitope, alpha=0.4, color="mediumseagreen", label="Epitope")
    ax2.set_ylabel("DiscoTope Score", fontsize=11, fontweight="bold")
    ax2.set_title(f"Conformational B-cell Epitopes (DiscoTope 3.0) — {target_name}", fontsize=13, fontweight="bold")
    if not has_graphbepi:
        ax2.set_xlabel("Residue Position", fontsize=11, fontweight="bold")
    ax2.legend()

    # GraphBepi panel (only if scores present)
    if has_graphbepi:
        ax3 = fig.add_subplot(gs[2], sharex=ax1)
        gb_epitope = df.get("is_epitope_graphbepi", False)
        ax3.plot(df["res_id"], df["graphbepi_score"], color="steelblue", linewidth=2, label="Score")
        ax3.axhline(y=0.1763, color="gray", linestyle="--", linewidth=1, label="Threshold (0.1763)")
        ax3.fill_between(df["res_id"], 0, df["graphbepi_score"],
                         where=gb_epitope, alpha=0.4, color="lightsteelblue", label="Epitope")
        ax3.set_ylabel("GraphBepi Score", fontsize=11, fontweight="bold")
        ax3.set_title(f"Conformational B-cell Epitopes (GraphBepi) — {target_name}", fontsize=13, fontweight="bold")
        ax3.set_xlabel("Residue Position", fontsize=11, fontweight="bold")
        ax3.legend()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
