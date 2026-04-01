"""Epitope score visualisation.

Produces a two-panel per-residue score plot (BepiPred + DiscoTope)
saved to the specified output path.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd


def plot_epitope_scores(
    df: pd.DataFrame,
    target_name: str,
    output_path: Path,
    interface_df: pd.DataFrame | None = None,
) -> None:
    """Save a multi-panel BepiPred + DiscoTope + (optional GraphBepi) + (optional Interface) plot."""
    has_graphbepi = "graphbepi_score" in df.columns and df["graphbepi_score"].any()
    has_interface = interface_df is not None
    base_panels = 3 if has_graphbepi else 2
    n_panels = base_panels + (1 if has_interface else 0)

    # Interface panel is shorter than predictor panels
    panel_heights = [4] * base_panels + ([1.5] if has_interface else [])

    sns.set_style("whitegrid")
    fig = plt.figure(figsize=(18, sum(panel_heights)))
    gs = fig.add_gridspec(n_panels, 1, hspace=0.5, height_ratios=panel_heights)

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
    if not has_graphbepi and not has_interface:
        ax2.set_xlabel("Residue Position", fontsize=11, fontweight="bold")
    ax2.legend()

    last_ax = ax2

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
        if not has_interface:
            ax3.set_xlabel("Residue Position", fontsize=11, fontweight="bold")
        ax3.legend()
        last_ax = ax3

    # Interface panel (only if interface_df provided)
    if has_interface:
        ax_iface = fig.add_subplot(gs[base_panels], sharex=ax1)
        contact_colors = {"direct": "#E24B4A", "adjacent": "#378ADD"}
        merged = df[["res_id"]].merge(
            interface_df[["res_id", "contact_type"]], on="res_id", how="left"
        )
        colors = merged["contact_type"].map(contact_colors).fillna("#D3D1C7")

        ax_iface.bar(merged["res_id"], 0.4, color=colors.tolist(), width=1.0, align="center")
        ax_iface.set_ylim(0, 1)
        ax_iface.set_yticks([])
        ax_iface.set_ylabel("Interface\nContact", fontsize=8, rotation=0, labelpad=45, va="center")
        ax_iface.set_xlabel("Residue Position", fontsize=11, fontweight="bold")
        legend_handles = [
            mpatches.Patch(color="#E24B4A", label="Direct"),
            mpatches.Patch(color="#378ADD", label="Adjacent"),
            mpatches.Patch(color="#D3D1C7", label="None"),
        ]
        ax_iface.legend(handles=legend_handles, fontsize=8)
        last_ax = ax_iface

    output_path.parent.mkdir(parents=True, exist_ok=True)

    if has_interface:
        save_path = output_path.parent / (output_path.stem + "_with_interface" + output_path.suffix)
    else:
        save_path = output_path

    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
