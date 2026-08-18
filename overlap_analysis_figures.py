import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

IN_PATH = "overlap_analysis_stats.xlsx"
TABLES_OUT = "overlap_summary_clean.xlsx"
FIGURES_DIR = "figures"
SIG_ALPHA = 0.05
LOG_FLOOR = -6

os.makedirs(FIGURES_DIR, exist_ok=True)

df = pd.read_excel(IN_PATH, sheet_name="Summary")

def label_mark_rep(row):
    rep = row.get("replicate")
    if pd.notna(rep):
        return f"{row['mark']} ({rep})"
    return str(row["mark"])


def comparison_kind(name):
    if "Peaks_vs_PI" in name:
        return "Peaks_vs_PI"
    if "Contacts_vs_Loops" in name:
        return "Contacts_vs_Loops"
    if "Loops_vs_PI" in name:
        return "Loops_vs_PI"
    if "Contacts_vs_PI" in name:
        return "Contacts_vs_PI"
    if "Loops_vs_Peaks" in name:
        return "Loops_vs_Peaks"
    if "Contacts_vs_Peaks" in name:
        return "Contacts_vs_Peaks"
    return "Other"


df["kind"] = df["comparison"].apply(comparison_kind)

WINDOW_ORDER = ["exact", "500bp", "1000bp"]
if "window" in df.columns:
    df["window"] = pd.Categorical(df["window"], categories=WINDOW_ORDER, ordered=True)

MARK_ORDER = ["K4", "K9", "K27", "K36", "K92"]
REPLICATE_ORDER = ["FL1", "FL2"]


def plot_sort_key(row):
    rep = row["replicate"]
    rep_rank = REPLICATE_ORDER.index(rep) if rep in REPLICATE_ORDER else 0
    mark_rank = MARK_ORDER.index(row["mark"]) if row["mark"] in MARK_ORDER else 99
    return (rep_rank, mark_rank)

tables = {}

sub = df[df["kind"].isin(["Contacts_vs_Peaks", "Loops_vs_Peaks"])].copy()
sub["label"] = sub.apply(label_mark_rep, axis=1)
sub["sig_both"] = sub["pvalue_both"] < SIG_ALPHA
sub["sig_either"] = sub["pvalue_either"] < SIG_ALPHA
cols = [
    "comparison", "label", "mark", "replicate", "both", "either", "total",
    "pct_both", "pct_either", "fold_enrichment_both", "pvalue_both", "sig_both",
    "fold_enrichment_either", "pvalue_either", "sig_either",
]
tables["ContactsLoops_vs_Peaks"] = sub[cols].sort_values(
    ["comparison", "pvalue_both"]
).reset_index(drop=True)

sub = df[df["kind"].isin(["Contacts_vs_PI", "Loops_vs_PI"])].copy()
sub["sig_both"] = sub["pvalue_both"] < SIG_ALPHA
sub["sig_either"] = sub["pvalue_either"] < SIG_ALPHA
cols = [
    "comparison", "window", "both", "either", "total", "pct_both", "pct_either",
    "fold_enrichment_both", "pvalue_both", "sig_both",
    "fold_enrichment_either", "pvalue_either", "sig_either",
]
tables["ContactsLoops_vs_PI"] = sub[cols].sort_values(
    ["comparison", "window"]
).reset_index(drop=True)

sub = df[df["kind"] == "Peaks_vs_PI"].copy()
sub["label"] = sub.apply(label_mark_rep, axis=1)
sub["sig"] = sub["pvalue"] < SIG_ALPHA
cols = [
    "comparison", "label", "mark", "replicate", "window", "overlapping", "total",
    "pct_overlapping", "fold_enrichment", "pvalue", "sig",
]
tables["Peaks_vs_PI"] = sub[cols].sort_values(
    ["comparison", "pvalue"]
).reset_index(drop=True)

sub = df[df["kind"] == "Contacts_vs_Loops"].copy()
cols = [
    "comparison", "contacts_with_loop_match", "total_contacts", "pct_contacts_matched",
    "loops_with_contact_match", "total_loops", "pct_loops_matched",
]
tables["Contacts_vs_Loops"] = sub[cols].reset_index(drop=True)

with pd.ExcelWriter(TABLES_OUT, engine="openpyxl") as writer:
    for name, t in tables.items():
        t.to_excel(writer, sheet_name=name[:31], index=False)

print(f"Saved {TABLES_OUT} with {len(tables)} cleaned tables.")

SIG_THRESHOLDS = [(0.001, "***"), (0.01, "**"), (0.05, "*")]

def sig_stars(p):
    if pd.isna(p):
        return ""
    for thresh, stars in SIG_THRESHOLDS:
        if p < thresh:
            return stars
    return ""

def transform_fold(values, log_scale):
    values = np.asarray(values, dtype=float)
    if not log_scale:
        return values, np.zeros(len(values), dtype=bool)
    floored = ~(values > 0)
    safe = np.where(values > 0, values, 1.0)
    out = np.where(values > 0, np.log2(safe), 0.0)
    return out, floored


def annotate_stars(ax, bar_info, extra_bottom_pad=0.0):
    ylim = ax.get_ylim()
    span = ylim[1] - ylim[0]
    if span < 0.5:
        ylim = (-1.0, 1.0)
        span = 2.0
        ax.set_ylim(*ylim)
    offset = 0.02 * span if span > 0 else 0.02
    new_top = ylim[1] + 0.12 * span if span > 0 else 1
    new_bottom = ylim[0] - extra_bottom_pad
    ax.set_ylim(new_bottom, new_top)
    for xpos, height, p in bar_info:
        stars = sig_stars(p)
        if stars:
            ax.text(xpos, height + offset, stars, ha="center", va="bottom",
                     fontsize=9, fontweight="bold")


def grouped_fold_bar(ax, sub, x_col, series, title, log_scale=False):
    x = np.arange(len(sub))
    n = len(series)
    width = 0.8 / max(n, 1)
    bar_info = []
    for i, (fold_col, p_col, disp) in enumerate(series):
        raw = sub[fold_col].fillna(0).to_numpy()
        heights, floored = transform_fold(raw, log_scale)
        pvals = sub[p_col].to_numpy() if p_col in sub.columns else [np.nan] * len(sub)
        xpos = x + (i - (n - 1) / 2) * width
        ax.bar(xpos, heights, width, label=disp)
        bar_info.extend(zip(xpos, heights, pvals))
        if log_scale:
            for xp, fl in zip(xpos, floored):
                if fl:
                    ax.text(xp, 0, "0", ha="center", va="top", fontsize=7,
                             style="italic", color="dimgray")
    baseline = 0 if log_scale else 1
    ax.axhline(baseline, color="black", linewidth=0.8, linestyle="--")
    ax.set_xticks(x)
    ax.set_xticklabels(sub[x_col], rotation=45, ha="right", fontsize=8)
    ax.set_title(title, fontsize=10)
    ax.set_ylabel("log2(Fold enrichment)" if log_scale else "Fold enrichment")
    annotate_stars(ax, bar_info)


def peaks_vs_pi_bar(ax, sub, title, log_scale=False):
    fold_pivot = sub.pivot_table(index="label", columns="window", values="fold_enrichment", observed=True)
    p_pivot = sub.pivot_table(index="label", columns="window", values="pvalue", observed=True)
    fold_pivot = fold_pivot[["exact", "500bp", "1000bp"]]
    p_pivot = p_pivot[["exact", "500bp", "1000bp"]]
    x = np.arange(len(fold_pivot))
    width = 0.8 / 3
    bar_info = []
    for i, wcol in enumerate(fold_pivot.columns):
        raw = fold_pivot[wcol].fillna(0).to_numpy()
        heights, floored = transform_fold(raw, log_scale)
        pvals = p_pivot[wcol].to_numpy()
        xpos = x + (i - 1) * width
        ax.bar(xpos, heights, width, label=wcol)
        bar_info.extend(zip(xpos, heights, pvals))
        if log_scale:
            for xp, fl in zip(xpos, floored):
                if fl:
                    ax.text(xp, 0, "0", ha="center", va="top", fontsize=7,
                             style="italic", color="dimgray")
    baseline = 0 if log_scale else 1
    ax.axhline(baseline, color="black", linewidth=0.8, linestyle="--")
    ax.set_xticks(x)
    ax.set_xticklabels(fold_pivot.index, rotation=45, ha="right", fontsize=8)
    ax.set_title(title, fontsize=10)
    ax.set_ylabel("log2(Fold enrichment)" if log_scale else "Fold enrichment")
    annotate_stars(ax, bar_info)


SCALE_SUFFIX = {False: "", True: "  [log2 scale]"}
SCALE_TAG = {False: "", True: "_log2"}


def make_anchor_pair_figure(kind, x_col, sort_for_peaks, name_stub, display_name, log_scale):
    sub_all = df[df["kind"] == kind].copy()
    if sub_all.empty:
        return
    sub_all["label"] = sub_all.apply(label_mark_rep, axis=1) if x_col == "label" else None
    comparisons = sub_all["comparison"].unique()
    n = len(comparisons)
    ncols = 2 if n == 4 else min(3, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4.5 * nrows), squeeze=False)
    for ax, comp in zip(axes.flat, comparisons):
        s = sub_all[sub_all["comparison"] == comp].copy()
        if sort_for_peaks:
            s["_sort_key"] = s.apply(plot_sort_key, axis=1)
            s = s.sort_values("_sort_key")
        grouped_fold_bar(
            ax, s, x_col,
            [
                ("fold_enrichment_both", "pvalue_both", "both anchors"),
                ("fold_enrichment_either", "pvalue_either", "either anchor"),
            ],
            comp, log_scale=log_scale,
        )
    axes.flat[0].legend(fontsize=8)
    for ax in axes.flat[n:]:
        ax.axis("off")
    fig.suptitle(f"Fold enrichment: {display_name}{SCALE_SUFFIX[log_scale]}  (* p<0.05, ** p<0.01, *** p<0.001)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(FIGURES_DIR, f"fold_enrichment_{name_stub}{SCALE_TAG[log_scale]}.png"), dpi=150)
    plt.close(fig)


for log_scale in (False, True):
    make_anchor_pair_figure("Contacts_vs_Peaks", "label", True, "contacts_vs_peaks", "Contacts vs. Peaks", log_scale)
    make_anchor_pair_figure("Loops_vs_Peaks", "label", True, "loops_vs_peaks", "Loops vs. Peaks", log_scale)
    make_anchor_pair_figure("Contacts_vs_PI", "window", False, "contacts_vs_pi", "Contacts vs. Parasitism Islands", log_scale)
    make_anchor_pair_figure("Loops_vs_PI", "window", False, "loops_vs_pi", "Loops vs. Parasitism Islands", log_scale)
    t = tables["Peaks_vs_PI"]
    comparisons = t["comparison"].unique()
    fig, axes = plt.subplots(1, len(comparisons), figsize=(7 * len(comparisons), 5), squeeze=False)
    for idx, comp in enumerate(comparisons):
        ax = axes.flat[idx]
        sub = t[t["comparison"] == comp].copy()
        peaks_vs_pi_bar(ax, sub, comp, log_scale=log_scale)
    axes.flat[0].legend(fontsize=8)
    fig.suptitle(f"Fold enrichment: Peaks vs. Parasitism Islands{SCALE_SUFFIX[log_scale]}  (* p<0.05, ** p<0.01, *** p<0.001)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(os.path.join(FIGURES_DIR, f"fold_enrichment_peaks_vs_pi{SCALE_TAG[log_scale]}.png"), dpi=150)
    plt.close(fig)

print(f"Saved figures to {FIGURES_DIR}/ (linear and log2 versions).")
