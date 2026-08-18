import os
import warnings

import matplotlib
import pandas as pd
from upsetplot import UpSet, from_memberships

matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=FutureWarning)

INPUT_XLSX = "functional_annotation/peaks_unfiltered_PFAM_by_dataset.xlsx"
SHEET_NAME = "PFAM_all"
OUT_DIR = "functional_annotation"
SIG_ALPHA_PLOT = 0.05
TITLE = "PFAM enrichment of genes overlapping peaks"

os.makedirs(OUT_DIR, exist_ok=True)

def build_upset_plot(pfam_all, sig_alpha, title, out_path):
    sig = pfam_all[pfam_all["padj_bh"] < sig_alpha]
    if sig.empty:
        print(f"No PFAM domains significant at padj_bh < {sig_alpha} -- skipping plot.")
        return

    memberships = []
    labels = []
    for domain, group in sig.groupby("PFAM"):
        memberships.append(sorted(group["dataset"].unique().tolist()))
        labels.append(domain)

    print(f"{len(memberships)} significant domains across {sig['dataset'].nunique()} datasets")

    data = from_memberships(memberships, data=labels)
    data = data.rename(columns={0: "domain"})

    n_sets = sig["dataset"].nunique()
    fig = plt.figure(figsize=(max(14, 1.1 * n_sets + 6), 8))
    upset = UpSet(data, subset_size="count", sort_by="cardinality", show_counts=True)
    upset.plot(fig=fig)
    fig.suptitle(f"{title} (padj_bh < {sig_alpha})", y=1.02)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_path}")

if __name__ == "__main__":
    pfam_all = pd.read_excel(INPUT_XLSX, sheet_name=SHEET_NAME)
    stem = os.path.splitext(os.path.basename(INPUT_XLSX))[0]
    out_path = os.path.join(OUT_DIR, f"{stem}_upset.png")
    build_upset_plot(pfam_all, SIG_ALPHA_PLOT, TITLE, out_path)
