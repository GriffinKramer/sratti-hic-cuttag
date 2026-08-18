import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd


def load_eigenvector(path):
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name", "score", "strand"],
    )
    df["compartment"] = df["score"].apply(
        lambda x: "A" if x > 0 else ("B" if x < 0 else "NA")
    )
    return dfs

#CONFIG: EDIT FOR EACH COMBO
# !!!
balanced_path = "fanc/Dovetail_FLF_eigenvector.bed"
sample_name = "FLF"
output_png = "fanc/Dovetail_FLF_compartments.png"
# !!!

balanced = load_eigenvector(balanced_path)
chromosomes = ["SRv7_Chr1", "SRv7_Chr2", "SRv7_ChrX1", "SRv7_ChrX2"]
colors = {"A": "#e74c3c", "B": "#3498db", "NA": "#cccccc"}

n_chrom = len(chromosomes)

fig, axes = plt.subplots(n_chrom, 1, figsize=(10, 7), sharey=False, squeeze=False)
fig.suptitle(
    f"{sample_name} A/B Compartments",
    fontsize=13,
    fontweight="bold",
)

for row, chrom in enumerate(chromosomes):
    ax = axes[row][0]
    sub = balanced[balanced["chrom"] == chrom]
    for _, r in sub.iterrows():
        ax.barh(
            0,
            r["end"] - r["start"],
            left=r["start"],
            color=colors.get(r["compartment"], "#cccccc"),
            height=0.8,
            edgecolor="none",
        )
    ax.set_yticks([])
    if not sub.empty:
        ax.set_xlim(0, sub["end"].max())
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_ylabel(
        chrom.replace("SRv7_", ""),
        fontsize=9,
        rotation=0,
        labelpad=45,
        va="center",
    )
    if row == n_chrom - 1:
        ax.set_xlabel("Genomic position (bp)", fontsize=9)

legend = [
    mpatches.Patch(color=colors["A"], label="A compartment"),
    mpatches.Patch(color=colors["B"], label="B compartment"),
]
fig.legend(handles=legend, loc="lower right", frameon=False)

plt.tight_layout()
plt.savefig(output_png, dpi=150, bbox_inches="tight")
print(f"Saved to {output_png}")
