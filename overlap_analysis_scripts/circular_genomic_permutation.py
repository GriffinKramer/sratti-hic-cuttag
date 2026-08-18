# imports
import os
from collections import defaultdict
import glob
import numpy as np
import pandas as pd

# config
IN_PATH = "overlap_analysis.xlsx"
OUT_PATH = "overlap_analysis_stats.xlsx"

ISLANDS_PATH = "genome/sratti_parasitism_islands.bed"
CHROM_SIZES_PATH = "genome/chrom.sizes"
WINDOWS = {"exact": 0, "500bp": 500, "1000bp": 1000}
Q_THRESHOLD = 0.01

N_PERMUTATIONS = 10000
SEED = 42

# fp
CONTACT_SOURCES = {
    "ArimaFLF": {
        "path": "fithic/ArimaFLF/ArimaFLF_10kb.spline_pass2.res10000.significances.txt.gz",
        "half_window": 5000,
    },
    "ArimaPF": {
        "path": "fithic/ArimaPF/ArimaPF_10kb.spline_pass2.res10000.significances.txt.gz",
        "half_window": 5000,
    },
    "DovetailFLF": {
        "path": "fithic/DovetailFLF/DovetailFLF_10kb.spline_pass2.res10000.significances.txt.gz",
        "half_window": 5000,
    },
    "DovetailPF": {
        "path": "fithic/DovetailPF/DovetailPF_10kb.spline_pass2.res10000.significances.txt.gz",
        "half_window": 5000,
    },
}

LOOP_SOURCES = {
    "DovetailFLF": "cooltools/DovetailFLF_dots.bedpe",
    "DovetailPF": "cooltools/DovetailPF_dots.bedpe",
}

SEACR_DIR = "SEACR"
PEAK_SAMPLES = {
    ("FL", "FL1", "K4"): "FLK4",
    ("FL", "FL1", "K9"): "FLK9",
    ("FL", "FL1", "K27"): "FLK27",
    ("FL", "FL1", "K36"): "FLK36",
    ("FL", "FL1", "K92"): "FLK92",
    ("FL", "FL2", "K4"): "FL2K4",
    ("FL", "FL2", "K9"): "FL2K9",
    ("FL", "FL2", "K27"): "FL2K27",
    ("FL", "FL2", "K36"): "FL2K36",
    ("FL", "FL2", "K92"): "FL2K92",
    ("PF", None, "K4"): "PFK4",
    ("PF", None, "K9"): "PFK9",
    ("PF", None, "K27"): "PFK27",
    ("PF", None, "K36"): "PFK36",
    ("PF", None, "K92"): "PFK92",
}

rng = np.random.default_rng(SEED)

# loading
def find_seacr_file(prefix):
    matches = glob.glob(f"{SEACR_DIR}/{prefix}_*_stringent.stringent.bed")
    if len(matches) != 1:
        raise ValueError(
            f"Expected 1 SEACR file for '{prefix}', found {len(matches)}: {matches}"
        )
    return matches[0]

def load_peaks(prefix):
    return pd.read_csv(
        find_seacr_file(prefix),
        sep="\t",
        header=None,
        names=["chr", "start", "end", "auc", "max", "region"],
    )

def load_islands():
    islands = pd.read_csv(
        ISLANDS_PATH,
        sep=r"\s+",
        header=None,
        names=["chr", "start", "end", "genes"],
        dtype={"start": int, "end": int},
    )
    print(f"Loaded {len(islands)} parasitism islands")
    return islands

def load_chrom_sizes():
    sizes = pd.read_csv(
        CHROM_SIZES_PATH, sep="\t", header=None, names=["chrom", "length"]
    )
    return dict(zip(sizes["chrom"], sizes["length"]))

def load_contacts(source_name):
    cfg = CONTACT_SOURCES[source_name]
    df = pd.read_csv(cfg["path"], sep="\t")
    sig = df[df["q-value"] < Q_THRESHOLD].copy()
    hw = cfg["half_window"]
    sig["start1"] = sig["fragmentMid1"] - hw
    sig["end1"] = sig["fragmentMid1"] + hw
    sig["start2"] = sig["fragmentMid2"] - hw
    sig["end2"] = sig["fragmentMid2"] + hw
    print(f"{source_name}: {len(sig)} significant contacts (q<{Q_THRESHOLD})")
    return sig

def load_loops(source_name):
    df = pd.read_csv(LOOP_SOURCES[source_name], sep="\t")
    df = df.rename(columns={"chrom1": "chr1", "chrom2": "chr2"})
    print(f"{source_name}: {len(df)} loops")
    return df

def to_arrays(df, chr_col="chr", start_col="start", end_col="end"):
    return (
        df[chr_col].to_numpy(),
        df[start_col].to_numpy(dtype=np.int64),
        df[end_col].to_numpy(dtype=np.int64),
    )

def to_two_anchor_arrays(df):
    return (
        df["chr1"].to_numpy(),
        df["start1"].to_numpy(dtype=np.int64),
        df["end1"].to_numpy(dtype=np.int64),
        df["chr2"].to_numpy(),
        df["start2"].to_numpy(dtype=np.int64),
        df["end2"].to_numpy(dtype=np.int64),
    )

def build_index_from_arrays(chrom_arr, start_arr, end_arr):
    index = {}
    for c in np.unique(chrom_arr):
        mask = chrom_arr == c
        s, e = start_arr[mask], end_arr[mask]
        order = np.argsort(s, kind="mergesort")
        s, e = s[order], e[order]
        index[c] = (s, e)
    return index

def overlap_batch(starts, ends, q_starts, q_ends):
    if len(starts) == 0:
        return np.zeros(len(q_starts), dtype=bool)
    pos = np.searchsorted(starts, q_ends, side="left")
    has_candidate = pos > 0
    candidate_idx = np.clip(pos - 1, 0, len(starts) - 1)
    return has_candidate & (ends[candidate_idx] > q_starts)

def batched_overlap(chrom_arr, start_arr, end_arr, index):
    out = np.zeros(len(chrom_arr), dtype=bool)
    for c in np.unique(chrom_arr):
        mask = chrom_arr == c
        entry = index.get(c)
        if entry is None:
            continue
        idx_s, idx_e = entry
        out[mask] = overlap_batch(idx_s, idx_e, start_arr[mask], end_arr[mask])
    return out

def circular_shuffle_arrays(chrom_arr, start_arr, end_arr, chrom_sizes, rng):
    if len(chrom_arr) == 0:
            return chrom_arr, start_arr, end_arr
    out_chrom, out_start, out_end = [], [], []
    for c in np.unique(chrom_arr):
        mask = chrom_arr == c
        s, e = start_arr[mask], end_arr[mask]
        w = e - s
        length = chrom_sizes.get(c)
        if length is None or length <= 0:
            out_chrom.append(np.full(mask.sum(), c))
            out_start.append(s)
            out_end.append(e)
            continue
        offset = rng.integers(0, length)
        shifted_s = (s + offset) % length
        shifted_e = shifted_s + w
        wraps = shifted_e > length
        out_chrom.append(np.full(mask.sum(), c))
        out_start.append(shifted_s)
        out_end.append(np.where(wraps, length, shifted_e))
        if wraps.any():
            out_chrom.append(np.full(int(wraps.sum()), c))
            out_start.append(np.zeros(int(wraps.sum()), dtype=np.int64))
            out_end.append(shifted_e[wraps] - length)
    return np.concatenate(out_chrom), np.concatenate(out_start), np.concatenate(out_end)

def extend_arrays(chrom_arr, start_arr, end_arr, window_bp, chrom_sizes):
    if window_bp <= 0:
        return start_arr, end_arr
    new_start = np.clip(start_arr - window_bp, 0, None)
    lengths = np.array([chrom_sizes.get(c, -1) for c in chrom_arr], dtype=np.int64)
    default_end = end_arr + window_bp
    new_end = np.where(lengths > 0, np.minimum(default_end, lengths), default_end)
    return new_start.astype(np.int64), new_end.astype(np.int64)

def moving_index_from_peaks(peaks_arrays, chrom_sizes, rng):
    chrom_arr, start_arr, end_arr = peaks_arrays
    c2, s2, e2 = circular_shuffle_arrays(chrom_arr, start_arr, end_arr, chrom_sizes, rng)
    return build_index_from_arrays(c2, s2, e2)

def moving_index_from_islands(island_arrays, window_bp, chrom_sizes, rng):
    chrom_arr, start_arr, end_arr = island_arrays
    c2, s2, e2 = circular_shuffle_arrays(chrom_arr, start_arr, end_arr, chrom_sizes, rng)
    s2, e2 = extend_arrays(c2, s2, e2, window_bp, chrom_sizes)
    return build_index_from_arrays(c2, s2, e2)

def permutation_stats(observed, null_values):
    null_values = np.array(null_values)
    expected = null_values.mean()
    std = null_values.std()
    fold = (
        observed / expected
        if expected > 0
        else float("inf")
        if observed > 0
        else float("nan")
    )
    z = (observed - expected) / std if std > 0 else float("nan")
    p_over = (np.sum(null_values >= observed) + 1) / (len(null_values) + 1)
    p_under = (np.sum(null_values <= observed) + 1) / (len(null_values) + 1)
    p = min(2 * min(p_over, p_under), 1.0)
    return (
        round(expected, 2),
        round(fold, 2) if np.isfinite(fold) else fold,
        round(z, 2) if np.isfinite(z) else z,
        round(p, 4),
    )

def permute_vs_peaks(twoanchor_arrays, peaks_arrays, chrom_sizes, n_perm, rng):
    c1, s1, e1, c2, s2, e2 = twoanchor_arrays
    both_null = np.empty(n_perm, dtype=np.int64)
    either_null = np.empty(n_perm, dtype=np.int64)
    for i in range(n_perm):
        index = moving_index_from_peaks(peaks_arrays, chrom_sizes, rng) #circular shuffle
        o1 = batched_overlap(c1, s1, e1, index) # check overlap
        o2 = batched_overlap(c2, s2, e2, index) # check overlap
        both_null[i] = np.count_nonzero(o1 & o2)
        either_null[i] = np.count_nonzero(o1 | o2)
    return both_null, either_null

def permute_vs_islands(twoanchor_arrays, island_arrays, window_bp, chrom_sizes, n_perm, rng):
    c1, s1, e1, c2, s2, e2 = twoanchor_arrays
    both_null = np.empty(n_perm, dtype=np.int64)
    either_null = np.empty(n_perm, dtype=np.int64)
    for i in range(n_perm):
        index = moving_index_from_islands(island_arrays, window_bp, chrom_sizes, rng) #circular shuffle
        o1 = batched_overlap(c1, s1, e1, index) # check overlap
        o2 = batched_overlap(c2, s2, e2, index) # check overlap
        both_null[i] = np.count_nonzero(o1 & o2)
        either_null[i] = np.count_nonzero(o1 | o2)
    return both_null, either_null

def permute_peaks_vs_islands(peaks_arrays, island_arrays, window_bp, chrom_sizes, n_perm, rng):
    chrom_arr, start_arr, end_arr = island_arrays
    ext_start, ext_end = extend_arrays(chrom_arr, start_arr, end_arr, window_bp, chrom_sizes)
    fixed_index = build_index_from_arrays(chrom_arr, ext_start, ext_end)

    null = np.empty(n_perm, dtype=np.int64)
    for i in range(n_perm):
        pc2, ps2, pe2 = circular_shuffle_arrays(*peaks_arrays, chrom_sizes, rng)
        null[i] = np.count_nonzero(batched_overlap(pc2, ps2, pe2, fixed_index))
    return null

print(f"Reading existing workbook: {IN_PATH}")
existing_sheets = pd.read_excel(IN_PATH, sheet_name=None)

chrom_sizes = load_chrom_sizes()
islands = load_islands()
peaks = {key: load_peaks(prefix) for key, prefix in PEAK_SAMPLES.items()}
contacts = {name: load_contacts(name) for name in CONTACT_SOURCES}
loops = {name: load_loops(name) for name in LOOP_SOURCES}
island_arrays = to_arrays(islands)
peaks_arrays_by_key = {key: to_arrays(df) for key, df in peaks.items()}
contact_arrays_by_source = {name: to_two_anchor_arrays(df) for name, df in contacts.items()}
loop_arrays_by_source = {name: to_two_anchor_arrays(df) for name, df in loops.items()}

out_sheets = {}

CONTACTS_VS_PEAKS_SHEETS = {
    "ArimaFLF_Contacts_vs_Peaks": ("ArimaFLF", "FL"),
    "ArimaPF_Contacts_vs_Peaks": ("ArimaPF", "PF"),
    "DoveFLF_Contacts_vs_Peaks": ("DovetailFLF", "FL"),
    "DovePF_Contacts_vs_Peaks": ("DovetailPF", "PF"),
}
LOOPS_VS_PEAKS_SHEETS = {
    "DoveFLF_Loops_vs_Peaks": ("DovetailFLF", "FL"),
    "DovePF_Loops_vs_Peaks": ("DovetailPF", "PF"),
}
CONTACTS_VS_PI_SHEETS = {
    "ArimaFLF_Contacts_vs_PI": "ArimaFLF",
    "ArimaPF_Contacts_vs_PI": "ArimaPF",
    "DoveFLF_Contacts_vs_PI": "DovetailFLF",
    "DovePF_Contacts_vs_PI": "DovetailPF",
}
LOOPS_VS_PI_SHEETS = {
    "DoveFLF_Loops_vs_PI": "DovetailFLF",
    "DovePF_Loops_vs_PI": "DovetailPF",
}
PEAKS_VS_PI_SHEETS = {"FLF_Peaks_vs_PI": "FL", "PF_Peaks_vs_PI": "PF"}

for sheet_name, (source_key, stage) in {
    **CONTACTS_VS_PEAKS_SHEETS,
    **LOOPS_VS_PEAKS_SHEETS,
}.items():
    df = existing_sheets[sheet_name].copy()
    is_contacts = sheet_name in CONTACTS_VS_PEAKS_SHEETS
    twoanchor_arrays = (
        contact_arrays_by_source[source_key] if is_contacts else loop_arrays_by_source[source_key]
    )
    print(f"\n{sheet_name}: permuting {N_PERMUTATIONS}x per mark...")
    for i, row in df.iterrows():
        key = (
            stage,
            row["replicate"] if pd.notna(row["replicate"]) else None,
            row["mark"],
        )
        both_null, either_null = permute_vs_peaks(
            twoanchor_arrays, peaks_arrays_by_key[key], chrom_sizes, N_PERMUTATIONS, rng
        )
        exp_b, fold_b, z_b, p_b = permutation_stats(row["both"], both_null)
        exp_e, fold_e, z_e, p_e = permutation_stats(row["either"], either_null)
        df.loc[
            i, ["expected_both", "fold_enrichment_both", "zscore_both", "pvalue_both"]
        ] = [exp_b, fold_b, z_b, p_b]
        df.loc[
            i,
            [
                "expected_either",
                "fold_enrichment_either",
                "zscore_either",
                "pvalue_either",
            ],
        ] = [exp_e, fold_e, z_e, p_e]
        print(
            f"  {row['mark']} ({row.get('replicate', '')}): fold_both={fold_b}, p_both={p_b}"
        )
    out_sheets[sheet_name] = df

for sheet_name, source_key in {**CONTACTS_VS_PI_SHEETS, **LOOPS_VS_PI_SHEETS}.items():
    df = existing_sheets[sheet_name].copy()
    is_contacts = sheet_name in CONTACTS_VS_PI_SHEETS
    twoanchor_arrays = (
        contact_arrays_by_source[source_key] if is_contacts else loop_arrays_by_source[source_key]
    )
    print(f"\n{sheet_name}: permuting {N_PERMUTATIONS}x per window...")
    for i, row in df.iterrows():
        window_bp = WINDOWS[row["window"]]
        both_null, either_null = permute_vs_islands(
            twoanchor_arrays, island_arrays, window_bp, chrom_sizes, N_PERMUTATIONS, rng
        )
        exp_b, fold_b, z_b, p_b = permutation_stats(row["both"], both_null)
        exp_e, fold_e, z_e, p_e = permutation_stats(row["either"], either_null)
        df.loc[
            i, ["expected_both", "fold_enrichment_both", "zscore_both", "pvalue_both"]
        ] = [exp_b, fold_b, z_b, p_b]
        df.loc[
            i,
            [
                "expected_either",
                "fold_enrichment_either",
                "zscore_either",
                "pvalue_either",
            ],
        ] = [exp_e, fold_e, z_e, p_e]
        print(f"  {row['window']}: fold_both={fold_b}, p_both={p_b}")
    out_sheets[sheet_name] = df

for sheet_name, stage in PEAKS_VS_PI_SHEETS.items():
    df = existing_sheets[sheet_name].copy()
    print(f"\n{sheet_name}: permuting {N_PERMUTATIONS}x per mark/window...")
    for i, row in df.iterrows():
        key = (
            stage,
            row["replicate"] if pd.notna(row["replicate"]) else None,
            row["mark"],
        ) # gives a key to look up a specific peak
        window_bp = WINDOWS[row["window"]]
        null = permute_peaks_vs_islands(
            peaks_arrays_by_key[key], island_arrays, window_bp, chrom_sizes, N_PERMUTATIONS, rng
        )
        exp, fold, z, p = permutation_stats(row["overlapping"], null)
        df.loc[i, ["expected_overlapping", "fold_enrichment", "zscore", "pvalue"]] = [
            exp,
            fold,
            z,
            p,
        ]
        print(
            f"  {row['mark']} ({row.get('replicate', '')}) {row['window']}: fold={fold}, p={p}"
        )
    out_sheets[sheet_name] = df

for sheet_name in ["DoveFLF_Contacts_vs_Loops", "DovePF_Contacts_vs_Loops"]:
    out_sheets[sheet_name] = existing_sheets[sheet_name].copy()

# output
all_rows = []
for sheet_name, df in out_sheets.items():
    tagged = df.copy()
    tagged["comparison"] = sheet_name
    all_rows.append(tagged)
summary_df = pd.concat(all_rows, ignore_index=True, sort=False)
cols = ["comparison"] + [c for c in summary_df.columns if c != "comparison"]
summary_df = summary_df[cols]

with pd.ExcelWriter(OUT_PATH, engine="openpyxl") as writer:
    summary_df.to_excel(writer, sheet_name="Summary", index=False)
    for sheet_name, df in out_sheets.items():
        df.to_excel(writer, sheet_name=sheet_name[:31], index=False)

print(f"\nSaved {OUT_PATH} with {len(out_sheets) + 1} sheets.")
