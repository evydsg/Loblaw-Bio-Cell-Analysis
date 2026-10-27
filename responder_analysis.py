"""Part 3: do cell population frequencies differ between responders and
non-responders?

Group: melanoma patients on miraclib, PBMC samples. The main analysis uses
baseline (day 0) only, so each patient contributes exactly one sample (samples
are independent). All timepoints pooled is also reported, as the spec does not
restrict time; note that pooling counts each patient up to three times.

Run with:  python responder_analysis.py
Writes outputs/part3_boxplots_<label>.png and outputs/part3_stats_<label>.csv
for label in day0 and all_timepoints.
"""
import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # draw to a file, no window needed
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import false_discovery_control, mannwhitneyu

from frequencies import get_frequency_table

ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "loblaw.db"
OUTPUT_DIR = ROOT / "outputs"

ALPHA = 0.05
TIMEPOINT = 0  # baseline: before treatment has had any effect; None = pooled

GROUP_SQL = """
SELECT s.sample_id AS sample, s.subject_id, u.response
FROM samples s
JOIN subjects u ON u.subject_id = s.subject_id
WHERE u.condition = 'melanoma'
  AND u.treatment = 'miraclib'
  AND s.sample_type = 'PBMC'
  AND u.response IN ('yes', 'no')
"""


def get_group_data(db_path=DB_PATH, timepoint=TIMEPOINT) -> pd.DataFrame:
    """Frequency rows for the melanoma / miraclib / PBMC group, with response.
    timepoint=None pools all timepoints."""
    sql, params = GROUP_SQL, ()
    if timepoint is not None:
        sql, params = sql + "  AND s.time_from_treatment_start = ?\n", (timepoint,)
    with sqlite3.connect(db_path) as conn:
        meta = pd.read_sql_query(sql, conn, params=params)
    freq = get_frequency_table(db_path)
    return freq.merge(meta, on="sample", how="inner")


def run_tests(data: pd.DataFrame) -> pd.DataFrame:
    """Mann-Whitney U per population, with Benjamini-Hochberg correction."""
    rows = []
    for pop, d in data.groupby("population"):
        yes = d.loc[d.response == "yes", "percentage"]
        no = d.loc[d.response == "no", "percentage"]
        u, p = mannwhitneyu(yes, no, alternative="two-sided")
        rows.append({
            "population": pop,
            "n_responders": len(yes),
            "n_non_responders": len(no),
            "median_responders": yes.median(),
            "median_non_responders": no.median(),
            "median_difference": yes.median() - no.median(),
            "U": u,
            "p_value": p,
            # rank-biserial correlation: effect size from -1 to +1
            "effect_size": 2 * u / (len(yes) * len(no)) - 1,
        })
    res = pd.DataFrame(rows)
    res["p_adjusted"] = false_discovery_control(res["p_value"], method="bh")
    res["significant"] = res["p_adjusted"] < ALPHA
    return res


def timepoint_label(timepoint) -> str:
    return "all timepoints" if timepoint is None else f"day {timepoint}"


def make_boxplots(data: pd.DataFrame, res: pd.DataFrame, timepoint=TIMEPOINT):
    """Return a figure with one responder vs non-responder boxplot per population."""
    pops = sorted(data["population"].unique())
    fig, axes = plt.subplots(1, len(pops), figsize=(3 * len(pops), 4.5), sharey=False)
    for ax, pop in zip(axes, pops):
        d = data[data.population == pop]
        ax.boxplot(
            [d.loc[d.response == "yes", "percentage"], d.loc[d.response == "no", "percentage"]],
            tick_labels=["Responders", "Non-responders"], widths=0.6,
        )
        padj = res.loc[res.population == pop, "p_adjusted"].iloc[0]
        ax.set_title(f"{pop}\nadj. p = {padj:.3f}")
        ax.set_ylabel("Relative frequency (%)")
        ax.tick_params(axis="x", labelsize=8)
    fig.suptitle(f"Melanoma, miraclib, PBMC, {timepoint_label(timepoint)}: "
                 "responders vs non-responders")
    fig.tight_layout()
    return fig


def main():
    OUTPUT_DIR.mkdir(exist_ok=True)
    for timepoint, name in [(TIMEPOINT, "day0"), (None, "all_timepoints")]:
        data = get_group_data(timepoint=timepoint)
        res = run_tests(data)
        res.to_csv(OUTPUT_DIR / f"part3_stats_{name}.csv", index=False)
        fig = make_boxplots(data, res, timepoint)
        fig.savefig(OUTPUT_DIR / f"part3_boxplots_{name}.png", dpi=150)
        plt.close(fig)
        print(f"=== {timepoint_label(timepoint)}: {data['sample'].nunique()} samples, "
              f"{data['subject_id'].nunique()} patients | alpha = {ALPHA}")
        print(res.round(4).to_string(index=False))
        sig = res.loc[res.significant, "population"].tolist()
        print("Significant after correction:", sig if sig else "none", "\n")


if __name__ == "__main__":
    main()