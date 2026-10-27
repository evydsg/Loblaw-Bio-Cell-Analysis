"""Part 4: data subset analysis.

Baseline (day 0) melanoma PBMC samples from miraclib-treated patients,
with breakdowns by project, response, and sex.

Run with:  python subset_analysis.py
"""
import sqlite3
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "loblaw.db"

BASELINE_SQL = """
SELECT
    s.sample_id,
    s.subject_id,
    u.project_id,
    u.response,
    u.sex
FROM samples s
JOIN subjects u ON u.subject_id = s.subject_id
WHERE u.condition = 'melanoma'
  AND s.sample_type = 'PBMC'
  AND s.time_from_treatment_start = 0
  AND u.treatment = 'miraclib'
ORDER BY s.sample_id
"""

AVG_B_CELL_SQL = """
SELECT AVG(cc.count) AS avg_b_cell, COUNT(*) AS n_samples,
       COUNT(DISTINCT s.subject_id) AS n_subjects
FROM samples s
JOIN subjects u ON u.subject_id = s.subject_id
JOIN cell_counts cc ON cc.sample_id = s.sample_id
WHERE u.condition = 'melanoma'
  AND u.sex = 'M'
  AND u.response = 'yes'
  AND s.time_from_treatment_start = 0
  AND cc.population = 'b_cell'
"""


def get_baseline_samples(db_path=DB_PATH) -> pd.DataFrame:
    """All melanoma PBMC baseline samples from miraclib-treated patients."""
    with sqlite3.connect(db_path) as conn:
        return pd.read_sql_query(BASELINE_SQL, conn)


def samples_per_project(df: pd.DataFrame) -> pd.DataFrame:
    """Number of samples from each project."""
    out = df.groupby("project_id")["sample_id"].nunique()
    return out.rename("n_samples").reset_index()


def subjects_by_response(df: pd.DataFrame) -> pd.DataFrame:
    """Number of distinct subjects who are responders / non-responders."""
    out = df.groupby("response")["subject_id"].nunique()
    return out.rename("n_subjects").reset_index()


def subjects_by_sex(df: pd.DataFrame) -> pd.DataFrame:
    """Number of distinct subjects who are male / female."""
    out = df.groupby("sex")["subject_id"].nunique()
    return out.rename("n_subjects").reset_index()


def avg_b_cell_male_responders(db_path=DB_PATH) -> float:
    """Average B cell count: melanoma males, responders, time 0, all sample
    types and all treatments. Returned rounded to two decimals."""
    with sqlite3.connect(db_path) as conn:
        avg = conn.execute(AVG_B_CELL_SQL).fetchone()[0]
    return round(avg, 2)


def main():
    df = get_baseline_samples()
    print(f"Baseline melanoma PBMC samples on miraclib: {len(df)}\n")
    print("Samples per project:")
    print(samples_per_project(df).to_string(index=False), "\n")
    print("Subjects by response:")
    print(subjects_by_response(df).to_string(index=False), "\n")
    print("Subjects by sex:")
    print(subjects_by_sex(df).to_string(index=False), "\n")
    print("Average B cell count, melanoma male responders at time 0 "
          f"(all sample/treatment types): {avg_b_cell_male_responders():.2f}")


if __name__ == "__main__":
    main()