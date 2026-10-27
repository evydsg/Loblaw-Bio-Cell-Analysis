"""Part 2: relative frequency of each immune cell population per sample.

Run with:  python frequencies.py
Reads loblaw.db (created by load_data.py) and writes outputs/frequency_table.csv
"""
import sqlite3
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "loblaw.db"
OUTPUT_DIR = ROOT / "outputs"

FREQUENCY_SQL = """
SELECT
    cc.sample_id                              AS sample,
    totals.total_count                        AS total_count,
    cc.population                             AS population,
    cc.count                                  AS count,
    100.0 * cc.count / totals.total_count     AS percentage
FROM cell_counts AS cc
JOIN (
    SELECT sample_id, SUM(count) AS total_count
    FROM cell_counts
    GROUP BY sample_id
) AS totals ON totals.sample_id = cc.sample_id
ORDER BY cc.sample_id, cc.population
"""


def get_frequency_table(db_path=DB_PATH) -> pd.DataFrame:
    """Return one row per (sample, population) with count and percentage."""
    with sqlite3.connect(db_path) as conn:
        return pd.read_sql_query(FREQUENCY_SQL, conn)


def main():
    df = get_frequency_table()
    OUTPUT_DIR.mkdir(exist_ok=True)
    out = OUTPUT_DIR / "frequency_table.csv"
    df.to_csv(out, index=False)
    print(f"{len(df)} rows written to {out.relative_to(ROOT)}")
    print(df.head(10).to_string(index=False, float_format=lambda x: f"{x:.2f}"))


if __name__ == "__main__":
    main()