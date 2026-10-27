"""Create the SQLite database and load cell-count.csv into it.

Run with:  python load_data.py
"""
import csv
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "loblaw.db"
CSV_PATH = ROOT / "cell-count.csv"
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

SCHEMA = """
PRAGMA foreign_keys = ON;

CREATE TABLE projects (
    project_id TEXT PRIMARY KEY
);

CREATE TABLE subjects (
    subject_id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL REFERENCES projects(project_id),
    condition  TEXT NOT NULL,
    age        INTEGER,
    sex        TEXT CHECK (sex IN ('M', 'F')),
    treatment  TEXT,
    response   TEXT CHECK (response IN ('yes', 'no'))   -- NULL for healthy controls
);

CREATE TABLE samples (
    sample_id   TEXT PRIMARY KEY,
    subject_id  TEXT NOT NULL REFERENCES subjects(subject_id),
    sample_type TEXT NOT NULL,
    time_from_treatment_start INTEGER NOT NULL
);

CREATE TABLE cell_counts (
    sample_id  TEXT NOT NULL REFERENCES samples(sample_id),
    population TEXT NOT NULL,
    count      INTEGER NOT NULL CHECK (count >= 0),
    PRIMARY KEY (sample_id, population)
);

CREATE INDEX idx_subjects_filters ON subjects(condition, treatment, response);
CREATE INDEX idx_samples_subject  ON samples(subject_id);
CREATE INDEX idx_samples_filters  ON samples(sample_type, time_from_treatment_start);
"""


def main():
    DB_PATH.unlink(missing_ok=True)  # start fresh so reruns are reproducible
    conn = sqlite3.connect(DB_PATH)
    conn.executescript(SCHEMA)

    with open(CSV_PATH, newline="") as f:
        for row in csv.DictReader(f):
            response = row["response"].strip() or None  # blank -> NULL
            conn.execute("INSERT OR IGNORE INTO projects VALUES (?)", (row["project"],))
            conn.execute(
                "INSERT OR IGNORE INTO subjects VALUES (?,?,?,?,?,?,?)",
                (row["subject"], row["project"], row["condition"],
                 int(row["age"]), row["sex"], row["treatment"], response),
            )
            conn.execute(
                "INSERT INTO samples VALUES (?,?,?,?)",
                (row["sample"], row["subject"], row["sample_type"],
                 int(row["time_from_treatment_start"])),
            )
            conn.executemany(
                "INSERT INTO cell_counts VALUES (?,?,?)",
                [(row["sample"], p, int(row[p])) for p in POPULATIONS],
            )

    conn.commit()
    for t in ["projects", "subjects", "samples", "cell_counts"]:
        n = conn.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
        print(f"{t:12s} {n:>6d} rows")
    conn.close()
    print(f"Database written to {DB_PATH.name}")


if __name__ == "__main__":
    main()