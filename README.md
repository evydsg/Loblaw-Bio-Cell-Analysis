# Loblaw Bio — Immune Cell Population Dashboard

A clinical trial analysis tool built for **Bob Loblaw** at Loblaw Bio. This project loads immune cell count data into a relational SQLite database, runs statistical analyses across five immune cell populations, and presents the results in an interactive dashboard.

---

## Quick Start (GitHub Codespaces)

Open this repository in GitHub Codespaces, then run the three commands below in order.

### 1. Install dependencies
```bash
make setup
```
Creates a Python virtual environment (`.venv`) and installs all required packages.

### 2. Run the pipeline
```bash
make pipeline
```
Runs all four analysis scripts in order:
- `load_data.py` — creates `loblaw.db` and loads `cell-count.csv`
- `frequencies.py` — computes relative cell population frequencies
- `responder_analysis.py` — compares responders vs non-responders (Part 3)
- `subset_analysis.py` — runs the baseline subset queries (Part 4)

### 3. Start the dashboard
```bash
make dashboard
```
Starts the Streamlit server on port 8501. In Codespaces, click the **Open in Browser** pop-up or go to the **Ports** tab and open port 8501.

---

## Live Dashboard

🔗 [View the hosted dashboard](<YOUR_STREAMLIT_CLOUD_URL_HERE>)

---

## Project Structure

```
.
├── load_data.py          # Part 1: creates loblaw.db from cell-count.csv
├── frequencies.py        # Part 2: relative frequency table
├── responder_analysis.py # Part 3: responder vs non-responder statistics
├── subset_analysis.py    # Part 4: baseline subset queries
├── app.py                # Streamlit dashboard
├── cell-count.csv        # Raw data (10,500 samples, 3,500 subjects)
├── requirements.txt      # Python dependencies
├── Makefile              # setup / pipeline / dashboard targets
└── outputs/              # Generated charts and CSV files (created by pipeline)
```

---

## Analysis Summary

### Part 1 — Database Design
The raw CSV was normalised into a SQLite database (`loblaw.db`) with four linked tables: `projects`, `subjects`, `samples`, and `cell_counts`. Storing each fact once avoids duplication and makes queries fast and reliable.

### Part 2 — Cell Population Frequencies
For each of the 10,500 samples, the five cell counts (b_cell, cd4_t_cell, cd8_t_cell, nk_cell, monocyte) are expressed as a percentage of that sample's total. This puts all samples on the same scale regardless of how many cells were collected.

### Part 3 — Responders vs Non-Responders
**Group:** melanoma patients treated with miraclib, PBMC samples, baseline (day 0) only. Each patient contributes one sample, keeping comparisons independent.

**Method:** two-sided Mann-Whitney U test per population, corrected for five simultaneous comparisons using the Benjamini-Hochberg method (α = 0.05). We used day 0 only because Bob's goal is to *predict* response before treatment takes effect. We also checked days 7 and 14; no population was significant after correction at either time point.

**Result:** no population showed a significant difference between responders (n = 331) and non-responders (n = 325). All adjusted p-values = 0.885. Median differences were at most 0.7 percentage points and effect sizes were negligible (all below 0.06). These five cell populations do not predict response to miraclib at baseline in this dataset.

### Part 4 — Subset Analysis
**Baseline group** (melanoma · miraclib · PBMC · day 0 · n = 656):

| Breakdown | Result |
|---|---|
| Samples per project | prj1: 384, prj3: 272 (prj2: none) |
| Subjects by response | Responders: 331, Non-responders: 325 |
| Subjects by sex | Male: 344, Female: 312 |

**Average B cell count** — melanoma males, responders, day 0, all sample and treatment types: **10,206.15**

---

## Dependencies

| Package | Version | Purpose |
|---|---|---|
| pandas | ≥ 2.0 | Data manipulation |
| scipy | ≥ 1.11 | Statistical tests |
| matplotlib | ≥ 3.9 | Charts and boxplots |
| streamlit | ≥ 1.30 | Interactive dashboard |

---

## Notes

- `.venv/` and `loblaw.db` are not committed to the repository. The database is rebuilt by `make pipeline` on every run.
- The dashboard reads directly from `loblaw.db`, not from saved CSV files.
- All analysis scripts can also be run individually: `python load_data.py`, `python frequencies.py`, etc.