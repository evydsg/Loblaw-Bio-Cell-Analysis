"""Dashboard for Bob Loblaw's clinical trial analysis.

Run with:  streamlit run app.py  (or: make dashboard)
"""
import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import streamlit as st

import load_data
from frequencies import get_frequency_table
from responder_analysis import ALPHA, get_group_data, make_boxplots, run_tests, timepoint_label
from subset_analysis import (
    avg_b_cell_male_responders_summary,
    get_baseline_samples,
    samples_per_project,
    subjects_by_response,
    subjects_by_sex,
)

ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "loblaw.db"

st.set_page_config(
    page_title="Loblaw Bio Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed",
)

@st.cache_resource
def ensure_database():
    if not DB_PATH.exists():
        load_data.main()

ensure_database()


def get_overview(db_path=DB_PATH) -> dict:
    """Row counts per table, project ids and sampling timepoints, from the database."""
    with sqlite3.connect(db_path) as conn:
        counts = {t: conn.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
                  for t in ["projects", "subjects", "samples", "cell_counts"]}
        projects = [r[0] for r in conn.execute("SELECT project_id FROM projects ORDER BY project_id")]
        timepoints = [r[0] for r in conn.execute(
            "SELECT DISTINCT time_from_treatment_start FROM samples ORDER BY 1")]
    return {"counts": counts, "projects": projects, "timepoints": timepoints}

# ── CSS ───────────────────────────────────────────────────────────────────────
st.markdown("""
<style>
    /* hide sidebar and its toggle */
    [data-testid="collapsedControl"] { display: none; }
    [data-testid="stSidebar"]        { display: none; }

    /* pull content up to remove Streamlit's default top gap */
    .block-container { padding-top: 0rem !important; margin-top: 0rem !important; }

    /* full-width navbar pinned to the top */
    .navbar {
        display: flex; align-items: center; gap: 0;
        background: #0d47a1;
        padding: 0 32px;
        height: 56px;
        margin-bottom: 28px;
        margin-left: -5rem; margin-right: -5rem;   /* bleed past Streamlit padding */
        width: calc(100% + 10rem);
    }
    .navbar-brand {
        font-size: 1.2em; font-weight: 700; color: white;
        margin-right: 32px; letter-spacing: 0.5px; white-space: nowrap;
    }
    .nav-link {
        color: rgba(255,255,255,0.75);
        padding: 6px 20px; border-radius: 20px;
        font-size: 0.95em; text-decoration: none;
        border: 2px solid transparent;
        transition: all 0.2s; white-space: nowrap;
    }
    .nav-link:hover { color: white; border-color: rgba(255,255,255,0.5); }
    .nav-link.active { color: #0d47a1; background: white; font-weight: 600; }

    /* section headers */
    .section-header {
        background: linear-gradient(90deg, #1565c0, #0d47a1);
        color: white; padding: 12px 20px; border-radius: 8px;
        margin: 24px 0 16px 0; font-size: 1.1em; font-weight: 600;
    }

    /* info cards */
    .info-card {
        background: #f0f4ff; border-left: 4px solid #1565c0;
        padding: 18px 22px; border-radius: 8px; margin-bottom: 12px;
    }
    .info-card h4 { margin: 0 0 6px 0; color: #0d47a1; }
    .info-card p  { margin: 0; color: #333; line-height: 1.6; }

    /* result boxes */
    .finding-box {
        background: #e8f5e9; border-left: 4px solid #2e7d32;
        padding: 16px; border-radius: 8px; margin-top: 16px;
    }
    .warning-box {
        background: #fff8e1; border-left: 4px solid #f9a825;
        padding: 16px; border-radius: 8px; margin-top: 16px;
    }

    /* style Streamlit buttons to look like nav pills */
    div[data-testid="stHorizontalBlock"] div[data-testid="stButton"] button {
        background: transparent;
        border: 2px solid rgba(255,255,255,0.45);
        color: white !important;
        border-radius: 20px;
        padding: 4px 16px;
        font-size: 0.9em;
        transition: all 0.2s;
    }
    div[data-testid="stHorizontalBlock"] div[data-testid="stButton"] button:hover {
        background: white;
        color: #0d47a1 !important;
        border-color: white;
        font-weight: 600;
    }
</style>
""", unsafe_allow_html=True)

# ── Navigation state ──────────────────────────────────────────────────────────
if "page" not in st.session_state:
    st.session_state.page = "Home"
if "font_size" not in st.session_state:
    st.session_state.font_size = 16
if "night_mode" not in st.session_state:
    st.session_state.night_mode = False

# ── Navbar ────────────────────────────────────────────────────────────────────
st.markdown("""
<div class="navbar">
    <span class="navbar-brand">🧬 Loblaw Bio</span>
</div>
""", unsafe_allow_html=True)

NAV_ITEMS = [
    ("Home",  "🏠", "Home",   "Overview of the trial and data"),
    ("Part2", "📊", "Part 2", "Cell population frequencies"),
    ("Part3", "🔬", "Part 3", "Responders vs non-responders"),
    ("Part4", "🔍", "Part 4", "Baseline subset counts"),
]

# Nav row: buttons + spacer + font slider + night mode
nav_col1, nav_col2, nav_col3, nav_col4, nav_spacer, font_col, night_col = \
    st.columns([1.1, 1.1, 1.1, 1.1, 2.5, 2, 1])

for col, (key, icon, label, tip) in zip(
    [nav_col1, nav_col2, nav_col3, nav_col4], NAV_ITEMS
):
    with col:
        is_active = st.session_state.page == key
        if st.button(
            f"{icon} {label}",
            width="stretch",
            help=tip,
            type="primary" if is_active else "secondary",
        ):
            st.session_state.page = key

with font_col:
    st.session_state.font_size = st.slider(
        "🔡 Text size", min_value=12, max_value=22,
        value=st.session_state.font_size, step=1,
        help="Adjust body text size",
    )

with night_col:
    st.markdown("<div style='height:8px'></div>", unsafe_allow_html=True)
    st.session_state.night_mode = st.toggle(
        "🌙 Night", value=st.session_state.night_mode,
        help="Switch between light and dark mode",
    )

# ── Dynamic theme ─────────────────────────────────────────────────────────────
if st.session_state.night_mode:
    bg        = "#0e1117"
    card_bg   = "#1e2130"
    card_border = "#3a4060"
    text      = "#e8eaf6"
    subtext   = "#9fa8da"
    header_bg = "linear-gradient(90deg, #1a237e, #283593)"
    find_bg   = "#1b2e1b"
    find_border = "#2e7d32"
    warn_bg   = "#2e2600"
    warn_border = "#f9a825"
    hr_color  = "#2a2f45"
else:
    bg        = "#ffffff"
    card_bg   = "#f0f4ff"
    card_border = "#1565c0"
    text      = "#111111"
    subtext   = "#555555"
    header_bg = "linear-gradient(90deg, #1565c0, #0d47a1)"
    find_bg   = "#e8f5e9"
    find_border = "#2e7d32"
    warn_bg   = "#fff8e1"
    warn_border = "#f9a825"
    hr_color  = "#e0e0e0"

st.markdown(f"""
<style>
    html, body, [class*="css"], p, li, span, div {{
        font-size: {st.session_state.font_size}px !important;
    }}
    .stApp {{ background-color: {bg} !important; color: {text} !important; }}
    .stApp p, .stApp li, .stApp span {{ color: {text} !important; }}
    .stApp h1, .stApp h2, .stApp h3, .stApp h4 {{ color: {text} !important; }}
    .section-header {{ background: {header_bg}; }}
    .info-card {{ background: {card_bg}; border-left-color: {card_border}; }}
    .info-card h4 {{ color: {card_border}; }}
    .info-card p  {{ color: {text}; }}
    .finding-box {{ background: {find_bg}; border-left-color: {find_border}; }}
    .warning-box {{ background: {warn_bg};  border-left-color: {warn_border}; }}
    hr {{ border-color: {hr_color} !important; }}
</style>
""", unsafe_allow_html=True)

st.markdown(f"<hr style='margin:4px 0 24px 0;border-color:{hr_color};'>", unsafe_allow_html=True)
page = st.session_state.page

# ── Home ──────────────────────────────────────────────────────────────────────
if page == "Home":
    st.markdown("# Loblaw Bio — Immune Cell Population Dashboard")
    st.markdown("### Clinical Trial Analysis: Miraclib in Melanoma & Carcinoma")
    st.write(
        "This dashboard presents the results of a clinical trial run by **Bob Loblaw** at Loblaw Bio. "
        "The trial investigates how the drug candidate **miraclib** affects immune cell populations "
        "in patients' blood, and whether those cell populations can help predict who will respond to treatment."
    )

    st.markdown('<div class="section-header">🗂️ About the Data</div>', unsafe_allow_html=True)

    overview = get_overview(DB_PATH)
    n = overview["counts"]
    days = ", ".join(str(t) for t in overview["timepoints"][:-1]) + f" and {overview['timepoints'][-1]}"

    c1, c2 = st.columns(2)
    with c1:
        st.markdown("""
        <div class="info-card">
            <h4>📦 Dataset</h4>
            <p>
                The data comes from <strong>cell-count.csv</strong>, a file with one row per biological sample.
                Each row records how many cells of each immune type were counted in that sample,
                along with details about the patient it came from.
            </p>
        </div>
        """, unsafe_allow_html=True)
        st.markdown("""
        <div class="info-card">
            <h4>🧪 Five Immune Cell Populations</h4>
            <p>
                <strong>B cells</strong> — produce antibodies to fight infection.<br>
                <strong>CD4 T cells</strong> — "helper" cells that coordinate the immune response.<br>
                <strong>CD8 T cells</strong> — "killer" cells that destroy infected or cancerous cells.<br>
                <strong>NK cells</strong> — natural killer cells that attack tumours without prior exposure.<br>
                <strong>Monocytes</strong> — patrol the blood and trigger inflammation.
            </p>
        </div>
        """, unsafe_allow_html=True)
    with c2:
        st.markdown(f"""
        <div class="info-card">
            <h4>📊 Scale of the Trial</h4>
            <p>
                <strong>{n['projects']:,}</strong> projects ({", ".join(overview["projects"])})<br>
                <strong>{n['subjects']:,}</strong> subjects (patients)<br>
                <strong>{n['samples']:,}</strong> samples — collected at day {days}<br>
                <strong>{n['cell_counts']:,}</strong> cell count measurements in total
            </p>
        </div>
        """, unsafe_allow_html=True)
        st.markdown("""
        <div class="info-card">
            <h4>💊 Treatments & Conditions</h4>
            <p>
                <strong>Conditions:</strong> Melanoma, Carcinoma, Healthy controls<br>
                <strong>Treatments:</strong> Miraclib, Phauximab, None (healthy controls)<br>
                <strong>Sample types:</strong> PBMC (peripheral blood), Whole Blood<br>
                <strong>Response:</strong> Yes / No (not recorded for healthy controls)
            </p>
        </div>
        """, unsafe_allow_html=True)

    st.markdown('<div class="section-header">🗄️ Part 1: Database Design</div>', unsafe_allow_html=True)
    st.write(
        "The raw CSV stored every patient detail on every row, repeating the same information thousands of times. "
        "We reorganised it into a **relational SQLite database** with four linked tables:"
    )
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Projects table", f"{n['projects']:,} rows", "One row per study")
    c2.metric("Subjects table", f"{n['subjects']:,} rows", "One row per patient")
    c3.metric("Samples table", f"{n['samples']:,} rows", "One row per blood draw")
    c4.metric("Cell counts table", f"{n['cell_counts']:,} rows", "One row per cell type per sample")
    st.write(
        "Storing each fact once avoids contradictions, makes queries fast, "
        "and means adding a sixth cell type later needs no redesign."
    )

    st.markdown('<div class="section-header">🧭 How to Use This Dashboard</div>', unsafe_allow_html=True)
    st.markdown("""
    | Section | What you will find |
    |---|---|
    | **Part 2** | The relative frequency of each cell type in every sample, with filters and charts |
    | **Part 3** | Statistical comparison of responders vs non-responders in melanoma patients on miraclib |
    | **Part 4** | Counts for the baseline melanoma PBMC group, and the average B cell count for male responders |
    """)

# ── Part 2 ────────────────────────────────────────────────────────────────────
elif page == "Part2":
    st.markdown("# Part 2: Relative Frequency of Immune Cell Populations")
    st.markdown("""
    **What we did:** Bob's first question was *"what is the frequency of each cell type in each sample?"*
    Raw cell counts aren't directly comparable between samples, because some blood draws simply contain more
    cells than others. To put every sample on the same scale, we calculated each population's count as a
    **percentage of that sample's total cell count**. A sample with 10,000 B cells out of 100,000 total
    has a B cell frequency of 10%, regardless of how large or small the sample was.
    """)

    df = get_frequency_table(DB_PATH)

    st.markdown('<div class="section-header">📊 Explore the Data</div>', unsafe_allow_html=True)
    c1, c2 = st.columns(2)
    with c1:
        chosen = st.selectbox("Filter by sample", ["All"] + sorted(df["sample"].unique().tolist()))
    with c2:
        chosen_pop = st.selectbox("Filter by population", ["All"] + sorted(df["population"].unique().tolist()))

    view = df.copy()
    if chosen != "All":
        view = view[view["sample"] == chosen]
    if chosen_pop != "All":
        view = view[view["population"] == chosen_pop]

    if chosen_pop != "All":
        st.markdown(f"#### Distribution of {chosen_pop} frequency across all samples")
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.hist(view["percentage"], bins=40, color="#1565c0", edgecolor="white", alpha=0.85, label=chosen_pop)
        ax.set_xlabel("Relative frequency (%)")
        ax.set_ylabel("Number of samples")
        ax.set_title(f"{chosen_pop} — frequency distribution")
        ax.legend(title="Population", loc="upper right")
        fig.tight_layout()
        st.pyplot(fig)
        plt.close(fig)
    elif chosen != "All":
        st.markdown(f"#### Cell population breakdown for {chosen}")
        pop_colors = {"b_cell": "#1565c0", "cd4_t_cell": "#0d47a1",
                      "cd8_t_cell": "#1976d2", "monocyte": "#42a5f5", "nk_cell": "#90caf9"}
        fig, ax = plt.subplots(figsize=(7, 3))
        bars = ax.barh(view["population"], view["percentage"],
                       color=[pop_colors.get(p, "#1565c0") for p in view["population"]], alpha=0.85)
        ax.set_xlabel("Relative frequency (%)")
        ax.set_title(f"Sample {chosen}")
        ax.set_xlim(0, max(view["percentage"]) * 1.2)
        ax.legend(bars, view["population"].tolist(), title="Population", loc="lower right")
        fig.tight_layout()
        st.pyplot(fig)
        plt.close(fig)
    else:
        st.info("Select a sample or population above to see a chart.")

    st.markdown("#### Data Table")
    if len(view) <= 10000:
        display = view.style.format({"percentage": "{:.2f}", "total_count": "{:,}", "count": "{:,}"})
    else:
        display = view
    st.dataframe(display, width="stretch", height=400)
    st.caption(f"Showing {len(view):,} of {len(df):,} rows.")

# ── Part 3 ────────────────────────────────────────────────────────────────────
elif page == "Part3":
    st.markdown("# Part 3: Responders vs Non-Responders")
    st.markdown("""
    **What we did:** Bob wants to know whether immune cell frequencies can predict who will respond
    to miraclib. We filtered the data to **melanoma patients on miraclib, PBMC samples**, and compared
    the five cell population frequencies between responders and non-responders using the
    **Mann-Whitney U test**, which does not assume normally distributed data. Because we ran five tests
    (one per population), we corrected for multiple comparisons using the **Benjamini-Hochberg method**
    (α = 0.05) to avoid false positives.

    **Recommended view: day 0 (baseline).** It shows the cells *before* treatment has had any effect,
    which is what matters for predicting response, and each patient contributes one sample, so the
    samples are independent as the test assumes. Pooling all timepoints counts each patient up to
    three times, which makes p-values look stronger than they should.
    """)

    tp_options = {"Day 0 (baseline, recommended)": 0, "Day 7": 7, "Day 14": 14, "All timepoints pooled": None}
    tp_choice = st.radio("Timepoint", list(tp_options), horizontal=True)
    timepoint = tp_options[tp_choice]

    data = get_group_data(DB_PATH, timepoint=timepoint)
    res  = run_tests(data)
    patients = data.drop_duplicates("subject_id")

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Samples analysed", f"{data['sample'].nunique():,}")
    c2.metric("Patients", f"{len(patients):,}")
    c3.metric("Responders", f"{(patients['response'] == 'yes').sum():,}")
    c4.metric("Non-responders", f"{(patients['response'] == 'no').sum():,}")

    st.markdown('<div class="section-header">📦 Boxplots: Cell Population Frequencies by Response</div>', unsafe_allow_html=True)
    st.write("Each box shows the spread of relative frequencies for one cell population. Responders are on the left, non-responders on the right.")
    fig = make_boxplots(data, res, timepoint)
    st.pyplot(fig)
    plt.close(fig)

    st.markdown('<div class="section-header">📋 Statistical Results</div>', unsafe_allow_html=True)
    st.write("Two-sided Mann-Whitney U test per population · Benjamini-Hochberg correction · α = 0.05 · "
             "effect size = rank-biserial correlation (−1 to +1)")
    display_res = res[[
        "population", "n_responders", "n_non_responders",
        "median_responders", "median_non_responders",
        "median_difference", "effect_size", "p_value", "p_adjusted", "significant"
    ]].copy()
    st.dataframe(
        display_res.style.format({
            "median_responders": "{:.2f}", "median_non_responders": "{:.2f}",
            "median_difference": "{:.2f}", "effect_size": "{:.3f}",
            "p_value": "{:.4f}", "p_adjusted": "{:.4f}",
        }),
        width="stretch",
    )

    sig = res.loc[res["significant"], "population"].tolist()
    nominal = res.loc[(res["p_value"] < ALPHA) & ~res["significant"], "population"].tolist()
    if sig:
        st.markdown(f'<div class="finding-box">✅ <strong>Significant populations after correction:</strong> {", ".join(sig)}</div>', unsafe_allow_html=True)
    else:
        nominal_note = (
            f" {', '.join(nominal)} reached p &lt; {ALPHA} before correction, but not after, "
            "so it should be treated as a lead to follow up rather than a finding."
            if nominal else ""
        )
        st.markdown(f"""
        <div class="warning-box">
        ⚠️ <strong>No significant differences after correction ({timepoint_label(timepoint)}).</strong><br><br>
        Smallest adjusted p = {res['p_adjusted'].min():.3f} (threshold {ALPHA}). Median differences
        are at most {res['median_difference'].abs().max():.2f} percentage points and the largest
        effect size is {res['effect_size'].abs().max():.3f}.{nominal_note}
        </div>
        """, unsafe_allow_html=True)

# ── Part 4 ────────────────────────────────────────────────────────────────────
elif page == "Part4":
    st.markdown("# Part 4: Subset Analysis")
    st.markdown("""
    **What we did:** Bob wanted to explore a specific slice of the data to understand early treatment effects.
    We queried the database for **melanoma patients treated with miraclib, PBMC samples only, at baseline
    (day 0)**, then broke this group down by project, treatment response, and sex. Finally, we answered
    a separate question about the average B cell count for a broader group: melanoma males who are
    responders at day 0, across all sample and treatment types.
    """)

    baseline = get_baseline_samples(DB_PATH)
    proj = samples_per_project(baseline)
    resp = subjects_by_response(baseline)
    sex  = subjects_by_sex(baseline)

    st.metric("Total samples in group (melanoma · miraclib · PBMC · day 0)", len(baseline))

    st.markdown('<div class="section-header">📊 Breakdowns</div>', unsafe_allow_html=True)
    c1, c2, c3 = st.columns(3)

    with c1:
        st.markdown("#### Samples per project")
        fig, ax = plt.subplots(figsize=(4, 3))
        bars = ax.bar(proj["project_id"], proj["n_samples"], color=["#1565c0","#0d47a1"], alpha=0.85)
        ax.set_ylabel("Samples")
        ax.set_title("Samples per project")
        ax.set_ylim(0, max(proj["n_samples"]) * 1.15)
        for bar, v in zip(bars, proj["n_samples"]):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 3, str(v), ha="center", fontweight="bold")
        ax.legend(bars, proj["project_id"].tolist(), title="Project", loc="upper right")
        fig.tight_layout()
        st.pyplot(fig)
        plt.close(fig)
        st.dataframe(proj, width="stretch", hide_index=True)
        missing = sorted(set(get_overview(DB_PATH)["projects"]) - set(proj["project_id"]))
        if missing:
            st.caption(f"Note: {', '.join(missing)} has no melanoma patients on miraclib "
                       "with PBMC samples at day 0.")

    with c2:
        st.markdown("#### Subjects by response")
        fig, ax = plt.subplots(figsize=(4, 3))
        colors = ["#c62828","#2e7d32"]
        bars = ax.bar(resp["response"], resp["n_subjects"], color=colors, alpha=0.85)
        ax.set_ylabel("Subjects")
        ax.set_title("Subjects by response")
        ax.set_ylim(0, max(resp["n_subjects"]) * 1.15)
        for bar, v in zip(bars, resp["n_subjects"]):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2, str(v), ha="center", fontweight="bold")
        ax.legend(bars, resp["response"].tolist(), title="Response", loc="upper right")
        fig.tight_layout()
        st.pyplot(fig)
        plt.close(fig)
        st.dataframe(resp, width="stretch", hide_index=True)

    with c3:
        st.markdown("#### Subjects by sex")
        fig, ax = plt.subplots(figsize=(4, 3))
        colors2 = ["#6a1b9a","#1565c0"]
        bars = ax.bar(sex["sex"], sex["n_subjects"], color=colors2, alpha=0.85)
        ax.set_ylabel("Subjects")
        ax.set_title("Subjects by sex")
        ax.set_ylim(0, max(sex["n_subjects"]) * 1.15)
        for bar, v in zip(bars, sex["n_subjects"]):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2, str(v), ha="center", fontweight="bold")
        ax.legend(bars, sex["sex"].tolist(), title="Sex", loc="upper right")
        fig.tight_layout()
        st.pyplot(fig)
        plt.close(fig)
        st.dataframe(sex, width="stretch", hide_index=True)

    b_cell = avg_b_cell_male_responders_summary(DB_PATH)
    st.markdown('<div class="section-header">🧫 Average B Cell Count — Melanoma Male Responders at Day 0</div>', unsafe_allow_html=True)
    st.write(
        "This question uses a broader filter: **melanoma males who are responders at time = 0**, "
        "across **all** sample types (PBMC and whole blood) and **all** treatments (miraclib and phauximab). "
        f"This group has {b_cell['n_samples']:,} samples from {b_cell['n_subjects']:,} distinct patients."
    )
    st.metric("Average B cell count", f"{b_cell['avg_b_cell']:.2f}")