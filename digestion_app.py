import streamlit as st
from Bio import Restriction, SeqIO
from Bio.Seq import Seq
from itertools import combinations
import numpy as np
import plotly.graph_objects as go
import re
import random
import io
import pandas as pd

# ── PAGE SETTINGS ──────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Restriction Digest Tool",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 Restriction Digest Tool")
st.markdown("Finds the best enzyme combinations for diagnostic restriction digests.")

# ── ENZYMES ────────────────────────────────────────────────────────────────────
DEFAULT_ENZYMES = [
    "EcoRI", "HindIII", "BamHI", "XbaI", "SalI", "PstI",
    "SphI", "KpnI", "SacI", "XhoI", "NcoI", "NheI",
    "NotI", "XmaI", "SmaI", "ClaI", "EcoRV", "NdeI",
    "ScaI", "PvuI", "PvuII", "BglII", "ApaI", "MluI",
    "SpeI", "NsiI", "AluI", "HaeIII", "TaqI", "RsaI",
    "DraI", "AscI", "PacI", "AgeI", "SacII", "MfeI",
    "AflII", "BsrGI", "PmeI", "SfiI", "AvrII", "SbfI",
    "BclI", "BssHII", "BstBI", "BstXI", "NarI", "BspHI",
]

# ── SEQUENCE LOADING ───────────────────────────────────────────────────────────
def read_sbd(raw):
    try:
        text = raw.decode("latin-1")
        matches = re.findall(r'[ATCGatcg]{100,}', text)
        if not matches:
            return None
        return max(matches, key=len).upper()
    except:
        return None

def load_sequence(uploaded_file):
    if uploaded_file is None:
        seq = "".join(random.choices(list("ATCG"), k=10000))
        return seq, True

    raw = uploaded_file.read()
    name = uploaded_file.name.lower()

    if name.endswith(".sbd"):
        seq = read_sbd(raw)
        if seq:
            return seq, False

    for fmt in ["fasta", "genbank"]:
        try:
            record = SeqIO.read(io.StringIO(raw.decode("utf-8")), fmt)
            return str(record.seq).upper(), False
        except:
            pass
        try:
            record = SeqIO.read(io.StringIO(raw.decode("utf-16")), fmt)
            return str(record.seq).upper(), False
        except:
            pass

    for enc in ["utf-8", "utf-16", "latin-1"]:
        try:
            text = raw.decode(enc)
            matches = re.findall(r'[ATCGatcg]{200,}', text)
            if matches:
                return max(matches, key=len).upper(), False
        except:
            pass

    return None, False

# ── ANALYSIS FUNCTIONS ─────────────────────────────────────────────────────────
def get_fragments(plasmid_seq, enzymes, plasmid_size):
    rb = Restriction.RestrictionBatch(enzymes)
    cut_sites = []
    for enz, sites in rb.search(plasmid_seq, linear=False).items():
        cut_sites.extend(sites)
    if not cut_sites:
        return []
    cut_sites = sorted(set(cut_sites))
    fragments = []
    for i in range(len(cut_sites)):
        start = cut_sites[i]
        end = cut_sites[(i + 1) % len(cut_sites)]
        fragments.append(end - start if end > start else plasmid_size - start + end)
    return sorted(fragments)

def score_combination(fragments, min_frag, max_frag, min_frags, max_frags, min_diff):
    n = len(fragments)
    if n < min_frags or n > max_frags:
        return None
    if any(f < min_frag for f in fragments):
        return None
    if max_frag and any(f > max_frag for f in fragments):
        return None
    for i in range(len(fragments) - 1):
        ratio = (fragments[i+1] - fragments[i]) / fragments[i+1]
        if ratio < min_diff:
            return None
    return np.std(fragments) / np.mean(fragments)

def find_best_digests(plasmid_sequence, selected_enzymes,
                      min_frag, max_frag, min_frags, max_frags,
                      min_diff, combo_size, top_n, prefer_short=False):
    plasmid_seq = Seq(plasmid_sequence.upper())
    plasmid_size = len(plasmid_seq)
    resolved = [getattr(Restriction, e) for e in selected_enzymes if hasattr(Restriction, e)]
    rb_all = Restriction.RestrictionBatch(resolved)
    search_results = rb_all.search(plasmid_seq, linear=False)
    cutting_enzymes = [enz for enz in resolved if search_results[enz]]

    results = []
    for size in range(combo_size[0], combo_size[1] + 1):
        for combo in combinations(cutting_enzymes, size):
            frags = get_fragments(plasmid_seq, list(combo), plasmid_size)
            if not frags:
                continue
            score = score_combination(frags, min_frag, max_frag,
                                      min_frags, max_frags, min_diff)
            if score is not None:
                results.append({
                    "enzymes": " + ".join(e.__name__ for e in combo),
                    "fragments": frags,
                    "n": len(frags),
                    "score": score
                })

    if prefer_short:
        results.sort(key=lambda x: max(x["fragments"]))
    else:
        results.sort(key=lambda x: x["score"])

    seen = []
    unique_results = []
    for r in results:
        if r["fragments"] not in seen:
            seen.append(r["fragments"])
            unique_results.append(r)
        if len(unique_results) >= top_n:
            break
    return unique_results, cutting_enzymes

# ── GEL VISUALIZATION ──────────────────────────────────────────────────────────
def draw_gel(results, plasmid_size, title_suffix=""):
    fig = go.Figure()
    marker_sizes = [10000, 8000, 6000, 5000, 4000, 3500, 3000,
                    2500, 2000, 1500, 1000, 750, 500, 250]
    marker_sizes = [m for m in marker_sizes if m <= plasmid_size * 1.5]
    thick_bands = {6000, 3000, 1000}
    y_min = np.log10(200)
    y_max = np.log10(plasmid_size * 1.5)

    def bp_to_y(bp):
        return np.log10(bp)

    lane_width = 0.25
    band_height = 0.006
    band_height_thick = 0.008

    colors = ["#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
              "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3", "#ff6b6b", "#4ecdc4"]

    # ── Marker ──
    max_marker = max(marker_sizes)
    for size in marker_sizes:
        y = bp_to_y(size)
        is_thick = size in thick_bands
        height = band_height_thick if is_thick else band_height
        intensity = 0.2 + 0.75 * (size / max_marker)
        fig.add_shape(type="rect",
                      x0=-lane_width, x1=lane_width,
                      y0=y - height, y1=y + height,
                      fillcolor="white", opacity=intensity, line=dict(width=0))
        fig.add_shape(type="rect",
                      x0=-lane_width + 0.04, x1=lane_width - 0.04,
                      y0=y - 0.002, y1=y + 0.002,
                      fillcolor="white", opacity=intensity * 0.2, line=dict(width=0))
        fig.add_annotation(x=-lane_width - 0.05, y=y,
                           text=f"<b>{size} bp</b>" if is_thick else f"{size} bp",
                           showarrow=False,
                           font=dict(color="white" if is_thick else "#aaaaaa",
                                     size=10 if is_thick else 9),
                           xanchor="right", yanchor="middle")

    fig.add_annotation(x=0, y=y_max + 0.05, text="GeneRuler 1kb",
                       showarrow=False,
                       font=dict(color="white", size=11, family="Arial Black"),
                       xanchor="center", yanchor="bottom")

    # ── Sample lanes ──
    for i, result in enumerate(results):
        lane_x = i + 1
        color = colors[i % len(colors)]
        fig.add_annotation(x=lane_x, y=y_max + 0.05,
                           text=result["enzymes"].replace(" + ", "<br>+ "),
                           showarrow=False,
                           font=dict(color=color, size=9),
                           xanchor="center", yanchor="bottom",
                           bgcolor="#2a2a4a", bordercolor=color, borderwidth=1)
        max_frag = max(result["fragments"])
        for frag in result["fragments"]:
            y = bp_to_y(frag)
            intensity = 0.2 + 0.75 * (frag / max_frag)
            fig.add_shape(type="rect",
                          x0=lane_x - lane_width, x1=lane_x + lane_width,
                          y0=y - band_height, y1=y + band_height,
                          fillcolor="white", opacity=intensity, line=dict(width=0))
            fig.add_shape(type="rect",
                          x0=lane_x - lane_width + 0.04, x1=lane_x + lane_width - 0.04,
                          y0=y - 0.002, y1=y + 0.002,
                          fillcolor="white", opacity=intensity * 0.3, line=dict(width=0))
            fig.add_trace(go.Scatter(
                x=[lane_x], y=[y], mode="markers",
                marker=dict(size=18, color="white", opacity=0),
                hovertemplate=f"<b>{result['enzymes']}</b><br>{frag} bp<extra></extra>",
                showlegend=False))

    n_lanes = len(results) + 1
    fig.update_layout(
        paper_bgcolor="#1a1a2e", plot_bgcolor="#1a1a2e",
        title=dict(
            text=f"Restriction Digest Simulation — Plasmid {plasmid_size} bp{title_suffix}",
            font=dict(color="white", size=14), x=0.5),
        xaxis=dict(showticklabels=False, showgrid=False, zeroline=False,
                   range=[-1.2, len(results) + 0.5]),
        yaxis=dict(showticklabels=False, showgrid=False, zeroline=False,
                   range=[y_min - 0.1, y_max + 0.3]),
        height=600,
        width=max(600, n_lanes * 75),
        margin=dict(t=100, b=20, l=80, r=20),
        hovermode="closest")
    return fig

# ── SIDEBAR ────────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("⚙️ Parameters")

    uploaded_file = st.file_uploader(
        "Upload plasmid sequence",
        type=["fasta", "fa", "fas", "gb", "gbk", "genbank", "sbd"],
        help="Supported formats: FASTA, GenBank, SeqBuilder (.sbd)")

    st.caption("💡 Without a file the tool runs in demo mode with a random 10kb plasmid.")

    run = st.button("▶  Run Analysis", type="primary", use_container_width=True)

    st.divider()
    st.subheader("🔧 Settings")

    min_frag = st.slider("Min fragment size (bp)", 100, 3000, 250, 50)
    max_frag = st.slider("Max fragment size (bp)", 1000, 50000, 8000, 500)
    min_frags = st.slider("Min number of bands (n)", 1, 8, 3)
    max_frags = st.slider("Max number of bands", 2, 10, 6)
    min_diff = st.slider("Min size difference", 0.05, 0.5, 0.15, 0.05)
    prefer_short = st.checkbox("Prefer short fragments (faster on gel)", value=False)
    combo_min = st.slider("Min enzymes per digest", 1, 3, 1)
    combo_max = st.slider("Max enzymes per digest", 1, 3, 2)
    top_n = st.slider("Top N results", 1, 20, 10)

    st.divider()
    st.subheader("🧪 Select enzymes")
    select_all = st.checkbox("Select all", value=True)
    if select_all:
        selected_enzymes = DEFAULT_ENZYMES
    else:
        selected_enzymes = st.multiselect(
            "Enzymes:", DEFAULT_ENZYMES, default=DEFAULT_ENZYMES[:10])

# ── MAIN AREA ──────────────────────────────────────────────────────────────────
if run:
    plasmid_seq, demo_modus = load_sequence(uploaded_file)

    if plasmid_seq is None:
        st.error("❌ Could not read file — please check the format!")
        st.stop()

    plasmid_size = len(plasmid_seq)

    if demo_modus:
        st.warning("⚠️ No file uploaded — running in demo mode with random 10kb plasmid")
    else:
        st.success(f"✅ Plasmid loaded: {plasmid_size} bp")

    for min_f in [min_frags, min_frags + 1]:
        st.divider()
        st.subheader(f"Analysis — minimum {min_f} bands")

        with st.spinner(f"Searching combinations with min {min_f} bands..."):
            best, cutting = find_best_digests(
                plasmid_seq, selected_enzymes,
                min_frag, max_frag, min_f, max_frags,
                min_diff, (combo_min, combo_max), top_n,
                prefer_short=prefer_short)

        st.caption(f"Enzymes that cut: {', '.join(e.__name__ for e in cutting)}")

        if not best:
            st.error("No combinations found — try adjusting the parameters!")
        else:
            df = pd.DataFrame([{
                "#": i+1,
                "Enzymes": r["enzymes"],
                "Fragments (bp)": ",  ".join(str(f) for f in r["fragments"]),
                "Number of bands": r["n"]
            } for i, r in enumerate(best)])
            st.dataframe(df, use_container_width=True, hide_index=True)

            fig = draw_gel(best, plasmid_size, title_suffix=f" (min {min_f} bands)")
            st.plotly_chart(fig, use_container_width=True)
else:
    st.info("👈 Set parameters and click **Run Analysis**!")
    st.markdown("""
    **Supported file formats:**
    - `.fasta` / `.fa` / `.fas`
    - `.gb` / `.gbk` / `.genbank`
    - `.sbd` (SeqBuilder Pro)

    **Without a file** the tool runs in demo mode with a random 10kb sequence.

    **Tip:** Start with max 2 enzymes per digest for faster results.
    """)
