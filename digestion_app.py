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
import math

# ── PAGE SETTINGS ──────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Plasmid Analysis Suite",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 Plasmid Analysis Suite")

st.markdown("""
<style>
    section[data-testid="stSidebar"] {
        width: 700px !important;
        min-width: 0px !important;
    }
</style>
""", unsafe_allow_html=True)

# ── TOOL SELECTION ─────────────────────────────────────────────────────────────
tool = st.sidebar.radio(
    "🔬 Select Tool",
    ["Restriction Digest Planner", "Multi-Plasmid Comparator", "Cut Site Map"],
    index=0
)

st.sidebar.divider()

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
            return None, None
        seq = max(matches, key=len).upper()
        # Try to extract name
        name_match = re.search(r'Written by SeqBuilder.*?([A-Za-z0-9_\-\.]+)\s*$',
                                text[:500], re.MULTILINE)
        return seq, None
    except:
        return None, None

def load_sequence(uploaded_file):
    if uploaded_file is None:
        seq = "".join(random.choices(list("ATCG"), k=10000))
        return seq, True, "Demo Plasmid"

    raw = uploaded_file.read()
    name = uploaded_file.name.lower()
    display_name = uploaded_file.name.rsplit(".", 1)[0]

    if name.endswith(".sbd"):
        seq, _ = read_sbd(raw)
        if seq:
            return seq, False, display_name

    for fmt in ["fasta", "genbank"]:
        try:
            record = SeqIO.read(io.StringIO(raw.decode("utf-8")), fmt)
            rec_name = record.id if record.id and record.id != "." else display_name
            return str(record.seq).upper(), False, rec_name
        except:
            pass
        try:
            record = SeqIO.read(io.StringIO(raw.decode("utf-16")), fmt)
            rec_name = record.id if record.id and record.id != "." else display_name
            return str(record.seq).upper(), False, rec_name
        except:
            pass

    for enc in ["utf-8", "utf-16", "latin-1"]:
        try:
            text = raw.decode(enc)
            matches = re.findall(r'[ATCGatcg]{200,}', text)
            if matches:
                return max(matches, key=len).upper(), False, display_name
        except:
            pass

    return None, False, display_name

# ── ANALYSIS FUNCTIONS ─────────────────────────────────────────────────────────
def get_cut_sites(plasmid_seq, enzymes):
    rb = Restriction.RestrictionBatch(enzymes)
    return rb.search(plasmid_seq, linear=False)

def get_fragments(plasmid_seq, enzymes, plasmid_size):
    search_results = get_cut_sites(plasmid_seq, enzymes)
    cut_sites = []
    for enz, sites in search_results.items():
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
                    "enzyme_list": [e.__name__ for e in combo],
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
def draw_gel(results, plasmid_size, title_suffix="", lane_labels=None):
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
                                     size=11 if is_thick else 10),
                           xanchor="right", yanchor="middle")

    fig.add_annotation(x=0, y=y_max + 0.08, text="GeneRuler 1kb",
                       showarrow=False,
                       font=dict(color="white", size=12, family="Arial Black"),
                       xanchor="center", yanchor="bottom")

    # ── Sample lanes ──
    for i, result in enumerate(results):
        lane_x = i + 1
        color = colors[i % len(colors)]

        top_label = lane_labels[i] if lane_labels else f"<b>#{i+1}</b>"
        fig.add_annotation(
            x=lane_x, y=y_max + 0.08,
            text=top_label,
            showarrow=False,
            font=dict(color=color, size=12, family="Arial Black"),
            xanchor="center", yanchor="bottom")

        fig.add_annotation(
            x=lane_x, y=y_max + 0.02,
            text=result["enzymes"].replace(" + ", "<br>+ "),
            showarrow=False,
            font=dict(color=color, size=11, family="Arial"),
            xanchor="center", yanchor="top",
            bgcolor="#2a2a4a",
            bordercolor=color,
            borderwidth=1.5,
            borderpad=5)

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
                marker=dict(size=20, color="white", opacity=0),
                hovertemplate=(
                    f"<b>{top_label} — {result['enzymes']}</b><br>"
                    f"Fragment: {frag} bp<br>"
                    f"Total bands: {result['n']}"
                    "<extra></extra>"),
                showlegend=False))

    n_lanes = len(results) + 1
    fig.update_layout(
        paper_bgcolor="#1a1a2e",
        plot_bgcolor="#1a1a2e",
        title=dict(
            text=f"Predicted Restriction Digest Pattern — {plasmid_size} bp plasmid{title_suffix}",
            font=dict(color="white", size=15, family="Arial Black"),
            x=0.5),
        xaxis=dict(showticklabels=False, showgrid=False, zeroline=False,
                   range=[-1.2, len(results) + 0.5]),
        yaxis=dict(showticklabels=False, showgrid=False, zeroline=False,
                   range=[y_min - 0.1, y_max + 0.35]),
        height=650,
        width=max(650, n_lanes * 85),
        margin=dict(t=120, b=20, l=90, r=20),
        hovermode="closest")
    return fig

# ── CUT SITE MAP ───────────────────────────────────────────────────────────────
def draw_cut_site_map(plasmid_seq, plasmid_size, plasmid_name, selected_enzymes):
    resolved = [getattr(Restriction, e) for e in selected_enzymes if hasattr(Restriction, e)]
    rb = Restriction.RestrictionBatch(resolved)
    search_results = rb.search(Seq(plasmid_seq), linear=False)
    cutting = {e.__name__: sites for e, sites in search_results.items() if sites}

    if not cutting:
        return None

    fig = go.Figure()
    colors = ["#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
              "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3", "#ff6b6b", "#4ecdc4"]

    # Draw plasmid circle
    theta = np.linspace(0, 2 * np.pi, 300)
    fig.add_trace(go.Scatter(
        x=np.cos(theta), y=np.sin(theta),
        mode="lines",
        line=dict(color="#aaaaaa", width=8),
        hoverinfo="skip", showlegend=False))

    # Draw cut sites
    for idx, (enz_name, sites) in enumerate(cutting.items()):
        color = colors[idx % len(colors)]
        for site in sites:
            angle = 2 * np.pi * site / plasmid_size - np.pi / 2
            inner = 0.85
            outer = 1.15
            x0 = inner * np.cos(angle)
            y0 = inner * np.sin(angle)
            x1 = outer * np.cos(angle)
            y1 = outer * np.sin(angle)

            fig.add_trace(go.Scatter(
                x=[x0, x1], y=[y0, y1],
                mode="lines",
                line=dict(color=color, width=2.5),
                hoverinfo="skip", showlegend=False))

            # Label
            label_r = 1.28
            fig.add_trace(go.Scatter(
                x=[label_r * np.cos(angle)],
                y=[label_r * np.sin(angle)],
                mode="markers+text",
                marker=dict(size=8, color=color),
                text=[f"{enz_name}<br>{site} bp"],
                textposition="middle center",
                textfont=dict(color=color, size=9),
                hovertemplate=f"<b>{enz_name}</b><br>Position: {site} bp<extra></extra>",
                showlegend=True,
                name=f"{enz_name} ({len(sites)}x)"))

    # Plasmid name in center
    fig.add_annotation(
        x=0, y=0.1,
        text=f"<b>{plasmid_name}</b>",
        showarrow=False,
        font=dict(color="white", size=14, family="Arial Black"),
        xanchor="center")
    fig.add_annotation(
        x=0, y=-0.1,
        text=f"{plasmid_size} bp",
        showarrow=False,
        font=dict(color="#aaaaaa", size=12),
        xanchor="center")

    fig.update_layout(
        paper_bgcolor="#1a1a2e",
        plot_bgcolor="#1a1a2e",
        title=dict(
            text=f"Restriction Site Map — {plasmid_name}",
            font=dict(color="white", size=15, family="Arial Black"),
            x=0.5),
        xaxis=dict(showticklabels=False, showgrid=False, zeroline=False,
                   range=[-1.7, 1.7]),
        yaxis=dict(showticklabels=False, showgrid=False, zeroline=False,
                   range=[-1.7, 1.7], scaleanchor="x"),
        height=650,
        width=700,
        showlegend=True,
        legend=dict(
            font=dict(color="white", size=10),
            bgcolor="#2a2a4a",
            bordercolor="#aaaaaa",
            borderwidth=1),
        margin=dict(t=80, b=20, l=20, r=20))
    return fig

# ══════════════════════════════════════════════════════════════════════════════
# TOOL 1 — RESTRICTION DIGEST PLANNER
# ══════════════════════════════════════════════════════════════════════════════
if tool == "Restriction Digest Planner":
    st.markdown("### 🔬 Restriction Digest Planner")
    st.markdown("Automatically identifies optimal enzyme combinations for diagnostic restriction analysis of circular plasmids.")

    with st.sidebar:
        st.header("⚙️ Parameters")

        uploaded_file = st.file_uploader(
            "Upload plasmid sequence",
            type=["fasta", "fa", "fas", "gb", "gbk", "genbank", "sbd"],
            help="Accepted formats: FASTA, GenBank, SeqBuilder Pro (.sbd)",
            key="uploader_1")

        st.caption("💡 No file required — without an upload the tool runs in demo mode using a randomly generated 10 kb plasmid.")

        run = st.button("▶  Run Analysis", type="primary", use_container_width=True, key="run_1")
        prefer_short = st.checkbox(
            "Prioritise short fragments",
            value=False,
            help="Ranks results by smallest largest fragment — minimises total gel run time")

        st.divider()
        st.subheader("🔧 Analysis Settings")
        min_frag = st.slider("Minimum fragment size (bp)", 100, 3000, 250, 50)
        max_frag = st.slider("Maximum fragment size (bp)", 1000, 50000, 8000, 500)
        min_frags = st.slider("Minimum number of bands (n)", 1, 8, 3,
                              help="Results will be displayed for n and n+1 bands simultaneously")
        max_frags = st.slider("Maximum number of bands", 2, 10, 6)
        min_diff = st.slider("Minimum relative size difference between adjacent bands",
                             0.05, 0.5, 0.15, 0.05)
        combo_min = st.slider("Minimum enzymes per digest", 1, 3, 1)
        combo_max = st.slider("Maximum enzymes per digest", 1, 3, 2,
                              help="⚠️ Setting this to 3 significantly increases computation time")
        top_n = st.slider("Number of results to display", 1, 20, 10)

        st.divider()
        st.subheader("🧪 Enzyme Selection")
        select_all = st.checkbox("Select all enzymes", value=True, key="sel_all_1")
        if select_all:
            selected_enzymes = DEFAULT_ENZYMES
        else:
            selected_enzymes = st.multiselect(
                "Select enzymes available in your laboratory:",
                DEFAULT_ENZYMES, default=DEFAULT_ENZYMES[:10], key="enz_1")

    if run:
        plasmid_seq, demo_modus, plasmid_name = load_sequence(uploaded_file)

        if plasmid_seq is None:
            st.error("❌ Unable to read the uploaded file. Please verify the file format and try again.")
            st.stop()

        plasmid_size = len(plasmid_seq)

        if demo_modus:
            st.warning("⚠️ No sequence file uploaded — running in demo mode with a randomly generated 10 kb plasmid.")
        else:
            st.success(f"✅ **{plasmid_name}** loaded successfully — {plasmid_size} bp")

        for min_f in [min_frags, min_frags + 1]:
            st.divider()
            st.subheader(f"Results — minimum {min_f} bands")

            with st.spinner(f"Searching for optimal enzyme combinations yielding at least {min_f} bands..."):
                best, cutting = find_best_digests(
                    plasmid_seq, selected_enzymes,
                    min_frag, max_frag, min_f, max_frags,
                    min_diff, (combo_min, combo_max), top_n,
                    prefer_short=prefer_short)

            st.caption(f"Enzymes with at least one recognition site: "
                       f"{', '.join(e.__name__ for e in cutting)}")

            if not best:
                st.error("No suitable enzyme combinations found. Consider relaxing the analysis parameters.")
            else:
                df = pd.DataFrame([{
                    "Rank": i+1,
                    "Enzyme(s)": r["enzymes"],
                    "Fragment sizes (bp)": ",  ".join(str(f) for f in r["fragments"]),
                    "No. of bands": r["n"],
                    "Largest fragment (bp)": max(r["fragments"])
                } for i, r in enumerate(best)])
                st.dataframe(df, use_container_width=True, hide_index=True)

                fig = draw_gel(best, plasmid_size,
                               title_suffix=f" — {plasmid_name} — min {min_f} bands")
                st.plotly_chart(fig, use_container_width=True)
    else:
        st.info("👈 Configure parameters in the sidebar and click **Run Analysis** to begin.")
        st.markdown("""
        **Accepted file formats:** FASTA, GenBank, SeqBuilder Pro (.sbd)

        **How it works:** Evaluates all possible enzyme combinations and ranks them by how
        evenly fragments are distributed across the detectable size range. Results are shown
        for n and n+1 minimum bands simultaneously.

        **Tip:** Limit to max 2 enzymes per digest for faster results.
        """)

# ══════════════════════════════════════════════════════════════════════════════
# TOOL 2 — MULTI-PLASMID COMPARATOR
# ══════════════════════════════════════════════════════════════════════════════
elif tool == "Multi-Plasmid Comparator":
    st.markdown("### 🔬 Multi-Plasmid Comparator")
    st.markdown("Compare restriction digest patterns across multiple plasmids — identifies enzymes that discriminate between constructs.")

    with st.sidebar:
        st.header("⚙️ Parameters")

        st.subheader("📂 Upload Plasmids")
        uploaded_files = st.file_uploader(
            "Upload 2 or more plasmid sequences",
            type=["fasta", "fa", "fas", "gb", "gbk", "genbank", "sbd"],
            accept_multiple_files=True,
            help="Upload at least 2 plasmids to compare",
            key="uploader_2")

        run2 = st.button("▶  Run Comparison", type="primary", use_container_width=True, key="run_2")

        st.divider()
        st.subheader("🔧 Analysis Settings")
        min_frag2 = st.slider("Minimum fragment size (bp)", 100, 3000, 250, 50, key="mf2")
        max_frag2 = st.slider("Maximum fragment size (bp)", 1000, 50000, 8000, 500, key="xf2")
        min_frags2 = st.slider("Minimum number of bands (n)", 1, 8, 3, key="mfr2")
        max_frags2 = st.slider("Maximum number of bands", 2, 10, 6, key="xfr2")
        min_diff2 = st.slider("Minimum relative size difference", 0.05, 0.5, 0.15, 0.05, key="md2")
        combo_max2 = st.slider("Maximum enzymes per digest", 1, 3, 2,
                               help="⚠️ Setting this to 3 significantly increases computation time",
                               key="cm2")
        top_n2 = st.slider("Number of results to display", 1, 20, 10, key="tn2")

        st.divider()
        st.subheader("🧪 Enzyme Selection")
        select_all2 = st.checkbox("Select all enzymes", value=True, key="sel_all_2")
        if select_all2:
            selected_enzymes2 = DEFAULT_ENZYMES
        else:
            selected_enzymes2 = st.multiselect(
                "Select enzymes:", DEFAULT_ENZYMES,
                default=DEFAULT_ENZYMES[:10], key="enz_2")

    if run2:
        if len(uploaded_files) < 2:
            st.error("❌ Please upload at least 2 plasmid sequences to compare.")
            st.stop()

        # Load all plasmids
        plasmids = []
        for f in uploaded_files:
            seq, demo, name = load_sequence(f)
            if seq:
                plasmids.append({"seq": seq, "name": name, "size": len(seq)})

        if len(plasmids) < 2:
            st.error("❌ Could not read at least 2 files. Please check the file formats.")
            st.stop()

        st.success(f"✅ Loaded {len(plasmids)} plasmids: {', '.join(p['name'] for p in plasmids)}")

        # Find cutting enzymes per plasmid
        resolved = [getattr(Restriction, e) for e in selected_enzymes2 if hasattr(Restriction, e)]
        cut_data = {}
        for p in plasmids:
            rb = Restriction.RestrictionBatch(resolved)
            results = rb.search(Seq(p["seq"]), linear=False)
            cut_data[p["name"]] = {e.__name__: sites for e, sites in results.items() if sites}

        # Find discriminating enzymes (cut in some but not all)
        all_cutting = set()
        for d in cut_data.values():
            all_cutting.update(d.keys())

        discriminating = []
        universal = []
        unique_per_plasmid = {p["name"]: [] for p in plasmids}

        for enz in all_cutting:
            cuts_in = [p["name"] for p in plasmids if enz in cut_data[p["name"]]]
            if len(cuts_in) == len(plasmids):
                universal.append(enz)
            else:
                discriminating.append(enz)
                for p in plasmids:
                    if enz in cut_data[p["name"]] and len(cuts_in) < len(plasmids):
                        unique_per_plasmid[p["name"]].append(enz)

        st.divider()
        st.subheader("🔍 Enzyme Overview")

        col1, col2 = st.columns(2)
        with col1:
            st.markdown("**✅ Cut in ALL plasmids (universal)**")
            if universal:
                st.write(", ".join(sorted(universal)))
            else:
                st.write("None")

        with col2:
            st.markdown("**⚡ Discriminating enzymes (cut in some but not all)**")
            if discriminating:
                st.write(", ".join(sorted(discriminating)))
            else:
                st.write("None — all plasmids have identical cut patterns")

        st.divider()
        st.subheader("📊 Cut Pattern Matrix")
        matrix_data = {}
        for enz in sorted(all_cutting):
            row = {}
            for p in plasmids:
                if enz in cut_data[p["name"]]:
                    row[p["name"]] = f"✅ {len(cut_data[p['name']][enz])}x"
                else:
                    row[p["name"]] = "❌"
            matrix_data[enz] = row

        matrix_df = pd.DataFrame(matrix_data).T
        matrix_df.index.name = "Enzyme"
        st.dataframe(matrix_df, use_container_width=True)

        # Find best discriminating digest combos
        st.divider()
        st.subheader("🏆 Best Discriminating Digest Combinations")
        st.markdown("Enzyme combinations that produce **clearly different band patterns** across plasmids:")

        # For each combination, compute fragments per plasmid and score difference
        disc_resolved = [getattr(Restriction, e) for e in discriminating if hasattr(Restriction, e)]
        also_use = [getattr(Restriction, e) for e in universal if hasattr(Restriction, e)]
        all_for_combo = disc_resolved + also_use

        combo_scores = []
        for size in range(1, combo_max2 + 1):
            for combo in combinations(all_for_combo, size):
                # Must include at least one discriminating enzyme
                combo_names = [e.__name__ for e in combo]
                if not any(e in discriminating for e in combo_names):
                    continue

                frag_sets = []
                valid = True
                for p in plasmids:
                    frags = get_fragments(p["seq"], list(combo), p["size"])
                    s = score_combination(frags, min_frag2, max_frag2,
                                          min_frags2, max_frags2, min_diff2)
                    if s is None:
                        valid = False
                        break
                    frag_sets.append(frags)

                if not valid:
                    continue

                # Score: how different are the patterns?
                # Compare all pairs of plasmids
                diff_score = 0
                for i in range(len(frag_sets)):
                    for j in range(i+1, len(frag_sets)):
                        # Different number of bands = good
                        if len(frag_sets[i]) != len(frag_sets[j]):
                            diff_score += 2
                        else:
                            # Same number — check size differences
                            diffs = [abs(a-b) / max(a,b)
                                     for a, b in zip(sorted(frag_sets[i]), sorted(frag_sets[j]))]
                            diff_score += sum(diffs)

                combo_scores.append({
                    "enzymes": " + ".join(combo_names),
                    "frag_sets": frag_sets,
                    "diff_score": diff_score
                })

        combo_scores.sort(key=lambda x: x["diff_score"], reverse=True)
        top_combos = combo_scores[:top_n2]

        if not top_combos:
            st.error("No discriminating combinations found with current settings.")
        else:
            # Show table
            rows = []
            for i, c in enumerate(top_combos):
                row = {"Rank": i+1, "Enzyme(s)": c["enzymes"]}
                for j, p in enumerate(plasmids):
                    row[f"{p['name']} fragments (bp)"] = ", ".join(str(f) for f in c["frag_sets"][j])
                rows.append(row)
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

            # Show gel for each combo — all plasmids side by side
            st.subheader("🧫 Predicted Gel — All Plasmids per Combination")
            max_size = max(p["size"] for p in plasmids)

            for i, combo in enumerate(top_combos[:5]):  # show top 5 gels
                st.markdown(f"**#{i+1} — {combo['enzymes']}**")
                gel_results = []
                lane_labels = []
                for j, p in enumerate(plasmids):
                    gel_results.append({
                        "enzymes": combo["enzymes"],
                        "fragments": combo["frag_sets"][j],
                        "n": len(combo["frag_sets"][j])
                    })
                    lane_labels.append(f"<b>{p['name']}</b>")

                fig = draw_gel(gel_results, max_size,
                               title_suffix=f" — {combo['enzymes']}",
                               lane_labels=lane_labels)
                st.plotly_chart(fig, use_container_width=True)
    else:
        st.info("👈 Upload at least 2 plasmid sequences and click **Run Comparison**.")
        st.markdown("""
        **What this tool does:**
        - Identifies enzymes that cut some plasmids but not others (discriminating enzymes)
        - Finds the best enzyme combinations to distinguish between your constructs
        - Shows a predicted gel with all plasmids side by side for each top combination

        **Use case:** Verify correct construct after cloning by comparing expected vs. actual digest pattern.
        """)

# ══════════════════════════════════════════════════════════════════════════════
# TOOL 3 — CUT SITE MAP
# ══════════════════════════════════════════════════════════════════════════════
elif tool == "Cut Site Map":
    st.markdown("### 🔬 Restriction Site Map")
    st.markdown("Visualise all restriction enzyme recognition sites on a circular plasmid map.")

    with st.sidebar:
        st.header("⚙️ Parameters")

        uploaded_file3 = st.file_uploader(
            "Upload plasmid sequence",
            type=["fasta", "fa", "fas", "gb", "gbk", "genbank", "sbd"],
            help="Accepted formats: FASTA, GenBank, SeqBuilder Pro (.sbd)",
            key="uploader_3")

        run3 = st.button("▶  Generate Map", type="primary", use_container_width=True, key="run_3")

        st.divider()
        st.subheader("🧪 Enzyme Selection")
        select_all3 = st.checkbox("Select all enzymes", value=True, key="sel_all_3")
        if select_all3:
            selected_enzymes3 = DEFAULT_ENZYMES
        else:
            selected_enzymes3 = st.multiselect(
                "Select enzymes to display:",
                DEFAULT_ENZYMES, default=DEFAULT_ENZYMES[:10], key="enz_3")

        show_non_cutters = st.checkbox(
            "Show enzymes with no recognition site",
            value=False)

    if run3:
        plasmid_seq3, demo3, plasmid_name3 = load_sequence(uploaded_file3)

        if plasmid_seq3 is None:
            st.error("❌ Unable to read file.")
            st.stop()

        plasmid_size3 = len(plasmid_seq3)

        if demo3:
            st.warning("⚠️ Running in demo mode with a randomly generated 10 kb plasmid.")
        else:
            st.success(f"✅ **{plasmid_name3}** — {plasmid_size3} bp")

        # Cut site table
        resolved3 = [getattr(Restriction, e) for e in selected_enzymes3 if hasattr(Restriction, e)]
        rb3 = Restriction.RestrictionBatch(resolved3)
        sr3 = rb3.search(Seq(plasmid_seq3), linear=False)

        cutting3 = {e.__name__: sites for e, sites in sr3.items() if sites}
        non_cutting3 = [e.__name__ for e, sites in sr3.items() if not sites]

        col1, col2 = st.columns([2, 1])

        with col1:
            fig_map = draw_cut_site_map(
                plasmid_seq3, plasmid_size3, plasmid_name3, selected_enzymes3)
            if fig_map:
                st.plotly_chart(fig_map, use_container_width=False)
            else:
                st.warning("No enzymes cut this plasmid with the current selection.")

        with col2:
            st.markdown("**Recognition sites found:**")
            if cutting3:
                site_df = pd.DataFrame([{
                    "Enzyme": enz,
                    "Cuts": len(sites),
                    "Positions (bp)": ", ".join(str(s) for s in sites)
                } for enz, sites in sorted(cutting3.items())])
                st.dataframe(site_df, use_container_width=True, hide_index=True)
            else:
                st.write("None")

            if show_non_cutters and non_cutting3:
                st.markdown("**No recognition site:**")
                st.caption(", ".join(sorted(non_cutting3)))
    else:
        st.info("👈 Upload a plasmid sequence and click **Generate Map**.")
        st.markdown("""
        **What this tool does:**
        - Draws a circular plasmid map with all restriction site positions
        - Lists all enzymes with their cut positions
        - Hover over a site to see enzyme name and exact position

        **Tip:** Use enzyme selection to focus on enzymes relevant to your cloning strategy.
        """)
