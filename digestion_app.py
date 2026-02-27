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
from scipy.optimize import linear_sum_assignment

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
st.sidebar.markdown("""
<div style="
    background: linear-gradient(135deg, #1e1e3f 0%, #2d1b69 100%);
    border: 1px solid #7c3aed;
    border-radius: 10px;
    padding: 1rem 1.1rem 0.6rem 1.1rem;
    margin-bottom: 0.8rem;
">
    <div style="font-size: 0.7rem; font-weight: 700; color: #a78bfa;
                letter-spacing: 0.12em; text-transform: uppercase; margin-bottom: 0.5rem;">
        Select Analysis Module
    </div>
    <div style="font-size: 1.25rem; font-weight: 800; color: #ffffff; line-height: 1.3;">
        🧬 Plasmid Analysis Suite
    </div>
</div>
""", unsafe_allow_html=True)

tool = st.sidebar.radio(
    label="",
    options=["Restriction Digest Planner", "Multi-Plasmid Comparator"],
    index=0,
    format_func=lambda x: (
        "🧪  Restriction Digest Planner" if x == "Restriction Digest Planner"
        else "🔀  Multi-Plasmid Comparator"
    )
)

st.sidebar.markdown("""<style>
div[role="radiogroup"] label {
    font-size: 1.05rem !important;
    font-weight: 600 !important;
    padding: 0.35rem 0 !important;
}
</style>""", unsafe_allow_html=True)

tool_descriptions = {
    "Restriction Digest Planner": "Identifies optimal enzyme combinations for diagnostic restriction analysis of circular plasmids. Results ranked by gel separation quality.",
    "Multi-Plasmid Comparator": "Compares restriction digest patterns across multiple constructs and identifies enzyme combinations that discriminate between plasmids.",
}

st.sidebar.caption(tool_descriptions[tool])
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
MAX_FILE_SIZE_BYTES = 5 * 1024 * 1024  # 5 MB

def read_sbd(raw, display_name="file"):
    try:
        text = raw.decode("latin-1")
        matches = re.findall(r'[ATCGatcg]{100,}', text)
        if not matches:
            st.warning(f"⚠️ **{display_name}.sbd**: No DNA sequence found. "
                       "The file may be corrupted or use an unsupported SeqBuilder version.")
            return None
        seq = max(matches, key=len).upper()
        if len(matches) > 1:
            total_found = sum(len(m) for m in matches)
            st.warning(f"⚠️ **{display_name}.sbd**: Found {len(matches)} sequence fragments "
                       f"({total_found:,} bp total). Using longest ({len(seq):,} bp). "
                       "If this seems wrong, try exporting as FASTA from SeqBuilder.")
        return seq
    except Exception as e:
        st.warning(f"⚠️ **{display_name}.sbd**: Failed to parse — {e}")
        return None

def load_sequence(uploaded_file):
    if uploaded_file is None:
        seq = "".join(random.choices(list("ATCG"), k=10000))
        return seq, True, "Demo Plasmid"

    raw = uploaded_file.read()
    name = uploaded_file.name.lower()
    display_name = uploaded_file.name.rsplit(".", 1)[0]

    # ── File size check (5 MB max) ──
    if len(raw) > MAX_FILE_SIZE_BYTES:
        st.error(f"❌ **{uploaded_file.name}** is too large "
                 f"({len(raw)/1024/1024:.1f} MB). Maximum allowed size is 5 MB. "
                 "For large sequences, export only the relevant region as FASTA.")
        return None, False, display_name

    if name.endswith(".sbd"):
        seq = read_sbd(raw, display_name)
        if seq:
            return seq, False, display_name
        return None, False, display_name

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
@st.cache_data(show_spinner=False)
def get_fragments(plasmid_seq, enzyme_names, plasmid_size):
    if not isinstance(plasmid_seq, Seq):
        plasmid_seq = Seq(str(plasmid_seq).upper())
    enzymes = [getattr(Restriction, e) for e in enzyme_names if hasattr(Restriction, e)]
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
    """Score how evenly bands are distributed — in log space (gel-physical)."""
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
    # Log-space CV: reflects visual band spacing on gel
    log_frags = np.log10(fragments)
    return np.std(log_frags) / np.mean(log_frags)

@st.cache_data(show_spinner=False)
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
            frags = get_fragments(str(plasmid_seq), tuple(e.__name__ for e in combo), plasmid_size)
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

# ── PLASMID MAP ────────────────────────────────────────────────────────────────
def draw_plasmid_map(plasmid_seq, plasmid_size, plasmid_name, enzyme_list):
    resolved = [getattr(Restriction, e) for e in enzyme_list if hasattr(Restriction, e)]
    if not resolved:
        return None
    rb = Restriction.RestrictionBatch(resolved)
    search_results = rb.search(Seq(plasmid_seq), linear=False)
    cutting = {e.__name__: sites for e, sites in search_results.items() if sites}
    if not cutting:
        return None

    fig = go.Figure()

    enz_colors = ["#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
                  "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3", "#ff6b6b", "#4ecdc4"]

    # Position 1 = top (12 o'clock), clockwise
    def pos_to_angle(pos):
        return np.pi / 2 - 2 * np.pi * (pos - 1) / plasmid_size

    # Backbone: clockwise from top
    theta = np.linspace(np.pi / 2, np.pi / 2 - 2 * np.pi, 500)
    fig.add_trace(go.Scatter(
        x=np.cos(theta), y=np.sin(theta),
        mode="lines",
        line=dict(color="#555577", width=12),
        hoverinfo="skip", showlegend=False))

    all_cuts = sorted(set(
        site for sites in cutting.values() for site in sites
    ))

    if all_cuts:
        n_cuts = len(all_cuts)
        # Distinct perceptually-spaced colors for fragments
        seg_palette = [
            "#4e79a7", "#f28e2b", "#e15759", "#76b7b2",
            "#59a14f", "#edc948", "#b07aa1", "#ff9da7",
            "#9c755f", "#bab0ac", "#d37295", "#a0cbe8",
        ]
        # Generate enough colors for any number of fragments
        seg_colors = [seg_palette[i % len(seg_palette)] for i in range(n_cuts)]
        for idx in range(n_cuts):
            start = all_cuts[idx]
            end = all_cuts[(idx + 1) % n_cuts]
            frag_size = end - start if end > start else plasmid_size - start + end

            a0 = pos_to_angle(start)
            a1 = pos_to_angle(end)

            if a1 >= a0:
                a1 -= 2 * np.pi

            n_pts = max(30, int(abs(a1 - a0) / (2 * np.pi) * 300) + 2)
            arc_angles = np.linspace(a0, a1, n_pts)

            r_inner, r_outer = 0.80, 1.00
            x_outer = r_outer * np.cos(arc_angles)
            y_outer = r_outer * np.sin(arc_angles)
            x_inner = r_inner * np.cos(arc_angles[::-1])
            y_inner = r_inner * np.sin(arc_angles[::-1])

            seg_color = seg_colors[idx]
            fig.add_trace(go.Scatter(
                x=np.concatenate([x_outer, x_inner]),
                y=np.concatenate([y_outer, y_inner]),
                fill="toself", fillcolor=seg_color,
                line=dict(width=0), mode="lines",
                text=f"<b>Fragment {idx+1}</b><br>Size: <b>{frag_size} bp</b><br>From: {start} bp → {all_cuts[(idx+1) % n_cuts]} bp",
                hoverinfo="text",
                hoveron="fills",
                showlegend=False, opacity=0.6))

    for enz_idx, (enz_name, sites) in enumerate(cutting.items()):
        color = enz_colors[enz_idx % len(enz_colors)]
        for site in sites:
            angle = pos_to_angle(site)
            fig.add_shape(type="line",
                          x0=0.76 * np.cos(angle), y0=0.76 * np.sin(angle),
                          x1=1.04 * np.cos(angle), y1=1.04 * np.sin(angle),
                          line=dict(color=color, width=2.5))
            # Hover marker placed outside the ring to avoid covering small fragments
            fig.add_trace(go.Scatter(
                x=[1.12 * np.cos(angle)], y=[1.12 * np.sin(angle)],
                mode="markers", marker=dict(size=6, color=color, opacity=0.01),
                hovertemplate=(
                    f"<b>{enz_name}</b><br>"
                    f"Position: {site} bp<br>"
                    f"({site/plasmid_size*100:.1f}% of plasmid)"
                    "<extra></extra>"),
                showlegend=False))
            fig.add_annotation(
                x=1.22 * np.cos(angle), y=1.22 * np.sin(angle),
                text=f"<b>{enz_name}</b><br><span style='font-size:9px;color:#aaa'>{site} bp</span>",
                showarrow=False, font=dict(color=color, size=10),
                xanchor="center", yanchor="middle")

        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode="markers",
            marker=dict(size=10, color=color, symbol="square"),
            name=f"{enz_name} ({len(sites)}×)", showlegend=True))

    fig.add_annotation(x=0, y=1.10, text="<b>1</b>",
                       showarrow=False, font=dict(color="#ffffff", size=11),
                       xanchor="center", yanchor="bottom")
    fig.add_shape(type="line", x0=0, y0=0.96, x1=0, y1=1.06,
                  line=dict(color="#ffffff", width=1.5, dash="dot"))
    fig.add_annotation(x=0, y=0.12, text=f"<b>{plasmid_name}</b>",
                       showarrow=False,
                       font=dict(color="white", size=13, family="Arial Black"),
                       xanchor="center")
    fig.add_annotation(x=0, y=-0.12, text=f"{plasmid_size:,} bp",
                       showarrow=False, font=dict(color="#aaaaaa", size=11),
                       xanchor="center")

    fig.update_layout(
        paper_bgcolor="#1a1a2e", plot_bgcolor="#1a1a2e",
        title=dict(text=f"Restriction Site Map — {plasmid_name}",
                   font=dict(color="white", size=14, family="Arial Black"), x=0.5),
        xaxis=dict(showticklabels=False, showgrid=False, zeroline=False, range=[-2.0, 2.0]),
        yaxis=dict(showticklabels=False, showgrid=False, zeroline=False,
                   range=[-2.0, 2.0], scaleanchor="x"),
        height=580, width=580, showlegend=True,
        legend=dict(font=dict(color="white", size=10),
                    bgcolor="#2a2a4a", bordercolor="#555", borderwidth=1),
        margin=dict(t=60, b=20, l=20, r=20), hovermode="closest")
    return fig



def parse_pasted_sequence(text, name="Pasted Sequence"):
    """Parse a raw sequence or FASTA string pasted by the user."""
    if not text or not text.strip():
        return None, name
    text = text.strip()
    # Try FASTA
    if text.startswith(">"):
        try:
            record = SeqIO.read(io.StringIO(text), "fasta")
            rec_name = record.id if record.id else name
            seq = str(record.seq).upper()
            seq = re.sub(r'[^ATCG]', '', seq)
            if len(seq) > 100:
                return seq, rec_name
        except:
            pass
    # Raw sequence — strip whitespace/numbers/non-DNA chars
    seq = re.sub(r'[^ATCGatcg]', '', text).upper()
    if len(seq) > 100:
        return seq, name
    return None, name

# ══════════════════════════════════════════════════════════════════════════════
# TOOL 1 — RESTRICTION DIGEST PLANNER
# ══════════════════════════════════════════════════════════════════════════════
if tool == "Restriction Digest Planner":
    st.markdown("### 🧪 Restriction Digest Planner")
    st.markdown("Automatically identifies optimal enzyme combinations for diagnostic restriction analysis of circular plasmids.")

    with st.sidebar:
        st.header("⚙️ Parameters")

        uploaded_file = st.file_uploader(
            "Upload plasmid sequence",
            type=["fasta", "fa", "fas", "gb", "gbk", "genbank", "sbd"],
            help="Accepted formats: FASTA, GenBank, SeqBuilder Pro (.sbd)",
            key="uploader_1")

        pasted_seq_1 = st.text_area(
            "— or paste sequence directly —",
            placeholder="Paste raw DNA sequence or FASTA here…",
            height=80,
            key="paste_1",
            help="Accepts raw sequence (ATCG) or FASTA format. Takes priority over file upload if both are provided."
        )

        st.caption("💡 No input required — without a sequence the tool runs in demo mode using a randomly generated 10 kb plasmid.")

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
        # Pasted sequence takes priority over file upload
        if pasted_seq_1 and pasted_seq_1.strip():
            _seq, _name = parse_pasted_sequence(pasted_seq_1, "Pasted Sequence")
            if _seq:
                plasmid_seq, demo_modus, plasmid_name = _seq, False, _name
            else:
                st.error("❌ Could not parse pasted sequence. Please check the input.")
                st.stop()
        else:
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
                    plasmid_seq, tuple(selected_enzymes),
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

                with st.expander("🗺️ Show plasmid restriction map for top result"):
                    top_enzymes = best[0]["enzyme_list"]
                    fig_map = draw_plasmid_map(
                        plasmid_seq, plasmid_size, plasmid_name, top_enzymes)
                    if fig_map:
                        st.plotly_chart(fig_map, use_container_width=False,
                                        key=f"map_t1_{min_f}")
                    else:
                        st.write("No cut sites to display.")
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
    st.markdown("### 🔀 Multi-Plasmid Comparator")
    st.markdown("Compare restriction digest patterns across multiple plasmids — identifies enzymes that discriminate between constructs.")

    with st.sidebar:
        st.header("⚙️ Parameters")

        uploaded_files = st.file_uploader(
            "Upload 2 or more plasmid sequences",
            type=["fasta", "fa", "fas", "gb", "gbk", "genbank", "sbd"],
            accept_multiple_files=True,
            help="Upload at least 2 plasmids to compare",
            key="uploader_2")

        st.caption("Optionally paste additional sequences below.")
        n_paste = st.number_input(
            "Number of pasted sequences", min_value=0, max_value=8, value=2, step=1,
            key="n_paste_2",
            help="Add up to 8 additional sequences by pasting directly.")
        pasted_seqs_2 = []
        labels = "ABCDEFGH"
        for _i in range(int(n_paste)):
            _label = labels[_i]
            _txt = st.text_area(
                f"Paste sequence {_label}",
                placeholder="Paste raw DNA or FASTA…",
                height=68, key=f"paste_2_{_label}",
                help="Name will be taken from FASTA header if present.")
            pasted_seqs_2.append((_txt, f"Sequence {_label}"))

        run2 = st.button("▶  Run Comparison", type="primary", use_container_width=True, key="run_2")

        st.divider()
        st.subheader("🔧 Analysis Settings")
        min_frag2 = st.slider("Minimum fragment size (bp)", 100, 3000, 250, 50, key="mf2")
        max_frag2 = st.slider("Maximum fragment size (bp)", 1000, 50000, 12000, 500, key="xf2")
        min_frags2 = st.slider("Minimum number of bands (n)", 1, 8, 1, key="mfr2")
        max_frags2 = st.slider("Maximum number of bands", 2, 10, 6, key="xfr2")
        min_diff2 = st.slider("Minimum relative size difference", 0.05, 0.5, 0.05, 0.05, key="md2")
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

        plasmids = []
        for f in uploaded_files:
            seq, demo, name = load_sequence(f)
            if seq:
                plasmids.append({"seq": seq, "name": name, "size": len(seq)})
        # Add pasted sequences
        for raw_paste, default_name in [(pasted_seq_2a, "Sequence A"), (pasted_seq_2b, "Sequence B")]:
            if raw_paste and raw_paste.strip():
                seq, name = parse_pasted_sequence(raw_paste, default_name)
                if seq:
                    plasmids.append({"seq": seq, "name": name, "size": len(seq)})
                else:
                    st.warning(f"⚠️ Could not parse pasted {default_name} — skipped.")

        if len(plasmids) < 2:
            st.error("❌ Could not read at least 2 files. Please check the file formats.")
            st.stop()

        st.success(f"✅ Loaded {len(plasmids)} plasmids: {', '.join(p['name'] for p in plasmids)}")

        resolved = [getattr(Restriction, e) for e in selected_enzymes2 if hasattr(Restriction, e)]
        cut_data = {}
        for p in plasmids:
            rb = Restriction.RestrictionBatch(resolved)
            results = rb.search(Seq(p["seq"]), linear=False)
            cut_data[p["name"]] = {e.__name__: sites for e, sites in results.items() if sites}

        all_cutting = set()
        for d in cut_data.values():
            all_cutting.update(d.keys())

        # Discriminating = absent in some plasmids OR different cut count
        discriminating = []
        universal_same = []

        for enz in all_cutting:
            counts = []
            for p in plasmids:
                if enz in cut_data[p["name"]]:
                    counts.append(len(cut_data[p["name"]][enz]))
                else:
                    counts.append(0)
            cuts_in_all = all(c > 0 for c in counts)
            all_same_count = len(set(counts)) == 1
            if not cuts_in_all or not all_same_count:
                discriminating.append(enz)
            else:
                universal_same.append(enz)

        st.divider()
        st.subheader("🔍 Enzyme Overview")
        col1, col2, col3 = st.columns(3)
        with col1:
            st.markdown("**✅ Identical in ALL plasmids**")
            st.write(", ".join(sorted(universal_same)) if universal_same else "None")
        with col2:
            st.markdown("**⚡ Discriminating (absent in some)**")
            absent_disc = [e for e in discriminating
                           if not all(e in cut_data[p["name"]] for p in plasmids)]
            st.write(", ".join(sorted(absent_disc)) if absent_disc else "None")
        with col3:
            st.markdown("**🔢 Discriminating (different cut count)**")
            count_disc = [e for e in discriminating
                          if all(e in cut_data[p["name"]] for p in plasmids)]
            st.write(", ".join(sorted(count_disc)) if count_disc else "None")

        if not discriminating:
            st.info("All enzymes produce identical patterns across all plasmids.")
        else:
            st.success(f"Found **{len(discriminating)} discriminating enzyme(s)** that produce different patterns.")

        st.divider()
        st.subheader("📊 Cut Pattern Matrix")
        matrix_data = {}
        for enz in sorted(all_cutting):
            row = {}
            for p in plasmids:
                if enz in cut_data[p["name"]]:
                    row[p["name"]] = f"✅ {len(cut_data[p['name']][enz])}×"
                else:
                    row[p["name"]] = "❌"
            matrix_data[enz] = row
        matrix_df = pd.DataFrame(matrix_data).T
        matrix_df.index.name = "Enzyme"

        def highlight_discriminating(row):
            enz = row.name
            if enz in discriminating:
                return ["background-color: #2a3a2a; font-weight: bold"] * len(row)
            return [""] * len(row)

        st.dataframe(
            matrix_df.style.apply(highlight_discriminating, axis=1),
            use_container_width=True
        )
        st.caption("🟢 Highlighted rows = discriminating enzymes (different pattern between plasmids)")

        st.divider()
        st.subheader("🏆 Best Discriminating Digest Combinations")

        disc_resolved = [getattr(Restriction, e) for e in discriminating if hasattr(Restriction, e)]
        also_use = [getattr(Restriction, e) for e in universal_same if hasattr(Restriction, e)]
        all_for_combo = disc_resolved + also_use

        combo_scores = []
        for size in range(1, combo_max2 + 1):
            for combo in combinations(all_for_combo, size):
                combo_names = [e.__name__ for e in combo]
                if not any(e in discriminating for e in combo_names):
                    continue

                frag_sets = []
                valid = True
                for p in plasmids:
                    frags = get_fragments(p["seq"], tuple(e.__name__ for e in combo), p["size"])
                    s = score_combination(frags, min_frag2, max_frag2,
                                          min_frags2, max_frags2, min_diff2)
                    if s is None:
                        valid = False
                        break
                    frag_sets.append(frags)

                if not valid:
                    continue

                # Gel-physical discrimination score: log-space Hungarian matching
                diff_score = 0
                for i in range(len(frag_sets)):
                    for j in range(i+1, len(frag_sets)):
                        fa = sorted(frag_sets[i])
                        fb = sorted(frag_sets[j])
                        na, nb = len(fa), len(fb)
                        if na != nb:
                            n = max(na, nb)
                            fa_pad = fa + [0] * (n - na)
                            fb_pad = fb + [0] * (n - nb)
                            cost = np.array([
                                [abs(np.log10(max(a, 1)) - np.log10(max(b, 1))) if (a > 0 and b > 0) else np.log10(max(max(fa), max(fb), 1))
                                 for b in fb_pad]
                                for a in fa_pad
                            ])
                            row_ind, col_ind = linear_sum_assignment(cost)
                            matched = cost[row_ind, col_ind].sum()
                            unmatched_penalty = abs(na - nb) * 1.5
                            diff_score += matched + unmatched_penalty
                        else:
                            cost = np.array([
                                [abs(np.log10(a) - np.log10(b)) for b in fb]
                                for a in fa
                            ])
                            row_ind, col_ind = linear_sum_assignment(cost)
                            diff_score += cost[row_ind, col_ind].sum()

                combo_scores.append({
                    "enzymes": " + ".join(combo_names),
                    "enzyme_list": combo_names,
                    "frag_sets": frag_sets,
                    "diff_score": diff_score
                })

        combo_scores.sort(key=lambda x: x["diff_score"], reverse=True)
        top_combos = combo_scores[:top_n2]

        if not top_combos:
            st.warning(
                "No discriminating combinations passed the current filter settings. "
                "Try relaxing: lower minimum band count, increase max fragment size, "
                "or reduce the minimum size difference between bands."
            )
            if discriminating:
                st.info(
                    f"Discriminating enzymes found: **{', '.join(sorted(discriminating))}** — "
                    "but their fragments don't meet the current size/band-count filters. "
                    "Try adjusting the Analysis Settings in the sidebar."
                )
        else:
            rows = []
            for i, c in enumerate(top_combos):
                row = {"Rank": i+1, "Enzyme(s)": c["enzymes"],
                       "Discrimination score": f"{c['diff_score']:.3f}"}
                for j, p in enumerate(plasmids):
                    row[f"{p['name']} fragments (bp)"] = ", ".join(str(f) for f in c["frag_sets"][j])
                rows.append(row)
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

            st.subheader("🧫 Predicted Gel — All Plasmids per Combination")
            max_size = max(p["size"] for p in plasmids)

            for i, combo in enumerate(top_combos[:5]):
                st.markdown(f"**#{i+1} — {combo['enzymes']}** (discrimination score: {combo['diff_score']:.3f})")
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

                with st.expander(f"🗺️ Show restriction maps for this digest combination"):
                    map_cols = st.columns(len(plasmids))
                    for j, p in enumerate(plasmids):
                        with map_cols[j]:
                            st.markdown(f"**{p['name']}**")
                            fig_map = draw_plasmid_map(
                                p["seq"], p["size"], p["name"], combo["enzyme_list"])
                            if fig_map:
                                st.plotly_chart(fig_map, use_container_width=True,
                                                key=f"map_t2_{i}_{j}")
                            else:
                                st.write("No cut sites.")
    else:
        st.info("👈 Upload at least 2 plasmid sequences and click **Run Comparison**.")
        st.markdown("""
        **What this tool does:**
        - Identifies enzymes that cut some plasmids but not others
        - Finds the best enzyme combinations to distinguish between constructs
        - Shows a predicted gel with all plasmids side by side

        **Use case:** Verify correct construct after cloning by comparing expected vs. actual digest pattern.
        """)
