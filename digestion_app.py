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

# PAGE SETTINGS
st.set_page_config(
    page_title="Plasmid Analysis Toolkit",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 Plasmid Analysis Toolkit")

st.markdown("""
<style>
    section[data-testid="stSidebar"] {
        width: 700px !important;
        min-width: 0px !important;
    }
</style>
""", unsafe_allow_html=True)

# MODULE SELECTION
# Module selection
if "active_tool" not in st.session_state:
    st.session_state.active_tool = None

st.sidebar.markdown("""
<style>
/* Card button base */
div[data-testid="stSidebar"] div[data-testid="stButton"] button {
    width: 100% !important;
    text-align: left !important;
    border-radius: 10px !important;
    padding: 0.75rem 1rem !important;
    margin-bottom: 0.5rem !important;
    font-size: 1.0rem !important;
    font-weight: 700 !important;
    line-height: 1.5 !important;
    white-space: normal !important;
    height: auto !important;
    transition: all 0.15s ease !important;
}
</style>
""", unsafe_allow_html=True)

st.sidebar.markdown(
    '<div style="font-size:2.5rem;font-weight:800;color:#ffffff;margin-bottom:0.2rem;line-height:1.1;">'
    '🧬</div>'
    '<div style="font-size:1.1rem;font-weight:800;color:#ffffff;margin-bottom:0.2rem;">'
    'Plasmid Analysis Toolkit</div>'
    '<div style="font-size:0.72rem;color:#a78bfa;font-weight:600;letter-spacing:0.1em;'
    'text-transform:uppercase;margin-bottom:1rem;">Select Analysis Module</div>',
    unsafe_allow_html=True)

_t1 = st.session_state.active_tool == "Restriction Digest Planner"
_t2 = st.session_state.active_tool == "Multi-Plasmid Comparator"
_t3 = st.session_state.active_tool == "Primer Analyzer"
_t4 = st.session_state.active_tool == "Feature Annotation Viewer"

if st.sidebar.button(
    ("✅ Restriction Digest Planner  ●  ACTIVE" if _t1
     else "🧪 Restriction Digest Planner"),
    key="btn_tool1",
    type="primary" if _t1 else "secondary",
    use_container_width=True,
    help="Identify optimal enzyme combinations for a single plasmid"):
    st.session_state.active_tool = "Restriction Digest Planner"
    st.rerun()

if st.sidebar.button(
    ("✅ Multi-Plasmid Comparator  ●  ACTIVE" if _t2
     else "🔀 Multi-Plasmid Comparator"),
    key="btn_tool2",
    type="primary" if _t2 else "secondary",
    use_container_width=True,
    help="Compare digest patterns and identify discriminating enzymes across constructs"):
    st.session_state.active_tool = "Multi-Plasmid Comparator"
    st.rerun()

if st.sidebar.button(
    ("✅ Primer Analyzer  ●  ACTIVE" if _t3
     else "🔬 Primer Analyzer"),
    key="btn_tool3",
    type="primary" if _t3 else "secondary",
    use_container_width=True,
    help="Calculate Tm, GC content, and check for hairpins and dimers"):
    st.session_state.active_tool = "Primer Analyzer"
    st.rerun()

if st.sidebar.button(
    ("✅ Feature Annotation Viewer  ●  ACTIVE" if _t4
     else "🗺️ Feature Annotation Viewer"),
    key="btn_tool4",
    type="primary" if _t4 else "secondary",
    use_container_width=True,
    help="Visualise annotated features on a plasmid map from GenBank files"):
    st.session_state.active_tool = "Feature Annotation Viewer"
    st.rerun()

tool = st.session_state.active_tool
st.sidebar.divider()

# ENZYMES
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

# DAM/DCM METHYLATION
# dam methylase blocks: GATC (EcoRV, ClaI, XbaI, MboI, BclI, BamHI overlap)
# dcm methylase blocks: CCWGG (EcoRII, SacII overlap)
DAM_BLOCKED = {"MboI", "BclI", "ClaI"}          # cut only in dam- strains
DCM_BLOCKED = {"SacII", "EcoRII"}                # cut only in dcm- strains
# Note: partial sensitivity — BamHI, EcoRI etc are NOT blocked by dam/dcm

# SEQUENCE LOADING
MAX_FILE_SIZE_BYTES = 5 * 1024 * 1024  # 5 MB

def read_sbd(raw, display_name="file"):
    """
    Parse SeqBuilder Pro (.sbd) binary files.
    SeqBuilder embeds the sequence as plain ASCII in the binary blob.
    Extract the longest continuous run of ATCG characters.
    """
    try:
        text = raw.decode("latin-1")
        matches = re.findall(r'[ATCGatcg]{100,}', text)
        if not matches:
            st.warning(
                f"⚠️ **{display_name}.sbd**: No DNA sequence found. "
                "The file may be corrupted or use an unsupported SeqBuilder version. "
                "Try exporting as FASTA or GenBank from SeqBuilder Pro.")
            return None
        seq = max(matches, key=len).upper()
        if len(matches) > 1:
            total_found = sum(len(m) for m in matches)
            st.warning(
                f"⚠️ **{display_name}.sbd**: Found {len(matches)} sequence regions "
                f"({total_found:,} bp total). Using longest ({len(seq):,} bp). "
                "If this seems wrong, export as FASTA from SeqBuilder.")
        return seq
    except Exception as e:
        st.warning(f"⚠️ **{display_name}.sbd**: Parse error: {e}")
        return None

def load_sequence(uploaded_file):
    if uploaded_file is None:
        seq = "".join(random.choices(list("ATCG"), k=10000))
        return seq, True, "Demo Plasmid"

    raw = uploaded_file.read()
    name = uploaded_file.name.lower()
    display_name = uploaded_file.name.rsplit(".", 1)[0]

    # File size check (5 MB max)
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

# ANALYSIS FUNCTIONS
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
    """Score how evenly bands are distributed across the gel (log space)."""
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

# GEL VISUALIZATION
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
        if is_thick:
            intensity = min(1.0, intensity * 1.35)
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

    # Global max fragment for consistent band intensity across all lanes
    global_max_frag = max(f for r in results for f in r["fragments"])

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

        for frag in result["fragments"]:
            y = bp_to_y(frag)
            intensity = 0.2 + 0.75 * (frag / global_max_frag)
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
                    f"<b>{top_label} | {result['enzymes']}</b><br>"
                    f"Fragment: {frag} bp<br>"
                    f"Total bands: {result['n']}"
                    "<extra></extra>"),
                showlegend=False))

    n_lanes = len(results) + 1
    fig.update_layout(
        paper_bgcolor="#1a1a2e",
        plot_bgcolor="#1a1a2e",
        title=dict(
            text=f"Predicted Restriction Digest Pattern  |  {plasmid_size} bp{title_suffix}",
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

# PLASMID MAP
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
        title=dict(text=f"Restriction Site Map  |  {plasmid_name}",
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
    # Raw sequence: strip whitespace, numbers, non-DNA chars
    seq = re.sub(r'[^ATCGatcg]', '', text).upper()
    if len(seq) > 100:
        return seq, name
    return None, name

# OVERVIEW
if tool is None:
    st.markdown("""
    <div style="text-align:center; padding: 2rem 0 1rem 0;">
        <div style="font-size:0.8rem;font-weight:700;color:#a78bfa;letter-spacing:0.14em;text-transform:uppercase;">
            🧬 Plasmid Analysis Toolkit
        </div>
        <div style="font-size:2.2rem;font-weight:800;color:#ffffff;margin-top:0.4rem;">
            Select an Analysis Module
        </div>
        <div style="font-size:1.05rem;color:#888888;margin-top:0.5rem;">
            Choose a module from the sidebar to get started.
        </div>
    </div>
    """, unsafe_allow_html=True)

    col1, col2 = st.columns(2, gap="large")

    with col1:
        st.markdown("""
        <div style="
            background: linear-gradient(135deg, rgba(124,58,237,0.15), rgba(124,58,237,0.05));
            border: 1.5px solid #7c3aed;
            border-radius: 14px;
            padding: 2rem;
            height: 100%;
        ">
            <div style="font-size:2.5rem;margin-bottom:0.75rem;">🧪</div>
            <div style="font-size:1.3rem;font-weight:800;color:#ffffff;margin-bottom:0.75rem;">
                Restriction Digest Planner
            </div>
            <div style="font-size:0.9rem;color:#c4b5fd;line-height:1.6;">
                Upload or paste a single plasmid sequence and automatically identify 
                the optimal enzyme combinations for a diagnostic digest.<br><br>
                Results are ranked by predicted gel separation quality and 
                visualised as a simulated agarose gel.
            </div>
            <div style="margin-top:1.5rem;">
                <div style="font-size:0.75rem;font-weight:700;color:#a78bfa;text-transform:uppercase;letter-spacing:0.08em;margin-bottom:0.5rem;">Accepts</div>
                <div style="font-size:0.85rem;color:#888;">FASTA · GenBank · SeqBuilder (.sbd) · Pasted sequence</div>
            </div>
        </div>
        """, unsafe_allow_html=True)
        st.markdown("<div style='height:0.75rem'></div>", unsafe_allow_html=True)
        if st.button("Open Restriction Digest Planner →", use_container_width=True, type="primary", key="ov_btn1"):
            st.session_state.active_tool = "Restriction Digest Planner"
            st.rerun()

    with col2:
        st.markdown("""
        <div style="
            background: linear-gradient(135deg, rgba(16,185,129,0.15), rgba(16,185,129,0.05));
            border: 1.5px solid #10b981;
            border-radius: 14px;
            padding: 2rem;
            height: 100%;
        ">
            <div style="font-size:2.5rem;margin-bottom:0.75rem;">🔀</div>
            <div style="font-size:1.3rem;font-weight:800;color:#ffffff;margin-bottom:0.75rem;">
                Multi-Plasmid Comparator
            </div>
            <div style="font-size:0.9rem;color:#6ee7b7;line-height:1.6;">
                Upload or paste two or more plasmid sequences and identify 
                enzyme combinations that produce distinct, discriminating 
                band patterns.<br><br>
                Ideal for colony screening and verification of multiple 
                constructs in parallel.
            </div>
            <div style="margin-top:1.5rem;">
                <div style="font-size:0.75rem;font-weight:700;color:#10b981;text-transform:uppercase;letter-spacing:0.08em;margin-bottom:0.5rem;">Accepts</div>
                <div style="font-size:0.85rem;color:#888;">FASTA · GenBank · SeqBuilder (.sbd) · Pasted sequences A–H</div>
            </div>
        </div>
        """, unsafe_allow_html=True)
        st.markdown("<div style='height:0.75rem'></div>", unsafe_allow_html=True)
        if st.button("Open Multi-Plasmid Comparator →", use_container_width=True, type="primary", key="ov_btn2"):
            st.session_state.active_tool = "Multi-Plasmid Comparator"
            st.rerun()

    st.markdown("<div style='height:1.5rem'></div>", unsafe_allow_html=True)
    col3, col4 = st.columns(2, gap="large")

    with col3:
        st.markdown("""
        <div style="
            background: linear-gradient(135deg, rgba(251,191,36,0.15), rgba(251,191,36,0.05));
            border: 1.5px solid #fbbf24;
            border-radius: 14px;
            padding: 2rem;
            height: 100%;
        ">
            <div style="font-size:2.5rem;margin-bottom:0.75rem;">🔬</div>
            <div style="font-size:1.3rem;font-weight:800;color:#ffffff;margin-bottom:0.75rem;">
                Primer Analyzer
            </div>
            <div style="font-size:0.9rem;color:#fde68a;line-height:1.6;">
                Paste one or more primer sequences to calculate melting temperature,
                GC content, and screen for potential hairpin structures and
                primer-dimer formation.<br><br>
                Supports both Nearest-Neighbor (SantaLucia 1998) and
                Wallace Rule Tm calculations.
            </div>
            <div style="margin-top:1.5rem;">
                <div style="font-size:0.75rem;font-weight:700;color:#fbbf24;text-transform:uppercase;letter-spacing:0.08em;margin-bottom:0.5rem;">Input</div>
                <div style="font-size:0.85rem;color:#888;">Paste raw primer sequence(s) directly</div>
            </div>
        </div>
        """, unsafe_allow_html=True)
        st.markdown("<div style='height:0.75rem'></div>", unsafe_allow_html=True)
        if st.button("Open Primer Analyzer →", use_container_width=True, type="primary", key="ov_btn3"):
            st.session_state.active_tool = "Primer Analyzer"
            st.rerun()

    with col4:
        st.markdown("""
        <div style="
            background: linear-gradient(135deg, rgba(239,68,68,0.15), rgba(239,68,68,0.05));
            border: 1.5px solid #ef4444;
            border-radius: 14px;
            padding: 2rem;
            height: 100%;
        ">
            <div style="font-size:2.5rem;margin-bottom:0.75rem;">🗺️</div>
            <div style="font-size:1.3rem;font-weight:800;color:#ffffff;margin-bottom:0.75rem;">
                Feature Annotation Viewer
            </div>
            <div style="font-size:0.9rem;color:#fca5a5;line-height:1.6;">
                Upload a GenBank file to visualise all annotated features —
                genes, promoters, terminators, primer sites, and restriction
                sites — on an interactive circular plasmid map.<br><br>
                Click any feature for details. Color-coded by feature type.
            </div>
            <div style="margin-top:1.5rem;">
                <div style="font-size:0.75rem;font-weight:700;color:#ef4444;text-transform:uppercase;letter-spacing:0.08em;margin-bottom:0.5rem;">Accepts</div>
                <div style="font-size:0.85rem;color:#888;">GenBank (.gb / .gbk) with annotations</div>
            </div>
        </div>
        """, unsafe_allow_html=True)
        st.markdown("<div style='height:0.75rem'></div>", unsafe_allow_html=True)
        if st.button("Open Feature Annotation Viewer →", use_container_width=True, type="primary", key="ov_btn4"):
            st.session_state.active_tool = "Feature Annotation Viewer"
            st.rerun()

# TOOL 1: RESTRICTION DIGEST PLANNER
elif tool == "Restriction Digest Planner":
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
            "Or paste sequence directly",
            placeholder="Paste raw DNA sequence or FASTA here…",
            height=80,
            key="paste_1",
            help="Accepts raw sequence (ATCG) or FASTA format. Takes priority over file upload if both are provided."
        )

        st.caption("💡 No input required. Without a sequence the tool runs in demo mode with a randomly generated 10 kb plasmid.")

        run = st.button("▶  Run Analysis", type="primary", use_container_width=True, key="run_1")
        prefer_short = st.checkbox(
            "Prioritise short fragments",
            value=False,
            help="Ranks results by smallest largest fragment, minimising total gel run time.")

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

        st.markdown("**Host strain methylation**")
        dam_col, dcm_col = st.columns(2)
        with dam_col:
            dam_plus = st.checkbox("dam+", value=True, key="dam_1",
                                   help="dam methylase blocks some enzymes at GATC sites")
        with dcm_col:
            dcm_plus = st.checkbox("dcm+", value=True, key="dcm_1",
                                   help="dcm methylase blocks some enzymes at CCWGG sites")

        select_all = st.checkbox("Select all enzymes", value=True, key="sel_all_1")
        if select_all:
            selected_enzymes = DEFAULT_ENZYMES
        else:
            selected_enzymes = st.multiselect(
                "Select enzymes available in your laboratory:",
                DEFAULT_ENZYMES, default=DEFAULT_ENZYMES[:10], key="enz_1")

        # Filter methylation-sensitive enzymes
        if dam_plus:
            selected_enzymes = [e for e in selected_enzymes if e not in DAM_BLOCKED]
        if dcm_plus:
            selected_enzymes = [e for e in selected_enzymes if e not in DCM_BLOCKED]

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
            st.warning("⚠️ No sequence uploaded. Running in demo mode with a randomly generated 10 kb plasmid.")
        else:
            st.success(f"✅ {plasmid_name} loaded ({plasmid_size} bp)")

        for min_f in [min_frags, min_frags + 1]:
            st.divider()
            st.subheader(f"Results: minimum {min_f} bands")

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
                               title_suffix=f" · {plasmid_name} · min {min_f} bands")
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

# TOOL 2: MULTI-PLASMID COMPARATOR
elif tool == "Multi-Plasmid Comparator":
    st.markdown("### 🔀 Multi-Plasmid Comparator")
    st.markdown("Compare restriction digest patterns across multiple plasmids to identify discriminating enzyme combinations.")

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

        st.markdown("**Host strain methylation**")
        dam_col2, dcm_col2 = st.columns(2)
        with dam_col2:
            dam_plus2 = st.checkbox("dam+", value=True, key="dam_2",
                                    help="dam methylase blocks some enzymes at GATC sites")
        with dcm_col2:
            dcm_plus2 = st.checkbox("dcm+", value=True, key="dcm_2",
                                    help="dcm methylase blocks some enzymes at CCWGG sites")

        select_all2 = st.checkbox("Select all enzymes", value=True, key="sel_all_2")
        if select_all2:
            selected_enzymes2 = DEFAULT_ENZYMES
        else:
            selected_enzymes2 = st.multiselect(
                "Select enzymes:", DEFAULT_ENZYMES,
                default=DEFAULT_ENZYMES[:10], key="enz_2")

        # Filter methylation-sensitive enzymes
        if dam_plus2:
            selected_enzymes2 = [e for e in selected_enzymes2 if e not in DAM_BLOCKED]
        if dcm_plus2:
            selected_enzymes2 = [e for e in selected_enzymes2 if e not in DCM_BLOCKED]

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
        for raw_paste, default_name in pasted_seqs_2:
            if raw_paste and raw_paste.strip():
                seq, name = parse_pasted_sequence(raw_paste, default_name)
                if seq:
                    plasmids.append({"seq": seq, "name": name, "size": len(seq)})
                else:
                    st.warning(f"⚠️ Could not parse pasted {default_name}. Skipped.")

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

                # Discrimination score: log-space optimal fragment matching
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
                    f"Discriminating enzymes found: **{', '.join(sorted(discriminating))}**. "
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

            st.subheader("🧫 Predicted Gel: All Plasmids per Combination")
            max_size = max(p["size"] for p in plasmids)

            for i, combo in enumerate(top_combos[:5]):
                st.markdown(f"**#{i+1}  {combo['enzymes']}**  (score: {combo['diff_score']:.3f})")
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
                               title_suffix=f" · {combo['enzymes']}",
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

# TOOL 3: PRIMER ANALYZER
elif tool == "Primer Analyzer":
    from Bio.SeqUtils.MeltingTemp import Tm_NN, Tm_Wallace
    from Bio.SeqUtils import gc_fraction

    st.markdown("### 🔬 Primer Analyzer")
    st.markdown("Calculate melting temperature, GC content, and screen for hairpin and primer-dimer formation.")

    with st.sidebar:
        st.header("⚙️ Parameters")
        st.subheader("Reaction Conditions")
        conc_primer = st.number_input("Primer concentration (nM)", value=250, min_value=1, max_value=10000, step=50)
        conc_na = st.number_input("Na⁺ concentration (mM)", value=50, min_value=1, max_value=1000, step=10)
        st.divider()
        st.subheader("Hairpin / Dimer Check")
        min_stem = st.slider("Minimum stem length (bp)", 3, 8, 4,
                             help="Minimum complementary run to flag as potential hairpin or dimer")

    st.markdown("#### Paste primer sequences")
    st.caption("One primer per line. Optionally prefix with a name: `FwdPrimer  ATCGATCG...`")
    raw_input = st.text_area("Primer sequences", height=160,
                              placeholder="Fwd_primer   ATCGATCGATCGATCG\nRev_primer   TAGCTAGCTAGCTAGC\nor just paste sequences without names")

    def parse_primers(text):
        primers = []
        for i, line in enumerate(text.strip().splitlines()):
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            # Check if first token looks like a name (no pure ATCG)
            if len(parts) >= 2 and re.search(r'[^ATCGatcg]', parts[0]):
                name = parts[0]
                seq = re.sub(r'[^ATCGatcg]', '', ''.join(parts[1:])).upper()
            else:
                seq = re.sub(r'[^ATCGatcg]', '', ''.join(parts)).upper()
                name = f"Primer {i+1}"
            if len(seq) >= 10:
                primers.append((name, seq))
        return primers

    def check_hairpin(seq, min_stem):
        """Check for potential hairpin: find complementary runs within the sequence."""
        comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        rev_comp = ''.join(comp.get(b, 'N') for b in reversed(seq))
        hits = []
        for i in range(len(seq) - min_stem * 2 - 3):
            for length in range(min_stem, min((len(seq) - i) // 2, 10)):
                stem = seq[i:i+length]
                rc_stem = ''.join(comp.get(b, 'N') for b in reversed(stem))
                if rc_stem in seq[i+length+3:]:
                    hits.append(f"stem: {stem} (pos {i+1}, {length} bp)")
                    break
        return hits[:3]  # return max 3 hits

    def check_dimer(seq1, seq2, min_stem):
        """Check for complementarity between two primers."""
        comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        rc2 = ''.join(comp.get(b, 'N') for b in reversed(seq2))
        hits = []
        for i in range(len(seq1) - min_stem + 1):
            for length in range(min_stem, min(len(seq1) - i, 12)):
                chunk = seq1[i:i+length]
                if chunk in rc2:
                    hits.append(f"{chunk} ({length} bp)")
                    break
        return hits[:3]

    if raw_input and raw_input.strip():
        primers = parse_primers(raw_input)
        if not primers:
            st.error("No valid primer sequences found. Sequences must be at least 10 bp.")
        else:
            st.markdown("---")
            st.subheader("Results")

            rows = []
            for name, seq in primers:
                seq_obj = Seq(seq)
                try:
                    tm_nn = Tm_NN(seq_obj,
                                  Na=conc_na,
                                  dnac1=conc_primer,
                                  dnac2=0)
                except Exception:
                    tm_nn = None
                try:
                    tm_w = Tm_Wallace(seq_obj)
                except Exception:
                    tm_w = None

                gc = gc_fraction(seq_obj) * 100
                hairpins = check_hairpin(seq, min_stem)

                rows.append({
                    "Name": name,
                    "Sequence (5'→3')": seq,
                    "Length (bp)": len(seq),
                    "Tm NN (°C)": f"{tm_nn:.1f}" if tm_nn else "n/a",
                    "Tm Wallace (°C)": f"{tm_w:.1f}" if tm_w else "n/a",
                    "GC (%)": f"{gc:.1f}",
                    "Hairpin risk": "⚠️ Possible" if hairpins else "✅ Low",
                })

            df = pd.DataFrame(rows)
            st.dataframe(df, use_container_width=True, hide_index=True)

            # Tm comparison chart
            if len(primers) > 1:
                st.subheader("Tm Comparison")
                fig = go.Figure()
                names = [r["Name"] for r in rows]
                tm_nn_vals = [float(r["Tm NN (°C)"]) if r["Tm NN (°C)"] != "n/a" else None for r in rows]
                tm_w_vals  = [float(r["Tm Wallace (°C)"]) if r["Tm Wallace (°C)"] != "n/a" else None for r in rows]

                fig.add_trace(go.Bar(name="Nearest-Neighbor", x=names, y=tm_nn_vals,
                                     marker_color="#7c3aed"))
                fig.add_trace(go.Bar(name="Wallace Rule", x=names, y=tm_w_vals,
                                     marker_color="#10b981"))
                fig.update_layout(
                    paper_bgcolor="#1a1a2e", plot_bgcolor="#1a1a2e",
                    font=dict(color="white"),
                    barmode="group",
                    yaxis=dict(title="Tm (°C)", gridcolor="#333"),
                    xaxis=dict(gridcolor="#333"),
                    legend=dict(bgcolor="#2a2a4a"),
                    height=350, margin=dict(t=30, b=20, l=40, r=20))
                st.plotly_chart(fig, use_container_width=True)

            # Dimer check
            if len(primers) >= 2:
                st.subheader("Primer-Dimer Screening")
                dimer_rows = []
                for i in range(len(primers)):
                    for j in range(i, len(primers)):
                        hits = check_dimer(primers[i][1], primers[j][1], min_stem)
                        label = "Self" if i == j else f"{primers[i][0]} + {primers[j][0]}"
                        dimer_rows.append({
                            "Pair": label,
                            "Risk": "⚠️ Possible" if hits else "✅ Low",
                            "Complementary regions": ", ".join(hits) if hits else "None detected"
                        })
                st.dataframe(pd.DataFrame(dimer_rows), use_container_width=True, hide_index=True)

            # Hairpin detail
            with st.expander("Hairpin detail"):
                for name, seq in primers:
                    hits = check_hairpin(seq, min_stem)
                    if hits:
                        st.markdown(f"**{name}**: {'; '.join(hits)}")
                    else:
                        st.markdown(f"**{name}**: No hairpin detected")
    else:
        st.info("👈 Paste one or more primer sequences to begin analysis.")


# TOOL 4: FEATURE ANNOTATION VIEWER
elif tool == "Feature Annotation Viewer":
    st.markdown("### 🗺️ Feature Annotation Viewer")
    st.markdown("Visualise annotated features on an interactive circular plasmid map from a GenBank file.")

    # Wong (2011) colorblind-safe palette — optimised for red-green colour vision deficiency
    FEATURE_COLORS = {
        "CDS":              "#0072B2",  # blue
        "gene":             "#777788",  # neutral grey — gene features are often just containers
        "promoter":         "#F0E442",  # yellow
        "terminator":       "#56B4E9",  # sky blue
        "rep_origin":       "#009E73",  # bluish green
        "primer_bind_fwd":  "#FFE500",  # bright yellow — forward primers
        "primer_bind_rev":  "#FF2D2D",  # bright red — reverse primers
        "primer_bind":      "#FFE500",  # fallback (no strand info)
        "misc_feature":     "#D55E00",  # vermillion
        "regulatory":       "#0096FF",  # bright blue
        "LTR":              "#F5C542",  # golden yellow
        "RBS":              "#56B4E9",  # sky blue
        "sig_peptide":      "#CC79A7",  # reddish purple
        "mat_peptide":      "#009E73",  # bluish green
    }
    DEFAULT_COLOR = "#aaaaaa"

    def get_label(feat):
        return (feat.qualifiers.get("product", [""])[0] or
                feat.qualifiers.get("gene",    [""])[0] or
                feat.qualifiers.get("label",   [""])[0] or
                feat.qualifiers.get("locus_tag",[""])[0] or
                feat.type)

    PRIMER_KEYWORDS = re.compile(
        r"primer|oligo|probe|oJW|oSB|oHH|oAB|oMB|oCW|^o[A-Z]{2}[0-9]|"
        r"fwd|rev|forward|reverse|^F_|^R_|_F$|_R$|_fwd|_rev",
        re.IGNORECASE)

    def primer_type_key(feat):
        """Return strand-aware primer type key, or None if not a primer."""
        is_p = False
        if feat.type == "primer_bind":
            is_p = True
        elif feat.type in ("misc_feature", "misc_binding"):
            if PRIMER_KEYWORDS.search(get_label(feat)):
                is_p = True
        if not is_p:
            return None
        strand = feat.location.strand
        if strand == 1:
            return "primer_bind_fwd"
        elif strand == -1:
            return "primer_bind_rev"
        return "primer_bind"  # unknown strand

    def is_primer(feat):
        return primer_type_key(feat) is not None

    def draw_annotation_map(record, show_labels=True, zoom_start=None, zoom_end=None):
        plasmid_size = len(record.seq)
        plasmid_name = record.id or record.name or "Plasmid"
        fig = go.Figure()

        def pos_to_angle(pos):
            return np.pi / 2 - 2 * np.pi * (pos - 1) / plasmid_size

        # Backbone ring
        theta = np.linspace(np.pi / 2, np.pi / 2 - 2 * np.pi, 600)
        fig.add_trace(go.Scatter(
            x=np.cos(theta), y=np.sin(theta),
            mode="lines", line=dict(color="#333355", width=18),
            hoverinfo="skip", showlegend=False))
        # Thin centre line of backbone
        fig.add_trace(go.Scatter(
            x=0.91 * np.cos(theta), y=0.91 * np.sin(theta),
            mode="lines", line=dict(color="#555577", width=1),
            hoverinfo="skip", showlegend=False))

        # Strand direction tick marks every ~5% of plasmid
        tick_every = max(1, plasmid_size // 20)
        for tick_pos in range(1, plasmid_size, tick_every):
            a = pos_to_angle(tick_pos)
            fig.add_shape(type="line",
                x0=0.88 * np.cos(a), y0=0.88 * np.sin(a),
                x1=0.94 * np.cos(a), y1=0.94 * np.sin(a),
                line=dict(color="#555577", width=1))

        # Features
        seen_types = {}
        occupied_tracks = []  # (start_bp, end_bp, track_index) for overlap avoidance
        for feat in record.features:
            if feat.type in ("source",):
                continue
            # Detect primers by type/name; use strand-aware key for color
            _pk = primer_type_key(feat)
            display_type = _pk if _pk else feat.type
            color = FEATURE_COLORS.get(display_type, DEFAULT_COLOR)
            start  = int(feat.location.start) + 1
            end    = int(feat.location.end)
            strand = feat.location.strand
            label  = get_label(feat)
            size_bp = end - start + 1

            a0 = pos_to_angle(start)
            a1 = pos_to_angle(end)

            if end < start:
                # Wrap-around feature (e.g. 9800..200): goes clockwise past origin
                size_bp = (plasmid_size - start) + end
                if a1 > a0:
                    a1 -= 2 * np.pi
            elif end == start:
                # Point feature (SNP, single base): force minimal visible arc
                a1 = a0 - np.radians(3)
            else:
                # Normal feature: a1 should be < a0 (clockwise)
                if a1 >= a0:
                    a1 -= 2 * np.pi

            # Hard cap: never draw more than one full circle
            if abs(a1 - a0) >= 2 * np.pi:
                a1 = a0 - (2 * np.pi - 0.01)

            # Minimum visible arc: 3° for very short features
            min_arc = np.radians(3)
            if abs(a1 - a0) < min_arc:
                mid_a_ext = (a0 + a1) / 2
                a0 = mid_a_ext + min_arc / 2
                a1 = mid_a_ext - min_arc / 2

            n_pts = max(40, int(abs(a1 - a0) / (2 * np.pi) * 500) + 2)
            arc = np.linspace(a0, a1, n_pts)

            # Fixed ring width — arc length is already proportional to bp position
            ring_w = 0.10

            # Two concentric rings: fwd (+) outer, rev (-) inner
            # Use track assignment to avoid overlapping features
            if strand == -1:
                base_r = 0.82
            else:
                base_r = 1.00

            # Find the first free track for this feature (greedy interval scheduling)
            def overlaps(s1, e1, s2, e2):
                # Handle wrap-around (circular)
                if s1 <= e1 and s2 <= e2:
                    return not (e1 < s2 or e2 < s1)
                return True  # conservative: assume overlap if wrap-around

            track = 0
            while True:
                blocked = any(
                    overlaps(start, end, us, ue) and ut == track
                    for us, ue, ut in occupied_tracks
                )
                if not blocked:
                    break
                track += 1

            occupied_tracks.append((start, end, track))

            track_offset = track * (ring_w + 0.02)
            if strand == -1:
                r_outer = base_r - track_offset
                r_inner = r_outer - ring_w
            else:
                r_outer = base_r + track_offset
                r_inner = r_outer - ring_w

            x_out = r_outer * np.cos(arc)
            y_out = r_outer * np.sin(arc)
            x_in  = r_inner * np.cos(arc[::-1])
            y_in  = r_inner * np.sin(arc[::-1])

            # Arrow head — width proportional to feature arc, capped at 15% of arc
            arc_span = abs(a1 - a0)  # total arc in radians
            # Arrow takes up at most 20% of the feature arc, min ~2°, max ~8°
            arrow_span = min(arc_span * 0.20, np.radians(8))
            arrow_span = max(arrow_span, np.radians(1.5))

            r_mid = (r_outer + r_inner) / 2
            tip_a  = a0 if strand != -1 else a1    # leading end (clockwise = decreasing angle)
            # For fwd strand tip is at a0 (start of arc), arrow points clockwise
            # For rev strand tip is at a1 (end of arc), arrow points counter-clockwise
            side_a = tip_a + (arrow_span if strand != -1 else -arrow_span)
            arrow_x = [r_outer * np.cos(tip_a),
                        r_mid   * np.cos(side_a),
                        r_inner * np.cos(tip_a),
                        r_outer * np.cos(tip_a)]
            arrow_y = [r_outer * np.sin(tip_a),
                        r_mid   * np.sin(side_a),
                        r_inner * np.sin(tip_a),
                        r_outer * np.sin(tip_a)]

            hover = (f"<b>{label}</b><br>"
                     f"Type: {display_type}<br>"
                     f"Position: {start:,} – {end:,} bp<br>"
                     f"Strand: {'(+) forward' if strand == 1 else '(-) reverse'}<br>"
                     f"Size: {size_bp:,} bp")

            # Feature arc body
            fig.add_trace(go.Scatter(
                x=np.concatenate([x_out, x_in]),
                y=np.concatenate([y_out, y_in]),
                fill="toself", fillcolor=color,
                line=dict(width=0, color=color), mode="lines",
                text=hover, hoverinfo="text", hoveron="fills",
                showlegend=display_type not in seen_types,
                name=display_type,
                legendgroup=display_type,
                opacity=0.88))
            seen_types[display_type] = True

            # Arrow head
            fig.add_trace(go.Scatter(
                x=arrow_x, y=arrow_y,
                fill="toself", fillcolor=color,
                line=dict(width=0), mode="lines",
                hoverinfo="skip", showlegend=False,
                legendgroup=feat.type, opacity=1.0))

            # Label — only if show_labels=True and feature is large enough to warrant one
            if show_labels:
                arc_fraction = abs(a1 - a0) / (2 * np.pi)
                if arc_fraction > 0.025 or size_bp > plasmid_size * 0.025:
                    mid_a = (a0 + a1) / 2
                    if strand != -1:
                        r_lbl = r_outer + 0.13
                    else:
                        r_lbl = r_inner - 0.10
                    # Scale font size and label length with feature size
                    lbl_size  = max(6, min(10, int(6 + arc_fraction * 40)))
                    lbl_chars = max(6, min(20, int(arc_fraction * 120)))
                    fig.add_annotation(
                        x=r_lbl * np.cos(mid_a),
                        y=r_lbl * np.sin(mid_a),
                        text=label[:lbl_chars],
                        showarrow=False,
                        font=dict(color=color, size=lbl_size, family="Arial"),
                        xanchor="center", yanchor="middle",
                        bgcolor="rgba(26,26,46,0.7)",
                        borderpad=1)

        # bp scale ticks (four cardinal points)
        for frac, pos_label in [(0, "1"), (0.25, f"{plasmid_size//4:,}"),
                                 (0.5, f"{plasmid_size//2:,}"), (0.75, f"{3*plasmid_size//4:,}")]:
            a = np.pi / 2 - 2 * np.pi * frac
            fig.add_annotation(
                x=1.28 * np.cos(a), y=1.28 * np.sin(a),
                text=f"<span style='font-size:9px;color:#666'>{pos_label}</span>",
                showarrow=False, font=dict(color="#777", size=9),
                xanchor="center", yanchor="middle")
            fig.add_shape(type="line",
                x0=1.01 * np.cos(a), y0=1.01 * np.sin(a),
                x1=1.06 * np.cos(a), y1=1.06 * np.sin(a),
                line=dict(color="#555577", width=1.5))

        # Center labels
        fig.add_annotation(x=0, y=0.10,
                           text=f"<b>{plasmid_name}</b>",
                           showarrow=False,
                           font=dict(color="white", size=13, family="Arial Black"),
                           xanchor="center")
        fig.add_annotation(x=0, y=-0.12,
                           text=f"{plasmid_size:,} bp",
                           showarrow=False,
                           font=dict(color="#aaaaaa", size=11),
                           xanchor="center")

        fig.update_layout(
            paper_bgcolor="#1a1a2e", plot_bgcolor="#1a1a2e",
            title=dict(text=f"Feature Map  |  {plasmid_name}",
                       font=dict(color="white", size=14, family="Arial Black"), x=0.5),
            xaxis=dict(showticklabels=False, showgrid=False, zeroline=False, range=[-1.85, 1.85]),
            yaxis=dict(showticklabels=False, showgrid=False, zeroline=False,
                       range=[-1.85, 1.85], scaleanchor="x"),
            height=700, showlegend=True,
            legend=dict(
                font=dict(color="white", size=10),
                bgcolor="#2a2a4a", bordercolor="#555", borderwidth=1,
                title=dict(text="Feature type  (click to toggle)", font=dict(color="#aaa", size=9)),
                itemclick="toggle", itemdoubleclick="toggleothers"),
            margin=dict(t=50, b=10, l=10, r=10),
            hovermode="closest",
            dragmode="pan")
        return fig, plasmid_size

    with st.sidebar:
        st.header("⚙️ Parameters")
        seq_file = st.file_uploader(
            "Upload sequence file",
            type=["gb", "gbk", "genbank", "fasta", "fa", "fas", "sbd"],
            help="GenBank with annotations (.gb/.gbk) or FASTA/SeqBuilder for sequence-only view",
            key="uploader_feat")

        if seq_file:
            st.divider()
            st.subheader("Feature Filter")
            all_types_placeholder = st.empty()

    # Alias for downstream code
    gb_file = seq_file if "seq_file" in dir() else None

    if seq_file:
        raw = seq_file.read()
        fname = seq_file.name.lower()
        record = None

        # Try GenBank first (has annotations)
        for enc in ("utf-8", "latin-1"):
            try:
                record = SeqIO.read(io.StringIO(raw.decode(enc)), "genbank")
                break
            except Exception:
                pass

        # Fall back to FASTA
        if record is None:
            for enc in ("utf-8", "latin-1"):
                try:
                    record = SeqIO.read(io.StringIO(raw.decode(enc)), "fasta")
                    break
                except Exception:
                    pass

        # Fall back to .sbd
        if record is None and fname.endswith(".sbd"):
            seq_str = read_sbd(raw, seq_file.name.rsplit(".", 1)[0])
            if seq_str:
                from Bio.SeqRecord import SeqRecord as _SR
                record = _SR(Seq(seq_str),
                             id=seq_file.name.rsplit(".", 1)[0],
                             name=seq_file.name.rsplit(".", 1)[0],
                             description="")

        if record is None:
            st.error("Could not parse the uploaded file. Please use GenBank, FASTA, or SeqBuilder format.")
            st.stop()

        if not record.features:
            st.info("ℹ️ No feature annotations found in this file. The map will show the plasmid backbone only. "
                    "For full annotation support, use a GenBank file exported from SnapGene, Benchling, or ApE.")

        plasmid_size = len(record.seq)
        plasmid_name = record.id or record.name or seq_file.name
        st.success(f"✅ {plasmid_name} loaded ({plasmid_size:,} bp)")

        feat_types = sorted(set(f.type for f in record.features if f.type != "source"))
        feat_counts = {t: sum(1 for f in record.features if f.type == t) for t in feat_types}

        with st.sidebar:
            selected_types = all_types_placeholder.multiselect(
                "Show feature types",
                options=feat_types,
                default=feat_types,
                format_func=lambda t: f"{t} ({feat_counts[t]})")
            st.divider()
            show_labels = st.checkbox("Show labels on map", value=True,
                help="Turn off for cleaner map when many features are present. Hover still works.")
            st.divider()
            st.markdown("**Linear zoom view**")
            st.caption("Inspect a specific region in detail.")
            zoom_enabled = st.checkbox("Enable zoom panel", value=False)
            if zoom_enabled:
                zoom_start = st.number_input("Start (bp)", min_value=1,
                    max_value=plasmid_size, value=1, step=100)
                zoom_end = st.number_input("End (bp)", min_value=1,
                    max_value=plasmid_size, value=min(plasmid_size, 2000), step=100)
            else:
                zoom_start, zoom_end = None, None

        # Filter record features
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature
        filtered = SeqRecord(record.seq, id=record.id, name=record.name,
                              description=record.description)
        filtered.features = [f for f in record.features
                              if f.type == "source" or f.type in selected_types]

        # Build feature table rows
        feat_rows = []
        for f in filtered.features:
            if f.type == "source":
                continue
            feat_rows.append({
                "Type": (primer_type_key(f) or f.type),
                "Label": get_label(f),
                "Start": int(f.location.start) + 1,
                "End": int(f.location.end),
                "Strand": "+" if f.location.strand == 1 else "-",
                "Size (bp)": int(f.location.end) - int(f.location.start),
            })
        feat_df = pd.DataFrame(feat_rows).sort_values("Start") if feat_rows else pd.DataFrame()

        col_map, col_table = st.columns([1.15, 1], gap="large")

        with col_map:
            with st.spinner("Rendering map..."):
                fig_annot, _ = draw_annotation_map(filtered, show_labels=show_labels)
                st.plotly_chart(fig_annot, use_container_width=True, key="feat_map")

            # Linear zoom panel
            if zoom_enabled and zoom_start and zoom_end and zoom_start < zoom_end:
                st.markdown(f"**Region {zoom_start:,} – {zoom_end:,} bp**")
                region_feats = [f for f in filtered.features
                                if f.type != "source"
                                and int(f.location.start) < zoom_end
                                and int(f.location.end) > zoom_start]
                if region_feats:
                    fig_linear = go.Figure()
                    fig_linear.add_shape(type="line",
                        x0=zoom_start, x1=zoom_end, y0=0.5, y1=0.5,
                        line=dict(color="#444466", width=6))
                    for fi, feat in enumerate(region_feats):
                        color = FEATURE_COLORS.get(feat.type, DEFAULT_COLOR)
                        fs = max(int(feat.location.start) + 1, zoom_start)
                        fe = min(int(feat.location.end), zoom_end)
                        strand = feat.location.strand
                        y_pos = 0.7 if strand != -1 else 0.3
                        lbl = get_label(feat)
                        fig_linear.add_shape(type="rect",
                            x0=fs, x1=fe, y0=y_pos - 0.12, y1=y_pos + 0.12,
                            fillcolor=color, opacity=0.85,
                            line=dict(width=0))
                        fig_linear.add_trace(go.Scatter(
                            x=[(fs + fe) / 2], y=[y_pos],
                            mode="markers", marker=dict(size=1, opacity=0),
                            text=f"<b>{lbl}</b><br>{feat.type}<br>{int(feat.location.start)+1:,}–{int(feat.location.end):,} bp",
                            hoverinfo="text", showlegend=False))
                        if fe - fs > (zoom_end - zoom_start) * 0.04:
                            fig_linear.add_annotation(
                                x=(fs + fe) / 2, y=y_pos,
                                text=lbl[:16], showarrow=False,
                                font=dict(color="white", size=9),
                                xanchor="center", yanchor="middle")
                    fig_linear.add_annotation(x=(zoom_start + zoom_end)/2, y=0.08,
                        text="← reverse strand", showarrow=False,
                        font=dict(color="#666", size=9))
                    fig_linear.add_annotation(x=(zoom_start + zoom_end)/2, y=0.92,
                        text="forward strand →", showarrow=False,
                        font=dict(color="#666", size=9))
                    fig_linear.update_layout(
                        paper_bgcolor="#1a1a2e", plot_bgcolor="#1a1a2e",
                        xaxis=dict(showgrid=False, color="#666",
                                   range=[zoom_start, zoom_end],
                                   tickformat=","),
                        yaxis=dict(showgrid=False, showticklabels=False,
                                   range=[0, 1]),
                        height=200, margin=dict(t=10, b=30, l=10, r=10),
                        hovermode="closest")
                    st.plotly_chart(fig_linear, use_container_width=True, key="linear_zoom")
                else:
                    st.caption("No features in selected region.")

        with col_table:
            n_feats = len(feat_df) if not feat_df.empty else 0
            st.subheader(f"Features ({n_feats})")

            # Search box
            search = st.text_input("Search features", placeholder="Gene name, type, position...",
                                   key="feat_search")
            if search and not feat_df.empty:
                mask = feat_df.apply(
                    lambda row: search.lower() in str(row).lower(), axis=1)
                feat_df = feat_df[mask]
                st.caption(f'{len(feat_df)} match(es) for "{search}"')

            if not feat_df.empty:
                def color_type(val):
                    c = FEATURE_COLORS.get(val, DEFAULT_COLOR)
                    return f"background-color: {c}22; color: {c}; font-weight: 600"

                st.dataframe(
                    feat_df.style.applymap(color_type, subset=["Type"]),
                    use_container_width=True, hide_index=True, height=540)
            else:
                st.caption("No features to display.")

            # GenBank download
            gb_out = io.StringIO()
            SeqIO.write(record, gb_out, "genbank")
            st.download_button(
                "Download GenBank",
                data=gb_out.getvalue(),
                file_name=f"{plasmid_name}.gb",
                mime="text/plain")
    else:
        st.info("👈 Upload a GenBank file (.gb or .gbk) with feature annotations.")
        st.markdown("""
        **What this module shows:**
        - All annotated features (CDS, promoters, terminators, primer sites, etc.) on a circular map
        - Color-coded by feature type, with hover details for each feature
        - Filterable feature table with start/end positions and strand orientation
        - GenBank download of the loaded file

        **Tip:** Most plasmid design tools (SnapGene, Benchling, ApE) can export GenBank files with annotations.
        """)
