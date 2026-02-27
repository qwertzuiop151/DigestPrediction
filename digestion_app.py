Skip to content
qwertzuiop151
DigestPrediction
Repository navigation
Code
Issues
Pull requests
Actions
Projects
Wiki
Security
Insights
Settings
Files
Go to file
t
digestion_app.py
requirements.txt
DigestPrediction
/
digestion_app.py
in
main

Edit

Preview
Indent mode

Spaces
Indent size

4
Line wrap mode

No wrap
Editing digestion_app.py file contents
  1
  2
  3
  4
  5
  6
  7
  8
  9
 10
 11
 12
 13
 14
 15
 16
 17
 18
 19
 20
 21
 22
 23
 24
 25
 26
 27
 28
 29
 30
 31
 32
 33
 34
 35
 36
 37
 38
 39
 40
 41
 42
 43
 44
 45
 46
 47
 48
 49
 50
 51
 52
 53
 54
 55
 56
 57
 58
 59
 60
 61
 62
 63
 64
 65
 66
 67
 68
 69
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
    ["Restriction Digest Planner", "Multi-Plasmid Comparator"],
    index=0
)

tool_descriptions = {
    "Restriction Digest Planner": "🧪 **Restriction Digest Planner** — Upload a plasmid sequence and automatically find the best enzyme combinations for a diagnostic digest. Ranked by band separation quality and visualised as a predicted agarose gel.",
    "Multi-Plasmid Comparator": "🔀 **Multi-Plasmid Comparator** — Upload 2 or more plasmids and identify which enzyme combinations produce distinct, distinguishable band patterns. Ideal for colony screening and construct verification.",
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
Use Control + Shift + m to toggle the tab key moving focus. Alternatively, use esc then tab to move to the next interactive element on the page.
