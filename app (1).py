"""MOLECULE — a Nothing OS inspired molecule explorer for Streamlit.

Single file. Run locally:   pip install streamlit rdkit && streamlit run app.py
Deploy: put this file in a GitHub repo as app.py, and add a requirements.txt
containing two lines:  streamlit  and  rdkit
(All styling is done in-code, so no .streamlit/config.toml is needed.)
"""
import random

import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

st.set_page_config(page_title="MOLECULE", page_icon="●", layout="wide")

# ----------------------------------------------------------------------------
# Example library
# ----------------------------------------------------------------------------
LIBRARY = {
    "Medicines": {
        "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "Paracetamol": "CC(=O)NC1=CC=C(O)C=C1",
        "Ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        "Naproxen": "COC1=CC2=CC(=CC=C2C=C1)C(C)C(=O)O",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "Metformin": "CN(C)C(=N)N=C(N)N",
        "Penicillin G": "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O",
        "Diazepam": "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc12",
        "Omeprazole": "COc1ccc2[nH]c(nc2c1)S(=O)Cc1ncc(C)c(OC)c1C",
        "Sildenafil": "CCCc1nn(C)c2c1nc([nH]c2=O)c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1",
        "Morphine": "CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O",
        "Quinine": "COc1ccc2nccc([C@@H](O)[C@H]3C[N@]4CC[C@H]3C[C@H]4C=C)c2c1",
    },
    "Vitamins": {
        "Vitamin C": "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O",
        "Niacin (B3)": "OC(=O)c1cccnc1",
        "Pyridoxine (B6)": "Cc1ncc(CO)c(CO)c1O",
        "Biotin (B7)": "O=C(O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12",
        "Folic acid (B9)": "Nc1nc2ncc(CNc3ccc(cc3)C(=O)N[C@@H](CCC(O)=O)C(O)=O)nc2c(=O)[nH]1",
        "Retinol (A)": "CC1=C(C(CCC1)(C)C)C=CC(C)=CC=CC(C)=CCO",
    },
    "Amino acids": {
        "Glycine": "NCC(=O)O",
        "Alanine": "C[C@H](N)C(=O)O",
        "Phenylalanine": "N[C@@H](Cc1ccccc1)C(=O)O",
        "Tyrosine": "N[C@@H](Cc1ccc(O)cc1)C(=O)O",
        "Tryptophan": "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O",
        "Glutamic acid": "N[C@@H](CCC(=O)O)C(=O)O",
    },
    "Brain & hormones": {
        "Dopamine": "NCCc1ccc(O)c(O)c1",
        "Serotonin": "NCCc1c[nH]c2ccc(O)cc12",
        "Melatonin": "COc1ccc2[nH]cc(CCNC(C)=O)c2c1",
        "Adrenaline": "CNC[C@H](O)c1ccc(O)c(O)c1",
        "GABA": "NCCCC(=O)O",
        "Acetylcholine": "CC(=O)OCC[N+](C)(C)C",
        "Testosterone": "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2O",
        "Estradiol": "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCc4cc(O)ccc34",
        "Cholesterol": "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
    },
    "Food & flavour": {
        "Glucose": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
        "Sucrose": "OC[C@H]1O[C@@](CO)(O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@@H]1O",
        "Citric acid": "OC(=O)CC(O)(CC(O)=O)C(O)=O",
        "Vanillin": "COc1cc(C=O)ccc1O",
        "Capsaicin": "COc1cc(CNC(=O)CCCCC=CC(C)C)ccc1O",
        "Menthol": "CC(C)[C@@H]1CC[C@@H](C)C[C@H]1O",
        "Limonene": "CC1=CCC(CC1)C(=C)C",
        "Nicotine": "CN1CCC[C@H]1c1cccnc1",
        "Resveratrol": "Oc1ccc(C=Cc2cc(O)cc(O)c2)cc1",
        "Curcumin": "COc1cc(C=CC(=O)CC(=O)C=Cc2ccc(O)c(OC)c2)ccc1O",
    },
    "Solvents & basics": {
        "Water": "O",
        "Methane": "C",
        "Ethanol": "CCO",
        "Acetone": "CC(C)=O",
        "Acetic acid": "CC(=O)O",
        "Urea": "NC(N)=O",
        "DMSO": "CS(C)=O",
        "THF": "C1CCOC1",
        "Chloroform": "ClC(Cl)Cl",
        "Benzene": "c1ccccc1",
        "Toluene": "Cc1ccccc1",
        "Naphthalene": "c1ccc2ccccc2c1",
    },
}
CUSTOM = "Custom SMILES"
ALL = {n: s for grp in LIBRARY.values() for n, s in grp.items()}

# ----------------------------------------------------------------------------
# Nothing OS styling
# ----------------------------------------------------------------------------
st.markdown(
    """
<style>
@import url('https://fonts.googleapis.com/css2?family=Doto:wght@400;700;900&family=Space+Mono:wght@400;700&display=swap');
:root { --red:#D71921; --fg:#FFFFFF; --dim:#8A8A8A; --line:#2A2A2A; --bg:#000000; }
html, body, [class*="css"], .stApp, .stMarkdown, label, input, textarea, button {
  font-family:'Space Mono', monospace !important;
}
.stApp {
  background-color:var(--bg);
  background-image:radial-gradient(#1c1c1c 1.2px, transparent 1.2px);
  background-size:18px 18px;
  color:var(--fg);
}
header[data-testid="stHeader"], footer, #MainMenu {display:none !important;}
.block-container {padding-top:2.2rem; max-width:1250px;}
.dot-title {font-family:'Doto', monospace !important; font-weight:900; font-size:4.4rem;
  letter-spacing:.06em; line-height:1; margin:0; color:var(--fg);}
.dot-title span {color:var(--red);}
.sub {color:var(--dim); font-size:.8rem; letter-spacing:.25em; text-transform:uppercase; margin:.4rem 0 1.6rem;}
.card {border:1px solid var(--line); border-radius:26px; padding:1.1rem 1.3rem; background:#050505;}
.lbl {color:var(--dim); font-size:.68rem; letter-spacing:.22em; text-transform:uppercase; margin-bottom:.3rem;}
.big {font-family:'Doto', monospace !important; font-weight:900; font-size:2.6rem; line-height:1.05; color:var(--fg);}
.big.red {color:var(--red);}
.unit {font-size:.75rem; color:var(--dim); font-family:'Space Mono', monospace;}
.name {font-family:'Doto', monospace !important; font-weight:900; font-size:2.2rem; letter-spacing:.04em;}
.smiles {color:var(--dim); font-size:.75rem; word-break:break-all;}
.dot {display:inline-block; width:11px; height:11px; border-radius:50%; margin-right:8px; border:1px solid var(--dim);}
.dot.on {background:var(--fg); border-color:var(--fg);} .dot.off {background:var(--red); border-color:var(--red);}
.rule {font-size:.82rem; margin:.25rem 0;}
div[data-baseweb="select"] > div, .stTextInput input {
  background:#000 !important; border:1px solid var(--line) !important; border-radius:999px !important; color:var(--fg) !important;}
.stButton > button, .stDownloadButton > button {
  background:transparent; color:var(--fg); border:1px solid var(--fg); border-radius:999px; padding:.4rem 1.4rem;
  letter-spacing:.12em; text-transform:uppercase; font-size:.75rem;}
.stButton > button:hover, .stDownloadButton > button:hover {background:var(--red); border-color:var(--red); color:#fff;}
.stRadio label, .stRadio div[role="radiogroup"] p {font-size:.8rem !important;}
svg text {fill:#fff;}
div[data-baseweb="popover"] ul, div[data-baseweb="popover"] li, div[data-baseweb="menu"] {
  background:#000 !important; color:#fff !important;}
div[data-baseweb="popover"] li:hover {background:var(--red) !important;}
.stTabs [data-baseweb="tab"] {color:var(--dim); letter-spacing:.15em;}
.stTabs [aria-selected="true"] {color:var(--fg) !important;}
.stTabs [data-baseweb="tab-highlight"] {background:var(--red) !important;}
.stTextInput input::placeholder {color:var(--dim);}
hr {border-color:var(--line);}
</style>
""",
    unsafe_allow_html=True,
)

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def parse(smiles: str):
    smiles = (smiles or "").strip()
    if not smiles:
        return None
    return Chem.MolFromSmiles(smiles)


def draw_2d(mol) -> str:
    d = rdMolDraw2D.MolDraw2DSVG(560, 420)
    o = d.drawOptions()
    o.setBackgroundColour((0, 0, 0, 0))
    o.bondLineWidth = 2
    o.updateAtomPalette({i: (1, 1, 1) for i in range(0, 120)})
    o.updateAtomPalette({8: (0.84, 0.1, 0.13)})  # oxygen = Nothing red
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    return d.GetDrawingText()


def mol_block_3d(mol):
    m = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(m, randomSeed=42) != 0:
        if AllChem.EmbedMolecule(m, useRandomCoords=True, randomSeed=42) != 0:
            return None
    try:
        AllChem.MMFFOptimizeMolecule(m, maxIters=500)
    except Exception:
        pass
    return Chem.MolToMolBlock(m)


def viewer_html(block: str, style: str, spin: bool) -> str:
    styles = {
        "Sticks": '{stick:{radius:0.14,colorscheme:"whiteCarbon"},sphere:{scale:0.22,colorscheme:"whiteCarbon"}}',
        "Spheres": '{sphere:{scale:0.35,colorscheme:"whiteCarbon"}}',
        "Lines": '{line:{colorscheme:"whiteCarbon"}}',
    }[style]
    return f"""
<div id="v" style="width:100%;height:430px;position:relative;border-radius:22px;overflow:hidden;"></div>
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"></script>
<script>
  const v = $3Dmol.createViewer("v", {{backgroundColor:"#000000"}});
  v.addModel(`{block}`, "mol");
  v.setStyle({{}}, {styles});
  v.zoomTo(); v.render();
  {"v.spin('y', 0.6);" if spin else ""}
</script>"""


def stat(label, value, unit="", red=False):
    st.markdown(
        f'<div class="card"><div class="lbl">{label}</div>'
        f'<div class="big{" red" if red else ""}">{value}<span class="unit"> {unit}</span></div></div>',
        unsafe_allow_html=True,
    )


# ----------------------------------------------------------------------------
# Header
# ----------------------------------------------------------------------------
st.markdown('<p class="dot-title">MOLECULE<span>.</span></p>', unsafe_allow_html=True)
st.markdown('<div class="sub">Structure · Properties · 3D viewer</div>', unsafe_allow_html=True)

if "pick" not in st.session_state:
    st.session_state.pick = "Aspirin"
if "cat" not in st.session_state:
    st.session_state.cat = "Medicines"

# ----------------------------------------------------------------------------
# Input column
# ----------------------------------------------------------------------------
left, right = st.columns([1, 2.2], gap="large")

with left:
    st.markdown('<div class="lbl">Source</div>', unsafe_allow_html=True)
    cats = list(LIBRARY.keys()) + [CUSTOM]
    cat = st.selectbox("Category", cats, key="cat", label_visibility="collapsed")

    if cat == CUSTOM:
        smiles = st.text_input("SMILES", value="CCO", placeholder="Paste a SMILES string")
        name = "Custom molecule"
    else:
        names = list(LIBRARY[cat].keys())
        if st.session_state.pick not in names:
            st.session_state.pick = names[0]
        name = st.selectbox("Molecule", names, key="pick", label_visibility="collapsed")
        smiles = LIBRARY[cat][name]

    if st.button("Surprise me"):
        c = random.choice(list(LIBRARY.keys()))
        st.session_state.cat = c
        st.session_state.pick = random.choice(list(LIBRARY[c].keys()))
        st.rerun()

    st.markdown(f'<div class="lbl" style="margin-top:1.4rem">Library · {len(ALL)} molecules</div>', unsafe_allow_html=True)
    search = st.text_input("Search", placeholder="search name…", label_visibility="collapsed")
    if search:
        hits = [n for n in ALL if search.lower() in n.lower()]
        if hits:
            st.caption("  ·  ".join(hits[:12]))
            hit = st.selectbox("Load", hits, label_visibility="collapsed")
            if st.button("Load molecule"):
                for c, grp in LIBRARY.items():
                    if hit in grp:
                        st.session_state.cat, st.session_state.pick = c, hit
                st.rerun()
        else:
            st.caption("no match")

mol = parse(smiles)

with right:
    if mol is None:
        st.markdown(
            '<div class="card"><div class="name" style="color:var(--red)">INVALID SMILES</div>'
            '<div class="smiles">Check the string and try again.</div></div>',
            unsafe_allow_html=True,
        )
        st.stop()

    formula = rdMolDescriptors.CalcMolFormula(mol)
    st.markdown(
        f'<div class="card"><div class="lbl">Selected</div><div class="name">{name.upper()}</div>'
        f'<div class="smiles">{formula} &nbsp;·&nbsp; {smiles}</div></div>',
        unsafe_allow_html=True,
    )

# ----------------------------------------------------------------------------
# Properties
# ----------------------------------------------------------------------------
mw = Descriptors.MolWt(mol)
logp = Crippen.MolLogP(mol)
tpsa = rdMolDescriptors.CalcTPSA(mol)
hbd = Lipinski.NumHDonors(mol)
hba = Lipinski.NumHAcceptors(mol)
rot = Lipinski.NumRotatableBonds(mol)
rings = rdMolDescriptors.CalcNumRings(mol)
heavy = mol.GetNumHeavyAtoms()

st.write("")
tabs = st.tabs(["2D", "3D", "PROPERTIES", "DRUG-LIKENESS"])

with tabs[0]:
    st.markdown(f'<div class="card" style="text-align:center">{draw_2d(mol)}</div>', unsafe_allow_html=True)

with tabs[1]:
    c1, c2 = st.columns([3, 1])
    with c2:
        style = st.radio("Style", ["Sticks", "Spheres", "Lines"], label_visibility="collapsed")
        spin = st.toggle("Spin", value=True)
    block = mol_block_3d(mol)
    with c1:
        if block:
            components.html(viewer_html(block, style, spin), height=450)
        else:
            st.warning("Could not generate a 3D conformer for this molecule.")
    if block:
        st.download_button("Download .mol", block, file_name=f"{name.replace(' ', '_')}.mol")

with tabs[2]:
    r1 = st.columns(4)
    r2 = st.columns(4)
    with r1[0]: stat("Molar mass", f"{mw:.1f}", "g/mol", red=True)
    with r1[1]: stat("LogP", f"{logp:.2f}")
    with r1[2]: stat("TPSA", f"{tpsa:.1f}", "Å²")
    with r1[3]: stat("Heavy atoms", heavy)
    with r2[0]: stat("H-donors", hbd)
    with r2[1]: stat("H-acceptors", hba)
    with r2[2]: stat("Rot. bonds", rot)
    with r2[3]: stat("Rings", rings)

with tabs[3]:
    rules = [
        ("Molar mass ≤ 500", mw <= 500),
        ("LogP ≤ 5", logp <= 5),
        ("H-donors ≤ 5", hbd <= 5),
        ("H-acceptors ≤ 10", hba <= 10),
    ]
    passed = sum(ok for _, ok in rules)
    html = '<div class="card"><div class="lbl">Lipinski rule of five</div>'
    for text, ok in rules:
        html += f'<div class="rule"><span class="dot {"on" if ok else "off"}"></span>{text}</div>'
    verdict = "ORAL DRUG-LIKE" if passed >= 3 else "UNLIKELY ORAL DRUG"
    html += f'<div class="big{"" if passed >= 3 else " red"}" style="margin-top:.9rem">{passed}/4</div>'
    html += f'<div class="lbl">{verdict}</div></div>'
    st.markdown(html, unsafe_allow_html=True)
    st.caption("Rough screening heuristic only — not medical or pharmaceutical advice.")
