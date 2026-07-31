#!/usr/bin/env python3
"""
protein_rna_contacts_full_pipeline_WWscale.py
=======================================
Single self-contained script. For 4 MDA5 structures (+/-MA2, crossed with
+/-ADP-AlF4), analyzed COMPLETELY INDEPENDENTLY -- no Kabsch superposition,
no RNA residue-registration/renumbering, no attempt to match a given RNA
nucleotide's identity across structures. Every residue number reported is
read directly from that structure's own native PDB numbering.

Produces, at a single fixed heavy-atom cutoff (default 4.5 A):

  1. FOUR CSVs, one per condition/state combination, each listing every
     INDIVIDUAL heavy-atom pair contact (not aggregated to residue level):
       atom_contacts_ADPAlF4_plusMA2.csv
       atom_contacts_ADPAlF4_minusMA2.csv
       atom_contacts_noADPAlF4_plusMA2.csv
       atom_contacts_noADPAlF4_minusMA2.csv
     Columns: partner_type, prot_chain/resnum/resname/atom,
     partner_chain/resnum/resname/atom, distance_A.
     Covers protein-protein, protein-RNA (both strands, both protomers), and
     protein-ligand contacts.

  2. FOUR plots, one per condition x subunit (A/B), showing per-protein-
     residue net RNA-contact change (summed across every RNA nucleotide that
     residue touches, collapsing nucleotide identity):
       protRNA_A_ADPAlF4_native.html   protRNA_B_ADPAlF4_native.html
       protRNA_A_noADPAlF4_native.html protRNA_B_noADPAlF4_native.html
     x = protein residue (native numbering), y = net contact change,
     size = total RNA contact count in whichever structure (+/-MA2) is
     stronger for that residue, color = Wimley-White hydrophobicity scale
     of the protein residue.

USAGE
-----
    python3 protein_rna_contacts_full_pipeline_WWscale.py

Edit CUTOFF / FILES / chain IDs / OUTDIR in the CONFIG section as needed.
REQUIRES: numpy, scipy, pandas, plotly
"""

import os
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import plotly.graph_objects as go


# =====================================================================================
# CONFIG
# =====================================================================================

CUTOFF = 4.5   # heavy-atom contact distance cutoff, Angstrom

CONDITIONS = {
    # condition_label: (plus_MA2_path, minus_MA2_path)
    "ADPAlF4":   ("MA2-ADPAlFx.pdb",  "noMA2-ADPAlF4.pdb"),
    "noADPAlF4": ("MA2.pdb",          "noMA2.pdb"),
}

PROT_A, PROT_B = "A", "B" #Specify chain names here. It should match the chain names in the PDB file
RNA_X, RNA_Y = "X", "Y"
LIGAND_RESNAME = "LIG"

OUTDIR = "outputs"


# =====================================================================================
# PDB PARSING (native numbering preserved throughout -- never modified)
# =====================================================================================

AA3 = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
RNA_RESNAMES = set("A C G U".split())


def parse_pdb(path, ligand_resname=LIGAND_RESNAME):
    atoms = []
    with open(path) as f:
        for line in f:
            rec = line[0:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue   
            if line[16] not in (" "):
                continue
            atomname = line[12:16].strip()
            element = line[76:78].strip()
            if element == "":
                element = atomname[0] if not atomname[0].isdigit() else atomname[1]
            if element == "H":
                continue   #skip over hydrogen atoms
            resname = line[17:20].strip()
            chain = line[21]
            resnum = int(line[22:26])   # native PDB residue number
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            if resname in AA3:
                mol_type = "protein"
            elif resname in RNA_RESNAMES:
                mol_type = "rna"
            elif resname == ligand_resname:
                mol_type = "ligand"
            else:
                mol_type = "other"
            atoms.append(dict(chain=chain, resname=resname, resnum=resnum, atomname=atomname,
                               element=element, x=x, y=y, z=z, mol_type=mol_type))
    return atoms


def select(atoms, chain, mol_type):
    return [a for a in atoms if a["chain"] == chain and a["mol_type"] == mol_type]


# =====================================================================================
# PER-ATOM HYDROPHOBICITY / ELECTROSTATIC ATOM-TYPE SCORING (+1 hydrophobic, 0 polar, -1 charged)
# =====================================================================================

CHARGED_PROTEIN_ATOMS = {
    ("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2"),
    ("LYS", "NZ"), ("ARG", "NH1"), ("ARG", "NH2"), ("ARG", "NE"),
    ("HIS", "ND1"), ("HIS", "NE2"),
}
BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}
RNA_PHOSPHATE_ATOMS = {"P", "OP1", "OP2", "O1P", "O2P"}
RNA_RIBOSE_ATOMS = {"O2'", "O3'", "O4'", "O5'", "C1'", "C2'", "C3'", "C4'", "C5'"}


def protein_atom_score(resname, atomname, element):
    if (resname, atomname) in CHARGED_PROTEIN_ATOMS:
        return -1.0
    if atomname in BACKBONE_ATOMS:
        return 0.0 if atomname in ("N", "O", "OXT") else 0.3
    if element == "S":
        return 0.5
    if element == "C":
        return 1.0
    return 0.0


def rna_atom_score(resname, atomname, element):
    if atomname in RNA_PHOSPHATE_ATOMS or element == "P":
        return -1.0
    if atomname in RNA_RIBOSE_ATOMS:
        return 0.3 if element == "C" else 0.0
    if element == "C":
        return 1.0
    return 0.0


def ligand_atom_score(resname, atomname, element):
    if element == "C":
        return 1.0
    if element == "S":
        return 0.5
    return 0.0


def atom_score(mol_type, resname, atomname, element):
    if mol_type == "protein":
        return protein_atom_score(resname, atomname, element)
    if mol_type == "rna":
        return rna_atom_score(resname, atomname, element)
    if mol_type == "ligand":
        return ligand_atom_score(resname, atomname, element)
    return None


# =====================================================================================
# WIMLEY-WHITE WHOLE-RESIDUE INTERFACIAL HYDROPHOBICITY SCALE
# Wimley, W.C. & White, S.H. (1996) "Experimentally determined hydrophobicity scale
# for proteins at membrane interfaces." Nat. Struct. Biol. 3, 842-848. Table 1.
#
# Values below are transcribed directly from Table 1 of the paper using the pH 8 
# whole-residue free energy of transfer, dG^residue, in kcal/mol. 
# For the four titratable side chains this is the IONIZED value
# (footnote 4: Asp-, Glu-, Lys+, Arg+, obtained via the pH 2 measurement per Eq. 1b),
# and for His this is the NEUTRAL/un-ionized value (footnote 5, His0, the
# physiologically dominant form at pH 7).
#
# SIGN CONVENTION:
# the paper defines dG as the free energy of transfer FROM BILAYER TO WATER (see
# Table 1 header and Fig. 1/2 axis labels), so:
#     MORE POSITIVE  = more costly to leave the bilayer = MORE HYDROPHOBIC
#     MORE NEGATIVE   = favors water = MORE HYDROPHILIC / charged
# E.g. Trp = +1.85 (aromatic/hydrophobic), Glu = -2.02 (charged).
# =====================================================================================

WIMLEY_WHITE_SCALE = {
    "TRP": 1.85, "PHE": 1.13, "TYR": 0.94, "LEU": 0.56, "ILE": 0.31,
    "MET": 0.23, "CYS": 0.24, "VAL": -0.07, "ALA": -0.17, "THR": -0.14,
    "GLY": -0.01, "SER": -0.13, "GLN": -0.58, "ASN": -0.42, "PRO": -0.45,
    "HIS": -0.17,          # neutral His0
    "ASP": -1.23,          # ionized Asp-
    "GLU": -2.02,          # ionized Glu-
    "LYS": -0.99,          # ionized Lys+
    "ARG": -0.81,          # ionized Arg+
}


def wimley_white_score(resname):
    """Return the Wimley-White interfacial free energy of transfer, bilayer-to-water
    (kcal/mol), for a residue name, or NaN if not recognized. More positive = more
    hydrophobic; more negative = more hydrophilic/charged (see convention note above)."""
    return WIMLEY_WHITE_SCALE.get(resname, np.nan)


# =====================================================================================
# PART 1: ATOM-LEVEL CONTACT ENUMERATION (one row per atom pair, native numbering)
# =====================================================================================

def raw_atom_pairs(atoms_1, atoms_2, partner_label, cutoff=CUTOFF):
    """Every individual atom-atom contact within cutoff -- no aggregation."""
    if len(atoms_1) == 0 or len(atoms_2) == 0:
        return []
    coords1 = np.array([[a["x"], a["y"], a["z"]] for a in atoms_1])
    coords2 = np.array([[a["x"], a["y"], a["z"]] for a in atoms_2])
    tree2 = cKDTree(coords2)
    neighbor_lists = tree2.query_ball_point(coords1, r=cutoff)
    rows = []
    for i, neighbors in enumerate(neighbor_lists):
        if not neighbors:
            continue
        a1 = atoms_1[i]
        for j in neighbors:
            a2 = atoms_2[j]
            d = float(np.linalg.norm(coords1[i] - coords2[j]))
            rows.append(dict(
                partner_type=partner_label,
                prot_chain=a1["chain"], prot_resnum=a1["resnum"], prot_resname=a1["resname"],
                prot_atom=a1["atomname"],
                partner_chain=a2["chain"], partner_resnum=a2["resnum"], partner_resname=a2["resname"],
                partner_atom=a2["atomname"],
                distance_A=round(d, 3),
            ))
    return rows


def build_atom_level_csv(path, cutoff=CUTOFF):
    """All individual atom-pair contacts for ONE structure: protein-protein,
    protein-RNA (both strands, both protomers), protein-ligand."""
    atoms = parse_pdb(path)
    A_prot = select(atoms, PROT_A, "protein")
    B_prot = select(atoms, PROT_B, "protein")
    X_rna  = select(atoms, RNA_X, "rna")
    Y_rna  = select(atoms, RNA_Y, "rna")
    A_lig  = select(atoms, PROT_A, "ligand")
    B_lig  = select(atoms, PROT_B, "ligand")

    rows = []
    rows += raw_atom_pairs(A_prot, B_prot, "protein-protein(A-B)", cutoff)
    rows += raw_atom_pairs(A_prot, X_rna,  "protein-RNA(A-X)", cutoff)
    rows += raw_atom_pairs(A_prot, Y_rna,  "protein-RNA(A-Y)", cutoff)
    rows += raw_atom_pairs(A_prot, A_lig,  "protein-ligand(A)", cutoff)
    rows += raw_atom_pairs(B_prot, X_rna,  "protein-RNA(B-X)", cutoff)
    rows += raw_atom_pairs(B_prot, Y_rna,  "protein-RNA(B-Y)", cutoff)
    rows += raw_atom_pairs(B_prot, B_lig,  "protein-ligand(B)", cutoff)

    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values(["prot_chain", "prot_resnum", "distance_A"]).reset_index(drop=True)
    return df


# =====================================================================================
# PART 2: RESIDUE-LEVEL RNA CONTACT SUMMARY (for the plots only -- collapses across
# RNA nucleotide identity within each structure independently)
# =====================================================================================

def residue_pair_contacts_rna(atoms_prot, atoms_rna, cutoff=CUTOFF):
    """Residue-pair level (protein residue <-> RNA residue): atom-pair count and
    mean hydrophobicity for that pair. Used only as an intermediate for the
    per-residue collapse below (the atom-level CSV above is independent of this)."""
    if len(atoms_prot) == 0 or len(atoms_rna) == 0:
        return pd.DataFrame(columns=["prot_resnum", "prot_resname", "n_atom_contacts", "hydrophobicity_mean"])
    coords1 = np.array([[a["x"], a["y"], a["z"]] for a in atoms_prot])
    coords2 = np.array([[a["x"], a["y"], a["z"]] for a in atoms_rna])
    tree2 = cKDTree(coords2)
    neighbor_lists = tree2.query_ball_point(coords1, r=cutoff)
    pairs = defaultdict(list)
    for i, neighbors in enumerate(neighbor_lists):
        if not neighbors:
            continue
        a1 = atoms_prot[i]
        s1 = atom_score("protein", a1["resname"], a1["atomname"], a1["element"])
        for j in neighbors:
            a2 = atoms_rna[j]
            s2 = atom_score("rna", a2["resname"], a2["atomname"], a2["element"])
            key = (a1["resnum"], a1["resname"], a2["resnum"])
            pairs[key].append((s1 + s2) / 2.0)
    rows = []
    for (rn1, rname1, rn2), scores in pairs.items():
        rows.append(dict(prot_resnum=rn1, prot_resname=rname1,
                          n_atom_contacts=len(scores), hydrophobicity_mean=float(np.mean(scores))))
    return pd.DataFrame(rows)


def collapse_per_residue(atoms_prot, atoms_X, atoms_Y, cutoff=CUTOFF):
    """Sum total RNA contact count per protein residue (both strands combined),
    plus unweighted mean hydrophobicity across the distinct RNA residues contacted."""
    dfx = residue_pair_contacts_rna(atoms_prot, atoms_X, cutoff)
    dfy = residue_pair_contacts_rna(atoms_prot, atoms_Y, cutoff)
    df = pd.concat([dfx, dfy], ignore_index=True)
    if len(df) == 0:
        return pd.DataFrame(columns=["prot_resnum", "prot_resname", "n_total", "hphob_unweighted"])
    agg = df.groupby(["prot_resnum", "prot_resname"]).agg(
        n_total=("n_atom_contacts", "sum"),
        hphob_unweighted=("hydrophobicity_mean", "mean"),
    ).reset_index()
    return agg


def build_residue_net_change(plus_path, minus_path, chain, cutoff=CUTOFF):
    atoms_plus = parse_pdb(plus_path)
    atoms_minus = parse_pdb(minus_path)
    prot_p = select(atoms_plus, chain, "protein")
    prot_m = select(atoms_minus, chain, "protein")
    X_p, Y_p = select(atoms_plus, RNA_X, "rna"), select(atoms_plus, RNA_Y, "rna")
    X_m, Y_m = select(atoms_minus, RNA_X, "rna"), select(atoms_minus, RNA_Y, "rna")

    cp = collapse_per_residue(prot_p, X_p, Y_p, cutoff).set_index("prot_resnum")
    cm = collapse_per_residue(prot_m, X_m, Y_m, cutoff).set_index("prot_resnum")

    resname_lookup = {a["resnum"]: a["resname"] for a in prot_p}
    resname_lookup.update({a["resnum"]: a["resname"] for a in prot_m})

    all_res = sorted(set(cp.index) | set(cm.index))
    rows = []
    for r in all_res:
        n_plus = float(cp["n_total"].get(r, 0.0))
        n_minus = float(cm["n_total"].get(r, 0.0))
        strongest_total = n_plus if n_plus >= n_minus else n_minus
        strongest_state = "+MA2" if n_plus >= n_minus else "-MA2"
        resname = resname_lookup.get(r, "?")
        ww_score = wimley_white_score(resname)   # color = property of the amino acid itself, not the contacts
        rows.append(dict(prot_resnum=r, prot_resname=resname,
                          n_plus=n_plus, n_minus=n_minus, net=n_plus - n_minus,
                          strongest_total=strongest_total, strongest_state=strongest_state,
                          ww_hydrophobicity=ww_score))
    return pd.DataFrame(rows)


# =====================================================================================
# PLOTTING
# =====================================================================================

def make_plot(coll, chain_label, condition_label, xrange, yrange, outpath, cutoff=CUTOFF):
    # Wimley-White convention (paper's own dG = bilayer-to-water transfer energy):
    # POSITIVE = hydrophobic, NEGATIVE = hydrophilic/charged (see comment above
    # WIMLEY_WHITE_SCALE). Plotted with raw values, no sign flip. "RdBu_r" maps
    # high/positive values to red and low/negative values to blue (verified
    # empirically), so red = hydrophobic falls out correctly with no inversion.
    color_vals = coll["ww_hydrophobicity"]
    hover = [f"{r.prot_resname}{r.prot_resnum}<br>net:{r.net:+.0f}<br>"
             f"strongest structure: {r.strongest_state} (total={r.strongest_total:.0f})<br>"
             f"Wimley-White dG (bilayer-to-water): {r.ww_hydrophobicity:.2f} kcal/mol "
             f"({'hydrophobic' if r.ww_hydrophobicity > 0 else 'hydrophilic/charged'})"
             for r in coll.itertuples()]
    fig = go.Figure(go.Scatter(
        x=coll["prot_resnum"], y=coll["net"], mode="markers",
        marker=dict(size=8 + 3.2 * np.sqrt(coll["strongest_total"]), color=color_vals,
                    colorscale="RdBu_r", cmin=-2.1, cmax=2.1, showscale=True,
                    colorbar=dict(title="&Delta;G (Wimley-White)<br>bilayer-to-water<br>"
                                        "kcal/mol<br>positive = hydrophobic"),
                    line=dict(width=0.8, color="rgba(0,0,0,0.4)")),
        hovertext=hover, hoverinfo="text"))
    fig.add_hline(y=0, line=dict(color="black", width=1, dash="dot"))
    fig.update_xaxes(title_text="Residue number (native PDB numbering)", range=xrange)
    fig.update_yaxes(title_text="net RNA-contact change<br>(+MA2 minus -MA2)", range=yrange)
    fig.update_layout(
        height=500, width=850, template="plotly_white",
        title=f"Protein({chain_label})-RNA contact remodeling, {condition_label}<br>"
              f"<sup>No RNA registration/realignment used -- counts independent per structure "
              f"(cutoff={cutoff} A) &middot; color = Wimley-White residue hydrophobicity "
              f"(verified vs. Table 1), independent of contact data</sup>",
    )
    fig.write_html(outpath, include_plotlyjs="cdn")


# =====================================================================================
# MAIN
# =====================================================================================

def main():
    os.makedirs(OUTDIR, exist_ok=True)

    # ---- PART 1: 4 atom-level CSVs ----
    csv_name_map = {
        ("ADPAlF4", "+MA2"):   "atom_contacts_ADPAlF4_plusMA2.csv",
        ("ADPAlF4", "-MA2"):   "atom_contacts_ADPAlF4_minusMA2.csv",
        ("noADPAlF4", "+MA2"): "atom_contacts_noADPAlF4_plusMA2.csv",
        ("noADPAlF4", "-MA2"): "atom_contacts_noADPAlF4_minusMA2.csv",
    }
    print("=== Writing 4 atom-level contact CSVs ===")
    for condition_label, (plus_path, minus_path) in CONDITIONS.items():
        for path, state_label in [(plus_path, "+MA2"), (minus_path, "-MA2")]:
            df = build_atom_level_csv(path)
            fname = csv_name_map[(condition_label, state_label)]
            outpath = os.path.join(OUTDIR, fname)
            df.to_csv(outpath, index=False)
            print(f"  {fname}: {len(df)} individual atom-pair rows  (file={path})")

    # ---- PART 2: 4 plots ----
    print("\n=== Building 4 net-change plots ===")
    per_condition_chain = {}
    for condition_label, (plus_path, minus_path) in CONDITIONS.items():
        per_condition_chain[condition_label] = {}
        for chain in [PROT_A, PROT_B]:
            coll = build_residue_net_change(plus_path, minus_path, chain)
            per_condition_chain[condition_label][chain] = coll
            print(f"  {condition_label} chain {chain}: {len(coll)} residues, "
                  f"total net = {coll['net'].sum():+.0f}" if len(coll) else
                  f"  {condition_label} chain {chain}: no RNA contacts found")

    # shared x-range (scan protein resnums across all 4 files)
    all_resnums = []
    for condition_label, (plus_path, minus_path) in CONDITIONS.items():
        for path in [plus_path, minus_path]:
            atoms = parse_pdb(path)
            all_resnums += [a["resnum"] for a in atoms if a["mol_type"] == "protein" and a["chain"] in (PROT_A, PROT_B)]
    xrange = [min(all_resnums) - 5, max(all_resnums) + 5]

    # shared y-range
    all_nets = []
    for chains in per_condition_chain.values():
        for coll in chains.values():
            if len(coll):
                all_nets.extend(coll["net"].tolist())
    y_abs_max = max(abs(min(all_nets)), abs(max(all_nets))) if all_nets else 1
    yrange = [-(y_abs_max * 1.05), (y_abs_max * 1.05)]

    for condition_label, chains in per_condition_chain.items():
        for chain, coll in chains.items():
            outpath = os.path.join(OUTDIR, f"protRNA_{chain}_{condition_label}_native.html")
            make_plot(coll, chain, condition_label, xrange, yrange, outpath)
            print(f"  wrote {outpath}")

    print("\nDone.")


if __name__ == "__main__":
    main()
