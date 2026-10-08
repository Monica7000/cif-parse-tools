#!/usr/bin/env python3
"""Extract ALL ligands / HETATM residues from 4LXZ.cif.

Outputs 3 variants for each of the 2 ligand classes:
  (A) Small-molecule ligands      → 4LXZ_ligands_sm.{txt,tex,csv}
        (non-polymer HETATM, excluding water/single-atom ions if requested)
  (B) All HETATM residues         → 4LXZ_ligands_all.{txt,tex,csv}
        (includes water HOH, monoatomic ions like CL/MG/ZN etc.)

For every ligand residue the script lists:
    Ligand / ResName, Chain, ResID [+icode], AtomName, X, Y, Z.
"""
from __future__ import annotations

import csv
from collections import defaultdict, OrderedDict
from pathlib import Path

import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

# =========================================================================
CIF_FILE       = "/Users/becke410/Downloads/4LXZ.cif"
OUT_SM_TXT     = "/Users/becke410/Eval_code/4LXZ_ligands_sm.txt"
OUT_SM_TEX     = "/Users/becke410/Eval_code/4LXZ_ligands_sm.tex"
OUT_SM_CSV     = "/Users/becke410/Eval_code/4LXZ_ligands_sm.csv"
OUT_ALL_TXT    = "/Users/becke410/Eval_code/4LXZ_ligands_all.txt"
OUT_ALL_TEX    = "/Users/becke410/Eval_code/4LXZ_ligands_all.tex"
OUT_ALL_CSV    = "/Users/becke410/Eval_code/4LXZ_ligands_all.csv"

# ---- classification rule -------------------------------------------------
# Default "small-molecule / true ligand" list excludes:
#   * water (HOH / SOL / TIP / WAT / H2O)
#   * monoatomic ions (any resname of length <=2 letters, e.g. MG / ZN / NA /
#     CL / CA / FE / CU / MN / K / etc.)
# If you want to keep everything in 'SM' too, set EXCLUDE_WATER_IONS_IN_SM = False.
EXCLUDE_WATER_IONS_IN_SM = True

WATER_NAMES = {"HOH", "SOL", "WAT", "TIP", "TIP3", "H2O", "DOD", "OH2"}


# =========================================================================
# Helpers
# =========================================================================
def _first_list(d, key, default=None):
    v = d.get(key, default)
    if v is None:
        return None
    if isinstance(v, list) and len(v) == 0:
        return None
    if isinstance(v, list):
        return v
    return [v]


def _fmt_res(seq_num: int, icode: str | None) -> str:
    icode = "" if not icode or icode in ("?", ".", "") else icode.strip()
    return f"{int(seq_num)}{icode}"


def _to_int(s) -> int:
    try:
        return int(str(s).split(".")[0])
    except Exception:
        return 0


def _is_monoatomic_ion(resname: str) -> bool:
    """Heuristic: single-letter or two-letter element-only resname = ion."""
    if not resname:
        return False
    rn = resname.strip().upper()
    if len(rn) > 3:
        return False
    # If the resname is just 1-2 alphabetic characters -> ion (MG, NA, CL ...).
    # 3-letter cases like ZN already 2-letter; typical 3-letter ions like 1HP/2HP
    # are handled as ligand, not ion.
    if len(rn) == 1:
        return True
    if len(rn) == 2 and rn.isalpha():
        return True
    if len(rn) == 3 and rn[0].isdigit() and rn[1:].isalpha() and len(rn[1:])==1:
        return True   # e.g. "1K" / "2M" ... very rare.
    return False


def _atom_sort_key(name: str):
    """Keep atom names in a roughly PDB-recognizable order."""
    name = (name or "").strip()
    # Backbone-ish first: just a loose preference (ligands rarely have these)
    order = {"N": 0, "CA": 1, "C": 2, "O": 3}
    if name in order:
        return (0, order[name], name)
    # Element letter priority: C < N < O < S < P < F < Cl < Br < I < rest
    prefix = ""
    if len(name) >= 2 and name[:2].upper() in {"CL", "BR", "FE", "CU", "MN",
                                                "MG", "NA", "CA", "ZN", "CO",
                                                "NI", "CD", "HG", "PB", "AL"}:
        prefix = name[:2].upper()
    elif name and name[0].isalpha():
        prefix = name[0].upper()
    prio_map = {"C":0, "N":1, "O":2, "S":3, "P":4, "F":5, "CL":6, "BR":7,
               "I":8, "MG":9, "ZN":9, "CA":9, "NA":9, "K":9, "MN":9, "FE":9,
               "CU":9, "CL":6}
    return (1, prio_map.get(prefix, 20), name)


# =========================================================================
# Parsing
# =========================================================================
def parse_ligand_atoms(d):
    """Return list of dicts, one per HETATM atom, with all required fields.

    Strategy: walk `_atom_site` once.  Keep atoms where
      - group_PDB == 'HETATM'     (direct PDB-style marker of non-standard residue)
      OR
      - the residue's entity.type in {non-polymer, branched, macrolide, water}
    We only use group_PDB because it is always populated in wwPDB mmCIF and
    corresponds 1:1 to the "ligand + ion + water" set that users expect.
    """
    N = len(d.get("_atom_site.id", []))
    if N == 0:
        raise ValueError("No _atom_site records in cif file.")

    g = _first_list(d, "_atom_site.group_PDB", [])

    auth_chain  = _first_list(d, "_atom_site.auth_asym_id", []) or _first_list(d, "_atom_site.label_asym_id", []) or ["?"]*N
    label_chain = _first_list(d, "_atom_site.label_asym_id", []) or auth_chain
    auth_seq    = _first_list(d, "_atom_site.auth_seq_id", []) or _first_list(d, "_atom_site.label_seq_id", []) or [0]*N
    label_seq   = _first_list(d, "_atom_site.label_seq_id", []) or auth_seq
    pdb_ins     = _first_list(d, "_atom_site.pdbx_PDB_ins_code", []) or ["?"]*N
    label_alt   = _first_list(d, "_atom_site.label_alt_id",  []) or ["."]*N
    label_comp  = _first_list(d, "_atom_site.label_comp_id", []) or _first_list(d, "_atom_site.auth_comp_id", []) or ["?"]*N
    auth_comp   = _first_list(d, "_atom_site.auth_comp_id",  []) or label_comp
    label_atom  = _first_list(d, "_atom_site.label_atom_id", []) or _first_list(d, "_atom_site.auth_atom_id", []) or ["?"]*N
    auth_atom   = _first_list(d, "_atom_site.auth_atom_id",  []) or label_atom

    Cartn_x = _first_list(d, "_atom_site.Cartn_x", [])
    Cartn_y = _first_list(d, "_atom_site.Cartn_y", [])
    Cartn_z = _first_list(d, "_atom_site.Cartn_z", [])

    # entity id per atom -> allows later classification
    atom_entity = _first_list(d, "_atom_site.label_entity_id", []) or ["?"]*N

    atoms = []
    for i in range(N):
        group = (g[i] if i < len(g) else "?") or "?"
        if str(group).strip().upper() != "HETATM":
            continue
        chain = auth_chain[i] if i < len(auth_chain) and auth_chain[i] not in (None,"?","") else (label_chain[i] if i < len(label_chain) else "?")
        seqid = auth_seq[i]   if i < len(auth_seq)  and auth_seq[i]  not in (None,"?","") else (label_seq[i] if i < len(label_seq) else 0)
        ins   = pdb_ins[i]    if i < len(pdb_ins)   else "?"
        alt   = label_alt[i]  if i < len(label_alt) else "."
        comp  = auth_comp[i]  if i < len(auth_comp) and auth_comp[i]  not in (None,"?","") else (label_comp[i] if i < len(label_comp) else "?")
        atom  = auth_atom[i]  if i < len(auth_atom) and auth_atom[i]  not in (None,"?","") else (label_atom[i] if i < len(label_atom) else "?")
        ent   = atom_entity[i] if i < len(atom_entity) else "?"
        try:
            x = float(Cartn_x[i])
            y = float(Cartn_y[i])
            z = float(Cartn_z[i])
        except Exception:
            continue
        # Altloc filter: prefer alt == "." (main) OR "A" if several present
        # We'll keep everything first; de-dup in a 2nd pass.
        atoms.append(dict(
            chain=str(chain).strip(),
            seqid=_to_int(seqid),
            ins="" if ins in ("?","",".",None) else str(ins).strip(),
            alt=""  if alt in ("?","",".",None) else str(alt).strip(),
            resname=str(comp).strip().upper(),
            atom=str(atom).strip(),
            x=x, y=y, z=z,
            entity_id=str(ent).strip(),
        ))

    # ---- alt-loc dedup ----
    #   group atoms by (chain, seqid, ins, resname, atom).  If alt == "" is
    #   present keep it; else keep the first alt alphabetically ("A" usually).
    gmap = OrderedDict()
    for a in atoms:
        key = (a["chain"], a["seqid"], a["ins"], a["resname"], a["atom"])
        bucket = gmap.setdefault(key, [])
        bucket.append(a)
    filtered = []
    for bucket in gmap.values():
        if len(bucket) == 1:
            filtered.append(bucket[0]); continue
        main = [b for b in bucket if b["alt"] == ""]
        if main:
            filtered.append(main[0]); continue
        bucket.sort(key=lambda b: (b["alt"], b["x"]))
        filtered.append(bucket[0])
    return filtered


def classify_ligand_residues(atoms):
    """Group atoms to residues, attach is_water / is_ion / is_small_mol flags.

    Returns: residues_ordered (list of dict with metadata + sorted atoms list).
    """
    resmap = OrderedDict()
    for a in atoms:
        key = (a["chain"], a["seqid"], a["ins"], a["resname"])
        if key not in resmap:
            resmap[key] = dict(
                chain=a["chain"], seqid=a["seqid"], ins=a["ins"],
                resname=a["resname"], entity_id=a["entity_id"], atoms=[],
                is_water=a["resname"] in WATER_NAMES,
                is_ion=_is_monoatomic_ion(a["resname"]),
            )
        resmap[key]["atoms"].append(a)

    residues = list(resmap.values())
    for r in residues:
        r["num_atoms"] = len(r["atoms"])
        r["atoms"].sort(key=lambda a: _atom_sort_key(a["atom"]))
        r["is_small_mol"] = (not r["is_water"]) and (not r["is_ion"])
    residues.sort(key=lambda r: (
        0 if r["is_small_mol"] else (1 if r["is_ion"] else 2),  # small-mol first
        r["resname"], r["chain"], r["seqid"], r["ins"],
    ))
    return residues


# =========================================================================
# Writers
# =========================================================================
def _ligand_header_for_TXT(residues_list, subset_name: str) -> str:
    n_res  = len(residues_list)
    n_atom = sum(r["num_atoms"] for r in residues_list)
    uniq = sorted({r["resname"] for r in residues_list})
    return (f"4LXZ ligands — {subset_name}.   "
            f"Residues: {n_res}, Atoms: {n_atom}, "
            f"unique resnames: {len(uniq)} ({', '.join(uniq)})")


def write_txt(residues, path: Path, subset_name: str):
    intro = _ligand_header_for_TXT(residues, subset_name)
    hdr = (f"{'ResName':<10}{'Chain':>6}{'ResID':>8}  "
           f"{'Atom':<8}{'X':>10}{'Y':>10}{'Z':>10}")
    lines = [intro, "", "=" * len(hdr), hdr, "-" * len(hdr)]
    for i, r in enumerate(residues):
        resid_str = _fmt_res(r["seqid"], r["ins"])
        for j, a in enumerate(r["atoms"]):
            if j == 0:
                lines.append(
                    f"{r['resname']:<10}{r['chain']:>6}{resid_str:>8}  "
                    f"{a['atom']:<8}{a['x']:>10.3f}{a['y']:>10.3f}{a['z']:>10.3f}"
                )
            else:
                # blank out residue name / chain / resid after first row
                lines.append(
                    f"{'':<10}{'':>6}{'':>8}  "
                    f"{a['atom']:<8}{a['x']:>10.3f}{a['y']:>10.3f}{a['z']:>10.3f}"
                )
        if i != len(residues) - 1:
            lines.append("")
    path.write_text("\n".join(lines) + "\n")


def write_tex(residues, path: Path, subset_name: str):
    caption_suffix = subset_name.replace("_", r"\_")
    label = subset_name.lower().replace(" ", "_")
    lines = [
        f"% 4LXZ ligand extraction — {subset_name}",
        r"% Requires: \usepackage{booktabs,longtable}",
        r"\begin{longtable}{@{}lcrr@{}}",
        r"  \toprule",
        r"  Residue (Chain, ID) & Atom & \multicolumn{1}{c}{X} & \multicolumn{1}{c}{Y} & \multicolumn{1}{c}{Z} \\",
        r"  \midrule",
        r"  \endfirsthead",
        r"  \multicolumn{5}{c}{\textit{Continued from previous page}} \\",
        r"  \toprule Residue (Chain, ID) & Atom & \multicolumn{1}{c}{X} & \multicolumn{1}{c}{Y} & \multicolumn{1}{c}{Z} \\ \midrule",
        r"  \endhead",
        r"  \bottomrule",
        r"  \endfoot",
        r"  \bottomrule",
    ]
    cap = (r"  \caption{All "
           + caption_suffix
           + r" residues of PDB 4LXZ with atom names and Cartesian coordinates.}")
    lines.append(cap)
    lines.append(f"  \\label{{tab:4lxz_lig_{label}}}")
    lines.append(r"  \endlastfoot")

    for i, r in enumerate(residues):
        resid_str = _fmt_res(r["seqid"], r["ins"])
        label_tex = (
            f"{r['resname']} ({r['chain']}, {resid_str})"
        ).replace("&", r"\&").replace("%", r"\%").replace("_", r"\_")
        for j, a in enumerate(r["atoms"]):
            if j == 0:
                first_col = label_tex
            else:
                first_col = ""
            lines.append(
                f"  {first_col} & {a['atom']} & "
                f"{a['x']:.3f} & {a['y']:.3f} & {a['z']:.3f} \\\\"
            )
        if i != len(residues) - 1:
            lines.append(r"  \cmidrule(r){1-5}")
    lines.append(r"\end{longtable}")
    path.write_text("\n".join(lines) + "\n")


def write_csv(residues, path: Path):
    fields = ["resname", "chain", "seqid", "ins", "resid",
              "atom", "x", "y", "z", "num_atoms_in_residue",
              "is_water", "is_ion", "is_small_mol", "entity_id"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in residues:
            resid = _fmt_res(r["seqid"], r["ins"])
            for a in r["atoms"]:
                w.writerow(dict(
                    resname=r["resname"], chain=r["chain"], seqid=r["seqid"],
                    ins=r["ins"], resid=resid, atom=a["atom"],
                    x=a["x"], y=a["y"], z=a["z"],
                    num_atoms_in_residue=r["num_atoms"],
                    is_water=r["is_water"], is_ion=r["is_ion"],
                    is_small_mol=r["is_small_mol"], entity_id=r["entity_id"],
                ))


# =========================================================================
# Main
# =========================================================================
def main():
    print(f"Parsing {CIF_FILE} ...")
    d = MMCIF2Dict(CIF_FILE)
    atoms = parse_ligand_atoms(d)
    residues = classify_ligand_residues(atoms)

    n_water = sum(1 for r in residues if r["is_water"])
    n_ion   = sum(1 for r in residues if r["is_ion"])
    n_sm    = sum(1 for r in residues if r["is_small_mol"])
    uniq_sm  = sorted({r["resname"] for r in residues if r["is_small_mol"]})
    uniq_ion = sorted({r["resname"] for r in residues if r["is_ion"]})
    print(f"Total HETATM residues parsed: {len(residues)}")
    print(f"  small-molecule ligands: {n_sm}  (unique names: {uniq_sm})")
    print(f"  water residues:         {n_water}")
    print(f"  monoatomic ions:        {n_ion}  (unique names: {uniq_ion})")

    res_sm  = [r for r in residues if r["is_small_mol"]]
    res_all = residues

    # ---- write small-molecules ----
    write_txt(res_sm,  Path(OUT_SM_TXT),  "small_molecule_ligands (no water/ions)")
    write_tex(res_sm,  Path(OUT_SM_TEX),  "small_molecule_ligands")
    write_csv(res_sm,  Path(OUT_SM_CSV))
    print(f"\nSmall-molecule ligands only → wrote:")
    print(f"  TXT : {OUT_SM_TXT}")
    print(f"  TEX : {OUT_SM_TEX}")
    print(f"  CSV : {OUT_SM_CSV}")

    # ---- write ALL HETATM (incl water/ions) ----
    write_txt(res_all, Path(OUT_ALL_TXT), "all_HETATM (ligands + ions + water)")
    write_tex(res_all, Path(OUT_ALL_TEX), "all_HETATM")
    write_csv(res_all, Path(OUT_ALL_CSV))
    print(f"\nAll HETATM residues → wrote:")
    print(f"  TXT : {OUT_ALL_TXT}")
    print(f"  TEX : {OUT_ALL_TEX}")
    print(f"  CSV : {OUT_ALL_CSV}")


if __name__ == "__main__":
    main()
