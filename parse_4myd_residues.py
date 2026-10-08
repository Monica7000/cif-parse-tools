#!/usr/bin/env python3
"""Parse 4MYD.cif: extract all atoms for Chain A Ser86, His232, Asp210
   and calculate all inter-residue N↔O hydrogen-bond-like distances
   between these three residues.

Outputs:
    - 4myd_selected_residues.txt  : per-atom coordinates
    - 4myd_selected_residues.tex  : inline LaTeX table
    - 4myd_NO_hbonds.txt          : inter-residue N↔O pair distances
"""

from pathlib import Path
from itertools import combinations
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

CIF_FILE = "/Users/becke410/Downloads/4MYD.cif"
TARGET_CHAIN = "A"
TARGET_RESIDUES = [
    (86,  "SER", "Ser86"),
    (232, "HIS", "His232"),
    (210, "ASP", "Asp210"),
]
OUT_TXT     = "/Users/becke410/Eval_code/4myd_selected_residues.txt"
OUT_TEX     = "/Users/becke410/Eval_code/4myd_selected_residues.tex"
OUT_HBONDS  = "/Users/becke410/Eval_code/4myd_NO_hbonds.txt"

# H-bond geometric cutoffs (standard: 2.5–3.5 Å is typical for N-O)
HBOND_STRONG_CUTOFF = 3.0   # ≤ 3.0 Å → strong
HBOND_CUTOFF        = 3.5   # ≤ 3.5 Å → putative H-bond (moderate/weak)


def parse_atoms(cif_path):
    mmcif_dict = MMCIF2Dict(cif_path)
    n = len(mmcif_dict.get("_atom_site.group_PDB", []))
    atoms = []
    for i in range(n):
        group    = mmcif_dict.get("_atom_site.group_PDB", [""])[i]
        chain    = mmcif_dict.get("_atom_site.auth_asym_id",
                   [mmcif_dict.get("_atom_site.label_asym_id", [""])[i]])[i]
        auth_seq = mmcif_dict.get("_atom_site.auth_seq_id", [""])[i]
        lab_seq  = mmcif_dict.get("_atom_site.label_seq_id", [""])[i]
        resname  = mmcif_dict.get("_atom_site.auth_comp_id",
                   [mmcif_dict.get("_atom_site.label_comp_id", [""])[i]])[i]
        atom     = mmcif_dict.get("_atom_site.auth_atom_id",
                   [mmcif_dict.get("_atom_site.label_atom_id", [""])[i]])[i]
        alt      = mmcif_dict.get("_atom_site.label_alt_id", [""])[i]
        try:
            seq_id = int(auth_seq) if auth_seq and auth_seq != "." else (
                     int(lab_seq)  if lab_seq  and lab_seq  != "." else None)
        except ValueError:
            seq_id = None
        try:
            x = float(mmcif_dict["_atom_site.Cartn_x"][i])
            y = float(mmcif_dict["_atom_site.Cartn_y"][i])
            z = float(mmcif_dict["_atom_site.Cartn_z"][i])
        except (ValueError, KeyError):
            continue
        if seq_id is None:
            continue
        # Only take primary altloc ('' or 'A') if multiple exist
        if alt and alt not in (" ", ".", "?") and alt != "A":
            continue
        atoms.append({
            "group":   group,
            "chain":   chain,
            "seq_id":  seq_id,
            "resname": resname,
            "atom":    atom,
            "x": x, "y": y, "z": z,
        })
    return atoms


def select_targets(atoms, chain, targets):
    """Return list of (display_name, [atom_dicts...]) preserving target order."""
    target_ids = {(r[0], r[1]): r[2] for r in targets}
    seq_order  = {r[0]: i for i, r in enumerate(targets)}
    selected   = {r[0]: [] for r in targets}
    for a in atoms:
        if a["chain"] != chain:
            continue
        key = (a["seq_id"], a["resname"].upper())
        if key not in target_ids:
            continue
        selected[a["seq_id"]].append(a)
    result = []
    for seq_id, display in sorted(((r[0], target_ids[(r[0], r[1])]) for r in targets),
                                  key=lambda kv: seq_order[kv[0]]):
        # Standard atom sort: backbone first, then sidechain by name
        def atom_key(atom_name):
            backbone = ["N", "CA", "C", "O", "OXT"]
            if atom_name in backbone:
                return (0, backbone.index(atom_name), atom_name)
            return (1, atom_name[0], atom_name)
        atoms_sorted = sorted(selected[seq_id], key=lambda a: atom_key(a["atom"]))
        result.append((display, seq_id, atoms_sorted))
    return result


def atom_element_ON(atom_name):
    """Return 'O'/'N' after stripping leading digits, else None."""
    stripped = atom_name.lstrip("0123456789")
    if not stripped:
        return None
    first = stripped[0]
    return first if first in ("O", "N") else None


def collect_ON(res_atoms):
    """Return list of (atom_dict, element) for O/N atoms in a residue."""
    out = []
    for a in res_atoms:
        e = atom_element_ON(a["atom"])
        if e is not None:
            out.append((a, e))
    return out


def classify(d):
    if d <= HBOND_STRONG_CUTOFF:
        return "strong"
    if d <= HBOND_CUTOFF:
        return "putative"
    return "long-range"


def compute_inter_res_NO_hbonds(blocks):
    """Compute all cross-residue N↔O pair distances for every residue pair.

    Skips pairs within the same residue. Returns sorted list of dicts:
      {donor_res, acceptor_res, donor_atom, acceptor_atom, dist, direction, type}
    """
    results = []
    for (name_i, seq_i, atoms_i), (name_j, seq_j, atoms_j) in combinations(blocks, 2):
        on_i = collect_ON(atoms_i)
        on_j = collect_ON(atoms_j)
        if not on_i or not on_j:
            continue
        for (ai, ei) in on_i:
            for (aj, ej) in on_j:
                if ei == ej:        # only N↔O, skip N↔N / O↔O
                    continue
                d = ((ai["x"]-aj["x"])**2 + (ai["y"]-aj["y"])**2 + (ai["z"]-aj["z"])**2) ** 0.5
                # Direction: use N → O as the "canonical" H-bond direction
                if ei == "N":
                    don_res, don_atom, acc_res, acc_atom = name_i, ai, name_j, aj
                else:
                    don_res, don_atom, acc_res, acc_atom = name_j, aj, name_i, ai
                direction = f"{don_atom['atom']}({don_res})→{acc_atom['atom']}({acc_res})"
                results.append({
                    "don_res":    don_res,
                    "acc_res":    acc_res,
                    "don_atom":   don_atom,
                    "acc_atom":   acc_atom,
                    "dist":       d,
                    "direction":  direction,
                    "type":       classify(d),
                })
    results.sort(key=lambda r: r["dist"])
    return results


def write_hbond_report(hbonds, path: Path):
    header = (
        f"{'D_res':<8}{'A_res':<8}{'D_atom':<12}{'A_atom':<12}"
        f"{'Dist(Å)':>10}  {'Type':<10}  Direction"
    )
    lines = [
        f"Inter-residue N↔O distances for {', '.join(d for *_, d in TARGET_RESIDUES)} (Chain A)",
        f"Cutoffs: strong ≤ {HBOND_STRONG_CUTOFF}Å, putative ≤ {HBOND_CUTOFF}Å, else 'long-range'.",
        "",
        f"{'='*len(header)}",
        header,
        f"{'-'*len(header)}",
    ]
    if not hbonds:
        lines.append("(no N↔O cross pairs found between these residues)")
    else:
        for h in hbonds:
            lines.append(
                f"{h['don_res']:<8}{h['acc_res']:<8}"
                f"{h['don_atom']['atom']:<12}{h['acc_atom']['atom']:<12}"
                f"{h['dist']:>10.3f}  {h['type']:<10}  {h['direction']}"
            )
        lines.append("")
        strong  = sum(1 for h in hbonds if h["type"] == "strong")
        putativ = sum(1 for h in hbonds if h["type"] == "putative")
        long_r  = sum(1 for h in hbonds if h["type"] == "long-range")
        lines.append(f"Total N↔O pairs: {len(hbonds)}")
        lines.append(f"  strong   (≤ {HBOND_STRONG_CUTOFF}Å): {strong}")
        lines.append(f"  putative (≤ {HBOND_CUTOFF}Å):        {putativ}")
        lines.append(f"  long-range (>{HBOND_CUTOFF}Å):       {long_r}")
    path.write_text("\n".join(lines) + "\n")


def write_txt(blocks, path: Path):
    lines = ["ResID   ResName   AtomName           X         Y         Z"]
    for _name, seq_id, atoms in blocks:
        if not atoms:
            continue
        resname = atoms[0]["resname"]
        for a in atoms:
            lines.append(
                f"{seq_id:<7d}{resname:<10s}{a['atom']:<12s}"
                f"{a['x']:>10.3f}{a['y']:>10.3f}{a['z']:>10.3f}"
            )
    path.write_text("\n".join(lines) + "\n")


def write_latex(blocks, path: Path):
    tex = []
    tex.append(r"% Inline LaTeX table — requires \usepackage{booktabs,graphicx} in preamble.")
    tex.append(r"\begin{table}[htbp]")
    tex.append(r"\centering")
    tex.append(r"\resizebox{0.7\textwidth}{!}{%")
    tex.append(r"\begin{tabular}{@{}llcrrr@{}}")
    tex.append(r"\toprule")
    tex.append(
        r"\textbf{ResID} & \textbf{ResName} & \textbf{AtomName} & "
        r"\textbf{X (\AA)} & \textbf{Y (\AA)} & \textbf{Z (\AA)} \\"
    )
    first = True
    all_rows = 0
    for display_name, seq_id, atoms in blocks:
        if not atoms:
            continue
        resname = atoms[0]["resname"]
        if first:
            tex.append(r"\midrule")
            first = False
        else:
            tex.append(r"\cmidrule(lr){1-6}")
        title = f"{display_name} (Chain A, Res {seq_id}, {resname})"
        tex.append("\\multicolumn{6}{@{}l}{\\textbf{" + title + "}} \\\\")
        tex.append(r"\cmidrule(r){1-6}")
        for a in atoms:
            tex.append(
                f"{seq_id:<4d} & {resname:<4s} & {a['atom']:<5s} & "
                f"{a['x']:7.3f} & {a['y']:7.3f} & {a['z']:7.3f} \\\\"
            )
            all_rows += 1
    tex.append(r"\bottomrule")
    tex.append(r"\end{tabular}%")
    tex.append(r"}")
    caption = (
        f"Atoms of 4MYD chain A {', '.join(d for d, _, _ in blocks)}; "
        f"total {all_rows} atoms. Coordinates in \\AA ngstr\\\"oms."
    )
    tex.append(f"\\caption{{{caption}}}")
    tex.append(r"\label{tab:4myd_ser86_his232_asp210}")
    tex.append(r"\end{table}")
    path.write_text("\n".join(tex) + "\n")


def main():
    atoms = parse_atoms(CIF_FILE)
    blocks = select_targets(atoms, TARGET_CHAIN, TARGET_RESIDUES)

    for display, seq_id, block in blocks:
        on_atoms = collect_ON(block)
        n_O = sum(1 for _, e in on_atoms if e == "O")
        n_N = sum(1 for _, e in on_atoms if e == "N")
        print(f"{display} (Res {seq_id}): {len(block)} atoms  "
              f"(O={n_O}, N={n_N})")

    write_txt(blocks, Path(OUT_TXT))
    print(f"Wrote coordinate text  → {OUT_TXT}")

    write_latex(blocks, Path(OUT_TEX))
    print(f"Wrote coordinate LaTeX → {OUT_TEX}")

    # ----- Inter-residue N↔O H-bonds -----
    hbonds = compute_inter_res_NO_hbonds(blocks)
    write_hbond_report(hbonds, Path(OUT_HBONDS))
    print(f"Wrote N↔O H-bond list  → {OUT_HBONDS}")

    # Print a compact summary to stdout
    print("\n--- N↔O distance summary (<=3.5Å only) ---")
    short = [h for h in hbonds if h["dist"] <= HBOND_CUTOFF]
    if not short:
        print("  (no N↔O pairs within the putative H-bond cutoff of "
              f"{HBOND_CUTOFF}Å)")
    else:
        for h in short:
            print(f"  {h['dist']:.3f} Å  {h['type']:<9s}  {h['direction']}")


if __name__ == "__main__":
    main()
