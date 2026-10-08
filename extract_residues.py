#!/usr/bin/env python3
import os
from collections import defaultdict
import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

CIF_FILE = "/Users/becke410/Downloads/4LXZ.cif"
TARGET_CHAIN = "A"
TARGET_RESID = 301
DISTANCE_CUTOFF = 3.5


def parse_atoms(cif_path):
    mmcif_dict = MMCIF2Dict(cif_path)

    cols = {
        "group_PDB": mmcif_dict.get("_atom_site.group_PDB", []),
        "id": mmcif_dict.get("_atom_site.id", []),
        "type_symbol": mmcif_dict.get("_atom_site.type_symbol", []),
        "label_atom_id": mmcif_dict.get("_atom_site.label_atom_id", []),
        "label_alt_id": mmcif_dict.get("_atom_site.label_alt_id", []),
        "label_comp_id": mmcif_dict.get("_atom_site.label_comp_id", []),
        "label_asym_id": mmcif_dict.get("_atom_site.label_asym_id", []),
        "label_entity_id": mmcif_dict.get("_atom_site.label_entity_id", []),
        "label_seq_id": mmcif_dict.get("_atom_site.label_seq_id", []),
        "pdbx_PDB_ins_code": mmcif_dict.get("_atom_site.pdbx_PDB_ins_code", []),
        "Cartn_x": mmcif_dict.get("_atom_site.Cartn_x", []),
        "Cartn_y": mmcif_dict.get("_atom_site.Cartn_y", []),
        "Cartn_z": mmcif_dict.get("_atom_site.Cartn_z", []),
        "occupancy": mmcif_dict.get("_atom_site.occupancy", []),
        "B_iso_or_equiv": mmcif_dict.get("_atom_site.B_iso_or_equiv", []),
        "auth_seq_id": mmcif_dict.get("_atom_site.auth_seq_id", []),
        "auth_comp_id": mmcif_dict.get("_atom_site.auth_comp_id", []),
        "auth_asym_id": mmcif_dict.get("_atom_site.auth_asym_id", []),
        "auth_atom_id": mmcif_dict.get("_atom_site.auth_atom_id", []),
        "pdbx_PDB_model_num": mmcif_dict.get("_atom_site.pdbx_PDB_model_num", []),
    }

    n = len(cols["group_PDB"])
    atoms = []
    for i in range(n):
        model_num = cols["pdbx_PDB_model_num"][i] if cols["pdbx_PDB_model_num"] else "1"
        if model_num not in ("1", 1, " 1") and str(model_num).strip() != "1":
            continue

        auth_chain = cols["auth_asym_id"][i].strip() if cols["auth_asym_id"] else ""
        label_chain = cols["label_asym_id"][i].strip() if cols["label_asym_id"] else ""
        chain_id = auth_chain if auth_chain else label_chain

        try:
            auth_seq = int(cols["auth_seq_id"][i]) if cols["auth_seq_id"][i] and cols["auth_seq_id"][i].strip() != "." else None
        except (ValueError, TypeError):
            auth_seq = None
        try:
            label_seq = int(float(cols["label_seq_id"][i])) if cols["label_seq_id"][i] and cols["label_seq_id"][i].strip() != "." and cols["label_seq_id"][i].strip() != "?" else None
        except (ValueError, TypeError):
            label_seq = None

        resid = auth_seq if auth_seq is not None else label_seq
        if resid is None:
            continue

        icode_raw = cols["pdbx_PDB_ins_code"][i] if cols["pdbx_PDB_ins_code"] else ""
        icode = icode_raw.strip() if icode_raw else ""
        if icode == "?" or icode == ".":
            icode = ""

        resname = cols["auth_comp_id"][i].strip() if cols["auth_comp_id"] else cols["label_comp_id"][i].strip()
        atom_name = cols["auth_atom_id"][i].strip() if cols["auth_atom_id"] else cols["label_atom_id"][i].strip()

        try:
            x = float(cols["Cartn_x"][i])
            y = float(cols["Cartn_y"][i])
            z = float(cols["Cartn_z"][i])
        except (ValueError, TypeError):
            continue

        alt_loc = cols["label_alt_id"][i].strip() if cols["label_alt_id"] else ""
        if alt_loc and alt_loc not in (".", "?", "A"):
            continue

        atoms.append({
            "chain": chain_id,
            "resid": resid,
            "icode": icode,
            "resname": resname,
            "atom_name": atom_name,
            "x": x,
            "y": y,
            "z": z,
        })

    return atoms


def is_hetero_ON(atom_name):
    """Heteroatoms O/N: atom_name starts with O or N (ignoring digit prefixes like 1H)."""
    stripped = atom_name.lstrip("0123456789")
    if not stripped:
        return False
    return stripped[0] in ("O", "N")


def atom_element_ON(atom_name):
    """Return 'O' / 'N' for an O/N atom (after stripping leading digits), else None."""
    stripped = atom_name.lstrip("0123456789")
    if not stripped:
        return None
    first = stripped[0]
    return first if first in ("O", "N") else None


def compute_NO_distances(target_atoms, surrounding_res_atoms):
    """Return list of (target_atom, res_atom, distance) filtered to **cross** N-O pairs only.

    i.e. (ligand N ↔ residue O)  OR  (ligand O ↔ residue N).  O-O or N-N pairs are dropped.
    """
    target_het = []
    for a in target_atoms:
        e = atom_element_ON(a["atom_name"])
        if e is not None:
            target_het.append((a, e))
    res_het = []
    for a in surrounding_res_atoms:
        e = atom_element_ON(a["atom_name"])
        if e is not None:
            res_het.append((a, e))

    pairs = []
    for (t, te) in target_het:
        for (r, re) in res_het:
            if te == re:          # skip O-O and N-N
                continue
            d = np.sqrt((t["x"] - r["x"]) ** 2 + (t["y"] - r["y"]) ** 2 + (t["z"] - r["z"]) ** 2)
            pairs.append((t, r, d))
    pairs.sort(key=lambda p: p[2])
    return pairs


def main():
    print("Parsing CIF file...")
    atoms = parse_atoms(CIF_FILE)
    print(f"Total atoms parsed: {len(atoms)}")

    target_atoms = [
        a for a in atoms
        if a["chain"] == TARGET_CHAIN and a["resid"] == TARGET_RESID
    ]

    if not target_atoms:
        print(f"Error: Could not find chain {TARGET_CHAIN} residue {TARGET_RESID}")
        print("Available chains and residue range:")
        chain_res = defaultdict(list)
        for a in atoms:
            chain_res[a["chain"]].append(a["resid"])
        for ch in sorted(chain_res):
            rs = sorted(set(chain_res[ch]))
            print(f"  Chain {ch}: residues {rs[0]} - {rs[-1]} ({len(rs)} residues)")
        return

    target_resname = target_atoms[0]["resname"]
    target_O = sum(1 for a in target_atoms if atom_element_ON(a["atom_name"]) == "O")
    target_N = sum(1 for a in target_atoms if atom_element_ON(a["atom_name"]) == "N")
    print(f"\nTarget: Chain {TARGET_CHAIN}, Residue {TARGET_RESID} ({target_resname}) "
          f"with {len(target_atoms)} atoms ({target_O} O, {target_N} N)")
    print(f"Searching for atoms within {DISTANCE_CUTOFF}Å...  "
          f"(heteroatom distance table will show **N↔O cross pairs only**)\n")

    coords_all = np.array([[a["x"], a["y"], a["z"]] for a in atoms], dtype=np.float64)
    coords_target = np.array([[a["x"], a["y"], a["z"]] for a in target_atoms], dtype=np.float64)

    diff = coords_all[:, np.newaxis, :] - coords_target[np.newaxis, :, :]
    dist_matrix = np.sqrt(np.sum(diff ** 2, axis=2))
    min_dist_per_atom = dist_matrix.min(axis=1)
    nearby_mask = min_dist_per_atom <= DISTANCE_CUTOFF

    selected_indices = np.where(nearby_mask)[0]

    residue_keys = set()
    res_key_to_atoms = defaultdict(list)
    target_key = None

    for idx in selected_indices:
        a = atoms[idx]
        key = (a["chain"], a["resid"], a["icode"], a["resname"])
        residue_keys.add(key)
        res_key_to_atoms[key].append(a)
        if a["chain"] == TARGET_CHAIN and a["resid"] == TARGET_RESID:
            target_key = key

    sorted_residue_keys = sorted(
        residue_keys,
        key=lambda k: (k[0], k[1], k[2], k[3]),
    )
    if target_key and target_key in sorted_residue_keys and sorted_residue_keys[0] != target_key:
        sorted_residue_keys.remove(target_key)
        sorted_residue_keys.insert(0, target_key)

    nearby_count = len(sorted_residue_keys) - 1
    print(f"Found {len(sorted_residue_keys)} residues total (1 target + {nearby_count} surrounding):")
    print("-" * 80)
    header = f"{'Chain':<6}{'ResID':<8}{'ResName':<10}{'AtomName':<10}{'X':>10}{'Y':>10}{'Z':>10}"
    print(header)
    print("-" * 80)

    output_lines = [header]
    summary_rows = []
    distance_block = []  # per-residue heteroatom O/N distances to ligand O/N

    for key in sorted_residue_keys:
        chain_id, resid, icode, resname = key
        resid_str = f"{resid}{icode}" if icode else str(resid)
        res_atoms = res_key_to_atoms[key]
        res_atoms_sorted = sorted(res_atoms, key=lambda a: a["atom_name"])
        for a in res_atoms_sorted:
            line = f"{chain_id:<6}{resid_str:<8}{resname:<10}{a['atom_name']:<10}{a['x']:>10.3f}{a['y']:>10.3f}{a['z']:>10.3f}"
            print(line)
            output_lines.append(line)
        summary_rows.append((key, len(res_atoms_sorted)))

        if key != target_key:
            pairs = compute_NO_distances(target_atoms, res_atoms_sorted)
            distance_block.append((key, pairs))

    print("-" * 80)

    output_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "extracted_atoms.txt")
    with open(output_file, "w") as f:
        f.write("\n".join(output_lines) + "\n")
    print(f"\nResults saved to: {output_file}")

    print("\nSummary of selected residues:")
    print(f"{'Chain':<6}{'ResID':<8}{'ResName':<10}{'#Atoms':<8}{'Role'}")
    print("-" * 50)
    for key, num_atoms in summary_rows:
        chain_id, resid, icode, resname = key
        resid_str = f"{resid}{icode}" if icode else str(resid)
        role = "TARGET" if key == target_key else "SURROUNDING"
        print(f"{chain_id:<6}{resid_str:<8}{resname:<10}{num_atoms:<8}{role}")

    # --- Heteroatom N↔O (cross pairs only) distance block ---
    lig_label = f"{TARGET_CHAIN}:{TARGET_RESID}({target_resname})"
    dist_header = (
        f"\n\n=== N↔O cross heteroatom distances to ligand {lig_label} ==="
        f"\n{'Chain':<6}{'ResID':<8}{'ResName':<10}"
        f"{'LigandAtom':<12}{'ResAtom':<12}{'Dist(Å)':>10}"
        f"{'Type':<10}"
    )
    print(dist_header)
    print("-" * 80)

    dist_lines = [dist_header.strip(), "-" * 80]
    total_pairs = 0
    for key, pairs in distance_block:
        chain_id, resid, icode, resname = key
        resid_str = f"{resid}{icode}" if icode else str(resid)
        if not pairs:
            line_no_pairs = (f"{chain_id:<6}{resid_str:<8}{resname:<10}"
                             f"{'—':<12}{'(no N↔O)':<12}{'—':>10}"
                             f"{'':<10}")
            print(line_no_pairs)
            dist_lines.append(line_no_pairs)
            continue
        for i, (t, r, d) in enumerate(pairs):
            total_pairs += 1
            te = atom_element_ON(t["atom_name"])
            re = atom_element_ON(r["atom_name"])
            pair_type = f"{te}→{re}"
            if i == 0:
                line = (f"{chain_id:<6}{resid_str:<8}{resname:<10}"
                        f"{t['atom_name']:<12}{r['atom_name']:<12}{d:>10.3f}"
                        f"{pair_type:<10}")
            else:
                line = (f"{'':<6}{'':<8}{'':<10}"
                        f"{t['atom_name']:<12}{r['atom_name']:<12}{d:>10.3f}"
                        f"{pair_type:<10}")
            print(line)
            dist_lines.append(line)

    dist_out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "hetero_distances.txt")
    with open(dist_out, "w") as f:
        f.write("\n".join(dist_lines) + "\n")
    print(f"\nN↔O cross heteroatom distances saved to: {dist_out}  "
          f"({total_pairs} pair{'s' if total_pairs != 1 else ''})")


if __name__ == "__main__":
    main()
