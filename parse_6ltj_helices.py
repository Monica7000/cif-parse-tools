#!/usr/bin/env python3
"""Parse helix information (alpha helix / pi helix / 3_10 helix) from 6LTJ.cif.

Source: mmCIF category `_struct_conf` (secondary structure annotations provided
by the PDB authors, not recomputed from coordinates).

Outputs:
    6ltj_helices.txt : human-readable table (chain | id | start | end | length | type)
    6ltj_helices.tex : inline LaTeX fragment (paste into a paper)
    6ltj_helices.csv : CSV (chain,id,beg,beg_ins,end,end_ins,length,type,details)
"""

from pathlib import Path
import csv
from collections import defaultdict
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

CIF_FILE = "/Users/becke410/Downloads/7C4J.cif"
OUT_TXT  = "/Users/becke410/Eval_code/7C4J_helices.txt"
OUT_TEX  = "/Users/becke410/Eval_code/7C4J_helices.tex"
OUT_CSV  = "/Users/becke410/Eval_code/7C4J_helices.csv"
OUT_CHAINS_TXT = "/Users/becke410/Eval_code/7C4J_chains.txt"
OUT_CHAINS_TEX = "/Users/becke410/Eval_code/7C4J_chains.tex"
OUT_CHAINS_CSV = "/Users/becke410/Eval_code/7C4J_chains.csv"

NON_HELIX_PREFIX = ("SHEET", "STRN", "TURN")


def _get_list(d, keys, default=None):
    for k in keys:
        if k in d and isinstance(d[k], list) and len(d[k]):
            return d[k]
    if default is not None:
        return default
    return None


def _first_list(d, keys, n):
    for k in keys:
        v = d.get(k)
        if isinstance(v, list) and len(v) == n:
            return v
        if isinstance(v, list) and len(v) > 0:
            return v
    return [""] * n


def _toint(x):
    try:  return int(x)
    except Exception:  return None


def parse_helices(cif_path: str):
    """Return list[dict(chain,id,beg,beg_ins,end,end_ins,length,type,details)]."""
    d = MMCIF2Dict(cif_path)
    return parse_helices_from_dict(d)


def parse_helices_from_dict(d):
    """Same as parse_helices, but from an already-parsed MMCIF2Dict."""
    helices = []

    ctype = d.get("_struct_conf.conf_type_id", [])
    if not ctype:
        return helices
    n = len(ctype)

    beg_chain = _first_list(d, ["_struct_conf.beg_auth_asym_id",
                                "_struct_conf.beg_label_asym_id"], n)
    end_chain = _first_list(d, ["_struct_conf.end_auth_asym_id",
                                "_struct_conf.end_label_asym_id"], n)
    beg_seq   = _first_list(d, ["_struct_conf.beg_auth_seq_id",
                                "_struct_conf.beg_label_seq_id"], n)
    end_seq   = _first_list(d, ["_struct_conf.end_auth_seq_id",
                                "_struct_conf.end_label_seq_id"], n)
    beg_ins   = d.get("_struct_conf.beg_pdbx_PDB_ins_code", [""] * n)
    end_ins   = d.get("_struct_conf.end_pdbx_PDB_ins_code", [""] * n)
    ids       = d.get("_struct_conf.id", [f"H{i+1}" for i in range(n)])
    details   = d.get("_struct_conf.details", [""] * n)

    for i in range(n):
        ct_up = (ctype[i] or "").upper()
        is_helix = ("HELX" in ct_up) or ("ALPHA" in ct_up) or \
                   ("PI" in ct_up.split("_")) or ("310" in ct_up) or ("3_10" in ct_up)
        if any(bad in ct_up for bad in NON_HELIX_PREFIX):
            is_helix = False
        if not is_helix:
            continue

        bseq = _toint(beg_seq[i]); eseq = _toint(end_seq[i])
        if bseq is None or eseq is None:
            continue
        c1 = beg_chain[i] if i < len(beg_chain) else ""
        c2 = end_chain[i] if i < len(end_chain) else c1
        chain = (c1 or c2) or "?"

        bi = (beg_ins[i] if i < len(beg_ins) else "").strip()
        ei = (end_ins[i] if i < len(end_ins) else "").strip()
        hid = ids[i] if i < len(ids) else f"H{i+1}"
        det = details[i] if i < len(details) else ""

        helices.append(dict(
            id=hid, chain=chain,
            beg=bseq, end=eseq,
            beg_ins=bi, end_ins=ei,
            length=eseq - bseq + 1,
            type=ctype[i] or "HELX_P",
            details=det or "",
        ))

    helices.sort(key=lambda h: (h["chain"], h["beg"], h["id"]))
    return helices


def _fmt_res(seq, ins):
    return f"{seq}{ins}" if ins else str(seq)


def write_txt(helices, path: Path):
    per_chain = defaultdict(list)
    for h in helices:
        per_chain[h["chain"]].append(h)

    header = (f"{'Chain':<7}{'ID':<11}{'Start':>7}{'End':>7}{'Length':>8}  "
              f"{'UniProt':<10}{'Protein'}")
    L = len(header) + 40
    lines = []
    lines.append(f"Helices extracted from {Path(CIF_FILE).name}")
    lines.append(f"Total helices: {len(helices)}   "
                 f"Chains with helices: {sorted(per_chain.keys())}")
    lines.append("")
    lines.append("=" * L)
    lines.append(header)
    lines.append("-" * L)
    for h in helices:
        prot_short = (h["protein"].split(",", 1)[0][:58]) if h["protein"] else ""
        lines.append(
            f"{h['chain']:<7}{h['id']:<11}"
            f"{_fmt_res(h['beg'], h['beg_ins']):>7}"
            f"{_fmt_res(h['end'], h['end_ins']):>7}"
            f"{h['length']:>8}  {h.get('uniprot',''):<10}{prot_short}"
        )
    lines.append("")
    lines.append("Per-chain summary:")
    lines.append(f"  {'Chain':<7}{'#Helices':>10}{'Residues':>12}"
                 f"{'AvgLen':>10}  {'UniProt':<10}  Protein")
    for c in sorted(per_chain):
        hs = per_chain[c]
        nres = sum(h["length"] for h in hs)
        avg = (nres / len(hs)) if hs else 0.0
        uni = hs[0].get("uniprot", "")
        prot_short = (hs[0]["protein"].split(",", 1)[0][:50]) if hs[0].get("protein") else ""
        lines.append(f"  {c:<7}{len(hs):>10}{nres:>12}{avg:>10.1f}  {uni:<10}  {prot_short}")
    path.write_text("\n".join(lines) + "\n")


def write_tex(helices, path: Path):
    lines = []
    lines.append(r"% Helix table for 6LTJ.cif with protein name & UniProt")
    lines.append(r"% \usepackage{booktabs} required")
    lines.append(r"\begin{longtable}{@{}llrrrrllp{5.2cm}@{}}")
    lines.append(r"  \toprule")
    lines.append(r"  Chain & ID & Start & End & Length & UniProt & db\_code & Protein \\")
    lines.append(r"  \midrule")
    lines.append(r"  \endfirsthead")
    lines.append(r"  \multicolumn{8}{c}{\textit{Continued from previous page}} \\")
    lines.append(r"  \toprule")
    lines.append(r"  Chain & ID & Start & End & Length & UniProt & db\_code & Protein \\")
    lines.append(r"  \midrule")
    lines.append(r"  \endhead")
    lines.append(r"  \bottomrule")
    lines.append(r"  \endfoot")
    lines.append(r"  \bottomrule")
    lines.append(r"  \caption{Helices in PDB 6LTJ with chain-level protein names and UniProt accessions.}")
    lines.append(r"  \label{tab:6ltj_helices_uniprot}")
    lines.append(r"  \endlastfoot")
    for h in helices:
        prot = (h.get("protein") or "").split(",", 1)[0]
        prot = prot.replace("&", r"\&").replace("%", r"\%").replace("_", r"\_")
        lines.append(
            f"  {h['chain']} & {h['id']} & "
            f"{_fmt_res(h['beg'], h['beg_ins'])} & "
            f"{_fmt_res(h['end'], h['end_ins'])} & "
            f"{h['length']} & "
            f"{h.get('uniprot','')} & {h.get('db_code','')} & "
            f"{prot} \\\\"
        )
    lines.append(r"\end{longtable}")
    path.write_text("\n".join(lines) + "\n")


def write_csv(helices, path: Path):
    fields = ["chain", "entity_id", "uniprot", "db_code", "organism", "protein",
              "id", "beg", "beg_ins", "end", "end_ins", "length", "type", "details"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for h in helices:
            row = {k: h.get(k, "") for k in fields}
            w.writerow(row)


def parse_chain_proteins(d):
    """Parse per-chain protein name + UniProt info.

    Mapping pipeline:
        chain_id (asym_id / pdbx_strand_id)
          ← _entity_poly.pdbx_strand_id (comma-separated, e.g. "A,E")
        ← entity_id
        → _entity_name_com.name        (protein name, comma separated aliases)
        → _struct_ref.pdbx_db_accession (UniProt accession; db_name='UNP')
        → _struct_ref.db_code          (UniProt mnemonic id, e.g. SMCA4_HUMAN)
        → _entity_src_gen.pdbx_gene_src_scientific_name (organism)

    Falls back to _struct_ref_seq.pdbx_strand_id for chains with per-chain
    alignment rows (covers chains that share an entity_id but have their own
    ref_seq lines).  Returns rows ordered by chain_id alphabetical.
    """
    info = {}           # entity_id -> dict
    chain_to_entity = {}  # chain -> entity_id (from _entity_poly)
    chain_meta_fallback = {}  # chain -> dict override from _struct_ref_seq

    # ---- _entity_poly: entity_id -> strand list -----------------------------
    ep_ent = d.get("_entity_poly.entity_id", []) or []
    ep_strands = d.get("_entity_poly.pdbx_strand_id", []) or []
    for ei, strands in zip(ep_ent, ep_strands):
        if strands in (None, "?", ""):
            continue
        for s in strands.split(","):
            s = s.strip()
            if not s:
                continue
            chain_to_entity[s] = str(ei)
            info.setdefault(str(ei), {"entity_id": str(ei),
                                      "name": "", "uniprot": "", "db_code": "", "organism": ""})

    # ---- _entity.pdbx_description (best first-pass name for ALL entities) ---
    for ent, name in zip(d.get("_entity.id", []) or [],
                         d.get("_entity.pdbx_description", []) or []):
        row = info.setdefault(str(ent), {"entity_id": str(ent),
                                         "name": "", "uniprot": "", "db_code": "", "organism": ""})
        if name and name != "?" and not row["name"]:
            row["name"] = name if isinstance(name, str) else "".join(name)

    # ---- _entity_name_com: entity -> protein name (preferred, may be longer) -
    for ent, name in zip(d.get("_entity_name_com.entity_id", []) or [],
                         d.get("_entity_name_com.name", []) or []):
        row = info.setdefault(str(ent), {"entity_id": str(ent),
                                         "name": "", "uniprot": "", "db_code": "", "organism": ""})
        if name and name != "?":
            row["name"] = name if isinstance(name, str) else "".join(name)

    # ---- _struct_ref: entity -> UniProt -------------------------------------
    for ent, db_name, acc, code in zip(
            d.get("_struct_ref.entity_id", []) or [],
            d.get("_struct_ref.db_name", []) or [],
            d.get("_struct_ref.pdbx_db_accession", []) or [],
            d.get("_struct_ref.db_code", []) or []):
        if (db_name or "").upper() not in ("UNP", "UNIPROT", "UNIPROTKB", "SP", "SWISSPROT", "TREMBL"):
            continue
        row = info.setdefault(str(ent), {"entity_id": str(ent),
                                         "name": "", "uniprot": "", "db_code": "", "organism": ""})
        if acc and acc != "?":
            row["uniprot"] = acc
        if code and code != "?":
            row["db_code"] = code

    # ---- _struct_ref_seq: chain-level override (per-strand UniProt) ---------
    for strand, ref_id, acc in zip(
            d.get("_struct_ref_seq.pdbx_strand_id", []) or [],
            d.get("_struct_ref_seq.ref_id", []) or [],
            d.get("_struct_ref_seq.pdbx_db_accession", []) or []):
        if not strand or strand == "?":
            continue
        # Look up entity_id: either from chain_to_entity or match ref_id
        chain_meta_fallback.setdefault(strand, {})
        if acc and acc != "?":
            chain_meta_fallback[strand]["uniprot"] = acc

    # ---- _entity_src_gen: organism ------------------------------------------
    for ent, org, gene in zip(
            d.get("_entity_src_gen.entity_id", []) or [],
            d.get("_entity_src_gen.pdbx_gene_src_scientific_name", []) or [],
            d.get("_entity_src_gen.pdbx_gene_src_common_name", []) or []):
        row = info.setdefault(str(ent), {"entity_id": str(ent),
                                         "name": "", "uniprot": "", "db_code": "", "organism": ""})
        if org and org != "?":
            row["organism"] = org
        elif gene and gene != "?" and not row.get("organism"):
            row["organism"] = gene

    # ---- _entity_src_syn: handle synthetically produced entities too --------
    for ent, org in zip(
            d.get("_entity_src_syn.entity_id", []) or [],
            d.get("_entity_src_syn.organism_scientific", []) or []):
        row = info.setdefault(str(ent), {"entity_id": str(ent),
                                         "name": "", "uniprot": "", "db_code": "", "organism": ""})
        if org and org != "?" and not row.get("organism"):
            row["organism"] = org

    # ---- Build per-chain result list ----------------------------------------
    all_chains = set(chain_to_entity.keys()) | set(chain_meta_fallback.keys())
    rows = []
    for chain in sorted(all_chains):
        ent_id = chain_to_entity.get(chain, "?")
        base = dict(info.get(ent_id, {"entity_id": ent_id,
                                      "name": "", "uniprot": "", "db_code": "", "organism": ""}))
        base["chain"] = chain
        base["entity_id"] = ent_id
        for k, v in chain_meta_fallback.get(chain, {}).items():
            if v:
                base[k] = v
        rows.append(base)

    # Sort: A-Z then numeric suffixes naturally (ASCII works for PDB chains)
    rows.sort(key=lambda r: r["chain"])
    return rows


def write_chains_txt(rows, path: Path):
    header = (f"{'Chain':<7}{'Entity':>6}  {'UniProt':<12}{'db_code':<18}"
              f"{'Organism':<30}  Protein name")
    L = len(header) + 30
    lines = [
        f"Chain / protein / UniProt map for {Path(CIF_FILE).name}",
        f"Total chains: {len(rows)}",
        "", "=" * L, header, "-" * L,
    ]
    for r in rows:
        org = r["organism"] or ""
        if len(org) > 30: org = org[:28] + ".."
        nm = r["name"] or ""
        lines.append(
            f"{r['chain']:<7}{r['entity_id']:>6}  "
            f"{(r['uniprot'] or ''):<12}{(r['db_code'] or ''):<18}"
            f"{org:<30}  {nm}"
        )
    path.write_text("\n".join(lines) + "\n")


def write_chains_tex(rows, path: Path):
    lines = []
    lines.append(r"% Per-chain protein name + UniProt for 6LTJ.cif")
    lines.append(r"% \usepackage{booktabs} required")
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"  \centering\small")
    lines.append(r"  \begin{tabular}{@{}llrllp{6.5cm}@{}}")
    lines.append(r"    \toprule")
    lines.append(r"    Chain & Entity & \multicolumn{1}{c}{UniProt} & "
                 r"db\_code & Organism & Protein name \\")
    lines.append(r"    \midrule")
    for r in rows:
        nm = (r["name"] or "").replace("&", r"\&").replace("%", r"\%").replace("_", r"\_")
        org = (r["organism"] or "").replace("&", r"\&").replace("%", r"\%").replace("_", r"\_")
        lines.append(
            f"    {r['chain']} & {r['entity_id']} & "
            f"{r.get('uniprot') or ''} & "
            f"{r.get('db_code') or ''} & "
            f"{org} & {nm} \\\\"
        )
    lines.append(r"    \bottomrule")
    lines.append(r"  \end{tabular}")
    lines.append(r"  \caption{Protein chains of PDB 6LTJ with their UniProt "
                 r"accessions and source organisms.}")
    lines.append(r"  \label{tab:6ltj_chains_uniprot}")
    lines.append(r"\end{table}")
    path.write_text("\n".join(lines) + "\n")


def write_chains_csv(rows, path: Path):
    fields = ["chain", "entity_id", "uniprot", "db_code", "organism", "name"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


def main():
    d = MMCIF2Dict(CIF_FILE)

    # --- Chain/Protein info (entity → chains + protein name + UniProt) ---
    chains_info = parse_chain_proteins(d)
    chains_by_asym = {c["chain"]: c for c in chains_info}

    per_chain_order = [c["chain"] for c in chains_info]
    print(f"\n=== Protein chains (with UniProt) from {Path(CIF_FILE).name} ===")
    print(f"{'Chain':<7}{'Entity':>6}  {'UniProt':<10}{'db_code':<16}{'Organism':<28}  Protein name")
    print("-" * 120)
    for c in chains_info:
        name_short = c["name"].split(",", 1)[0] if c["name"] else c["name"]
        print(f"{c['chain']:<7}{c['entity_id']:>6}  {c['uniprot']:<10}{c['db_code']:<16}{c['organism']:<28}  {name_short}")
    print()

    # --- Helices ---
    helices = parse_helices_from_dict(d)
    if not helices:
        print("(no helix records found in _struct_conf)")
    else:
        # attach protein name + uniprot to each helix row
        for h in helices:
            c = chains_by_asym.get(h["chain"])
            h["entity_id"] = c["entity_id"] if c else "?"
            h["protein"]   = c["name"]      if c else "?"
            h["uniprot"]   = c["uniprot"]   if c else "?"
            h["db_code"]   = c["db_code"]   if c else "?"
            h["organism"]  = c["organism"]  if c else "?"

        per_chain = defaultdict(list)
        for h in helices:
            per_chain[h["chain"]].append(h)

        print(f"Found {len(helices)} helices across chains {sorted(per_chain.keys())}:")
        print()
        print(f"{'Chain':<7}{'ID':<10}{'Start':>7}{'End':>7}{'Len':>5}  "
              f"{'UniProt':<10}{'Protein (short)'}")
        print("-" * 110)
        for h in helices:
            short = (h["protein"].split(",", 1)[0][:52]) if h["protein"] else ""
            print(f"{h['chain']:<7}{h['id']:<10}"
                  f"{_fmt_res(h['beg'], h['beg_ins']):>7}"
                  f"{_fmt_res(h['end'], h['end_ins']):>7}{h['length']:>5}  "
                  f"{h['uniprot']:<10}{short}")
        print()
        for c in per_chain_order:
            if c not in per_chain:
                continue
            hs = per_chain[c]
            nres = sum(h["length"] for h in hs)
            avg = (nres / len(hs)) if hs else 0.0
            name = chains_by_asym[c]["name"].split(",", 1)[0][:60] if chains_by_asym.get(c) else ""
            print(f"  Chain {c}: {len(hs):>3} helices, {nres:>4} res, avg {avg:>5.1f} aa  | {name}")

    # --- Writers ---
    if helices:
        write_txt(helices, Path(OUT_TXT))
        write_tex(helices, Path(OUT_TEX))
        write_csv(helices, Path(OUT_CSV))
        print()
        print(f"Wrote helix TXT → {OUT_TXT}")
        print(f"Wrote helix TEX → {OUT_TEX}")
        print(f"Wrote helix CSV → {OUT_CSV}")

    write_chains_txt(chains_info, Path(OUT_CHAINS_TXT))
    write_chains_tex(chains_info, Path(OUT_CHAINS_TEX))
    write_chains_csv(chains_info, Path(OUT_CHAINS_CSV))
    print(f"Wrote chains TXT → {OUT_CHAINS_TXT}")
    print(f"Wrote chains TEX → {OUT_CHAINS_TEX}")
    print(f"Wrote chains CSV → {OUT_CHAINS_CSV}")


if __name__ == "__main__":
    main()
