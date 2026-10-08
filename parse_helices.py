#!/usr/bin/env python3
"""Parse helices from ANY mmCIF file and produce TWO table variants.

Variant 1 (simple):  Protein name + Helix Start + Helix End
Variant 2 (detailed): Protein name + Helix Start + Helix End + Helix Length
                      + footer: per-subunit total helix length

Outputs (txt + inline LaTeX + csv):
    {prefix}_helices_v1.txt / .tex       # protein name | start | end
    {prefix}_helices_v2.txt / .tex       # protein name | start | end | length + per-subunit totals
    {prefix}_helices.csv                 # combined raw data (csv, chain/entity/uniprot/protein/beg/end/length/type)
    {prefix}_protein_aa_summary.txt / .tex / .csv
"""

from pathlib import Path
import argparse
import csv
from collections import defaultdict, OrderedDict
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

NON_HELIX_PREFIX = ("SHEET", "STRN", "TURN")


def _first_list(d, keys, n):
    for k in keys:
        v = d.get(k)
        if isinstance(v, list) and len(v) >= n:
            return v
        if isinstance(v, list) and len(v) > 0:
            return v
    return [""] * n


def _toint(x):
    try:  return int(x)
    except Exception:  return None


def _fmt_res(seq, ins):
    return f"{seq}{ins}" if ins else str(seq)


# ---------------------------------------------------------------------------
# Parse helix & chain/protein metadata
# ---------------------------------------------------------------------------
def parse_mmcif(cif_path):
    return MMCIF2Dict(cif_path)


def parse_helices(d):
    ctype = d.get("_struct_conf.conf_type_id", [])
    helices = []
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
    conf_type = d.get("_struct_conf.conf_type_id", [""] * n)

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
        helices.append(dict(
            id       = ids[i] if i < len(ids) else f"H{i+1}",
            chain    = chain,
            beg      = bseq, end = eseq,
            beg_ins  = bi,   end_ins = ei,
            length   = eseq - bseq + 1,
            type     = conf_type[i] or "HELX_P",
        ))

    helices.sort(key=lambda h: (h["chain"], h["beg"], h["id"]))
    return helices


def parse_chain_proteins(d):
    """chain → entity_id, protein name, uniprot, db_code, organism, total_aa."""
    info = {}
    chain_to_entity = {}

    ep_ent     = d.get("_entity_poly.entity_id", []) or []
    ep_strands = d.get("_entity_poly.pdbx_strand_id", []) or []
    ep_seq     = d.get("_entity_poly.pdbx_seq_one_letter_code_can",
                       d.get("_entity_poly.pdbx_seq_one_letter_code", [])) or []
    # Entity-level total aa (length of canonical sequence or entity_poly_seq)
    epoly_seq_ent = d.get("_entity_poly_seq.entity_id", []) or []
    epoly_seq_num = d.get("_entity_poly_seq.num", []) or []
    entity_seq_len_from_poly = {}
    if epoly_seq_ent and epoly_seq_num:
        mx = {}
        for eid, num in zip(epoly_seq_ent, epoly_seq_num):
            try:
                n = int(num)
                mx[eid] = max(mx.get(eid, 0), n)
            except Exception:
                pass
        entity_seq_len_from_poly = mx
    entity_seq_len_from_seq = {}
    for eid, seq in zip(ep_ent, ep_seq):
        if seq and seq != "?":
            clean = "".join(c for c in seq if c.isalpha())
            entity_seq_len_from_seq[str(eid)] = len(clean)

    for ei, strands in zip(ep_ent, ep_strands):
        if strands in (None, "?", ""):
            continue
        for s in strands.split(","):
            s = s.strip()
            if not s:
                continue
            chain_to_entity[s] = str(ei)
            info.setdefault(str(ei), {"entity_id": str(ei), "name": "",
                                      "uniprot": "", "db_code": "", "organism": "",
                                      "total_aa": 0})
    # Prefer the entity-level sequence length (pdbx_seq_one_letter_code_can) then fallback entity_poly_seq max num
    for eid, row in info.items():
        row["total_aa"] = entity_seq_len_from_seq.get(eid, entity_seq_len_from_poly.get(eid, 0))

    for ent, name in zip(d.get("_entity.id", []) or [],
                         d.get("_entity.pdbx_description", []) or []):
        row = info.setdefault(str(ent), {"entity_id": str(ent), "name": "",
                                          "uniprot": "", "db_code": "", "organism": ""})
        if name and name != "?" and not row["name"]:
            row["name"] = name if isinstance(name, str) else "".join(name)

    for ent, name in zip(d.get("_entity_name_com.entity_id", []) or [],
                         d.get("_entity_name_com.name", []) or []):
        row = info.setdefault(str(ent), {"entity_id": str(ent), "name": "",
                                          "uniprot": "", "db_code": "", "organism": ""})
        if name and name != "?":
            row["name"] = name if isinstance(name, str) else "".join(name)

    for ent, db_name, acc, code in zip(
            d.get("_struct_ref.entity_id", []) or [],
            d.get("_struct_ref.db_name", []) or [],
            d.get("_struct_ref.pdbx_db_accession", []) or [],
            d.get("_struct_ref.db_code", []) or []):
        if (db_name or "").upper() not in ("UNP","UNIPROT","UNIPROTKB","SP","SWISSPROT","TREMBL"):
            continue
        row = info.setdefault(str(ent), {"entity_id": str(ent), "name": "",
                                          "uniprot": "", "db_code": "", "organism": ""})
        if acc and acc != "?":  row["uniprot"] = acc
        if code and code != "?": row["db_code"] = code

    for ent, org, gene in zip(
            d.get("_entity_src_gen.entity_id", []) or [],
            d.get("_entity_src_gen.pdbx_gene_src_scientific_name", []) or [],
            d.get("_entity_src_gen.pdbx_gene_src_common_name", []) or []):
        row = info.setdefault(str(ent), {"entity_id": str(ent), "name": "",
                                          "uniprot": "", "db_code": "", "organism": ""})
        if org and org != "?":          row["organism"] = org
        elif gene and gene != "?" and not row.get("organism"):
            row["organism"] = gene

    for strand, acc in zip(
            d.get("_struct_ref_seq.pdbx_strand_id", []) or [],
            d.get("_struct_ref_seq.pdbx_db_accession", []) or []):
        if not strand or strand == "?":
            continue
        eid = chain_to_entity.get(strand)
        if eid and acc and acc != "?":
            info[eid]["uniprot"] = info[eid].get("uniprot") or acc

    all_chains = sorted(set(chain_to_entity.keys()))
    rows = []
    for ch in all_chains:
        eid = chain_to_entity.get(ch, "?")
        base = dict(info.get(eid, {"entity_id": eid, "name": "", "uniprot": "",
                                    "db_code": "", "organism": "", "total_aa": 0}))
        base["chain"] = ch
        base["entity_id"] = eid
        base["total_aa"] = base.get("total_aa", 0) or 0
        rows.append(base)
    return rows


def annotate_helices(helices, chain_info):
    by_chain = {c["chain"]: c for c in chain_info}
    for h in helices:
        c = by_chain.get(h["chain"], {})
        h["entity_id"] = c.get("entity_id", "?")
        h["protein"]   = c.get("name") or "(unnamed)"
        h["uniprot"]   = c.get("uniprot") or ""
        h["db_code"]   = c.get("db_code") or ""
        h["organism"]  = c.get("organism") or ""
        h["total_aa"]  = c.get("total_aa") or 0
    return helices


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------
def _short_protein(name):
    """Take first comma-separated alias, trim to a sensible max length."""
    n = (name or "")
    n = n.split(",", 1)[0].strip()
    return n


def _short_protein_code(name):
    """Return the short protein code (e.g. 'SWP73', 'SNF2', 'ARP7', 'TYE2', 'GAM3').

    Strategy (ordered, longest match wins so 'SWP73' beats 'SWI'):
      1) Prefer tokens of the form 'letters+digits' where letters are uppercase
         and the numeric part has >=1 digit (SWP73, SNF6, SNF5, SNF2, ARP7, ARP9,
         SWP82, TYE2, GAM3).  Pick the LONGEST match (ties: alphabetical last =>
         prefer code-style over generic 'SWI'-style tokens).
      2) Otherwise use longest all-uppercase 3-10 letter identifier.
      3) Fallback: last titlecase-ish token of length 3-12; or truncate alias0.
    """
    alias0 = (name or "").split(",", 1)[0]
    tokens = alias0.replace("-", " ").replace("/", " ").replace("(", " ").replace(")", " ").split()
    import re
    code_re = re.compile(r"^([A-Z]{2,})(\d+)$")
    letter_digit = []
    for t in tokens:
        m = code_re.match(t)
        if m and m.group(1) and m.group(2):
            letter_digit.append(t)
    if letter_digit:
        letter_digit.sort(key=lambda x: (len(x), x), reverse=True)
        return letter_digit[0]
    # --- 处理 'Regulator of Ty1 transposition protein 102' → 'RTT102' 风格的 fallback ---
    # 提取所有 "首字母大写字母串 + 数字" 的组合（跨 token，如 "Ty1" + "protein 102" → TYP102 / Ty102...），
    # 以及所有独立 token 末尾带数字的 (Ty1 / GAM3 / ... )
    mixed_digit = []
    for t in tokens:
        if len(t) >= 3 and t[0].isalpha() and (t[-1].isdigit() or any(c.isdigit() for c in t)):
            mixed_digit.append(t)
    if mixed_digit:
        # 选取含 digit 最多 + 最长的
        mixed_digit.sort(key=lambda x: (sum(c.isdigit() for c in x), len(x), x), reverse=True)
        t0 = mixed_digit[0]
        # 把前面所有字母字符拼起来 + 保留末尾数字，全大写（Ty1 → TY1 / TYP102 样式）
        head = "".join(c for c in t0 if c.isalpha())
        tail = "".join(c for c in t0 if c.isdigit())
        if head and tail:
            return (head + tail).upper()
    all_upper = [t for t in tokens if t.isalpha() and t.isupper() and 3 <= len(t) <= 10]
    if all_upper:
        all_upper.sort(key=len, reverse=True)
        return all_upper[0]
    short = [t for t in tokens if 3 <= len(t) <= 12 and t[0].isupper()]
    if short:
        return short[-1]
    return _short_protein(name)[:12]


def write_v1_txt(helices, path: Path, prefix: str=""):
    header = f"{'Protein code':<16}{'Helix start':>14}{'Helix end':>12}"
    lines = [
        f"{prefix.upper()} helices — variant 1 (protein CODE only: SWP73/SNF2 style).  Total helices: {len(helices)}",
        "", "=" * len(header), header, "-" * len(header),
    ]
    for h in helices:
        code = _short_protein_code(h["protein"])
        lines.append(
            f"{code:<16}"
            f"{_fmt_res(h['beg'], h['beg_ins']):>14}"
            f"{_fmt_res(h['end'], h['end_ins']):>12}"
        )
    path.write_text("\n".join(lines) + "\n")


def write_v1_tex(helices, path: Path, prefix: str=""):
    lines = [
        rf"% Variant 1: short protein CODE (SWP73 / SNF2 / ... style) + start + end ({prefix.lower()})",
        r"% Requires: \usepackage{booktabs,longtable}",
        r"\begin{longtable}{@{}lrr@{}}",
        r"  \toprule",
        r"  Protein code & Helix start & Helix end \\",
        r"  \midrule",
        r"  \endfirsthead",
        r"  \multicolumn{3}{c}{\textit{Continued from previous page}} \\",
        r"  \toprule Protein code & Helix start & Helix end \\ \midrule",
        r"  \endhead",
        r"  \bottomrule",
        r"  \endfoot",
        r"  \bottomrule",
        rf"  \caption{{Helices of PDB {prefix.upper()} using short protein codes (e.g. SWP73, SNF2, ARP7).}}",
        rf"  \label{{tab:{prefix.lower()}_helices_v1}}",
        r"  \endlastfoot",
    ]
    for h in helices:
        code = _short_protein_code(h["protein"])
        code = code.replace("&", r"\&").replace("%", r"\%").replace("_", r"\_")
        lines.append(
            f"  {code} & {_fmt_res(h['beg'], h['beg_ins'])} & "
            f"{_fmt_res(h['end'], h['end_ins'])} \\\\"
        )
    lines.append(r"\end{longtable}")
    path.write_text("\n".join(lines) + "\n")


def write_v2_txt(helices, path: Path, prefix: str=""):
    subunit_len = defaultdict(int)
    subunit_cnt = defaultdict(int)
    subunit_aa  = {}
    ordered_keys = []
    for h in helices:
        k = (h["chain"], h["protein"])
        if k not in subunit_len:
            ordered_keys.append(k)
            subunit_aa[k] = h.get("total_aa") or 0
        subunit_len[k] += h["length"]
        subunit_cnt[k] += 1

    header = (f"{'Subunit (Chain)':<40}{'Helix start':>13}"
              f"{'Helix end':>12}{'Length':>9}")
    lines = [
        f"{prefix.upper()} helices — variant 2 (subunit name + chain | start | end | length).  "
        f"Total helices: {len(helices)}",
        "", "=" * len(header), header, "-" * len(header),
    ]
    for h in helices:
        name = _short_protein(h["protein"])
        label = f"{name} ({h['chain']})"
        lines.append(
            f"{label:<40}"
            f"{_fmt_res(h['beg'], h['beg_ins']):>13}"
            f"{_fmt_res(h['end'], h['end_ins']):>12}"
            f"{h['length']:>9}"
        )

    lines.append("")
    lines.append("=" * len(header))
    lines.append("Per-subunit summary:")
    lines.append(
        f"  {'Subunit (Chain)':<42}"
        f"{'#hel':>6}{'HelLen':>10}{'TotalAA':>10}{'%Hel':>8}"
    )
    lines.append("  " + "-" * 78)
    grand_hel = 0; grand_cnt = 0; grand_aa = 0
    for (ch, prot) in ordered_keys:
        cnt = subunit_cnt[(ch, prot)]
        hel_len = subunit_len[(ch, prot)]
        total_aa = subunit_aa.get((ch, prot), 0)
        pct = (hel_len / total_aa * 100.0) if total_aa else float("nan")
        name = _short_protein(prot)
        label = f"{name} ({ch})"
        pct_s = f"{pct:>7.1f}%" if pct == pct else "     n/a"   # nan check
        lines.append(
            f"  {label:<42}{cnt:>6}{hel_len:>10}{total_aa:>10}{pct_s}"
        )
        grand_cnt += cnt; grand_hel += hel_len; grand_aa += total_aa
    lines.append("  " + "-" * 78)
    grand_pct = (grand_hel / grand_aa * 100.0) if grand_aa else float("nan")
    grand_pct_s = f"{grand_pct:>7.1f}%" if grand_pct == grand_pct else "     n/a"
    lines.append(
        f"  {'GRAND TOTAL':<42}"
        f"{grand_cnt:>6}{grand_hel:>10}{grand_aa:>10}{grand_pct_s}"
    )
    path.write_text("\n".join(lines) + "\n")


def write_v2_tex(helices, path: Path, prefix: str=""):
    subunit_len = defaultdict(int)
    subunit_cnt = defaultdict(int)
    subunit_aa  = {}
    ordered_keys = []
    for h in helices:
        k = (h["chain"], h["protein"])
        if k not in subunit_len:
            ordered_keys.append(k)
            subunit_aa[k] = h.get("total_aa") or 0
        subunit_len[k] += h["length"]
        subunit_cnt[k] += 1

    lines = [
        r"% Variant 2: subunit (chain) + start + end + length + per-subunit total length / total aa / % helical",
        r"% Requires: \usepackage{booktabs,longtable}",
        r"\begin{longtable}{@{}p{7.5cm}rrrl@{}}",
        r"  \toprule",
        r"  Subunit (Chain) & Helix start & Helix end & \multicolumn{1}{c}{Length} \\",
        r"  \midrule",
        r"  \endfirsthead",
        r"  \multicolumn{5}{c}{\textit{Continued from previous page}} \\",
        r"  \toprule Subunit (Chain) & Helix start & Helix end & \multicolumn{1}{c}{Length} \\ \midrule",
        r"  \endhead",
        r"  \bottomrule",
        r"  \endfoot",
        r"  \bottomrule",
        rf"  \caption{{Helices of PDB {prefix.upper()} with per-subunit total helix length, total amino acids, and \% helical content.}}",
        rf"  \label{{tab:{prefix.lower()}_helices_v2}}",
        r"  \endlastfoot",
    ]
    for h in helices:
        name = _short_protein(h["protein"])
        name = name.replace("&", r"\&").replace("%", r"\%").replace("_", r"\_")
        chain = h["chain"].replace("&", r"\&")
        lines.append(
            f"  {name} ({chain}) & "
            f"{_fmt_res(h['beg'], h['beg_ins'])} & "
            f"{_fmt_res(h['end'], h['end_ins'])} & "
            f"\\multicolumn{{1}}{{c}}{{{h['length']}}} \\\\"
        )

    # Footer: per-subunit totals + grand total
    lines.append(r"  \midrule")
    lines.append(r"  \multicolumn{5}{@{}l}{\textbf{Per-subunit summary:}} \\")
    lines.append(r"  \cmidrule(r){1-5}")
    lines.append(
        r"  \textbf{Subunit (Chain)} & \multicolumn{1}{r}{\textbf{\#hel}} & "
        r"\multicolumn{1}{r}{\textbf{HelLen}} & "
        r"\multicolumn{1}{r}{\textbf{TotalAA}} & "
        r"\multicolumn{1}{c}{\textbf{\%Hel}} \\"
    )
    lines.append(r"  \cmidrule(r){1-5}")
    grand_hel = 0; grand_cnt = 0; grand_aa = 0
    for (ch, prot) in ordered_keys:
        cnt = subunit_cnt[(ch, prot)]
        hel_len = subunit_len[(ch, prot)]
        total_aa = subunit_aa.get((ch, prot), 0)
        pct = (hel_len / total_aa * 100.0) if total_aa else float("nan")
        name = _short_protein(prot)
        name = name.replace("&", r"\&").replace("%", r"\%").replace("_", r"\_")
        chain = ch.replace("&", r"\&")
        pct_s = f"{pct:.1f}\%" if pct == pct else "n/a"
        lines.append(
            f"  {name} ({chain}) & "
            f"\\multicolumn{{1}}{{r}}{{{cnt}}} & "
            f"\\multicolumn{{1}}{{r}}{{{hel_len}}} & "
            f"\\multicolumn{{1}}{{r}}{{{total_aa}}} & "
            f"\\multicolumn{{1}}{{c}}{{{pct_s}}} \\\\"
        )
        grand_cnt += cnt; grand_hel += hel_len; grand_aa += total_aa
    lines.append(r"  \cmidrule(r){1-5}")
    grand_pct = (grand_hel / grand_aa * 100.0) if grand_aa else float("nan")
    grand_pct_s = f"{grand_pct:.1f}\%" if grand_pct == grand_pct else "n/a"
    lines.append(
        f"  \\textbf{{GRAND TOTAL}} & "
        f"\\multicolumn{{1}}{{r}}{{{grand_cnt}}} & "
        f"\\multicolumn{{1}}{{r}}{{{grand_hel}}} & "
        f"\\multicolumn{{1}}{{r}}{{{grand_aa}}} & "
        f"\\multicolumn{{1}}{{c}}{{{grand_pct_s}}} \\\\"
    )

    lines.append(r"\end{longtable}")
    path.write_text("\n".join(lines) + "\n")


def write_csv(helices, path: Path):
    fields = ["chain", "entity_id", "uniprot", "db_code", "organism", "protein",
              "id", "beg", "beg_ins", "end", "end_ins", "length", "type", "total_aa"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for h in helices:
            row = {k: h.get(k, "") for k in fields}
            w.writerow(row)


def build_aa_summary(chain_info):
    """Return deduplicated list of {code, name, chains_str, total_aa} per protein.

    Dedup rule: two chains are "same protein" if they have the same entity_id
    OR same protein name AND same total_aa length.  Each unique protein gets
    its short protein code (from _short_protein_code) plus a comma-joined list
    of its chains (e.g. "B,D").
    """
    from collections import OrderedDict as _OD
    keyed = _OD()                 # unique key -> row dict
    for c in chain_info:
        name = c.get("name") or ""
        ent  = c.get("entity_id") or "?"
        aa   = c.get("total_aa") or 0
        chain = c.get("chain")
        code = _short_protein_code(name) if name else chain
        # Fallback for very short/unrecognizable names (e.g. "Unkown", "Unkn"):
        #   - strip common typo variants "unkown" / "unknown" / "Unkn" / "?", etc.
        #   - if name contains a real UniProt ID (e.g. chain G = '7C4J' accession style),
        #     keep it as is; else collapse all UNK variants to UNK.
        def _looks_like_unknown(s):
            if not s: return True
            sl = s.lower()
            if s in ("?", "UNKNOWN", "UNK", "UNKN", "Unkn"): return True
            if sl in ("unknown", "unkown", "unk", "unkn"): return True
            if sl.startswith("unkn") and len(sl) <= 7: return True
            if sl.startswith("unk") and len(sl) <= 5: return True
            return False
        if _looks_like_unknown(code) or _looks_like_unknown(name):
            # preserve a real-looking uniprot style short id like 7C4J only if name == code
            if code and code.isalnum() and 3 <= len(code) <= 6 and not code.upper().startswith("UNK"):
                pass
            else:
                code = f"UNK{chain}" if chain else "UNK"
        if not code or len(code) < 3 or not any(ch.isalpha() for ch in code):
            letters = "".join(ch for ch in (name or chain or "") if ch.isalpha())
            if letters:
                code = letters[:6].upper()
            else:
                code = f"UNK{chain}" if chain else "UNK"
        # Uniqueness key: prefer entity_id; fall back (name, length) for robustness
        key = ("e", ent) if ent != "?" else ("n", name, aa)
        if key not in keyed:
            keyed[key] = dict(code=code, name=name, chains=[], total_aa=aa,
                              uniprot=c.get("uniprot",""), db_code=c.get("db_code",""),
                              organism=c.get("organism",""))
        if chain and chain not in keyed[key]["chains"]:
            keyed[key]["chains"].append(chain)
    rows = []
    for r in keyed.values():
        r2 = dict(r)
        r2["chains"] = ",".join(r2["chains"])
        rows.append(r2)
    # Sort by protein code then by (name length desc => longest first)
    rows.sort(key=lambda r: (r["code"], -len(r["name"])))
    return rows


def write_aa_txt(rows, path: Path, prefix: str=""):
    header = (f"{'Code':<10}{'Chains':<10}{'TotalAA':>9}  "
              f"{'UniProt':<12}{'db_code':<18}Protein name")
    lines = [
        f"{prefix.upper()} per-protein total amino-acid summary (unique protein codes).",
        f"Unique proteins: {len(rows)}",
        "",
        "=" * len(header),
        header,
        "-" * len(header),
    ]
    grand = 0
    for r in rows:
        lines.append(
            f"{r['code']:<10}{r['chains']:<10}{r['total_aa']:>9}  "
            f"{r.get('uniprot',''):<12}{r.get('db_code',''):<18}"
            f"{_short_protein(r['name'])}"
        )
        grand += r["total_aa"]
    lines.append("-" * len(header))
    lines.append(f"{'TOTAL':<10}{'':<10}{grand:>9}")
    path.write_text("\n".join(lines) + "\n")


def write_aa_tex(rows, path: Path, prefix: str=""):
    lines = [
        rf"% {prefix.upper()} unique protein code => total amino acids summary.",
        r"% Requires: \usepackage{booktabs}",
        r"\begin{table}[htbp]",
        r"  \centering\small",
        r"  \begin{tabular}{@{}llrllp{6.8cm}@{}}",
        r"    \toprule",
        r"    Code & Chains & \multicolumn{1}{c}{Total AA} & UniProt & db\_code & Protein name \\",
        r"    \midrule",
    ]
    grand = 0
    for r in rows:
        nm = _short_protein(r["name"]).replace("&", r"\&").replace("%", r"\%").replace("_", r"\_")
        lines.append(
            f"    {r['code']} & {r['chains']} & "
            f"\\multicolumn{{1}}{{c}}{{{r['total_aa']}}} & "
            f"{r.get('uniprot','')} & {r.get('db_code','')} & {nm} \\\\"
        )
        grand += r["total_aa"]
    lines.append(r"    \midrule")
    lines.append(
        f"    \\textbf{{TOTAL}} & & "
        f"\\multicolumn{{1}}{{c}}{{\\textbf{{{grand}}}}} & & & \\\\"
    )
    lines.append(r"    \bottomrule")
    lines.append(r"  \end{tabular}")
    lines.append(rf"  \caption{{Unique protein chains of PDB {prefix.upper()} with short protein codes "
                 r"and total number of amino acids. Chains sharing an entity (e.g. B \& D = TYE2) "
                 r"are collapsed into one row.}}")
    lines.append(rf"  \label{{tab:{prefix.lower()}_protein_aa}}")
    lines.append(r"\end{table}")
    path.write_text("\n".join(lines) + "\n")


def write_aa_csv(rows, path: Path):
    fields = ["code", "chains", "total_aa", "uniprot", "db_code", "organism", "name"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


# ---------------------------------------------------------------------------
def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Parse helices from any mmCIF file and produce summary tables (txt/LaTeX/csv)."
    )
    parser.add_argument(
        "cif",
        type=Path,
        help="Path to the input mmCIF file (e.g. structure.cif)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory for generated files. Defaults to current working directory.",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default=None,
        help="Prefix used for output filenames and headers (e.g. 7C4J). Defaults to cif file stem.",
    )
    return parser


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    cif_path: Path = args.cif
    if not cif_path.exists():
        raise SystemExit(f"Error: CIF file not found: {cif_path}")

    prefix: str = args.prefix or cif_path.stem
    out_dir: Path = args.out_dir or Path.cwd()
    out_dir.mkdir(parents=True, exist_ok=True)

    OUT_V1_TXT = out_dir / f"{prefix}_helices_v1.txt"
    OUT_V1_TEX = out_dir / f"{prefix}_helices_v1.tex"
    OUT_V2_TXT = out_dir / f"{prefix}_helices_v2.txt"
    OUT_V2_TEX = out_dir / f"{prefix}_helices_v2.tex"
    OUT_CSV    = out_dir / f"{prefix}_helices.csv"
    OUT_AA_TXT = out_dir / f"{prefix}_protein_aa_summary.txt"
    OUT_AA_TEX = out_dir / f"{prefix}_protein_aa_summary.tex"
    OUT_AA_CSV = out_dir / f"{prefix}_protein_aa_summary.csv"

    d = parse_mmcif(str(cif_path))
    chain_info = parse_chain_proteins(d)
    helices    = annotate_helices(parse_helices(d), chain_info)

    # Only keep helices whose chain was matched to a protein (skip DNA/ion if any)
    helices_prot = [h for h in helices if h.get("protein") and h["protein"] != "(unnamed)"]
    print(f"Parsed {len(helices_prot)} protein helices across "
          f"{len({h['chain'] for h in helices_prot})} chains")

    # --- Screen summary ---
    subunit_len = defaultdict(int); subunit_cnt = defaultdict(int); keys_ordered = []
    for h in helices_prot:
        k = (h["chain"], h["protein"])
        if k not in subunit_len: keys_ordered.append(k)
        subunit_len[k] += h["length"]; subunit_cnt[k] += 1
    print("\nPer-subunit summary:")
    grand = 0
    for (ch, prot) in keys_ordered:
        print(f"  Chain {ch}: {_short_protein(prot):<40s}  "
              f"#helices={subunit_cnt[(ch,prot)]:>3}  total_len={subunit_len[(ch,prot)]:>4}")
        grand += subunit_len[(ch, prot)]
    print(f"  {'GRAND TOTAL':<55}  {grand} aa total in helices")

    write_v1_txt(helices_prot, Path(OUT_V1_TXT), prefix=prefix)
    write_v1_tex(helices_prot, Path(OUT_V1_TEX), prefix=prefix)
    print(f"\nWrote v1 (name|start|end)         → {OUT_V1_TXT}  / {OUT_V1_TEX}")

    write_v2_txt(helices_prot, Path(OUT_V2_TXT), prefix=prefix)
    write_v2_tex(helices_prot, Path(OUT_V2_TEX), prefix=prefix)
    print(f"Wrote v2 (name|start|end|length)  → {OUT_V2_TXT}  / {OUT_V2_TEX}")

    write_csv(helices_prot, Path(OUT_CSV))
    print(f"\nWrote CSV (incl. total_aa)      → {OUT_CSV}")

    # --- NEW: Unique protein code => total amino acids table ---
    aa_rows = build_aa_summary(chain_info)
    write_aa_txt(aa_rows, Path(OUT_AA_TXT), prefix=prefix)
    write_aa_tex(aa_rows, Path(OUT_AA_TEX), prefix=prefix)
    write_aa_csv(aa_rows, Path(OUT_AA_CSV))
    print(f"\nWrote protein code / TotalAA table:")
    print(f"  TXT  → {OUT_AA_TXT}")
    print(f"  TEX  → {OUT_AA_TEX}")
    print(f"  CSV  → {OUT_AA_CSV}")
    print(f"\nProtein code summary (inline):")
    print(f"  {'Code':<10}{'Chains':<10}{'TotalAA':>9}  Protein")
    grand = 0
    for r in aa_rows:
        print(f"  {r['code']:<10}{r['chains']:<10}{r['total_aa']:>9}  {_short_protein(r['name'])}")
        grand += r["total_aa"]
    print(f"  {'TOTAL':<10}{'':<10}{grand:>9}")


if __name__ == "__main__":
    main()
