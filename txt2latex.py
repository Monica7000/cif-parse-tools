#!/usr/bin/env python3
"""Convert extracted_atoms.txt → inline LaTeX table (no Chain column).

Usage:
    python txt2latex.py extracted_atoms.txt
    python txt2latex.py extracted_atoms.txt -o custom.tex
"""

import argparse
import sys
from pathlib import Path

HEADER_KEYS = {"Chain", "ResID", "ResName", "AtomName", "X", "Y", "Z"}


def parse_txt(path: Path):
    """Return (target_label, rows) where each row is (resid, resname, atom, x, y, z)."""
    rows = []
    target_resid = None
    first_resname = None
    with path.open() as f:
        for line in f:
            line = line.rstrip("\n")
            parts = line.split()
            if not parts:
                continue
            if parts[0] in HEADER_KEYS:
                continue
            if len(parts) < 7:
                continue
            chain, resid, resname, atom = parts[0], parts[1], parts[2], parts[3]
            try:
                x, y, z = float(parts[4]), float(parts[5]), float(parts[6])
            except ValueError:
                continue
            if target_resid is None:
                target_resid = resid
                first_resname = resname
            rows.append((resid, resname, atom, x, y, z))
    target_label = (target_resid, first_resname) if target_resid else None
    return target_label, rows


def build_latex(target_label, rows, out_path: Path):
    target_resid, target_resname = target_label or (None, None)

    # Insert \cmidrule after the target block (first contiguous ResID).
    split_idx = len(rows)
    for i, (resid, _rn, _a, _x, _y, _z) in enumerate(rows):
        if resid != target_resid:
            split_idx = i
            break

    lines = []
    lines.append(r"% Inline LaTeX table — copy this snippet directly into your LaTeX document.")
    lines.append(r"% Requires in preamble:")
    lines.append(r"%   \usepackage{booktabs}")
    lines.append(r"%   \usepackage{graphicx}   % for \resizebox")
    lines.append("")
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\resizebox{0.8\textwidth}{!}{%")
    lines.append(r"\begin{tabular}{@{}llcrrr@{}}")
    lines.append(r"\toprule")
    lines.append(
        r"\textbf{ResID} & \textbf{ResName} & \textbf{AtomName} & "
        r"\textbf{X (\AA)} & \textbf{Y (\AA)} & \textbf{Z (\AA)} \\"
    )
    lines.append(r"\midrule")

    def fmt_row(r):
        resid, resname, atom, x, y, z = r
        return f"{resid:<4s} & {resname:<4s} & {atom:<5s} & {x:7.3f} & {y:7.3f} & {z:7.3f} \\\\"

    for r in rows[:split_idx]:
        lines.append(fmt_row(r))

    if 0 < split_idx < len(rows):
        lines.append(r"\cmidrule(lr){1-6}")

    for r in rows[split_idx:]:
        lines.append(fmt_row(r))

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}%")
    lines.append(r"}")
    caption_parts = []
    if target_label:
        caption_parts.append(
            f"Target residue {target_resid} ({target_resname}): {split_idx} atoms"
        )
    n_sur = len(rows) - split_idx
    if n_sur:
        caption_parts.append(f"surrounding: {n_sur} atoms")
    caption = (
        "; ".join(caption_parts)
        + ". Coordinates in \\AA ngstr\\\"oms."
    )
    lines.append(f"\\caption{{{caption}}}")
    label = (
        f"tab:atoms_{target_resid}_{target_resname}"
        if target_label
        else "tab:atoms"
    )
    lines.append(f"\\label{{{label}}}")
    lines.append(r"\end{table}")
    lines.append("")

    out_path.write_text("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(description="Convert extracted_atoms.txt to inline LaTeX table (no Chain column).")
    parser.add_argument("input", type=Path, help="Path to the input extracted_atoms.txt file")
    parser.add_argument("-o", "--output", type=Path, default=None, help="Output .tex file path (default: replace input suffix with .tex)")
    args = parser.parse_args()

    in_path: Path = args.input
    out_path: Path = args.output if args.output is not None else in_path.with_suffix(".tex")

    if not in_path.exists():
        print(f"ERROR: input file not found: {in_path}", file=sys.stderr)
        sys.exit(1)

    target_label, rows = parse_txt(in_path)
    if not rows:
        print("ERROR: no data rows parsed from input.", file=sys.stderr)
        sys.exit(1)

    build_latex(target_label, rows, out_path)
    print(f"Wrote {len(rows)} rows → {out_path}")
    n_target = sum(1 for r in rows if r[0] == (target_label[0] if target_label else None))
    n_sur = len(rows) - n_target
    print(f"  Target    : {n_target} atoms  (ResID {target_label[0] if target_label else '?'}, "
          f"{target_label[1] if target_label else '?'})")
    print(f"  Surrounding: {n_sur} atoms")


if __name__ == "__main__":
    main()
