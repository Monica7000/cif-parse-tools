# mmCIF / PDB Structure Parsing Scripts

A small collection of Python scripts that extract structural details from **mmCIF (`.cif`)**
files downloaded from the PDB or AlphaFold DB — helices, ligands, residue-level
coordinates, chain → protein name / UniProt maps, and H-bond-like N↔O contacts.

All outputs are written in **three formats** so you can plug them straight into
analysis (`.csv`), a manuscript (`.tex` — inline LaTeX with `booktabs`/`longtable`),
or quick interactive review (`.txt`).

---

## Quick start

```bash
# 1. Clone & install
git clone <your-repo-url>
cd cif-parse-tools
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt

# 2. Parse helices from any .cif
python parse_helices_2.py ~/Downloads/6LTJ.cif --out-dir out/ --prefix 6LTJ

# 3. Parse ligands
python parse_ligands.py ~/Downloads/4LXZ.cif --out-dir out/

# 4. Extract atoms within N Å of a target residue + compute N↔O hetero-distances
python extract_residues.py ~/Downloads/4LXZ.cif \
       --chain A --resid 301 --cutoff 3.5 --out-dir out/

# 5. Turn the extracted-atom TXT into a LaTeX table
python txt2latex.py out/4LXZ_extracted_atoms.txt -o out/4LXZ_extracted_atoms.tex
```

---

## Requirements

| Package   | Version | Why                                                    |
|-----------|---------|--------------------------------------------------------|
| biopython | 1.85    | `Bio.PDB.MMCIF2Dict` — fast low-level mmCIF key→list parser |
| numpy     | 1.26.4  | vectorised distance matrix in `extract_residues.py`    |

Install with:

```bash
pip install -r requirements.txt
```

Python 3.9+ recommended (tested on 3.9).

---

## Scripts

Every script shares the same CLI shape where it makes sense:

```
script.py INPUT.cif [--out-dir DIR] [--prefix STR] [extra options…]
```

- `--out-dir` — where to write `.txt` / `.tex` / `.csv` (default: the CIF's own
  directory, or current working dir). Always auto-created if missing.
- `--prefix` — filename prefix for outputs (default: the CIF filename stem,
  e.g. `7C4J.cif` → prefix `7C4J` → `7C4J_helices.txt` etc).

### `parse_helices_2.py` — helices + per-chain protein / UniProt map

Parses the `_struct_conf` mmCIF category (author-deposited secondary structure,
*not* recomputed from coordinates).  Recognises α-helices, π-helices, and 3₁₀
helices; also walks the `_entity_poly` / `_struct_ref` / `_entity_name_com`
tables so each helix and each chain row is annotated with:

- protein name (first comma-separated alias)
- **UniProt accession** (`_struct_ref.pdbx_db_accession` when `db_name = UNP`)
- UniProt mnemonic code (`db_code`, e.g. `SMCA4_HUMAN`)
- source organism (`_entity_src_gen` / `_entity_src_syn`)

**Outputs** (`PREFIX` default = CIF stem):

| File                              | Description                                       |
|-----------------------------------|---------------------------------------------------|
| `PREFIX_helices.txt/.tex/.csv`    | One row per helix: Chain, ID, Start, End, Length, UniProt, Protein |
| `PREFIX_chains.txt/.tex/.csv`     | One row per chain: Chain, Entity, UniProt, db_code, Organism, Protein |

**Usage:**
```bash
python parse_helices_2.py path/to/6LTJ.cif
python parse_helices_2.py path/to/7EGM.cif --out-dir results --prefix 7EGM
```

---

### `parse_helices.py` — helices with short protein codes & % helical content

Same helix parser as above but with **three extra views** tailored for tables in
a paper:

1. **Variant 1 (compact)** — short protein *code* (`SWP73`, `SNF2`, `ARP7`, … style
   using a longest-match heuristic over the protein name aliases) + start + end.
   Great for tight supplementary tables.
2. **Variant 2 (detailed)** — subunit name + chain, start, end, **length**, and
   a `\cmidrule` footer row with *#helices / total helix residues / total aa /
   % helical content* per subunit, plus a GRAND TOTAL row.
3. **Protein AA summary** — unique protein codes, chains collapsed by entity-id,
   total amino acids, UniProt, full name.

**Outputs:**

| File                                      | Description                                       |
|-------------------------------------------|---------------------------------------------------|
| `PREFIX_helices_v1.txt/.tex`              | Short code | start | end                          |
| `PREFIX_helices_v2.txt/.tex`              | Subunit | start | end | length + footer summary      |
| `PREFIX_helices.csv`                      | Combined CSV with `total_aa` per helix row        |
| `PREFIX_protein_aa_summary.txt/.tex/.csv` | Unique protein codes → total aa table             |

**Usage:**
```bash
python parse_helices.py path/to/7c4j.cif --out-dir out --prefix 7C4J
```

---

### `parse_ligands.py` — ligand & HETATM atom coordinates

Dumps every **HETATM** residue (all `_atom_site` rows with `group_PDB == HETATM`)
in the file, one atom per line, split into *two classes*:

- **small-molecule / true ligands** (`_sm`) — non-polymer HETATM, excludes
  water (`HOH/WAT/SOL/…`) and 1–2 letter monoatomic ions (`MG`, `ZN`, `NA`,
  `CL`, `CA`, …) by default. Toggle with `--keep-water-ions-in-sm`.
- **all HETATM** (`_all`) — includes water + monoatomic ions.

Every atom row lists: `ResName, Chain, ResID[+icode], AtomName, X, Y, Z`.
Alternative locations are de-duplicated (prefer alt `""`/`"."` else first
alphabetical altloc).

**Outputs:**

| File                              | Description                                       |
|-----------------------------------|---------------------------------------------------|
| `PREFIX_ligands_sm.txt/.tex/.csv`  | True small-molecule ligands only                  |
| `PREFIX_ligands_all.txt/.tex/.csv` | Every HETATM (incl. water, monoatomic ions)       |

CSV columns include boolean flags `is_water / is_ion / is_small_mol` so you can
re-filter downstream without re-parsing the CIF.

**Usage:**
```bash
# Just the "real" ligands (no water / ions)
python parse_ligands.py 4LXZ.cif --out-dir out

# Keep water + ions inside the `_sm` output too (all HETATM is always written)
python parse_ligands.py 4LXZ.cif --keep-water-ions-in-sm
```

---

### `extract_residues.py` — local environment around a target residue

Given a **chain ID + residue number**, extracts:

1. All atoms of the target residue.
2. All atoms in *any* residue within `--cutoff` Å (default 3.5 Å) of *any*
   target atom.  Distance matrix is computed with numpy for speed.
3. A second block: **cross N↔O heteroatom distances** between the target and
   each surrounding residue (O-O and N-N pairs are filtered out; only N↔O are
   shown, tagged `N→O` or `O→N`).  Useful for eyeballing H-bond /
   hydrogen-bond-like contacts between a ligand and nearby protein residues.

**Outputs:**

| File                               | Description                                         |
|------------------------------------|-----------------------------------------------------|
| `PREFIX_extracted_atoms.txt`       | Chain / ResID / ResName / AtomName / X / Y / Z      |
| `PREFIX_hetero_distances.txt`      | Target-labelled N↔O cross distances + pair type     |

(Feed the first `.txt` to `txt2latex.py` to get a LaTeX table.)

**Usage:**
```bash
python extract_residues.py 4LXZ.cif \
       --chain A --resid 301 --cutoff 3.5 \
       --out-dir out --prefix 4LXZ
```

---

### `parse_residues.py` — selected residues + inter-residue N↔O H-bond table

Focused tool: for a list of *specific* residues (default: Chain A Ser86, His232,
Asp210 — the classic catalytic triad), extracts:

- every atom (standard order: N, CA, C, O first, then side-chain by name)
- **every inter-residue N↔O pairwise distance**, classified:
  - `strong`   ≤ `--strong-cutoff` Å (default 3.0)
  - `putative` ≤ `--hbond-cutoff`  Å (default 3.5)
  - `long-range` otherwise

**Three ways to choose the target residues:**

1. **Interactively** — run the script with no `--residue` and it prompts for a
   chain + residues one-by-one.  Blank first entry → uses the default triad.
2. **On the command line** — repeatable `--residue SEQ:RESNAME:LABEL`,
   e.g. `--residue 86:SER:Ser86 --residue 232:HIS:His232` (skips the prompt).
3. **Non-interactive default** — add `--no-prompt` to skip prompts entirely and
   fall back to the built-in catalytic Ser-His-Asp triad.

**Outputs:**

| File                                 | Description                                        |
|--------------------------------------|----------------------------------------------------|
| `PREFIX_selected_residues.txt/.tex`  | Per-atom coordinates table                         |
| `PREFIX_NO_hbonds.txt`               | Sorted N↔O cross-pair distances + classification + summary counts |

**Usage:**
```bash
# 1. Interactive mode: prompts for chain + residues (good for exploration)
python parse_residues.py 4MYD.cif --out-dir out

# 2. Non-interactive: same default triad, guaranteed no TTY prompts (CI / scripts)
python parse_residues.py 4MYD.cif --out-dir out --no-prompt

# 3. Fully explicit command-line set (catalytic Ser-His-Asp in another PDB)
python parse_residues.py 1GCI.cif \
       --chain A \
       --residue 195:SER:Ser195 \
       --residue 57:HIS:His57 \
       --residue 102:ASP:Asp102 \
       --out-dir out --prefix 1GCI_triad
```

---

### `txt2latex.py` — TXT → LaTeX table helper

Converts the columnar `.txt` produced by `extract_residues.py` (and similar
whitespace-separated tables with header
`Chain ResID ResName AtomName X Y Z`) into an inline LaTeX `booktabs` table
snippet.  Automatically:

- drops the `Chain` column (kept in the caption context)
- inserts a `\cmidrule` break between the target residue block and the
  surrounding-residue block
- writes a `\caption` citing counts + the target residue label

**Usage:**
```bash
# Default: reads extracted_atoms.txt → writes extracted_atoms.tex (adjacent)
python txt2latex.py extracted_atoms.txt

# Explicit input/output
python txt2latex.py out/4LXZ_extracted_atoms.txt -o out/4LXZ_extracted_atoms.tex
```

Preamble requirements for the generated `.tex`:
```latex
\usepackage{booktabs}
\usepackage{graphicx}   % for \resizebox (used in a few scripts)
\usepackage{longtable}  % for multi-page helix / ligand tables
```

---

## Example outputs bundled in the repo

All PDB identifiers used in examples / tests:

| PDB  | Files                                                                 |
|------|-----------------------------------------------------------------------|
| 6LTJ | `6ltj_helices.*` / `6ltj_chains.*`                                    |
| 7C4J | `7C4J_helices*.txt,tex` / `7C4J_protein_aa_summary.*` / chains + csv  |
| 7EGM | `7EGM_helices.*` / `7EGM_chains.*`                                    |
| 4LXZ | `4LXZ_ligands_sm.*` / `4LXZ_ligands_all.*` / `extracted_atoms.*` / `hetero_distances.txt` |
| 4MYD | `4myd_selected_residues.*` / `4myd_NO_hbonds.txt`                    |

These are *sample outputs* — regenerate them any time with the commands shown
above to verify against the latest wwPDB mmCIF.

---

## Notes

- All scripts use **Biopython's `MMCIF2Dict`** directly rather than the full
  `MMCIFParser` — this avoids building an atom-structure object when we only
  want tabular fields, and is significantly faster on large assemblies.
- Chain IDs are read from `auth_asym_id` (PDB author-facing IDs) and fall back
  to `label_asym_id` — matching the numbering you see on the wwPDB website.
- Hardcoded absolute paths present in the original working copies have all been
  replaced with the CLI interface shown above.  The only place your filesystem
  appears is the arguments you pass on the command line.

---

## License

MIT — use freely for academic or commercial work.  Citation of the underlying
PDB entries (and Biopython) is always appreciated when you publish figures or
tables produced with these scripts.
