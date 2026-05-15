# restretto-py

This directory contains a Python reimplementation of the core REstretto workflow
described in `restretto_implementation.md`.

The implementation uses OpenBabel for chemistry file IO, bond perception,
rotor detection, canonical smiles for ligand-sized molecules, and SDF output.
It currently supports:

- REstretto-style config parsing and validation.
- X-Score atom constants and Vina-like energy terms.
- Basic vector, atom, molecule, graph-distance, RMSD, and intra-energy helpers.
- Binary `.grid` read/write compatible with the documented layout.
- `.mol2`, `.pdbqt`, `.pdb`, and `.sdf` reading through OpenBabel.
- CLI commands: `atomgrid-gen`, `score-only`, `intraenergy-only`, `decompose`,
  `conformer-docking`, `atom-docking`, and `easytest-docking`.

The docking command is a deterministic simplified pipeline: it reads receptor
and ligand files with OpenBabel, scores input conformers, sorts them by score,
and writes SDF output with OpenBabel. Full C++ parity for fragment reuse search
and local optimization remains a future extension point.

## Environment

Use the local virtual environment only:

```powershell
py -3 -m venv .venv
.venv\Scripts\python -m pip install -e .
```

If pip networking is blocked, download/install wheels inside `.venv`; the only
required runtime package is `openbabel-wheel`.

## Test

```powershell
.venv\Scripts\python -m unittest discover -s tests -v
```

## Example

```powershell
.venv\Scripts\python -m restretto.cli atomgrid-gen references\restretto\testdata\testgrid.in
.venv\Scripts\python -m restretto.cli conformer-docking references\restretto\testdata\testgrid.in
```
