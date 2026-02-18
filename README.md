# The Cageinator
A tool for automatically constructing and screening metal-organic cages from linker and node building blocks

---

## ⚠️ Attention Notice ⚠️
This is a preliminary release and is currently in active development. Release versions of the style **0.0.x** indicate a development stage. The package itself, file formats, and behavior may change at any time without prior notice, and backwards compatibility is not guaranteed. Further, current/newly implemented functionalities may show unexpected behaviour until bug fixes have been applied. Use at your own risk.

Since this tool is primarily developed for UNIX-based operating systems such as Linux and macOS, running Cageinator on Windows may produce unexpected behavior. While the core functionality should work, the use with external components (e.g., xTB or Open Babel) is not guaranteed in the current release. Therefore, using WSL is recommended until a future release provides better Windows support.

---

## Setting up the Cageinator

### Requirements

- Python 3.10+ recommended
- External tools required if optimizations are required:
  - xTB (S. Grimme, C. Bannwarth, P. Shushkov, J. Chem. Theory Comput. 2017, 13, 1989–2009. https://doi.org/10.1021/acs.jctc.7b00118)
  - Open Babel (N. M. O’Boyle, M. Banck, C. A. James, C. Morley, T. Vandermeersch, G. R. Hutchison, J. Cheminf. 2011, 3, 33. https://doi.org/10.1186/1758-2946-3-33)

### Installation

At a directory of choice:

    git clone https://github.com/JanL80/cageinator
    cd ./cageinator
    pip install .

---

## Usage

### Overview

The cageinator builds metal-organic cages by combining node JSON files with linker input files, assembling one or more cage shapes, and optionally optimizing the resulting structures with Open Babel (force-field) and/or xTB

The program (currently) supports two modes:

- **CLI batch mode (default)**: runs a full batch pipeline from folders
- **Interactive menu mode (deprecated; will be removed in future versions)**: launches an interactive interface

### Input expectations:

- Nodes: `*.json` files (node building units)
- Linkers: `*.json`, `*.xyz`, or `*.mol` files
  - `*.xyz`/`*.mol` are converted to JSON internally
  - currently, xyz files are recommended as input

### Commandline call

```bash
cageinator --nodes NODES_DIR --linkers LINKERS_DIR --out OUT_DIR [OPTIONS]
```

To view the help section of the cageinator:
```bash
cageinator --help
```

---

## Example Inputs

### Example 1: Basic assembly (build all shapes)

```bash
cageinator --nodes ./nodes --linkers ./linkers --out ./project_out
```

This will:
- copy/prepare node JSONs into `OUT_DIR/nodes_json/`
- copy/convert linker files into `OUT_DIR/linkers_json/`
- build all supported cage shapes
- place results under `OUT_DIR/assemblies/` (each assembly into its own subdirectory)

### Example 2: Select specific shapes

Shapes can be provided as space-separated values, and each token may also contain comma-separated values

Allowed shapes:
- `cuboctaeder`
- `octaeder`
- `cuboctaeder_d3h`
- `lantern`
  
The shapes can be entered space separated (with annotations) or comma-separated:

```bash
### Space-separated
cageinator --nodes ./nodes --linkers ./linkers --out ./out --shape "cuboctaeder octaeder"

### Comma-separated
cageinator --nodes ./nodes --linkers ./linkers --out ./out --shape cuboctaeder,lantern
```

### Example 3: xTB optimization

Enable xTB optimization of all cages after assembly:

```bash
cageinator --nodes ./nodes --linkers ./linkers --out ./out --xtb-opt
```

Pass extra flags to the xTB driver using `--xtb-flags`:

```bash
cageinator --nodes ./nodes --linkers ./linkers --out ./out --xtb-opt --xtb-flags "--method gfn2 --threads 8"
```

Show xTB specific help section:

```bash
cageinator --xtb-help
```

Notes:
- If `xtb_opt.xyz` already exists in an assembly subfolder and is non-empty, that folder is skipped

### Example 4: Open Babel optimization (force-field)

Run Open Babel optimization (force-field) on all built cages:

```bash
cageinator --nodes ./nodes --linkers ./linkers --out ./out --obabel-opt
```

Pass extra flags to the Open Babel optimizer module via `--obabel-flags`:

```bash
cageinator --nodes ./nodes --linkers ./linkers --out ./out --obabel-opt --obabel-flags "--ff mmff94 --steps 500"
```

Notes:
- In **FF-only mode** (Open Babel without xTB), folders with existing `ff_opt.mol` or `xtb_opt.xyz` are skipped

### Example 5: Combined FF pre-optimization + xTB

If both `--obabel-opt` and `--xtb-opt` are enabled, force-field optimization is performed as a **pre-step** inside the xTB driver (per assembly folder), and then xTB is run on the FF-optimized structure:

```bash
cageinator --nodes ./nodes --linkers ./linkers --out ./out --obabel-opt --xtb-opt --obabel-flags "--ff gaff --steps 500" --xtb-flags "--method gfn2 --threads 8"
```

### Interactive menu mode (deprecated; will be removed in future versions)

Launch the interactive interface:

```bash
cageinator --menu
```

---

## Output structure

Typical output directories under `OUT_DIR`:

- `OUT_DIR/nodes_json/`  
  Node JSONs copied from `--nodes`

- `OUT_DIR/linkers_json/`  
  Linker JSONs copied or generated from `--linkers`

- `OUT_DIR/assemblies/`  
  Assembly outputs, reorganized into one subfolder per produced file/stem

Within each assembly subfolder, optimization outputs may include:
- `ff_opt.mol`, `ff_failed_last.mol`, `ff_opt.log`
- `xtb_opt.xyz`, `xtb_failed_last.xyz`, `xtb_opt.log`, and `tb_opt.trj`

---

## Exit codes

- `0`: success
- `2`: usage error / missing inputs / invalid shape / filesystem error
- `130`: interrupted (Ctrl+C)

---

## To-Do
- Robust error handling across the pipeline
- Resume/checkpointing to avoid rerunning failed or completed steps
- Output handling for generated ".out" file(s) and general results management
- Expand chemical/structural coverage: more linker types, more GBU geometries + additional assembly models
- Adding support for more input formats (e.g., CIF, PDB, Gaussian/Turbomole outputs)
- Adding an Encapsulation workflow for host-guest description
- Conformer handling: multiple orientations and a conformer search workflow (e.g., CREST)
- Evaluation/energy comparison to rank most stable/likely cages
- Adding an in-depth documentation
- Adding a node creation interface
- Bug fixes and minor polishing

---

## Citation

If you use this software in academic work, please cite in this style:


**Cageinator**: J. Leodolter, *cageinator* (Version X.Y.Z), GitHub: GitHub: https://github.com/JanL80/cageinator, manuscript in preparation (2026).  

<br/>

Additionally, please cite the external tools used by this workflow as applicable:

- S. Grimme, C. Bannwarth, P. Shushkov, J. Chem. Theory Comput. 2017, 13, 1989–2009. https://doi.org/10.1021/acs.jctc.7b00118  
- (If using GFN2-xTB) <br/>
C. Bannwarth, S. Ehlert, S. Grimme, J. Chem. Theory Comput. 2019, 15, 1652–1671. https://doi.org/10.1021/acs.jctc.8b01176  
- N. M. O’Boyle, M. Banck, C. A. James, C. Morley, T. Vandermeersch, G. R. Hutchison, J. Cheminf. 2011, 3, 33. https://doi.org/10.1186/1758-2946-3-33

