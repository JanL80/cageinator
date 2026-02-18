
import sys, re, json, math
from pathlib import Path
from collections import Counter
import numpy as np



COVALENT_R = {
    "H": 0.31, "B": 0.85, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
    "Si": 1.11
}



CONNECTOR_ELEMENTS = {
    "N","O","S","P","Cl","Br","I"
}



def bond_cutoff(el1, el2, scale=1.25):
    r = COVALENT_R.get(el1, 0.77) + COVALENT_R.get(el2, 0.77)
    return scale * r



def u(v):
    v = np.asarray(v, float)
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else v



def dumps_inline(obj, indent=2, inline_keys={"xyz","vector"}):
    ### Pretty JSON with single-line arrays for xyz/vector
    def fmt(v, level=0, parent_key=None):
        sp = " " * (indent * level)
        if isinstance(v, dict):
            if not v: return "{}"
            items = list(v.items())
            lines = ["{"]
            for i, (k, vv) in enumerate(items):
                comma = "," if i < len(items) - 1 else ""
                lines.append(" "*(indent*(level+1)) + json.dumps(k) + ": " + fmt(vv, level+1, k) + comma)
            lines.append(sp + "}")
            return "\n".join(lines)
        if isinstance(v, list):
            if parent_key in inline_keys and all(not isinstance(x, (list, dict)) for x in v):
                return "[" + ", ".join(json.dumps(x) for x in v) + "]"
            if not v: return "[]"
            lines = ["["]
            for i, x in enumerate(v):
                comma = "," if i < len(v) - 1 else ""
                lines.append(" "*(indent*(level+1)) + fmt(x, level+1) + comma)
            lines.append(sp + "]")
            return "\n".join(lines)
        return json.dumps(v)
    return fmt(obj) + "\n"



### Parsing different filetypes section

def parse_xyz(path: Path):
    lines = [ln.rstrip() for ln in path.read_text().splitlines() if ln.strip()]
    atoms = []
    comment = ""
    if lines and re.match(r"^\s*\d+\s*$", lines[0]):
        n = int(lines[0].strip())
        comment = lines[1] if len(lines) > 1 else ""
        data = lines[2:2+n]
    else:
        data = lines
    for i, ln in enumerate(data):
        parts = ln.split()
        if len(parts) < 4: continue
        el = parts[0]
        try:
            x, y, z = map(float, parts[1:4])
        except ValueError:
            continue
        atoms.append({"idx": len(atoms), "el": el, "xyz": [x, y, z]})
    return atoms, comment

def parse_mol_obabel(path: Path):
    try:
        from openbabel import openbabel as ob
    except Exception:
        import openbabel as ob

    ext = path.suffix.lower()
    infmt = "sdf" if ext in {".sdf", ".sd"} else "mdl"
    conv = ob.OBConversion()
    if not conv.SetInFormat(infmt):
        raise ValueError(f"OpenBabel: cannot set input format '{infmt}'")

    mol = ob.OBMol()
    if not conv.ReadFile(mol, str(path)):
        raise ValueError(f"OpenBabel: failed to read '{path}'")

    atoms = []
    for a in ob.OBMolAtomIter(mol):
        atoms.append({
            "idx": a.GetIdx() - 1,
            "el": ob.GetSymbol(a.GetAtomicNum()),
            "xyz": [float(a.GetX()), float(a.GetY()), float(a.GetZ())],
        })

    bonds = []
    for b in ob.OBMolBondIter(mol):
        order = int(b.GetBondOrder())
        bonds.append({
            "a": b.GetBeginAtomIdx() - 1,
            "b": b.GetEndAtomIdx() - 1,
            "type": str(order)
        })

    comment = ""
    return atoms, bonds, comment

def parse_mol2(path: Path):
    lines = path.read_text().splitlines()
    atoms = []
    bonds = []
    in_atom = False
    in_bond = False

    ### Temporary storage with original (1-based) mol2 atom indices
    tmp_atoms = {}

    for ln in lines:
        line = ln.strip()
        if not line:
            continue

        if line.startswith("@<TRIPOS>ATOM"):
            in_atom = True
            in_bond = False
            continue
        if line.startswith("@<TRIPOS>BOND"):
            in_atom = False
            in_bond = True
            continue
        if line.startswith("@<TRIPOS>"):
            in_atom = False
            in_bond = False
            continue

        if in_atom:
            parts = line.split()
            ### id, name, x, y, z, atom_type, ...
            if len(parts) < 6:
                continue
            atom_id = int(parts[0])
            x, y, z = map(float, parts[2:5])
            atom_type = parts[5]
            ### crude element guess, leading letters before digits or dot
            m = re.match(r"([A-Za-z]+)", atom_type)
            el = m.group(1) if m else atom_type
            tmp_atoms[atom_id] = {"orig_idx": atom_id, "el": el, "xyz": [x, y, z]}

        elif in_bond:
            parts = line.split()
            ### id, a1, a2, type
            if len(parts) < 4:
                continue
            a1 = int(parts[1])
            a2 = int(parts[2])
            btype = parts[3]
            bonds.append({"a_orig": a1, "b_orig": a2, "type": btype})

    ### Reindex atoms to 0..N-1 in ascending original index
    sorted_ids = sorted(tmp_atoms.keys())
    idx_map = {orig_id: new_idx for new_idx, orig_id in enumerate(sorted_ids)}

    atoms = []
    for orig_id in sorted_ids:
        a = tmp_atoms[orig_id]
        atoms.append({
            "idx": idx_map[orig_id],
            "el": a["el"],
            "xyz": a["xyz"]
        })

    ### Remap bonds to new 0-based indices
    remapped_bonds = []
    for b in bonds:
        if b["a_orig"] not in idx_map or b["b_orig"] not in idx_map:
            continue
        remapped_bonds.append({
            "a": idx_map[b["a_orig"]],
            "b": idx_map[b["b_orig"]],
            "type": b["type"]
        })

    comment = ""  ### no charge/mult parsing from mol2 here
    return atoms, remapped_bonds, comment



def build_bonds(atoms):
    bonds = []
    coords = np.array([a["xyz"] for a in atoms], float)
    elems  = [a["el"] for a in atoms]
    n = len(atoms)
    for i in range(n):
        for j in range(i+1, n):
            d = np.linalg.norm(coords[i]-coords[j])
            if d <= bond_cutoff(elems[i], elems[j]):
                bonds.append({"a": i, "b": j})
    return bonds



def neighbor_map(bonds, n_atoms):
    nbrs = {i: [] for i in range(n_atoms)}
    for b in bonds:
        i, j = b["a"], b["b"]
        nbrs[i].append(j); nbrs[j].append(i)
    return nbrs



def guess_connectors(atoms, bonds):
    coords = np.array([a["xyz"] for a in atoms], float)
    elems  = [a["el"] for a in atoms]
    nbrs = neighbor_map(bonds, len(atoms))
    conns = []
    for i, a in enumerate(atoms):
        el = elems[i]
        if el not in CONNECTOR_ELEMENTS:
            continue
        heavy = [j for j in nbrs[i] if elems[j] != "H"]                 ### heavy neighbors
        if not heavy:
            v = u(coords[i] - coords.mean(axis=0))                      ### point away from centroid if isolated
        else:
            vecs = [u(coords[j] - coords[i]) for j in heavy[:3]]        ### donor vector; opposite of average bond directions
            v = u(-np.sum(vecs, axis=0))
        conns.append({
            "id": f"{el}{i}",
            "atom_index": i,
            "element": el,
            "vector": [float(v[0]), float(v[1]), float(v[2])]
        })
    return conns



def hill_formula(atoms):
    counts = Counter(a["el"] for a in atoms)
    parts = []
    if "C" in counts:
        c = counts.pop("C")
        parts.append(f"C{'' if c==1 else c}")
    if "H" in counts:
        h = counts.pop("H")
        parts.append(f"H{'' if h==1 else h}")
    for el in sorted(counts.keys()):
        n = counts[el]
        parts.append(f"{el}{'' if n==1 else n}")
    return "".join(parts)



def parse_charge_mult(comment):
    charge = 0
    mult = 1
    if comment:
        m = re.search(r"charge\s*=\s*([+-]?\d+)", comment, re.I)
        if m: charge = int(m.group(1))
        m = re.search(r"(mult|multiplicity)\s*=\s*(\d+)", comment, re.I)
        if m: mult = int(m.group(2))
    return charge, mult



def main():
    if len(sys.argv) < 2:
        print("Usage: python coordinate_adder.py input.(xyz|mol2) [output.json]", file=sys.stderr)
        sys.exit(1)

    in_geom = Path(sys.argv[1])
    out_json = Path(sys.argv[2]) if len(sys.argv) > 2 else in_geom.with_suffix(".json")

    ext = in_geom.suffix.lower()

    if ext == ".mol2":
        atoms, bonds, comment = parse_mol2(in_geom)
    elif ext in {".mol", ".sdf", ".sd"}:
        atoms, bonds, comment = parse_mol_obabel(in_geom)
    else:
        atoms, comment = parse_xyz(in_geom)
        if not atoms:
            print("No atoms parsed.", file=sys.stderr); sys.exit(2)
        bonds = build_bonds(atoms)

    if not atoms:
        print("No atoms parsed.", file=sys.stderr); sys.exit(2)

    connectors = guess_connectors(atoms, bonds)

    data = {
        "unit_id": in_geom.stem,
        "class": "",
        "gbu_type": "",
        "gbu_subtype": "",
        "composition": {
            "formula": hill_formula(atoms),
            "charge": None,
            "spin_mult": None
        },
        "atoms": atoms,
        "bonds": bonds,  ### may include "type" field if from MOL2
        "connectors": connectors,
        "coordination_count": len(connectors),
        "coordination_atoms": [c["atom_index"] for c in connectors]
    }

    ch, mult = parse_charge_mult(comment)  ### fill charge/multiplicity if present in comment, else defaults (0,1)
    data["composition"]["charge"] = ch
    data["composition"]["spin_mult"] = mult

    out_json.write_text(dumps_inline(data, inline_keys={"xyz","vector"}))
    print("\n")
    print(str(out_json))



if __name__ == "__main__":
    main()
    