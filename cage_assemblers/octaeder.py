
import sys, json, re, math
from pathlib import Path
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict
import numpy as np
import itertools
import unicodedata
import copy



shape_name = "octaeder"                     ### part 1/2 that should be changed



### Geometry Definition (part 2/2 to be changed per shape)
def shape_vertices():
    return np.array([
        [ 1, 0, 0],
        [-1, 0, 0],
        [ 0, 1, 0],
        [ 0,-1, 0],
        [ 0, 0, 1],
        [ 0, 0,-1],
    ], float)



def shape_edges(V):
    E=[]
    for i in range(len(V)):
        for j in range(i+1, len(V)):
            if np.isclose(np.linalg.norm(V[j]-V[i]), math.sqrt(2.0)):
                E.append((i,j))
    return E



def best_rotation_and_assignment(site_dirs_local: List[np.ndarray], target_dirs: List[np.ndarray]):
    P = np.stack(site_dirs_local, axis=0)
    best = (1e9, None, None)
    for perm in itertools.permutations(range(4)):
        Q = np.stack([target_dirs[k] for k in perm], axis=0)
        R, t = kabsch(P, Q)
        D = (R @ P.T).T + t
        rms = float(np.sqrt(np.mean(np.sum((D-Q)**2, axis=1))))
        if rms < best[0]:
            best = (rms, R, perm)
    return best[1], best[2]



def _shape_stoich():                          ### read out number of nodes for total formal charge determination
    V0 = shape_vertices()
    return len(V0), len(shape_edges(V0))



M_MULTIPLICITY, L_MULTIPLICITY = _shape_stoich()



PERIODIC = {
    "H","He","Li","Be","B","C","N","O","F","Ne",
    "Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
    "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
    "Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
    "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
    "Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
    "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
    "Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
    "Tl","Pb","Bi","Po","At","Rn"
}



def sanitize_element(tok: str) -> str:
    tok = re.sub(r"[^A-Za-z]", "", tok or "")
    if not tok: return "X"
    if len(tok)==1: el = tok.upper()
    else: el = tok[0].upper()+tok[1:].lower()
    return el if el in PERIODIC else "X"



def ascii_only(s: str) -> str:
    return unicodedata.normalize("NFKD", s).encode("ascii","ignore").decode("ascii")



### currently unused; old xyz writer, just for reference (may be removed in the future)
### replaced with mol writer (see write_mol_strict below)
def write_xyz_strict(path: Path, elems: List[str], xyz: np.ndarray, comment: str="") -> None:

    E = [sanitize_element(e) for e in elems]                                                            ### sanitize elements and ensure finite coordinates
    if any(e=="X" for e in E):
        bad = {e for e in E if e=="X"}
        raise ValueError(f"Invalid element symbol(s) in XYZ: {bad}")
    if not np.isfinite(xyz).all():
        raise ValueError("Non-finite coordinates in XYZ")

    com = xyz.mean(axis=0)                                                                              ### center structure to improve viewer heuristics
    Xc = xyz - com

    cmt = ascii_only(comment).replace("\n"," ").strip()                                                 ### ASCII-only, single-line comment

    with open(path, "w", encoding="ascii", newline="\n") as f:                                          ### write
        f.write(f"{len(E)}\n{cmt}\n")
        for e,(x,y,z) in zip(E, Xc):
            f.write(f"{e:<2s} {x: .6f} {y: .6f} {z: .6f}\n")

    with open(path, "r", encoding="ascii") as f:                                                        ### read-back validation
        lines = f.read().splitlines()
    try:
        n = int(lines[0].strip())
    except Exception as e:
        raise ValueError("XYZ header count not integer") from e
    if len(lines) != 2 + n:
        raise ValueError(f"XYZ line count mismatch: header {n}, file has {len(lines)-2}")
    for k in range(2, 2+n):
        parts = lines[k].split()
        if len(parts)!=4:
            raise ValueError(f"Bad XYZ line {k+1}: {lines[k]!r}")
        _el, x, y, z = parts
        float(x); float(y); float(z)                                                                    ### will raise on failure



def write_mol_strict(path: Path,
                     elems: List[str],
                     xyz: np.ndarray,
                     bonds: Optional[List[Dict]] = None,
                     name: str = "assembled",
                     comment: str = "") -> None:
    bonds = bonds or []

    ### sanitize elements and ensure finite coordinates
    E = [sanitize_element(e) for e in elems]
    if any(e == "X" for e in E):
        bad = {e for e in E if e == "X"}
        raise ValueError(f"Invalid element symbol(s) in mol: {bad}")
    if not np.isfinite(xyz).all():
        raise ValueError("Non-finite coordinates in mol")

    ### center structure
    com = xyz.mean(axis=0)
    Xc = xyz - com

    ### ASCII-only name and comment
    mol_name = ascii_only(name).replace("\n", " ").strip() or "assembled"
    cmt = ascii_only(comment).replace("\n", " ").strip()

    natoms = int(len(E))
    nbonds = int(len(bonds))

    def bond_order_from_type(t) -> int:
        if t is None:
            return 1
        s = str(t).strip().lower()
        if s in {"1", "single"}:
            return 1
        if s in {"2", "double"}:
            return 2
        if s in {"3", "triple"}:
            return 3
        if s in {"4", "ar", "aro", "aromatic"}:
            return 4  # aromatic in MDL bond order
        if s in {"am"}:
            return 1
        try:
            v = int(float(s))
            return v if v in (1, 2, 3, 4) else 1
        except Exception:
            return 1

    ### Validate bond indices
    for b in bonds:
        a = int(b["a"])
        c = int(b["b"])
        if not (0 <= a < natoms and 0 <= c < natoms):
            raise ValueError(f"Bond atom index out of range: a={a}, b={c}, natoms={natoms}")

    with open(path, "w", encoding="ascii", newline="\n") as f:
        ### MDL molfile header (3 lines)
        f.write(f"{mol_name}\n")
        f.write(f"{cmt}\n" if cmt else "\n")
        f.write("Generated by assembly script\n")

        ### V3000 marker line (counts are in the V3000 COUNTS record)
        ### Keep standard spacing
        f.write("  0  0  0  0  0  0            999 V3000\n")

        ### CTAB begin
        f.write("M  V30 BEGIN CTAB\n")
        f.write(f"M  V30 COUNTS {natoms} {nbonds} 0 0 0\n")

        ### Atoms
        f.write("M  V30 BEGIN ATOM\n")
        for i, (el, (x, y, z)) in enumerate(zip(E, Xc), start=1):
            f.write(f"M  V30 {i} {el} {x:.6f} {y:.6f} {z:.6f} 0\n")
        f.write("M  V30 END ATOM\n")

        ### Bonds
        f.write("M  V30 BEGIN BOND\n")
        for k, b in enumerate(bonds, start=1):
            a = int(b["a"]) + 1
            c = int(b["b"]) + 1
            order = bond_order_from_type(b.get("type", "1"))
            f.write(f"M  V30 {k} {order} {a} {c}\n")
        f.write("M  V30 END BOND\n")

        f.write("M  V30 END CTAB\n")
        f.write("M  END\n")



### json I/O
def load_json_lenient(path: Path):
    txt = path.read_text(encoding="utf-8")
    try:
        return json.loads(txt)
    except json.JSONDecodeError:
        return json.loads(re.sub(r",(\s*[}\]])", r"\1", txt))



### Linear Algebra section

def nrm(v):
    v = np.array(v, float)
    n = np.linalg.norm(v)
    if n == 0: raise ValueError("zero vector")
    return v / n

def kabsch(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1,:] *= -1
        R = Vt.T @ U.T
    t = Q.mean(axis=0) - R @ P.mean(axis=0)
    return R, t

def golden_section(f, a, b, iters=48):
    phi=(1+5**0.5)/2; inv=1.0/phi
    c=b-inv*(b-a); d=a+inv*(b-a)
    fc=f(c); fd=f(d)
    for _ in range(iters):
        if fc<fd: b,d,fd=d,c,fc; c=b-inv*(b-a); fc=f(c)
        else: a,c,fc=c,d,fd; d=a+inv*(b-a); fd=f(d)
    x=0.5*(a+b); return x, f(x)



### Parse section

def unit_atoms(u):
    elems = [(a.get("el") or a.get("element") or "").capitalize() for a in u["atoms"]]
    xyz   = np.array([a["xyz"] for a in u["atoms"]], float)
    return elems, xyz

def get_connectors(u, role=None):
    return [c for c in u.get("connectors", []) if role is None or c.get("role")==role]

def _make_nbrs(atoms, bonds):
    nbrs = [[] for _ in range(len(atoms))]
    for b in bonds or []:
        i, j = int(b["a"]), int(b["b"])
        nbrs[i].append(j); nbrs[j].append(i)
    return nbrs

def _remove_atoms_inplace(linker, remove_idxs):
    if not remove_idxs:
        return
    atoms = linker.get("atoms", []) or []
    bonds = linker.get("bonds", []) or []
    conns = linker.get("connectors", []) or []
    ### only allow removing H
    rm = {int(i) for i in remove_idxs
          if 0 <= int(i) < len(atoms) and atoms[int(i)]["el"] == "H"}
    if not rm:
        return
    idx_map = {}
    new_atoms = []
    for old_i, a in enumerate(atoms):
        if old_i in rm:
            idx_map[old_i] = None
        else:
            idx_map[old_i] = len(new_atoms)
            new_atoms.append(a)
    new_bonds = []
    for b in bonds or []:
        a = idx_map.get(int(b["a"]))
        c = idx_map.get(int(b["b"]))
        if a is None or c is None:
            continue
        new_bonds.append({"a": a, "b": c})
    new_conns = []
    for c in conns or []:
        ai = c.get("atom_index")
        if ai is None:
            continue
        m = idx_map.get(int(ai))
        if m is None:
            continue
        c2 = dict(c); c2["atom_index"] = m
        new_conns.append(c2)
    linker["atoms"] = new_atoms
    linker["bonds"] = new_bonds
    linker["connectors"] = new_conns

def _find_carboxyl_groups(atoms, bonds):
    nbrs = _make_nbrs(atoms, bonds)
    out = []
    for ci, a in enumerate(atoms):
        if a["el"] != "C":
            continue
        onbrs = [j for j in nbrs[ci] if atoms[j]["el"] == "O"]
        if len(onbrs) != 2:
            continue
        oh = None
        h_on_oh = None
        ### pick the oxygen that actually bears a terminal H (degree-1 H)
        for o in onbrs:
            hnei = [j for j in nbrs[o]
                    if atoms[j]["el"] == "H" and len(nbrs[j]) == 1]
            if hnei:
                oh = o
                h_on_oh = hnei[0]
                break
        out.append({"C": ci, "O": onbrs, "OH": oh, "H": h_on_oh})
    return out

def _is_carboxylate_oxygen(i, atoms, bonds):
    if atoms[i]["el"] != "O":
        return False
    nbrs = _make_nbrs(atoms, bonds)
    if any(atoms[j]["el"] == "H" for j in nbrs[i]):                 ### OH excluded
        return False
    cnbrs = [j for j in nbrs[i] if atoms[j]["el"] == "C"]
    if len(cnbrs) != 1:
        return False
    c = cnbrs[0]
    onbrs = [j for j in nbrs[c] if atoms[j]["el"] == "O"]
    return len(onbrs) == 2

def _is_N_donor(i, atoms, bonds):
    if atoms[i]["el"] != "N":
        return False
    nbrs = _make_nbrs(atoms, bonds)
    heavy_deg = sum(1 for j in nbrs[i] if atoms[j]["el"] != "H")
    return heavy_deg <= 3                                           ### not quaternary

def _carboxylate_O_parent_map(atoms, bonds):
    nbrs = _make_nbrs(atoms, bonds)
    m = {}
    for c in range(len(atoms)):
        if atoms[c]["el"] != "C":
            continue
        onbrs = [j for j in nbrs[c] if atoms[j]["el"] == "O"]
        if len(onbrs) != 2:
            continue
        for o in onbrs:
            if any(atoms[h]["el"] == "H" for h in nbrs[o]):
                continue
            m[o] = c
    return m

def _remove_atoms_with_map_inplace(linker, remove_idxs):
    if not remove_idxs:
        n = len(linker.get("atoms", []) or [])
        return {i: i for i in range(n)}
    atoms = linker.get("atoms", []) or []
    bonds = linker.get("bonds", []) or []
    conns = linker.get("connectors", []) or []
    ### only allow removing H
    rm = {int(i) for i in remove_idxs
          if 0 <= int(i) < len(atoms) and atoms[int(i)]["el"] == "H"}
    idx_map = {}
    new_atoms = []
    for old_i, a in enumerate(atoms):
        if old_i in rm:
            idx_map[old_i] = None
        else:
            idx_map[old_i] = len(new_atoms)
            new_atoms.append(a)
    new_bonds = []
    for b in bonds or []:
        a = idx_map.get(int(b["a"]))
        c = idx_map.get(int(b["b"]))
        if a is None or c is None:
            continue
        new_bonds.append({"a": a, "b": c})
    new_conns = []
    for c in conns or []:
        ai = c.get("atom_index")
        if ai is None:
            continue
        m = idx_map.get(int(ai))
        if m is None:
            continue
        c2 = dict(c); c2["atom_index"] = m
        new_conns.append(c2)
    linker["atoms"] = new_atoms
    linker["bonds"] = new_bonds
    linker["connectors"] = new_conns
    return idx_map



### charge logic section

def _sanitize_atoms_inplace(u):
    for a in (u.get("atoms") or []):
        a["el"] = sanitize_element(
            a.get("el") or a.get("element") or ""
        )

def _count_preexisting_carboxylate_deprot(u) -> int:
    atoms = (u.get("atoms") or [])
    bonds = (u.get("bonds") or [])
    if not atoms:
        return 0
    _sanitize_atoms_inplace(u)
    groups = _find_carboxyl_groups(atoms, bonds)
    ### COO− if no OH hydrogen present
    return sum(1 for g in groups if g.get("H") is None)



### Donor inference section

def _choose_donor_pair_geom(linker, donors):

    atoms = linker.get("atoms", []) or []
    donors = list(dict.fromkeys(int(i) for i in donors))

    if len(donors) < 2:
        return None

    ### If the JSON provides coordination_atoms, restrict to those
    coord_set = set(linker.get("coordination_atoms") or [])
    if coord_set:
        coord_donors = [i for i in donors if i in coord_set]
        if len(coord_donors) >= 2:
            donors = coord_donors

    if len(donors) == 2:
        return int(donors[0]), int(donors[1])

    ### Max distance pair among donors
    coords = np.array([atoms[i]["xyz"] for i in donors], float)
    best_d2 = -1.0
    best_pair = (donors[0], donors[1])

    for ii in range(len(donors)):
        for jj in range(ii + 1, len(donors)):
            d2 = float(np.sum((coords[jj] - coords[ii])**2))
            if d2 > best_d2:
                best_d2 = d2
                best_pair = (donors[ii], donors[jj])

    return int(best_pair[0]), int(best_pair[1])

def infer_donors(linker):
    _sanitize_atoms_inplace(linker)
    atoms = linker.get("atoms", []) or []
    bonds = linker.get("bonds", []) or []
    if not atoms:
        raise ValueError("no atoms in linker")

    ### detect carboxyl groups and map O->group
    groups = _find_carboxyl_groups(atoms, bonds)
    o_to_group = {}
    for g in groups:
        for o in g["O"]:
            o_to_group[o] = g

    ### normalized group records, deterministic order by parent carbon index
    carb_groups = []
    for g in groups:
        carb_groups.append({
            "C": g["C"],
            "O": list(g["O"]),
            "has_H": (g.get("H") is not None),
            "H": g.get("H")
        })
    carb_groups.sort(key=lambda x: x["C"])

    donors: list[int] = []
    seen_groups: set[int] = set()

    ### explicit donor connectors
    for c in (linker.get("connectors") or []):
        if c.get("role") != "donor":
            continue
        ai = c.get("atom_index")
        if ai is None:
            continue
        ai = int(ai)
        if not (0 <= ai < len(atoms)):
            continue

        if atoms[ai]["el"] == "N":
            if _is_N_donor(ai, atoms, bonds):
                donors.append(ai)
        else:
            g = o_to_group.get(ai)
            if g is not None:
                pc = g["C"]
                if pc in seen_groups:
                    continue
                donors.append(ai)
                seen_groups.add(pc)

    ### COO− groups
    for cg in carb_groups:
        if not cg["has_H"]:
            oidx = int(min(cg["O"]))
            if cg["C"] not in seen_groups and oidx not in donors:
                donors.append(oidx)
                seen_groups.add(cg["C"])

    ### N donors
    Ns = [i for i in range(len(atoms)) if _is_N_donor(i, atoms, bonds)]
    for n in sorted(Ns):
        if n not in donors:
            donors.append(int(n))

    ### fallback COOH groups
    for cg in carb_groups:
        if cg["has_H"]:
            oidx = int(min(cg["O"]))
            if cg["C"] not in seen_groups and oidx not in donors:
                donors.append(oidx)
                seen_groups.add(cg["C"])

    ### dedupe for safety
    donors = list(dict.fromkeys(donors))

    if len(donors) < 2:
        name = linker.get("name") or linker.get("id") or "linker"
        print(f"Skipping {name}: only {len(donors)} eligible donors")
        return None

    ### geometric choice of the two least sterically blocked donors
    chosen = _choose_donor_pair_geom(linker, donors)
    if chosen is None:
        name = linker.get("name") or linker.get("id") or "linker"
        print(f"Skipping {name}: unable to choose donor pair geometrically")
        return None

    ### deprotonation logic
    rb_atoms = copy.deepcopy(linker.get("atoms", []) or [])
    rb_bonds = copy.deepcopy(linker.get("bonds", []) or [])
    rb_conns = copy.deepcopy(linker.get("connectors", []) or [])

    remove: List[int] = []
    nbrs = _make_nbrs(atoms, bonds)
    for d in chosen:
        g = o_to_group.get(d)
        if g is None:
            continue
        h_idx = g.get("H")
        if h_idx is None:
            continue
        if 0 <= h_idx < len(atoms) and atoms[h_idx]["el"] == "H" and len(nbrs[h_idx]) == 1:
            remove.append(int(h_idx))

    if remove:
        def _elt_counts(A):
            t = {}
            for a in A:
                t[a["el"]] = t.get(a["el"], 0) + 1
            return t

        before = _elt_counts(linker.get("atoms", []) or [])
        idx_map = _remove_atoms_with_map_inplace(linker, remove)
        after = _elt_counts(linker.get("atoms", []) or [])

        nonH_bad = any(e != "H" and after.get(e, 0) < before.get(e, 0) for e in before)
        total_bad = (sum(after.values()) != sum(before.values()) - len(remove))
        if nonH_bad or total_bad:
            linker["atoms"] = rb_atoms
            linker["bonds"] = rb_bonds
            linker["connectors"] = rb_conns
            remove = []
            idx_map = {i: i for i in range(len(linker["atoms"]))}

        if remove:
            d0 = idx_map.get(chosen[0]); d1 = idx_map.get(chosen[1])
            if d0 is None or d1 is None:
                linker["atoms"] = rb_atoms
                linker["bonds"] = rb_bonds
                linker["connectors"] = rb_conns
                return None
            chosen = (int(d0), int(d1))

    linker["_deprot_count"] = int(len(remove))
    return chosen



### dataclass section

@dataclass(eq=False)
class NodeType:
    elems: List[str]
    xyz_local: np.ndarray
    site_dirs_local: List[np.ndarray]
    bonds: List[Dict]
    metal_idx: int
    S: int

@dataclass(eq=False)
class LinkerType:
    elems: List[str]
    xyz: np.ndarray
    bonds: List[Dict]
    donor_idx: Tuple[int,int]
    d12: float
    deprot_per_linker: int = 0



def make_node_type(node_json) -> NodeType:
    _sanitize_atoms_inplace(node_json)
    ne, nxyz = unit_atoms(node_json)
    bonds = node_json.get("bonds", []) or []

    if "metal_atom_index" in node_json:
        mi = int(node_json["metal_atom_index"])
    else:
        masses = {
            "H":1,"C":12,"N":14,"O":16,"F":19,
            "P":31,"S":32,"Cl":35,"Pd":106,"Cu":64,
            "Zn":65,"Ag":108,"Pt":195,"Au":197
        }
        mi = int(np.argmax([masses.get(e, 0) for e in ne]))
    mpos = nxyz[mi].copy()
    sites = get_connectors(node_json, "metal_site")
    if len(sites)!=4:
        raise ValueError(f"node must have 4 metal_site connectors for {shape_name}")
    dirs = [ nrm(np.array(c["vector"], float)) for c in sites ]
    return NodeType(ne, nxyz - mpos, dirs, bonds, mi, 4)



def make_linker_type(link_json) -> Optional[LinkerType]:
    _sanitize_atoms_inplace(link_json)
    pair = infer_donors(link_json)
    if pair is None:
        name = link_json.get("name") or link_json.get("id") or "linker"
        print(f"Skipping {name}: insufficient eligible donors")
        return None

    le, lxyz = unit_atoms(link_json)
    bonds = link_json.get("bonds", []) or []

    i1, i2 = pair
    d12 = float(np.linalg.norm(lxyz[i2] - lxyz[i1]))
    if d12 < 1e-6:
        name = link_json.get("name") or link_json.get("id") or "linker"
        print(f"Skipping {name}: donor1 and donor2 coincide")
        return None
    deprot = int(link_json.get("_deprot_count", 0))
    for k in ("_rollback_atoms", "_rollback_bonds", "_rollback_connectors"):
        if k in link_json:
            link_json.pop(k, None)
    return LinkerType(le, lxyz, bonds, (i1, i2), d12, deprot_per_linker=deprot)



### Build section

def _int_or_zero(x):
    try:
        return int(round(float(x)))
    except Exception:
        return 0

def assemble_shape(node_json, linker_json, out_dir: Path):

    node_name   = Path(sys.argv[1]).stem
    linker_name = Path(sys.argv[2]).stem

    nt = make_node_type(node_json)

    ### count pre-existing COO− BEFORE mutating linker
    pre_dep0 = _count_preexisting_carboxylate_deprot(copy.deepcopy(linker_json))

    lt = make_linker_type(linker_json)
    if lt is None:
        print("No build: linker ineligible")
        return None
    M_N = float(node_json.get("M_N", 2.05))

    ### ignore declared linker charge; derive from geometry
    M = M_MULTIPLICITY
    L = L_MULTIPLICITY
    q_node = _int_or_zero((node_json.get("composition") or {}).get("charge"))
    added_dep = int(getattr(lt, "deprot_per_linker", 0))
    q_linker_final = - (pre_dep0 + added_dep)
    q_total = M * q_node + L * q_linker_final

    V0 = shape_vertices()
    E  = shape_edges(V0)
    nbrs = {i: [] for i in range(len(V0))}
    for i, j in E:
        nbrs[i].append(j); nbrs[j].append(i)

    target_dirs_by_vertex = {
        i: [nrm(V0[j] - V0[i]) for j in nbrs[i]]
        for i in range(len(V0))
    }

    Rnode = {}
    site_for_neighbor = {}
    for i in range(len(V0)):
        R, perm = best_rotation_and_assignment(nt.site_dirs_local, target_dirs_by_vertex[i])
        if R is None:
            raise RuntimeError("node-site fit failed")
        Rnode[i] = R
        inv = [None] * 4
        for local_idx, k in enumerate(perm):
            inv[k] = local_idx
        for k, j in enumerate(nbrs[i]):
            site_for_neighbor[(i, j)] = inv[k]

    def anchor(i, j, s):
        R = Rnode[i]
        site = site_for_neighbor[(i, j)]
        return s * V0[i] + (R @ nt.site_dirs_local[site]) * M_N

    def golden_objective(s):
        err = 0.0
        for (i, j) in E:
            A = anchor(i, j, s); B = anchor(j, i, s)
            err += abs(np.linalg.norm(B - A) - lt.d12)
        return err / len(E)

    s_opt, _ = golden_section(golden_objective, 0.1, 30.0)

    node_coords = []
    for i in range(len(V0)):
        R = Rnode[i]; t = s_opt * V0[i]
        nX = (R @ nt.xyz_local.T).T + t
        node_coords.append(nX)

    mass_lut = {
        "H":1.0079,"B":10.81,"C":12.0107,"N":14.0067,"O":15.999,
        "F":18.998,"Na":22.99,"Mg":24.305,"Al":26.982,"Si":28.085,
        "P":30.973,"S":32.065,"Cl":35.453,"K":39.098,"Ca":40.078,
        "Sc":44.956,"Ti":47.867,"V":50.942,"Cr":51.996,"Mn":54.938,
        "Fe":55.845,"Co":58.933,"Ni":58.693,"Cu":63.546,"Zn":65.38,
        "Ga":69.723,"Ge":72.63,"As":74.922,"Se":78.971,"Br":79.904,
        "Kr":83.798,"Rb":85.468,"Sr":87.62,"Y":88.906,"Zr":91.224,
        "Nb":92.906,"Mo":95.95,"Pd":106.42,"Ag":107.8682,"Cd":112.414,
        "In":114.818,"Sn":118.71,"Sb":121.76,"Te":127.60,"I":126.90447,
        "Xe":131.293,"Pt":195.084,"Au":196.966569,"Hg":200.592,"Tl":204.38,
        "Pb":207.2,"Bi":208.9804
    }

    def amass(el: str) -> float:
        return float(mass_lut.get(el, 12.0))

    masses_nodes = []
    coords_nodes = []
    for nX in node_coords:
        for e, r in zip(nt.elems, nX):
            masses_nodes.append(amass(e)); coords_nodes.append(r)
    masses_nodes = np.asarray(masses_nodes, float)
    coords_nodes = np.asarray(coords_nodes, float)
    COM_nodes = (masses_nodes[:, None] * coords_nodes).sum(axis=0) / masses_nodes.sum()

    def rotate_about_axis(X: np.ndarray, origin: np.ndarray, axis_unit: np.ndarray, theta: float) -> np.ndarray:
        u = nrm(axis_unit)
        c = math.cos(theta); s = math.sin(theta)
        ux, uy, uz = u
        K = np.array([[0, -uz, uy],[uz, 0, -ux],[-uy, ux, 0]], float)
        R = c*np.eye(3) + s*K + (1.0-c)*np.outer(u,u)
        return (R @ (X - origin).T).T + origin

    donor1, donor2 = lt.donor_idx
    all_elems: List[str] = []
    all_xyz:   List[np.ndarray] = []
    all_bonds: List[Dict] = []
    metal_global_idx: Dict[int, int] = {}

    ### place nodes once, collect internal node bonds and metal positions
    for v_idx, nX in enumerate(node_coords):
        base = len(all_elems)
        all_elems.extend(nt.elems)
        all_xyz.extend(list(nX))

        ### internal node bonds
        for b in (nt.bonds or []):
            ai = base + int(b["a"])
            bi = base + int(b["b"])
            rec = {"a": ai, "b": bi}
            if "type" in b:
                rec["type"] = b["type"]
            all_bonds.append(rec)

        ### remember where this nodes metal ended up globally
        metal_global_idx[v_idx] = base + nt.metal_idx

    ### place one linker per edge, collect linker bonds and metal-donor bonds
    for (i, j) in E:
        A = anchor(i, j, s_opt)
        B = anchor(j, i, s_opt)
        mid = 0.5 * (A + B)
        u_ab = nrm(B - A)

        ### orient linker so donor1 -> A, donor2 -> B
        Q = np.vstack([lt.xyz[donor1], lt.xyz[donor2]])     ### original donor positions
        P = np.vstack([A, B])                               ### target positions
        Rl, tl = kabsch(Q, P)
        lX = (Rl @ lt.xyz.T).T + tl

        ### COM-based “inside/outside” orientation
        lCOM = lX.mean(axis=0)
        v = lCOM - mid
        r = mid - COM_nodes

        r_perp = r - np.dot(r, u_ab) * u_ab
        if np.linalg.norm(r_perp) < 1e-10:
            r_perp = r
        v_perp = v - np.dot(v, u_ab) * u_ab

        nr = np.linalg.norm(r_perp)
        nv = np.linalg.norm(v_perp)
        if nr > 0 and nv > 0:
            r_hat = r_perp / nr
            v_hat = v_perp / nv
            cos_t = float(np.clip(np.dot(v_hat, r_hat), -1.0, 1.0))
            sin_t = float(np.dot(u_ab, np.cross(v_hat, r_hat)))
            theta = math.atan2(sin_t, cos_t)
            lX = rotate_about_axis(lX, origin=mid, axis_unit=u_ab, theta=theta)
            lCOM = lX.mean(axis=0)
            v = lCOM - mid
            v_perp = v - np.dot(v, u_ab) * u_ab

        if np.dot(v_perp, r_perp) <= 0:
            lX = rotate_about_axis(lX, origin=mid, axis_unit=u_ab, theta=math.pi)
            lCOM = lX.mean(axis=0)
            v = lCOM - mid
            v_perp = v - np.dot(v, u_ab) * u_ab

        ### add linker atoms
        base = len(all_elems)
        all_elems.extend(lt.elems)
        all_xyz.extend(list(lX))

        ### internal linker bonds
        for b in (lt.bonds or []):
            ai = base + int(b["a"])
            bi = base + int(b["b"])
            rec = {"a": ai, "b": bi}
            if "type" in b:
                rec["type"] = b["type"]
            all_bonds.append(rec)

        ### metal–donor bonds: donor1 is bound to node i, donor2 to node j
        metal_i = metal_global_idx[i]
        metal_j = metal_global_idx[j]
        donor1_global = base + donor1
        donor2_global = base + donor2

        all_bonds.append({"a": metal_i, "b": donor1_global, "type": "1"})
        all_bonds.append({"a": metal_j, "b": donor2_global, "type": "1"})

    all_xyz_np = np.array(all_xyz, float)
    stoich_tag = f"M{M}L{L}_{shape_name}"
    q_tag = f"Q{q_total:+d}"
    mol_name = f"{node_name}__{linker_name}__{stoich_tag}__{q_tag}"
    out = out_dir / f"{node_name}__{linker_name}__{stoich_tag}__{q_tag}.mol"

    comment = (
        f"{shape_name}; s={s_opt:.3f}; d12={lt.d12:.3f} A; M-N={M_N:.2f} A; "
        f"linkers aligned by COM-perp to node-COM; Q_total={q_total}"
    )

    write_mol_strict(out, all_elems, all_xyz_np,
                      bonds=all_bonds,
                      name=mol_name,
                      comment=comment)
    return str(out)



### Command-Line Interaction
def main():
    if len(sys.argv)!=4:
        print("Usage: python octaeder.py NODE.json LINKER.json OUT_DIR"); sys.exit(1)
    node = load_json_lenient(Path(sys.argv[1]))
    linker = load_json_lenient(Path(sys.argv[2]))

    out_dir = Path(sys.argv[3]); out_dir.mkdir(parents=True, exist_ok=True)
    try:
        fn = assemble_shape(node, linker, out_dir)
        print( "\n", "Wrote", fn)
    except Exception as e:
        print("Error:", e); sys.exit(2)



if __name__=="__main__":
    main()
