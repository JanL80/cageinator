import sys, json, re, math
from pathlib import Path
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict
import numpy as np
import itertools
import unicodedata
import copy



shape_name = "cuboctaeder"                     ### part 1/2 that should be changed



### Geometry Definition (part 2/2 to be changed per shape)
def shape_vertices():
    return np.array([
        [ 1, 1, 0],[ 1,-1, 0],
        [-1, 1, 0],[-1,-1, 0],
        [ 1, 0, 1],[ 1, 0,-1],
        [-1, 0, 1],[-1, 0,-1],
        [ 0, 1, 1],[ 0, 1,-1],
        [ 0,-1, 1],[ 0,-1,-1]
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



def _sanitize_atoms_inplace(u):
    for a in (u.get("atoms") or []):
        a["el"] = sanitize_element(a.get("el") or a.get("element") or "")



def _count_preexisting_carboxylate_deprot(u) -> int:
    atoms = (u.get("atoms") or [])
    bonds = (u.get("bonds") or [])
    if not atoms:
        return 0
    _sanitize_atoms_inplace(u)
    groups = _find_carboxyl_groups(atoms, bonds)
    ### COO− if no OH hydrogen present
    return sum(1 for g in groups if g.get("H") is None)



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
            return 4
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

def _is_true(v):
    """Safely parse booleans from JSON that might be typed as strings"""
    if isinstance(v, str):
        return v.strip().lower() in {"true", "1", "yes"}
    return bool(v)

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
        if 0 <= i < len(atoms) and 0 <= j < len(atoms):
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
        is_true_carboxylate = True
        
        ### Verify these are actually terminal carboxylate oxygens (not esters/anhydrides)
        for o in onbrs:
            heavy_nbrs = [j for j in nbrs[o] if atoms[j]["el"] != "H"]
            if len(heavy_nbrs) > 1:
                is_true_carboxylate = False
                break
                
            hnei = [j for j in nbrs[o]
                    if atoms[j]["el"] == "H" and len(nbrs[j]) == 1]
            if hnei:
                oh = o
                h_on_oh = hnei[0]
                
        if not is_true_carboxylate:
            continue
            
        out.append({"C": ci, "O": onbrs, "OH": oh, "H": h_on_oh})
    return out

def _is_carboxylate_oxygen(i, atoms, bonds):
    if atoms[i]["el"] != "O":
        return False
    nbrs = _make_nbrs(atoms, bonds)
    if any(atoms[j]["el"] == "H" for j in nbrs[i]):             ### OH excluded
        return False
        
    ### Prevent ester misidentification
    heavy_nbrs = [j for j in nbrs[i] if atoms[j]["el"] != "H"]
    if len(heavy_nbrs) > 1:
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
    heavy_nbrs = [j for j in nbrs[i] if atoms[j]["el"] != "H"]
    
    ### 1. Pyridine N has exactly 2 heavy neighbors (both C) and NO hydrogens
    if len(heavy_nbrs) != 2 or len(nbrs[i]) != 2:
        return False
    if atoms[heavy_nbrs[0]]["el"] != "C" or atoms[heavy_nbrs[1]]["el"] != "C":
        return False
        
    ### 2. Trace 6-membered rings containing this N using DFS
    def find_6_rings(start):
        rings = []
        def dfs(curr, path):
            if len(path) == 6:
                if start in nbrs[curr]:
                    rings.append(path)
                return
            for nxt in nbrs[curr]:
                if atoms[nxt]["el"] != "H" and nxt not in path:
                    dfs(nxt, path + [nxt])
        dfs(start, [start])
        return rings

    rings = find_6_rings(i)
    if not rings:
        return False
        
    ### 3. Ensure the ring is strictly a Pyridine (1 N, 5 C, sp2-hybridized)
    for r in rings:
        elements = [atoms[idx]["el"] for idx in r]
        
        if elements.count("N") == 1 and elements.count("C") == 5:
            ### Check for aromatic/sp2 nature: no ring atom can be sp3 (max 3 heavy neighbors)
            is_valid_aromatic = True
            for idx in r:
                h_deg = sum(1 for j in nbrs[idx] if atoms[j]["el"] != "H")
                if h_deg > 3:
                    is_valid_aromatic = False
                    break
            
            if is_valid_aromatic:
                return True
                
    return False

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

def find_template_matches(linker, template):
    """
    Finds all sub-graph matches of `template` within `linker`.
    Now features STRICT evaluation using on-demand cached DFS cycles.
    """
    l_atoms = linker.get("atoms", [])
    l_bonds = linker.get("bonds", [])
    l_nbrs = _make_nbrs(l_atoms, l_bonds)

    t_atoms = template.get("atoms", [])
    t_bonds = template.get("bonds", [])
    t_nbrs = _make_nbrs(t_atoms, t_bonds)

    node_rings_cache = {}

    def get_rings(node_idx, req_size):
        if (node_idx, req_size) not in node_rings_cache:
            rings = []
            def dfs(curr, path):
                if len(path) == req_size:
                    if node_idx in l_nbrs[curr]:
                        r = tuple(sorted(path))
                        if r not in rings:
                            rings.append(r)
                    return
                for nxt in l_nbrs[curr]:
                    if nxt not in path:
                        dfs(nxt, path + [nxt])
            dfs(node_idx, [node_idx])
            node_rings_cache[(node_idx, req_size)] = rings
        return node_rings_cache[(node_idx, req_size)]

    matches = []

    def backtrack(t_idx, current_mapping):
        if t_idx == len(t_atoms):
            if _is_true(template.get("exact_subgraph")) or _is_true(template.get("induced_subgraph")):
                for ta1 in range(len(t_atoms)):
                    for ta2 in range(ta1 + 1, len(t_atoms)):
                        t_bonded = ta2 in t_nbrs[ta1]
                        l_bonded = current_mapping[ta2] in l_nbrs[current_mapping[ta1]]
                        if t_bonded != l_bonded:
                            return
            
            matches.append(current_mapping.copy())
            return

        t_atom = t_atoms[t_idx]
        t_el = t_atom.get("el") or t_atom.get("element")
        t_deg = len(t_nbrs[t_idx])

        for l_idx, l_atom in enumerate(l_atoms):
            l_el = l_atom.get("el") or l_atom.get("element")
            
            ### Constraints for a valid match
            if sanitize_element(l_el) != sanitize_element(t_el): 
                continue
            if len(l_nbrs[l_idx]) < t_deg: 
                continue
            if l_idx in current_mapping.values(): 
                continue

            if "heavy_degree" in t_atom:
                l_heavy = sum(1 for n in l_nbrs[l_idx] if sanitize_element(l_atoms[n].get("el", "")) != "H")
                if l_heavy != int(t_atom["heavy_degree"]):
                    continue
                    
            if "total_degree" in t_atom:
                if len(l_nbrs[l_idx]) != int(t_atom["total_degree"]):
                    continue

            if "ring_size" in t_atom or "ring_heteroatoms" in t_atom:
                req_size = int(t_atom.get("ring_size", 6))
                my_rings = get_rings(l_idx, req_size)
                
                if "ring_size" in t_atom and not my_rings:
                    continue
                    
                if "ring_heteroatoms" in t_atom:
                    req_het = int(t_atom["ring_heteroatoms"])
                    valid_het = False
                    for r in my_rings:
                        het_count = sum(1 for idx in r if sanitize_element(l_atoms[idx].get("el", "")) not in ("C", "H"))
                        if het_count == req_het:
                            valid_het = True
                            break
                    if not valid_het:
                        continue

            valid_topology = True
            for prev_t_idx in t_nbrs[t_idx]:
                if prev_t_idx < t_idx:
                    prev_l_idx = current_mapping[prev_t_idx]
                    if prev_l_idx not in l_nbrs[l_idx]:
                        valid_topology = False
                        break
            
            if valid_topology:
                current_mapping[t_idx] = l_idx
                backtrack(t_idx + 1, current_mapping)
                del current_mapping[t_idx]

    backtrack(0, {})
    return matches



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

def infer_donors(linker, templates=None):
    _sanitize_atoms_inplace(linker)
    atoms = linker.get("atoms", []) or []
    bonds = linker.get("bonds", []) or []
    if not atoms:
        raise ValueError("no atoms in linker")

    ### detect carboxyl groups and map O->group for deprotonation and deduplication
    groups = _find_carboxyl_groups(atoms, bonds)
    o_to_group = {}
    for g in groups:
        for o in g["O"]:
            o_to_group[o] = g

    donors_set = set()
    seen_groups = set()  ### tracks parent carbons to prevent picking both O's in a carboxylate
    template_deprot_flags = {} ### maps matched donor atom index -> should it deprotonate?

    ### 1. Explicit donor connectors
    for c in (linker.get("connectors") or []):
        if c.get("role") != "donor":
            continue
        ai = c.get("atom_index")
        if ai is None:
            continue
        ai = int(ai)
        if not (0 <= ai < len(atoms)):
            continue
        
        g = o_to_group.get(ai)
        if g is not None:
            pc = g["C"]
            if pc in seen_groups:
                continue
            seen_groups.add(pc)
            
        donors_set.add(ai)

    ### 2. STRICT Mode vs Fallback Mode
    if templates is not None:
        for template in templates:
            donor_t_indices = [
                i for i, a in enumerate(template.get("atoms", [])) 
                if _is_true(a.get("is_donor", False))
            ]
            
            matches = find_template_matches(linker, template)
            for match in matches:
                for d_t_idx in donor_t_indices:
                    matched_idx = match[d_t_idx]
                    
                    g = o_to_group.get(matched_idx)
                    if g is not None:
                        pc = g["C"]
                        if pc in seen_groups:
                            continue
                        seen_groups.add(pc)
                        
                    donors_set.add(matched_idx)

                    if _is_true(template.get("atoms", [])[d_t_idx].get("deprotonate", False)):
                        template_deprot_flags[matched_idx] = True

        if not donors_set:
            name = linker.get("name") or linker.get("id") or "linker"
            print(f"Skipping {name}: Strict mode active, but no template matches found.")
            return None
            
    else:
        ### 3. Fallback Heuristics (only triggers if templates == None)
        carb_groups = []
        for g in groups:
            carb_groups.append({
                "C": g["C"],
                "O": list(g["O"]),
                "has_H": (g.get("H") is not None),
                "H": g.get("H")
            })
        carb_groups.sort(key=lambda x: x["C"])

        ### COO− groups
        for cg in carb_groups:
            if not cg["has_H"]:
                oidx = int(min(cg["O"]))
                if cg["C"] not in seen_groups and oidx not in donors_set:
                    donors_set.add(oidx)
                    seen_groups.add(cg["C"])

        ### N donors
        Ns = [i for i in range(len(atoms)) if _is_N_donor(i, atoms, bonds)]
        for n in sorted(Ns):
            donors_set.add(int(n))

        ### fallback COOH groups
        for cg in carb_groups:
            if cg["has_H"]:
                oidx = int(min(cg["O"]))
                if cg["C"] not in seen_groups and oidx not in donors_set:
                    donors_set.add(oidx)
                    seen_groups.add(cg["C"])

    ### Dedupe and convert to list
    donors = list(donors_set)

    if len(donors) < 2:
        name = linker.get("name") or linker.get("id") or "linker"
        print(f"Skipping {name}: only {len(donors)} eligible donors found.")
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

    remove = []
    nbrs = _make_nbrs(atoms, bonds)
    for d in chosen:
        g = o_to_group.get(d)
        if g is not None:
            h_idx = g.get("H")
            if h_idx is not None and 0 <= h_idx < len(atoms) and atoms[h_idx]["el"] == "H" and len(nbrs[h_idx]) == 1:
                remove.append(int(h_idx))
                
        elif template_deprot_flags.get(d, False):
            terminal_hs = [nbr for nbr in nbrs[d] if atoms[nbr]["el"] == "H" and len(nbrs[nbr]) == 1]
            if terminal_hs:
                remove.append(int(terminal_hs[0]))

    remove = list(set(remove))

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



def make_linker_type(link_json, templates=None) -> Optional[LinkerType]:
    _sanitize_atoms_inplace(link_json)
    pair = infer_donors(link_json, templates)
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

def assemble_shape(node_json, linker_json, out_dir: Path, templates=None):

    node_name   = Path(sys.argv[1]).stem
    linker_name = Path(sys.argv[2]).stem

    nt = make_node_type(node_json)

    ### count pre-existing COO− before mutating linker
    pre_dep0 = _count_preexisting_carboxylate_deprot(copy.deepcopy(linker_json))

    lt = make_linker_type(linker_json, templates)
    if lt is None:
        print("No build: linker ineligible")
        return None

    M_N = float(node_json.get("M_N", 2.05))

    ### ignore declared linker charge; derive from geometry
    M = M_MULTIPLICITY
    L = L_MULTIPLICITY
    
    ### Safely extract node charge from either the composition block or top-level
    q_comp = (node_json.get("composition") or {}).get("charge")
    q_top = node_json.get("charge")
    if q_comp is not None:
        q_node = _int_or_zero(q_comp)
    elif q_top is not None:
        q_node = _int_or_zero(q_top)
    else:
        q_node = 0
        
    added_dep = int(getattr(lt, "deprot_per_linker", 0))
    q_linker_final = - (pre_dep0 + added_dep)
    q_total = M * q_node + L * q_linker_final

    ### Determine names and create the cage directory early
    stoich_tag = f"M{M}L{L}_{shape_name}"
    q_tag = f"Q{q_total:+d}"
    mol_name = f"{node_name}__{linker_name}__{stoich_tag}__{q_tag}"
    
    cage_dir = out_dir / mol_name
    cage_dir.mkdir(parents=True, exist_ok=True)

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

    ### Backbone Isolation
    adj = {idx: [] for idx in range(len(lt.elems))}
    for b in (lt.bonds or []):
        adj[int(b["a"])].append(int(b["b"]))
        adj[int(b["b"])].append(int(b["a"]))
        
    queue_bb = [(donor1, [donor1])]
    visited_bb = {donor1}
    backbone_path = []
    while queue_bb:
        curr, path = queue_bb.pop(0)
        if curr == donor2:
            backbone_path = path
            break
        for nxt in adj[curr]:
            if nxt not in visited_bb:
                visited_bb.add(nxt)
                queue_bb.append((nxt, path + [nxt]))
                
    if not backbone_path:
        backbone_path = list(range(len(lt.elems)))

    ### Base function to generate oriented linker coordinates
    def _compute_linker_coords(edge_idx, flipped, twist_deg):
        i, j = E[edge_idx]
        A = anchor(i, j, s_opt)
        B = anchor(j, i, s_opt)
        
        ### Mirrored flip end-to-end
        if flipped:
            A, B = B, A
            
        mid = 0.5 * (A + B)
        u_ab = nrm(B - A)

        ### initial alignment of donor–donor axis to A–B
        Q = np.vstack([lt.xyz[donor1], lt.xyz[donor2]])
        P = np.vstack([A, B])
        Rl, tl = kabsch(Q, P)
        lX = (Rl @ lt.xyz.T).T + tl

        lCOM_all = lX.mean(axis=0)
        lCOM_backbone = lX[backbone_path].mean(axis=0)
        
        v_bb = lCOM_backbone - mid
        v_bb_perp = v_bb - np.dot(v_bb, u_ab) * u_ab
        
        if np.linalg.norm(v_bb_perp) > 0.2:
            v = v_bb
        else:
            v = lCOM_all - mid

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
            
            if np.linalg.norm(v_bb_perp) > 0.2:
                v = lX[backbone_path].mean(axis=0) - mid
            else:
                v = lX.mean(axis=0) - mid
            v_perp = v - np.dot(v, u_ab) * u_ab

        ### ensure linker points outward
        if np.dot(v_perp, r_perp) <= 0:
            lX = rotate_about_axis(lX, origin=mid, axis_unit=u_ab, theta=math.pi)

        ### Apply the corrective twist around the A-B axis
        twist_rad = math.radians(twist_deg)
        if abs(twist_rad) > 1e-6:
            lX = rotate_about_axis(lX, origin=mid, axis_unit=u_ab, theta=twist_rad)

        return lX

    ### Dynamic caching system
    coords_cache = {}
    def get_coords(edge_idx, flipped, twist_deg):
        key = (edge_idx, flipped, twist_deg)
        if key not in coords_cache:
            coords_cache[key] = _compute_linker_coords(edge_idx, flipped, twist_deg)
        return coords_cache[key]

    ### Build dynamic collision threshold matrix
    is_H = np.array([e == 'H' for e in lt.elems])
    thresh_matrix = np.full((len(lt.elems), len(lt.elems)), 1.2)   ### Heavy-Heavy clashes
    thresh_matrix[is_H, :] = 0.8                                   ### Heavy-Hydrogen
    thresh_matrix[:, is_H] = 0.8                                   ### Hydrogen-Heavy
    thresh_matrix[is_H[:, None] & is_H[None, :]] = 0.0             ### Ignore H-H overlaps completely
    thresh_sq = thresh_matrix ** 2

    def evaluate_clashes(coords_list):
        clashing_edges = []
        total_clashes = 0
        for idx1 in range(len(coords_list)):
            for idx2 in range(idx1 + 1, len(coords_list)):
                c1 = coords_list[idx1]
                c2 = coords_list[idx2]
                
                dists_sq = np.sum((c1[:, None, :] - c2[None, :, :])**2, axis=2)
                clashes_here = int(np.sum(dists_sq < thresh_sq))
                if clashes_here > 0:
                    total_clashes += clashes_here
                    clashing_edges.extend([idx1, idx2])
        return total_clashes, clashing_edges

    ### Organize targeted search sequence: small twists first, flipping as absolute last resort
    test_sequence = []
    for deg in [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8, -8, 9, -9, 10, -10, 11, -11, 12, -12, 13, -13, 14, -14, 15, -15]:
        test_sequence.append((False, deg))
    for deg in [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8, -8, 9, -9, 10, -10, 11, -11, 12, -12, 13, -13, 14, -14, 15, -15]:
        test_sequence.append((True, deg))


    ### Independent Cavity Depth Safeguard Check
    current_states = [(False, 0)] * len(E)
    current_coords = [get_coords(i, f, t) for i, (f, t) in enumerate(current_states)]
    
    ### Hard limit for minimum distance of ANY single atom to the cage center
    hard_depth_threshold = 5.0 
    
    initial_atomic_dists = [float(np.min(np.linalg.norm(c - COM_nodes, axis=1))) for c in current_coords]
    initial_min_dist = min(initial_atomic_dists) if initial_atomic_dists else float('inf')

    if initial_min_dist < hard_depth_threshold:
        print(f"Construction failed for {mol_name}: Atoms reaching too far into the cavity (Min atomic dist: {initial_min_dist:.2f} A).")
        return None

    ### Iterative Relaxation Loop (Greedy Hill-Climbing)
    loop_depth_threshold = 0.1
    max_allowed_clashes = 0
    max_iters = 50
    success = False

    for iteration in range(max_iters):
        clashes, clashing_edges = evaluate_clashes(current_coords)
        
        atomic_dists = [float(np.min(np.linalg.norm(c - COM_nodes, axis=1))) for c in current_coords]
        min_atomic_dist = min(atomic_dists) if atomic_dists else float('inf')

        ### Exit if requirements met
        if clashes <= max_allowed_clashes and min_atomic_dist >= loop_depth_threshold:
            success = True
            break
            
        improvement_found = False
        edges_to_check = list(set(clashing_edges)) if clashes > max_allowed_clashes else [int(np.argmin(atomic_dists))]

        for edge_idx in edges_to_check:
            original_f, original_t = current_states[edge_idx]
            
            for f, t in test_sequence:
                if (f, t) == (original_f, original_t):
                    continue
                
                ### Substitute the coordinate for this specific edge and test
                test_coords = current_coords.copy()
                test_coords[edge_idx] = get_coords(edge_idx, f, t)
                
                t_clashes, _ = evaluate_clashes(test_coords)
                t_min_atomic = min([float(np.min(np.linalg.norm(c - COM_nodes, axis=1))) for c in test_coords])
                
                ### Is this the smallest rotation that improves the structure?
                if (t_clashes < clashes) or (t_clashes == clashes and t_min_atomic > min_atomic_dist and min_atomic_dist < loop_depth_threshold):
                    current_states[edge_idx] = (f, t)
                    current_coords[edge_idx] = test_coords[edge_idx]
                    clashes = t_clashes
                    min_atomic_dist = t_min_atomic
                    improvement_found = True
                    break ### Apply improvement instantly and start next global iteration
            
            if improvement_found:
                break
                
        if not improvement_found:
            break ### Stuck in local minimum, could not resolve

    if not success:
        print(f"Construction failed for {mol_name}: Could not resolve steric clashes or depth limits after iterative relaxation (Remaining clashes: {clashes}, Min atomic dist: {min_atomic_dist:.2f} A).")
        return None

    chosen_linker_coords = current_coords

    ### Final Linker Placement
    for edge_idx, (i, j) in enumerate(E):
        lX = chosen_linker_coords[edge_idx]
        is_flipped, final_twist = current_states[edge_idx]
        
        ### place linker once and add its internal bonds
        offset = len(all_elems)
        all_elems.extend(lt.elems)
        all_xyz.extend(list(lX))

        for b in (lt.bonds or []):
            rec = {
                "a": offset + int(b["a"]),
                "b": offset + int(b["b"]),
            }
            if "type" in b:
                rec["type"] = b["type"]
            all_bonds.append(rec)

        ### connect donors to metals on nodes i and j
        metal_i = metal_global_idx[i]
        metal_j = metal_global_idx[j]

        donor1_global = offset + donor1  
        donor2_global = offset + donor2   

        if is_flipped:
            ### donor1 mapped to node j, donor2 mapped to node i
            all_bonds.append({"a": metal_j, "b": donor1_global, "type": "1"})
            all_bonds.append({"a": metal_i, "b": donor2_global, "type": "1"})
        else:
            ### donor1 mapped to node i, donor2 mapped to node j
            all_bonds.append({"a": metal_i, "b": donor1_global, "type": "1"})
            all_bonds.append({"a": metal_j, "b": donor2_global, "type": "1"})

    all_xyz_np = np.array(all_xyz, float)
    
    out = cage_dir / f"{mol_name}.mol"

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
    assemblers_dir = Path(__file__).resolve().parent
    template_dir = assemblers_dir / "coordination_templates"
    
    templates = None 
    
    if template_dir.exists() and template_dir.is_dir():
        templates = []
        for p in template_dir.glob("*.json"):
            templates.append(load_json_lenient(p))
        print(f"Strict Template Mode Active: Loaded {len(templates)} templates from '{template_dir}'")
    else:
        print(f"Warning: '{template_dir}' directory not found. Operating in Fallback Heuristic mode.")

    if len(sys.argv)!=4:
        print("Usage: python cuboctaeder.py NODE.json LINKER.json OUT_DIR"); sys.exit(1)
    node = load_json_lenient(Path(sys.argv[1]))
    linker = load_json_lenient(Path(sys.argv[2]))

    out_dir = Path(sys.argv[3]); out_dir.mkdir(parents=True, exist_ok=True)
    try:
        fn = assemble_shape(node, linker, out_dir, templates)
        if fn:
            print( "\n", "Wrote", fn)
    except Exception as e:
        print("Error:", e); sys.exit(2)



if __name__=="__main__":
    main()