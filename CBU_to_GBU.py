
import json, argparse
from pathlib import Path
import numpy as np



### math section

def u(v):
    v = np.asarray(v, float)
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else v

def angle_deg(a, b):
    c = np.clip(np.dot(u(a), u(b)), -1.0, 1.0)
    return float(np.degrees(np.arccos(c)))

def fit_plane(points):
    P = np.asarray(points, float)
    c = P.mean(axis=0)
    U, S, Vt = np.linalg.svd(P - c)
    n = Vt[-1]
    rms = float(np.sqrt(np.mean(((P - c) @ n) ** 2)))
    e1, e2 = Vt[0], Vt[1]
    return c, n, rms, e1, e2

def vec_is_numeric3(v):
    try:
        a = np.asarray(v, float)
        return a.shape == (3,) and np.all(np.isfinite(a))
    except Exception:
        return False

def _fibonacci_sphere(n=256):                           ### exterior bonding partner detection
    i = np.arange(n, dtype=float)
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    z = 1.0 - 2.0 * (i + 0.5) / n
    r = np.sqrt(np.clip(1.0 - z * z, 0.0, 1.0))
    theta = 2.0 * np.pi * i / phi
    x = r * np.cos(theta); y = r * np.sin(theta)
    V = np.stack([x, y, z], axis=1)
    return V / (np.linalg.norm(V, axis=1, keepdims=True) + 1e-12)

def exterior_atom_mask(coords, n_dirs=256, tol=1e-9):
    C = np.asarray(coords, float)
    if C.size == 0:
        return np.zeros(0, dtype=bool), np.zeros(3)
    c = C.mean(axis=0)
    V = _fibonacci_sphere(n_dirs)
    proj = (C - c) @ V.T                                ### [N, n_dirs]
    m = proj.max(axis=0, keepdims=True)                 ### [1, n_dirs]
    outside = np.any(proj >= m - tol, axis=1)
    return outside, c



### check if Amine or Carboxylate section

def _make_nbrs(atoms, bonds):
    nbrs = [[] for _ in range(len(atoms))]
    for b in bonds:
        i, j = int(b["a"]), int(b["b"])
        nbrs[i].append(j); nbrs[j].append(i)
    return nbrs

def _is_carboxylate_oxygen(ai, atoms, nbrs):
    ### O with no H neighbors, bonded to a carbonyl carbon that has exactly two O neighbors
    if atoms[ai]["el"] != "O":
        return False
    if any(atoms[j]["el"] == "H" for j in nbrs[ai]):
        return False
    cnbrs = [j for j in nbrs[ai] if atoms[j]["el"] == "C"]
    if len(cnbrs) != 1:
        return False
    c = cnbrs[0]
    onbrs = [j for j in nbrs[c] if atoms[j]["el"] == "O"]
    if len(onbrs) != 2:
        return False
    ### at least one O neighbor of the carbon is heavy-degree 1 (only bonded to that carbon)
    def heavy_deg(k): return sum(1 for m in nbrs[k] if atoms[m]["el"] != "H")
    if not any(heavy_deg(o) == 1 for o in onbrs):
        return False
    return True

def _is_nitrogen_donor(ai, atoms, nbrs):
    ### Any N that is not quaternary (≤3 heavy neighbors).
    if atoms[ai]["el"] != "N":
        return False
    heavy_deg = sum(1 for j in nbrs[ai] if atoms[j]["el"] != "H")
    return heavy_deg <= 3

def bonding_atom_mask_carboxylate_N(atoms, bonds):
    ### Boolean mask over atoms that may contribute to bonding:
    ### - All eligible nitrogens (non-quaternary).
    ### - Carboxylates: at most one O per -COO group. Choose the more exterior O
    N = len(atoms)
    mask = np.zeros(N, dtype=bool)
    if N == 0:
        return mask

    nbrs = _make_nbrs(atoms, bonds)
    coords = np.array([a["xyz"] for a in atoms], float)
    centroid = coords.mean(axis=0)

    ### Nitrogen donors
    for i in range(N):
        if _is_nitrogen_donor(i, atoms, nbrs):
            mask[i] = True

    ### Carboxylate O donors: keep at most one O per carbonyl carbon
    carb_to_Os = {}
    for i in range(N):
        if not _is_carboxylate_oxygen(i, atoms, nbrs):
            continue
        cnbrs = [j for j in nbrs[i] if atoms[j]["el"] == "C"]
        if len(cnbrs) != 1:
            continue
        c = cnbrs[0]
        carb_to_Os.setdefault(c, []).append(i)

    for c, o_list in carb_to_Os.items():
        if not o_list:
            continue
        ### choose the more exterior oxygen (farther from centroid)
        best_o = max(o_list, key=lambda k: float(np.linalg.norm(coords[k] - centroid)))
        mask[best_o] = True

    return mask



### connector handling section

def connector_vectors(obj):

    atoms = obj.get("atoms", [])
    bonds = obj.get("bonds", [])
    conns = obj.get("connectors", [])
    if not conns:
        name = obj.get("name") or obj.get("id") or obj.get("title") or "object"
        print(f"No connectors found in {name}")
        return []

    if all(("vector" in c and vec_is_numeric3(c["vector"])) for c in conns):                ### Use explicit vectors only if all are numeric 3-vectors
        return [u(c["vector"]) for c in conns]

    nbrs = {i: [] for i in range(len(atoms))}                                               ### Build neighbor list
    for b in bonds:
        i, j = int(b["a"]), int(b["b"])
        nbrs[i].append(j); nbrs[j].append(i)
    coords = np.array([a["xyz"] for a in atoms], float)
    elements = [a["el"] for a in atoms]

    vecs = []
    for c in conns:
        ai = int(c["atom_index"])
        heavy = [j for j in nbrs.get(ai, []) if elements[j] != "H"]                         ### heavy neighbors
        if heavy:
            heavy = sorted(heavy, key=lambda j: np.linalg.norm(coords[j]-coords[ai]))[:3]   ### opposite of average bond directions
            v = -np.sum([u(coords[j]-coords[ai]) for j in heavy], axis=0)
            vecs.append(u(v))
        else:
            cent = coords.mean(axis=0)                                                      ### point away from centroid
            vecs.append(u(coords[ai] - cent))
    return vecs

def connector_vectors_filtered(obj, outward_dot_min=0.2, n_dirs=256):
    atoms = obj.get("atoms", [])
    bonds = obj.get("bonds", [])
    conns = obj.get("connectors", [])
    if not conns:
        name = obj.get("name") or obj.get("id") or obj.get("title") or "object"
        print(f"No connectors found in {name}")
        return []

    coords = np.array([a["xyz"] for a in atoms], float) if atoms else np.zeros((0, 3))
    bonding_mask = bonding_atom_mask_carboxylate_N(atoms, bonds)

    ### If a connector lacks atom_index, exclude it here because donor type is unknown.
    explicit = all(("vector" in c and vec_is_numeric3(c["vector"])) for c in conns)
    vecs_raw, ais = [], []
    if explicit:
        for c in conns:
            if "atom_index" not in c:
                continue
            ai = int(c["atom_index"])
            if not bonding_mask[ai]:
                continue
            vecs_raw.append(u(c["vector"]))
            ais.append(ai)
    else:
        nbrs = _make_nbrs(atoms, bonds)
        elements = [a["el"] for a in atoms]
        for c in conns:
            ai = int(c["atom_index"])
            if not bonding_mask[ai]:
                continue
            heavy = [j for j in nbrs.get(ai, []) if elements[j] != "H"]
            if heavy:
                heavy = sorted(heavy, key=lambda j: np.linalg.norm(coords[j]-coords[ai]))[:3]
                v = -np.sum([u(coords[j]-coords[ai]) for j in heavy], axis=0)
                vecs_raw.append(u(v)); ais.append(ai)
            else:
                cent = coords.mean(axis=0) if len(coords) else np.zeros(3)
                vecs_raw.append(u(coords[ai] - cent)); ais.append(ai)

    if not vecs_raw:
        ### Nothing usable after filtering
        name = obj.get("name") or obj.get("id") or obj.get("title") or "object"
        print(f"No usable connectors after filtering in {name}")
        return []

    outside_mask, c = exterior_atom_mask(coords, n_dirs=n_dirs)

    vecs_f = []
    for v, ai in zip(vecs_raw, ais):
        if not outside_mask[ai]:
            continue
        r = coords[ai] - c
        if np.dot(v, r) <= outward_dot_min:
            continue
        vecs_f.append(u(v))
    if not vecs_f:
        name = obj.get("name") or obj.get("id") or obj.get("title") or "object"
        print(f"*** No outward exterior connectors in {name} ***")
    return vecs_f



### 5-planar test
def pentagon_planarity_and_spacing(vecs, planarity_tol=0.15, gap_tol_deg=15.0):
    tips = np.asarray([u(v) for v in vecs], float)
    _, _, rms, e1, e2 = fit_plane(tips)
    planar = rms < planarity_tol
    x = tips @ e1; y = tips @ e2
    ang = np.mod(np.degrees(np.arctan2(y, x)), 360.0)
    ang.sort()
    gaps = np.diff(np.r_[ang, ang[0] + 360.0])
    gap_err = float(np.std(gaps - 72.0))
    spacing_ok = gap_err <= gap_tol_deg
    return planar, spacing_ok, gaps.tolist(), rms, gap_err



### classification
def classify(obj,
             linear_threshold_deg=160.0,
             planarity_tol_3=0.15,
             planarity_tol_4=0.15,
             planarity_tol_5=0.15,
             angle_tol_deg=20.0,
             outward_dot_min=0.2,
             exterior_dirs=256):
    conns = obj.get("connectors", [])
    n = len(conns)
    vecs = connector_vectors_filtered(obj, outward_dot_min=outward_dot_min, n_dirs=exterior_dirs)
    n = len(vecs)
    cls = (obj.get("class") or "").lower()

    if n == 2:                                                                                  ### Treat any 2-connector as 2-ditopic; subtype by angle
        a = angle_deg(vecs[0], vecs[1])
        subtype = "2-linear" if a >= linear_threshold_deg else "2-bent"
        return {"gbu_type": "2-ditopic", "gbu_subtype": subtype, "angle_deg": a}

    if n == 3:
        tips = np.asarray([u(v) for v in vecs], float)
        _, _, rms, *_ = fit_plane(tips)
        triple = abs(np.dot(tips[0], np.cross(tips[1], tips[2])))
        if rms < planarity_tol_3 and triple < 0.10:
            return {"gbu_type": "3-planar", "details": {"plane_rms": rms, "triple": triple}}
        return {"gbu_type": "3-pyramidal", "details": {"plane_rms": rms, "triple": triple}}

    if n == 4:
        tips = np.asarray([u(v) for v in vecs], float)
        _, _, rms, *_ = fit_plane(tips)
        ang = sorted(angle_deg(tips[i], tips[j]) for i in range(4) for j in range(i+1,4))
        near_180 = sum(1 for a in ang if abs(a - 180) <= angle_tol_deg)
        near_90  = sum(1 for a in ang if abs(a -  90) <= angle_tol_deg)
        if rms < planarity_tol_4 and near_180 >= 2 and near_90 >= 4:
            return {"gbu_type": "4-planar", "details": {"plane_rms": rms, "angles_deg": ang}}
        return {"gbu_type": "4-pyramidal", "details": {"plane_rms": rms, "angles_deg": ang}}

    if n == 5:
        planar, spacing_ok, gaps, rms, gap_err = pentagon_planarity_and_spacing(
            vecs, planarity_tol=planarity_tol_5, gap_tol_deg=15.0
        )
        if planar and spacing_ok:
            return {"gbu_type": "5-planar",
                    "details": {"plane_rms": rms, "gaps_deg": gaps, "gap_err_deg": gap_err}}
        return {"gbu_type": "5-pyramidal",
                "details": {"plane_rms": rms, "gaps_deg": gaps, "gap_err_deg": gap_err}}

    return {"gbu_type": f"{n}-connector", "note": "no rule for this count or all filtered."}



### CLI Input
def main():
    ap = argparse.ArgumentParser(description="Infer GBU from a CBU JSON and optionally write it back.")
    ap.add_argument("json_path", type=Path, help="Path to CBU JSON")
    ap.add_argument("--write", action="store_true",
                    help="Only overwrite existing keys 'gbu_type' and 'gbu_subtype'.")
    args = ap.parse_args()

    data = json.loads(args.json_path.read_text())
    res = classify(data)
    print(json.dumps(res, indent=2))

    if args.write:
        changed = []
        if "gbu_type" in data and "gbu_type" in res:
            data["gbu_type"] = res["gbu_type"]; changed.append("gbu_type")
        if "gbu_subtype" in data and "gbu_subtype" in res:
            data["gbu_subtype"] = res["gbu_subtype"]; changed.append("gbu_subtype")
        if changed:
            args.json_path.write_text(json.dumps(data, indent=2))
            print(f"Updated keys: {', '.join(changed)}")
        else:
            print("No keys updated (keys absent or not applicable).")



if __name__ == "__main__":
    main()
