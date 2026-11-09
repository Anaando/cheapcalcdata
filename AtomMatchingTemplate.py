import argparse
import itertools
from collections import Counter

import numpy as np
from scipy.optimize import linear_sum_assignment


def read_xyz(filename):
    lines = open(filename, "r", encoding="utf-8").read().splitlines()
    if len(lines) < 2:
        raise ValueError(f"{filename}: not a valid XYZ file")
    n = int(lines[0].strip())
    if len(lines) < n + 2:
        raise ValueError(f"{filename}: expected {n} atom lines")
    symbols, xyz = [], []
    for line_no, line in enumerate(lines[2:2+n], start=3):
        fields = line.split()
        if len(fields) < 4:
            raise ValueError(f"{filename}:{line_no}: invalid atom line")
        symbols.append(fields[0])
        xyz.append([float(x) for x in fields[1:4]])
    return symbols, np.asarray(xyz, dtype=float), lines[1]


def kabsch(moving, target):
    """Return rotation and translation mapping moving onto target."""
    cm = moving.mean(axis=0)
    ct = target.mean(axis=0)
    a = moving - cm
    b = target - ct
    u, _, vt = np.linalg.svd(a.T @ b)
    d = np.linalg.det(u @ vt)
    correction = np.diag([1.0, 1.0, d])
    rotation = u @ correction @ vt
    translation = ct - cm @ rotation
    return rotation, translation


def pca_frame(coords):
    centered = coords - coords.mean(axis=0)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    return vt.T


def assignment(reference_symbols, reference_xyz, mobile_symbols, mobile_xyz):
    """Element-constrained Hungarian assignment for a given aligned geometry."""
    mapping = np.empty(len(reference_symbols), dtype=int)
    total = 0.0
    for element in sorted(set(reference_symbols)):
        ref_ids = np.array([i for i, s in enumerate(reference_symbols) if s == element])
        mob_ids = np.array([i for i, s in enumerate(mobile_symbols) if s == element])
        delta = reference_xyz[ref_ids, None, :] - mobile_xyz[None, mob_ids, :]
        cost = np.sum(delta * delta, axis=2)
        rows, cols = linear_sum_assignment(cost)
        mapping[ref_ids[rows]] = mob_ids[cols]
        total += cost[rows, cols].sum()
    return mapping, total


def fit_and_assign(ref_symbols, ref_xyz, mob_symbols, mob_xyz, max_iter=100):
    """Try PCA axis orientations, then alternate assignment and rigid alignment."""
    ref_frame = pca_frame(ref_xyz)
    mob_frame = pca_frame(mob_xyz)
    best = None

    # Enumerate axis permutations and sign flips to avoid depending on PCA signs.
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((-1.0, 1.0), repeat=3):
            transform = mob_frame[:, perm] @ np.diag(signs) @ ref_frame.T
            if np.linalg.det(transform) < 0:
                continue
            trial = (mob_xyz - mob_xyz.mean(axis=0)) @ transform + ref_xyz.mean(axis=0)

            previous = None
            for _ in range(max_iter):
                mapping, _ = assignment(ref_symbols, ref_xyz, mob_symbols, trial)
                # Correspondences: mobile atom mapping[i] should sit on reference atom i.
                rotation, translation = kabsch(mob_xyz[mapping], ref_xyz)
                aligned = mob_xyz @ rotation + translation
                if previous is not None and np.array_equal(mapping, previous):
                    trial = aligned
                    break
                previous = mapping.copy()
                trial = aligned

            mapping, cost = assignment(ref_symbols, ref_xyz, mob_symbols, trial)
            # Refit once to final assignment and calculate final RMSD.
            rotation, translation = kabsch(mob_xyz[mapping], ref_xyz)
            aligned = mob_xyz @ rotation + translation
            rmsd = np.sqrt(np.mean(np.sum((aligned[mapping] - ref_xyz) ** 2, axis=1)))
            if best is None or rmsd < best[0]:
                best = (rmsd, mapping.copy(), aligned.copy())

    return best


def main():
    parser = argparse.ArgumentParser(description="Reorder XYZ B to match atom order in XYZ A.")
    parser.add_argument("reference", help="Reference XYZ; its atom order is retained")
    parser.add_argument("mobile", help="XYZ to reorder")
    parser.add_argument("--output", default="mobile_reordered.xyz")
    parser.add_argument("--mapping", default="atom_mapping.csv")
    parser.add_argument("--warning-rmsd", type=float, default=1.0,
                        help="Warn if aligned RMSD exceeds this value in Angstrom (default: 1.0)")
    args = parser.parse_args()

    ref_s, ref_x, ref_comment = read_xyz(args.reference)
    mob_s, mob_x, mob_comment = read_xyz(args.mobile)

    if len(ref_s) != len(mob_s):
        raise SystemExit(
            f"STOP: atom counts differ ({len(ref_s)} vs {len(mob_s)}). "
            "A standard NEB endpoint pair needs the same atoms."
        )
    if Counter(ref_s) != Counter(mob_s):
        raise SystemExit("STOP: elemental compositions differ; cannot make a one-to-one mapping.")

    rmsd, mapping, aligned_mobile = fit_and_assign(ref_s, ref_x, mob_s, mob_x)
    reordered_xyz = mob_x[mapping]

    with open(args.output, "w", encoding="utf-8") as f:
        f.write(f"{len(mob_s)}\n")
        f.write(f"Reordered from {args.mobile}; order follows {args.reference}\n")
        for symbol, xyz in zip(ref_s, reordered_xyz):
            f.write(f"{symbol:2s} {xyz[0]: .10f} {xyz[1]: .10f} {xyz[2]: .10f}\n")

    with open(args.mapping, "w", encoding="utf-8") as f:
        f.write("reference_atom,mobile_original_atom,element,aligned_distance_A\n")
        for i, j in enumerate(mapping):
            distance = np.linalg.norm(aligned_mobile[j] - ref_x[i])
            f.write(f"{i+1},{j+1},{ref_s[i]},{distance:.6f}\n")

    print(f"Reference atoms: {len(ref_s)}")
    print("Elemental compositions: MATCH")
    print(f"Best geometry-alignment RMSD: {rmsd:.4f} Å")
    print(f"Reordered XYZ written to: {args.output}")
    print(f"Atom mapping written to: {args.mapping}")
    if rmsd > args.warning_rmsd:
        print("WARNING: RMSD is high. Do NOT assume this mapping is chemically correct.")
        print("Inspect the mapping CSV and structures; manual anchors/chemical checks may be needed.")
    else:
        print("NOTE: A low RMSD is encouraging, but still verify chemically before NEB.")


if __name__ == "__main__":
    main()
