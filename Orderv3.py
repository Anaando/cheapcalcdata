import argparse
import csv
import itertools
from collections import Counter

import numpy as np
from scipy.optimize import linear_sum_assignment


# Approximate covalent radii in Angstrom.
RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66,
    "F": 0.57, "P": 1.07, "S": 1.05, "Cl": 1.02,
    "Br": 1.20, "I": 1.39, "Zn": 1.22, "Mg": 1.41,
    "B": 0.84, "Si": 1.11, "Fe": 1.24, "Cu": 1.32,
    "Ni": 1.24, "Co": 1.26, "Mn": 1.39,
    "Pd": 1.39, "Pt": 1.36, "Ru": 1.46
}


def read_xyz(filename):

    with open(filename, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()

    if len(lines) < 2:
        raise ValueError(f"Invalid XYZ file: {filename}")

    n = int(lines[0].strip())

    if len(lines) < n + 2:
        raise ValueError(f"Incorrect atom count in {filename}")

    symbols = []
    coords = []

    for line in lines[2:2+n]:

        fields = line.split()

        if len(fields) < 4:
            raise ValueError("Invalid XYZ atom line")

        symbols.append(fields[0])
        coords.append([float(x) for x in fields[1:4]])

    return symbols, np.asarray(coords), lines[1]


def rotate_x_180(coords):
    """
    Rotate all coordinates by 180 degrees around the X axis.
    """

    rotation = np.array([
        [1.0,  0.0,  0.0],
        [0.0, -1.0,  0.0],
        [0.0,  0.0, -1.0]
    ])

    return coords @ rotation


def kabsch(moving, target):
    """
    Rigidly align moving coordinates onto target coordinates.
    """

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

    _, _, vt = np.linalg.svd(
        centered, full_matrices=False
    )

    return vt.T


def build_connectivity(symbols, coords, scale=1.20):
    """
    Infer approximate bonds using covalent radii and distances.
    """

    n = len(symbols)
    adjacency = np.zeros((n, n), dtype=bool)

    for i in range(n):

        for j in range(i + 1, n):

            ri = RADII.get(symbols[i], 0.77)
            rj = RADII.get(symbols[j], 0.77)

            cutoff = scale * (ri + rj)

            distance = np.linalg.norm(
                coords[i] - coords[j]
            )

            if 0.1 < distance <= cutoff:

                adjacency[i, j] = True
                adjacency[j, i] = True

    return adjacency


def graph_descriptors(symbols, adjacency):
    """
    Describe each atom by the elements and number
    of its approximate neighbouring atoms.
    """

    elements = sorted(set(symbols))
    descriptors = []

    for i in range(len(symbols)):

        neighbours = np.where(adjacency[i])[0]

        counts = [
            sum(symbols[j] == element for j in neighbours)
            for element in elements
        ]

        descriptors.append(counts + [len(neighbours)])

    return np.asarray(descriptors, dtype=float)


def read_anchors(filename, ref_symbols, mob_symbols):
    """
    Read manual atom correspondences.
    Atom numbers in CSV are one-based.
    """

    anchors = {}

    if filename is None:
        return anchors

    with open(filename, newline="", encoding="utf-8") as f:

        reader = csv.DictReader(f)

        required = {"reference_atom", "mobile_atom"}

        if not required.issubset(reader.fieldnames or []):
            raise ValueError(
                "Anchor CSV needs reference_atom and mobile_atom columns"
            )

        for row in reader:

            i = int(row["reference_atom"]) - 1
            j = int(row["mobile_atom"]) - 1

            if not (0 <= i < len(ref_symbols)):
                raise ValueError(
                    f"Invalid reference anchor: {i+1}"
                )

            if not (0 <= j < len(mob_symbols)):
                raise ValueError(
                    f"Invalid mobile anchor: {j+1}"
                )

            if ref_symbols[i] != mob_symbols[j]:
                raise ValueError(
                    f"Anchor element mismatch: reference {i+1}, mobile {j+1}"
                )

            if i in anchors.values() or j in anchors:
                raise ValueError("Duplicate anchor assignment")

            anchors[j] = i

    return anchors


def assignment(ref_symbols, mob_symbols, ref_graph, mob_graph,
               ref_coords, aligned_mobile, anchors, graph_weight):

    n = len(ref_symbols)
    mapping = np.full(n, -1, dtype=int)

    # Apply fixed manual correspondences.
    for mobile_id, reference_id in anchors.items():
        mapping[reference_id] = mobile_id

    # Fixed reference IDs are anchor VALUES.
    fixed_ref = set(anchors.values())
    fixed_mob = set(anchors.keys())

    total_cost = 0.0

    for element in sorted(set(ref_symbols)):

        ref_ids = [
            i for i, s in enumerate(ref_symbols)
            if s == element and i not in fixed_ref
        ]

        mob_ids = [
            j for j, s in enumerate(mob_symbols)
            if s == element and j not in fixed_mob
        ]

        if len(ref_ids) != len(mob_ids):
            raise ValueError(
                "Unequal number of atoms of an element"
            )

        if not ref_ids:
            continue

        # Geometrical matching cost.
        delta = (
            ref_coords[ref_ids, None, :]
            - aligned_mobile[None, mob_ids, :]
        )

        geometric_cost = np.sum(delta ** 2, axis=2)

        # Chemical-environment matching cost.
        graph_delta = (
            ref_graph[ref_ids, None, :]
            - mob_graph[None, mob_ids, :]
        )

        graph_cost = np.sum(np.abs(graph_delta), axis=2)

        cost = geometric_cost + graph_weight * graph_cost

        rows, cols = linear_sum_assignment(cost)

        for r, c in zip(rows, cols):
            mapping[ref_ids[r]] = mob_ids[c]

        total_cost += cost[rows, cols].sum()

    if np.any(mapping < 0):
        raise ValueError("Some reference atoms were not matched")

    return mapping, total_cost


def fit_and_assign(ref_symbols, ref_xyz,
                   mob_symbols, mob_xyz,
                   ref_graph, mob_graph,
                   anchors, graph_weight,
                   max_iter=100):

    best = None

    # Generate initial PCA orientations.
    ref_frame = pca_frame(ref_xyz)
    mob_frame = pca_frame(mob_xyz)

    initial_alignments = []

    for perm in itertools.permutations(range(3)):

        for signs in itertools.product((-1.0, 1.0), repeat=3):

            transform = (
                mob_frame[:, perm]
                @ np.diag(signs)
                @ ref_frame.T
            )

            if np.linalg.det(transform) < 0:
                continue

            aligned = (
                (mob_xyz - mob_xyz.mean(axis=0))
                @ transform
                + ref_xyz.mean(axis=0)
            )

            initial_alignments.append(aligned)

    # Iterative assignment and alignment.
    for initial in initial_alignments:

        trial = initial.copy()
        previous = None

        for _ in range(max_iter):

            mapping, _ = assignment(
                ref_symbols, mob_symbols,
                ref_graph, mob_graph,
                ref_xyz, trial,
                anchors, graph_weight
            )

            # Exclude manually anchored atoms from alignment fitting.
            fit_ids = np.array([
                i for i in range(len(ref_symbols))
                if i not in set(anchors.values())
            ], dtype=int)

            if len(fit_ids) < 3:
                raise ValueError(
                    "Not enough unanchored atoms for alignment"
                )

            rotation, translation = kabsch(
                mob_xyz[mapping[fit_ids]],
                ref_xyz[fit_ids]
            )

            aligned = mob_xyz @ rotation + translation

            if previous is not None and np.array_equal(
                mapping, previous
            ):
                trial = aligned
                break

            previous = mapping.copy()
            trial = aligned

        mapping, cost = assignment(
            ref_symbols, mob_symbols,
            ref_graph, mob_graph,
            ref_xyz, trial,
            anchors, graph_weight
        )

        fit_ids = np.array([
            i for i in range(len(ref_symbols))
            if i not in set(anchors.values())
        ], dtype=int)

        rotation, translation = kabsch(
            mob_xyz[mapping[fit_ids]],
            ref_xyz[fit_ids]
        )

        aligned = mob_xyz @ rotation + translation

        distances = np.linalg.norm(
            aligned[mapping] - ref_xyz, axis=1
        )

        rmsd = np.sqrt(np.mean(distances ** 2))

        if best is None or rmsd < best[0]:

            best = (
                rmsd,
                mapping.copy(),
                aligned.copy(),
                distances.copy()
            )

    return best


def write_outputs(args, ref_symbols, mapping,
                  aligned_mobile, distances, rmsd):

    # Write reordered and aligned coordinates.
    reordered = aligned_mobile[mapping]

    with open(args.output, "w", encoding="utf-8") as f:

        f.write(f"{len(mapping)}\n")
        f.write(
            f"Reordered from {args.mobile}; "
            f"order follows {args.reference}; "
            f"aligned coordinates\n"
        )

        for symbol, xyz in zip(ref_symbols, reordered):

            f.write(
                f"{symbol:2s} "
                f"{xyz[0]: .10f} "
                f"{xyz[1]: .10f} "
                f"{xyz[2]: .10f}\n"
            )

    # Write mapping CSV.
    with open(args.mapping, "w", newline="", encoding="utf-8") as f:

        writer = csv.writer(f)

        writer.writerow([
            "reference_atom",
            "mobile_original_atom",
            "element",
            "aligned_distance_A"
        ])

        for i, j in enumerate(mapping):

            writer.writerow([
                i + 1,
                j + 1,
                ref_symbols[i],
                f"{distances[i]:.6f}"
            ])

    # Write validation report.
    with open(args.validation, "w", encoding="utf-8") as f:

        f.write("ATOM MATCHING VALIDATION\n")
        f.write("========================\n\n")

        f.write(f"Reference: {args.reference}\n")
        f.write(f"Mobile: {args.mobile}\n")
        f.write(f"Number of atoms: {len(mapping)}\n")
        f.write(f"Overall RMSD: {rmsd:.4f} Angstrom\n")

        f.write(
            f"Maximum individual distance: "
            f"{distances.max():.4f} Angstrom\n\n"
        )

        f.write("Atoms with large individual distances:\n")

        flagged = np.where(
            distances > args.warning_distance
        )[0]

        if len(flagged) == 0:
            f.write("None\n")

        else:

            for i in flagged:

                f.write(
                    f"Reference atom {i+1}: "
                    f"mobile atom {mapping[i]+1}, "
                    f"distance {distances[i]:.4f} Angstrom\n"
                )

        if rmsd > args.warning_rmsd:
            f.write("\nWARNING: Overall RMSD is high.\n")

        if len(flagged) > 0:
            f.write(
                "\nWARNING: Some individual atom matches "
                "have large distances.\n"
            )

        f.write(
            "\nIMPORTANT: Geometry and inferred connectivity "
            "do not prove chemical identity. Inspect the mapping "
            "and molecular structures before NEB.\n"
        )


def main():

    parser = argparse.ArgumentParser(
        description="XYZ atom matcher with optional fixed rotation"
    )

    parser.add_argument("reference")
    parser.add_argument("mobile")

    parser.add_argument(
        "--output", default="mobile_reordered.xyz"
    )

    parser.add_argument(
        "--mapping", default="atom_mapping.csv"
    )

    parser.add_argument(
        "--validation", default="matching_validation.txt"
    )

    parser.add_argument("--anchors", default=None)

    parser.add_argument(
        "--graph-weight", type=float, default=0.5
    )

    parser.add_argument(
        "--connectivity-scale", type=float, default=1.20
    )

    parser.add_argument(
        "--warning-rmsd", type=float, default=1.0
    )

    parser.add_argument(
        "--warning-distance", type=float, default=2.0
    )

    parser.add_argument(
        "--rotate-x-180",
        action="store_true",
        help="Rotate the entire mobile structure 180 degrees around X"
    )

    args = parser.parse_args()

    ref_s, ref_x, _ = read_xyz(args.reference)
    mob_s, mob_x, _ = read_xyz(args.mobile)

    if len(ref_s) != len(mob_s):
        raise SystemExit(
            f"STOP: Different atom counts: "
            f"{len(ref_s)} and {len(mob_s)}"
        )

    if Counter(ref_s) != Counter(mob_s):
        raise SystemExit("STOP: Elemental compositions differ")

    anchors = read_anchors(
        args.anchors, ref_s, mob_s
    )

    print(f"Reference atoms: {len(ref_s)}")
    print(f"Mobile atoms: {len(mob_s)}")
    print(f"Manual anchors: {len(anchors)}")

    # Optional fixed 180 degree X rotation.
    if args.rotate_x_180:

        mob_x = rotate_x_180(mob_x)

        print(
            "Applied fixed 180 degree rotation around X axis."
        )

    else:

        print("No fixed rotation applied.")

    print("Building approximate molecular connectivity...")

    ref_adj = build_connectivity(
        ref_s, ref_x, args.connectivity_scale
    )

    mob_adj = build_connectivity(
        mob_s, mob_x, args.connectivity_scale
    )

    ref_graph = graph_descriptors(ref_s, ref_adj)
    mob_graph = graph_descriptors(mob_s, mob_adj)

    print("Aligning structures and matching atoms...")

    rmsd, mapping, aligned, distances = fit_and_assign(
        ref_s, ref_x,
        mob_s, mob_x,
        ref_graph, mob_graph,
        anchors, args.graph_weight
    )

    if len(set(mapping.tolist())) != len(mapping):
        raise SystemExit(
            "STOP: Duplicate atom assignments detected"
        )

    write_outputs(
        args, ref_s,
        mapping, aligned, distances, rmsd
    )

    print("\nMATCHING COMPLETE")
    print(f"Overall RMSD: {rmsd:.4f} Angstrom")
    print(
        f"Maximum atom distance: "
        f"{distances.max():.4f} Angstrom"
    )

    print(f"Reordered XYZ: {args.output}")
    print(f"Mapping CSV: {args.mapping}")
    print(f"Validation: {args.validation}")

    if rmsd > args.warning_rmsd or np.any(
        distances > args.warning_distance
    ):

        print("\nWARNING: Mapping requires manual inspection.")

    else:

        print("\nNOTE: Mapping passed numerical warning thresholds.")
        print("Chemical verification is still required.")


if __name__ == "__main__":
    main()
