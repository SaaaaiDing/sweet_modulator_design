#!/usr/bin/env python3

import csv
import numpy as np

from schrodinger.application.desmond.packages import analysis, traj_util

CMS = "binder1_200ns-out.cms"

# Binder1 = A
# TAS1R2 = B
ASL_A = '(chain.name "A")'
ASL_B = '(chain.name "B")'
ASL_AB = '((chain.name "A") OR (chain.name "B"))'

# For BSA, 1-ns sampling is enough for the trajectory trend:
# original trajectory = 0.1 ns/frame
STRIDE = 10

# Only use this as "post-equilibration" if the same criterion
# is applied consistently to all binders.
EQUIL_NS = 20.0


print("Reading trajectory...")
msys_model, cms_model, tr = traj_util.read_cms_and_traj(CMS)

tr_sub = tr[::STRIDE]

print(f"Original frames: {len(tr)}")
print(f"Analyzed frames: {len(tr_sub)}")
print(
    f"Time range: {tr_sub[0].time/1000:.2f} - "
    f"{tr_sub[-1].time/1000:.2f} ns"
)

# ------------------------------------------------------------
# Bound SASA:
# partner is present and therefore occludes the surface
# ------------------------------------------------------------

A_bound_analyzer = analysis.SolventAccessibleSurfaceArea(
    msys_model,
    cms_model,
    ASL_A
)

B_bound_analyzer = analysis.SolventAccessibleSurfaceArea(
    msys_model,
    cms_model,
    ASL_B
)

# ------------------------------------------------------------
# Free SASA:
# remove the opposite binding partner from the occlusion set
# while keeping exactly the same MD coordinates
# ------------------------------------------------------------

A_free_analyzer = analysis.SolventAccessibleSurfaceArea(
    msys_model,
    cms_model,
    ASL_A,
    exclude_asl=ASL_B
)

B_free_analyzer = analysis.SolventAccessibleSurfaceArea(
    msys_model,
    cms_model,
    ASL_B,
    exclude_asl=ASL_A
)

# Independent complex SASA for sanity checking
AB_analyzer = analysis.SolventAccessibleSurfaceArea(
    msys_model,
    cms_model,
    ASL_AB
)


print("Calculating Binder bound SASA...")
A_bound = np.asarray(
    analysis.analyze(tr_sub, A_bound_analyzer),
    dtype=float
)

print("Calculating TAS1R2 bound SASA...")
B_bound = np.asarray(
    analysis.analyze(tr_sub, B_bound_analyzer),
    dtype=float
)

print("Calculating Binder free SASA...")
A_free = np.asarray(
    analysis.analyze(tr_sub, A_free_analyzer),
    dtype=float
)

print("Calculating TAS1R2 free SASA...")
B_free = np.asarray(
    analysis.analyze(tr_sub, B_free_analyzer),
    dtype=float
)

print("Calculating complex SASA...")
AB = np.asarray(
    analysis.analyze(tr_sub, AB_analyzer),
    dtype=float
)


# ------------------------------------------------------------
# Buried surface area
# ------------------------------------------------------------

buried_A = A_free - A_bound
buried_B = B_free - B_bound

BSA_total = buried_A + buried_B

# Some papers call this interface area
interface_area_half = BSA_total / 2.0

# Independent formula:
# SASA_A_free + SASA_B_free - SASA_AB
BSA_from_complex = A_free + B_free - AB

times_ns = np.asarray(
    [fr.time / 1000.0 for fr in tr_sub],
    dtype=float
)


# ------------------------------------------------------------
# Save trajectory data
# ------------------------------------------------------------

outfile = "binder1_interface_BSA_vs_time.csv"

with open(outfile, "w", newline="") as f:

    w = csv.writer(f)

    w.writerow([
        "time_ns",
        "binder_free_SASA_A2",
        "binder_bound_SASA_A2",
        "receptor_free_SASA_A2",
        "receptor_bound_SASA_A2",
        "complex_SASA_A2",
        "binder_buried_SASA_A2",
        "receptor_buried_SASA_A2",
        "total_BSA_A2",
        "interface_area_half_BSA_A2",
        "BSA_complex_formula_A2"
    ])

    for i in range(len(times_ns)):
        w.writerow([
            times_ns[i],
            A_free[i],
            A_bound[i],
            B_free[i],
            B_bound[i],
            AB[i],
            buried_A[i],
            buried_B[i],
            BSA_total[i],
            interface_area_half[i],
            BSA_from_complex[i]
        ])


def summarize(label, mask):

    values = BSA_total[mask]
    half = interface_area_half[mask]

    print(f"\n{label}")
    print(
        f"Total buried SASA = "
        f"{np.mean(values):.1f} ? {np.std(values, ddof=1):.1f} ??"
    )
    print(
        f"Median total buried SASA = "
        f"{np.median(values):.1f} ??"
    )
    print(
        f"Range = "
        f"{np.min(values):.1f} - {np.max(values):.1f} ??"
    )

    print(
        f"Interface area (= BSA/2) = "
        f"{np.mean(half):.1f} ? {np.std(half, ddof=1):.1f} ??"
    )


summarize(
    "Full trajectory",
    np.ones(len(times_ns), dtype=bool)
)

posteq_mask = times_ns >= EQUIL_NS

if np.any(posteq_mask):
    summarize(
        f"{EQUIL_NS:.0f}-{times_ns[-1]:.0f} ns",
        posteq_mask
    )


# Sanity check:
diff = BSA_total - BSA_from_complex

print("\nConsistency check")
print(
    "Mean difference between two BSA formulations: "
    f"{np.mean(diff):.2f} ??"
)

print(f"\nDONE")
print(f"Output: {outfile}")

