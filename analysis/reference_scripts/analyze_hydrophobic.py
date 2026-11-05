#!/usr/bin/env python3

import csv
from collections import Counter

from schrodinger.application.desmond.packages import analysis, traj_util

CMS = "binder1_200ns-out.cms"
BINDER_CHAIN = "A"
RECEPTOR_CHAIN = "B"

print("Reading CMS and trajectory...")
msys_model, cms_model, tr = traj_util.read_cms_and_traj(CMS)

binder_aids = list(cms_model.chain[BINDER_CHAIN].getAtomIndices())
receptor_aids = list(cms_model.chain[RECEPTOR_CHAIN].getAtomIndices())

print(f"Frames: {len(tr)}")
print(f"Time range: {tr[0].time/1000:.2f} - {tr[-1].time/1000:.2f} ns")
print(f"Binder chain {BINDER_CHAIN}: {len(binder_aids)} atoms")
print(f"Receptor chain {RECEPTOR_CHAIN}: {len(receptor_aids)} atoms")

# General hydrophobic interaction finder for two atom sets
hyd = analysis.HydrophobicInterFinder(
    msys_model,
    cms_model,
    binder_aids,
    receptor_aids
)

print("Calculating hydrophobic interactions...")
results = analysis.analyze(tr, hyd)

def atom_info(aid):
    a = cms_model.atom[aid]
    return {
        "aid": aid,
        "chain": a.chain.strip(),
        "resname": a.pdbres.strip(),
        "resnum": a.resnum,
        "atom": a.pdbname.strip(),
    }

atom_pair_frames = Counter()
res_pair_frames = Counter()

# ------------------------------------------------------
# 1. Hydrophobic interaction count vs time
# ------------------------------------------------------
with open("binder1_hydrophobic_count_vs_time.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["frame", "time_ns", "hydrophobic_contact_count"])

    for i, contacts in enumerate(results):
        w.writerow([
            i,
            tr[i].time / 1000.0,
            len(contacts)
        ])

# ------------------------------------------------------
# 2. All atom-level hydrophobic contact events
# ------------------------------------------------------
with open("binder1_hydrophobic_all_events.csv", "w", newline="") as f:

    w = csv.writer(f)

    w.writerow([
        "frame",
        "time_ns",
        "binder_chain",
        "binder_resname",
        "binder_resnum",
        "binder_atom",
        "receptor_chain",
        "receptor_resname",
        "receptor_resnum",
        "receptor_atom"
    ])

    for iframe, contacts in enumerate(results):

        atom_pairs_this_frame = set()
        residue_pairs_this_frame = set()

        for pair in contacts:

            aid1 = int(pair[0])
            aid2 = int(pair[1])

            x1 = atom_info(aid1)
            x2 = atom_info(aid2)

            # Normalize output:
            # Binder always first, receptor always second
            if x1["chain"] == BINDER_CHAIN:
                binder = x1
                receptor = x2
            else:
                binder = x2
                receptor = x1

            w.writerow([
                iframe,
                tr[iframe].time / 1000.0,

                binder["chain"],
                binder["resname"],
                binder["resnum"],
                binder["atom"],

                receptor["chain"],
                receptor["resname"],
                receptor["resnum"],
                receptor["atom"]
            ])

            atom_pair = (
                binder["chain"],
                binder["resname"],
                binder["resnum"],
                binder["atom"],

                receptor["chain"],
                receptor["resname"],
                receptor["resnum"],
                receptor["atom"]
            )

            residue_pair = (
                (
                    binder["chain"],
                    binder["resname"],
                    binder["resnum"]
                ),
                (
                    receptor["chain"],
                    receptor["resname"],
                    receptor["resnum"]
                )
            )

            atom_pairs_this_frame.add(atom_pair)
            residue_pairs_this_frame.add(residue_pair)

        # Count at most once per frame for occupancy
        for p in atom_pairs_this_frame:
            atom_pair_frames[p] += 1

        for p in residue_pairs_this_frame:
            res_pair_frames[p] += 1


nframes = len(tr)

# ------------------------------------------------------
# 3. Atom-pair occupancy
# ------------------------------------------------------
with open(
    "binder1_hydrophobic_atom_pair_occupancy.csv",
    "w",
    newline=""
) as f:

    w = csv.writer(f)

    w.writerow([
        "binder_chain",
        "binder_resname",
        "binder_resnum",
        "binder_atom",

        "receptor_chain",
        "receptor_resname",
        "receptor_resnum",
        "receptor_atom",

        "frames_present",
        "occupancy_percent"
    ])

    for pair, n in atom_pair_frames.most_common():

        w.writerow(
            list(pair)
            +
            [
                n,
                100.0 * n / nframes
            ]
        )

# ------------------------------------------------------
# 4. Residue-pair occupancy
# ------------------------------------------------------
with open(
    "binder1_hydrophobic_residue_pair_occupancy.csv",
    "w",
    newline=""
) as f:

    w = csv.writer(f)

    w.writerow([
        "binder_chain",
        "binder_resname",
        "binder_resnum",

        "receptor_chain",
        "receptor_resname",
        "receptor_resnum",

        "frames_present",
        "occupancy_percent"
    ])

    for (bres, rres), n in res_pair_frames.most_common():

        w.writerow([
            bres[0],
            bres[1],
            bres[2],

            rres[0],
            rres[1],
            rres[2],

            n,
            100.0 * n / nframes
        ])


counts = [len(x) for x in results]

print("\nDONE")
print(
    f"Mean intermolecular hydrophobic atom contacts/frame: "
    f"{sum(counts)/len(counts):.2f}"
)
print(f"Minimum: {min(counts)}")
print(f"Maximum: {max(counts)}")

print("\nTop residue-pair hydrophobic contacts:")

for (bres, rres), n in res_pair_frames.most_common(20):

    print(
        f"{bres[0]}:{bres[1]}{bres[2]}"
        f"  <->  "
        f"{rres[0]}:{rres[1]}{rres[2]}"
        f" : {100*n/nframes:.1f}%"
    )

