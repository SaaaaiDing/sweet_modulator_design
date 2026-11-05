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
print(f"Binder chain {BINDER_CHAIN}: {len(binder_aids)} atoms")
print(f"Receptor chain {RECEPTOR_CHAIN}: {len(receptor_aids)} atoms")

sb = analysis.SaltBridgeFinder(
    msys_model,
    cms_model,
    binder_aids,
    receptor_aids
)

print("Calculating salt bridges...")
results = analysis.analyze(tr, sb)

def atom_info(aid):
    a = cms_model.atom[aid]
    return {
        "chain": a.chain.strip(),
        "resname": a.pdbres.strip(),
        "resnum": a.resnum,
        "atom": a.pdbname.strip(),
    }

res_pair_frames = Counter()
atom_pair_frames = Counter()

# 1. Count vs time
with open("binder1_saltbridge_count_vs_time.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["frame", "time_ns", "salt_bridge_count"])
    for i, bridges in enumerate(results):
        w.writerow([i, tr[i].time / 1000.0, len(bridges)])

# 2. All events
with open("binder1_saltbridge_all_events.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow([
        "frame", "time_ns",
        "anion_chain", "anion_resname", "anion_resnum", "anion_atom",
        "cation_chain", "cation_resname", "cation_resnum", "cation_atom"
    ])

    for iframe, bridges in enumerate(results):

        residue_pairs_this_frame = set()
        atom_pairs_this_frame = set()

        for anion_aid, cation_aid in bridges:
            an = atom_info(anion_aid)
            cat = atom_info(cation_aid)

            w.writerow([
                iframe, tr[iframe].time / 1000.0,
                an["chain"], an["resname"], an["resnum"], an["atom"],
                cat["chain"], cat["resname"], cat["resnum"], cat["atom"]
            ])

            atom_pairs_this_frame.add((
                an["chain"], an["resname"], an["resnum"], an["atom"],
                cat["chain"], cat["resname"], cat["resnum"], cat["atom"]
            ))

            # Normalize output as Binder residue -> Receptor residue
            if an["chain"] == BINDER_CHAIN:
                binder_res = (an["chain"], an["resname"], an["resnum"])
                receptor_res = (cat["chain"], cat["resname"], cat["resnum"])
            else:
                binder_res = (cat["chain"], cat["resname"], cat["resnum"])
                receptor_res = (an["chain"], an["resname"], an["resnum"])

            residue_pairs_this_frame.add((binder_res, receptor_res))

        for p in atom_pairs_this_frame:
            atom_pair_frames[p] += 1

        for p in residue_pairs_this_frame:
            res_pair_frames[p] += 1

nframes = len(tr)

# 3. Residue-pair occupancy
with open("binder1_saltbridge_residue_pair_occupancy.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow([
        "binder_chain", "binder_resname", "binder_resnum",
        "receptor_chain", "receptor_resname", "receptor_resnum",
        "frames_present", "occupancy_percent"
    ])

    for (bres, rres), n in res_pair_frames.most_common():
        w.writerow([
            bres[0], bres[1], bres[2],
            rres[0], rres[1], rres[2],
            n, 100.0 * n / nframes
        ])

counts = [len(x) for x in results]

print("\nDONE")
print(f"Mean intermolecular salt bridges/frame: {sum(counts)/len(counts):.2f}")
print(f"Minimum: {min(counts)}")
print(f"Maximum: {max(counts)}")

print("\nTop residue-pair salt bridges:")
for (bres, rres), n in res_pair_frames.most_common(15):
    print(
        f"{bres[0]}:{bres[1]}{bres[2]}  <->  "
        f"{rres[0]}:{rres[1]}{rres[2]} : "
        f"{100*n/nframes:.1f}%"
    )
