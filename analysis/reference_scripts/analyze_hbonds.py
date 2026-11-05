#!/usr/bin/env python3

import csv
from collections import Counter, defaultdict

from schrodinger.application.desmond.packages import analysis, traj_util

CMS = "binder1_200ns-out.cms"
BINDER_CHAIN = "A"
RECEPTOR_CHAIN = "B"

print("Reading CMS and trajectory...")
msys_model, cms_model, tr = traj_util.read_cms_and_traj(CMS)

print(f"Frames: {len(tr)}")
print(f"Time range: {tr[0].time/1000:.2f} - {tr[-1].time/1000:.2f} ns")

# Chain atom IDs (AIDs)
binder_aids = list(cms_model.chain[BINDER_CHAIN].getAtomIndices())
receptor_aids = list(cms_model.chain[RECEPTOR_CHAIN].getAtomIndices())

print(f"Binder chain {BINDER_CHAIN}: {len(binder_aids)} atoms")
print(f"Receptor chain {RECEPTOR_CHAIN}: {len(receptor_aids)} atoms")

# Intermolecular H-bonds only: Binder A <-> Receptor B
hb = analysis.HydrogenBondFinder(
    msys_model,
    cms_model,
    binder_aids,
    receptor_aids
)

print("Calculating H-bonds...")
results = analysis.analyze(tr, hb)

# ---------- helper ----------
def atom_info(aid):
    a = cms_model.atom[aid]
    return {
        "aid": aid,
        "chain": a.chain.strip(),
        "resname": a.pdbres.strip(),
        "resnum": a.resnum,
        "atom": a.pdbname.strip(),
    }

def resid_key(info):
    return (info["chain"], info["resname"], info["resnum"])

def resid_label(info):
    return f'{info["chain"]}:{info["resname"]}{info["resnum"]}'

# ---------- 1. H-bond count vs time ----------
with open("binder1_hbond_count_vs_time.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["frame", "time_ns", "hbond_count"])
    for i, bonds in enumerate(results):
        w.writerow([i, tr[i].time / 1000.0, len(bonds)])

# ---------- occupancies ----------
atom_pair_frames = Counter()
res_pair_frames = Counter()

# Save detailed frame-by-frame bonds as well
with open("binder1_hbond_all_events.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow([
        "frame", "time_ns",
        "donor_chain", "donor_resname", "donor_resnum", "donor_atom",
        "acceptor_chain", "acceptor_resname", "acceptor_resnum", "acceptor_atom"
    ])

    for iframe, bonds in enumerate(results):

        # use sets so one pair contributes at most once per frame to occupancy
        atom_pairs_this_frame = set()
        residue_pairs_this_frame = set()

        for acceptor_aid, donor_aid in bonds:
            acc = atom_info(acceptor_aid)
            don = atom_info(donor_aid)

            w.writerow([
                iframe, tr[iframe].time / 1000.0,
                don["chain"], don["resname"], don["resnum"], don["atom"],
                acc["chain"], acc["resname"], acc["resnum"], acc["atom"]
            ])

            atom_pair = (
                don["chain"], don["resname"], don["resnum"], don["atom"],
                acc["chain"], acc["resname"], acc["resnum"], acc["atom"]
            )
            atom_pairs_this_frame.add(atom_pair)

            # Normalize residue pair order as Binder -> Receptor
            if don["chain"] == BINDER_CHAIN:
                binder_res = resid_key(don)
                receptor_res = resid_key(acc)
            else:
                binder_res = resid_key(acc)
                receptor_res = resid_key(don)

            residue_pairs_this_frame.add((binder_res, receptor_res))

        for p in atom_pairs_this_frame:
            atom_pair_frames[p] += 1

        for p in residue_pairs_this_frame:
            res_pair_frames[p] += 1

nframes = len(tr)

# ---------- 2. atom-pair occupancy ----------
with open("binder1_hbond_atom_pair_occupancy.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow([
        "donor_chain", "donor_resname", "donor_resnum", "donor_atom",
        "acceptor_chain", "acceptor_resname", "acceptor_resnum", "acceptor_atom",
        "frames_present", "occupancy_percent"
    ])

    for pair, n in atom_pair_frames.most_common():
        w.writerow(list(pair) + [n, 100.0 * n / nframes])

# ---------- 3. residue-pair occupancy ----------
with open("binder1_hbond_residue_pair_occupancy.csv", "w", newline="") as f:
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
print(f"Mean intermolecular H-bonds/frame: {sum(counts)/len(counts):.2f}")
print(f"Minimum: {min(counts)}")
print(f"Maximum: {max(counts)}")

print("\nTop residue-pair H-bonds:")
for (bres, rres), n in res_pair_frames.most_common(15):
    print(
        f"{bres[0]}:{bres[1]}{bres[2]}  <->  "
        f"{rres[0]}:{rres[1]}{rres[2]} : "
        f"{100*n/nframes:.1f}%"
    )
