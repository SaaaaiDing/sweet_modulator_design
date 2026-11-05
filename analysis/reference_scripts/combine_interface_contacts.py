#!/usr/bin/env python3

import csv
from collections import defaultdict, Counter

NFRAMES = 2001
BINDER_CHAIN = "A"
RECEPTOR_CHAIN = "B"

HBOND_FILE = "binder1_hbond_all_events.csv"
SALT_FILE = "binder1_saltbridge_all_events.csv"
HYDRO_FILE = "binder1_hydrophobic_all_events.csv"

# key:
# ((binder_resname, binder_resnum),
#  (receptor_resname, receptor_resnum))
#
# value:
# set(frames)

hbond_frames = defaultdict(set)
salt_frames = defaultdict(set)
hydro_frames = defaultdict(set)
all_frames = defaultdict(set)


def normalize_pair(chain1, resname1, resnum1,
                   chain2, resname2, resnum2):

    resnum1 = int(resnum1)
    resnum2 = int(resnum2)

    if chain1 == BINDER_CHAIN and chain2 == RECEPTOR_CHAIN:
        return (
            (resname1.strip(), resnum1),
            (resname2.strip(), resnum2)
        )

    if chain2 == BINDER_CHAIN and chain1 == RECEPTOR_CHAIN:
        return (
            (resname2.strip(), resnum2),
            (resname1.strip(), resnum1)
        )

    return None


# ==========================================================
# H-BONDS
# ==========================================================

with open(HBOND_FILE) as f:
    reader = csv.DictReader(f)

    for row in reader:
        frame = int(row["frame"])

        pair = normalize_pair(
            row["donor_chain"],
            row["donor_resname"],
            row["donor_resnum"],
            row["acceptor_chain"],
            row["acceptor_resname"],
            row["acceptor_resnum"]
        )

        if pair:
            hbond_frames[pair].add(frame)
            all_frames[pair].add(frame)


# ==========================================================
# SALT BRIDGES
# ==========================================================

with open(SALT_FILE) as f:
    reader = csv.DictReader(f)

    for row in reader:
        frame = int(row["frame"])

        pair = normalize_pair(
            row["anion_chain"],
            row["anion_resname"],
            row["anion_resnum"],
            row["cation_chain"],
            row["cation_resname"],
            row["cation_resnum"]
        )

        if pair:
            salt_frames[pair].add(frame)
            all_frames[pair].add(frame)


# ==========================================================
# HYDROPHOBIC CONTACTS
# ==========================================================

with open(HYDRO_FILE) as f:
    reader = csv.DictReader(f)

    for row in reader:
        frame = int(row["frame"])

        pair = normalize_pair(
            row["binder_chain"],
            row["binder_resname"],
            row["binder_resnum"],
            row["receptor_chain"],
            row["receptor_resname"],
            row["receptor_resnum"]
        )

        if pair:
            hydro_frames[pair].add(frame)
            all_frames[pair].add(frame)


# ==========================================================
# WRITE COMBINED TABLE
# ==========================================================

outfile = "binder1_combined_interface_occupancy.csv"

rows = []

for pair, frames in all_frames.items():

    binder, receptor = pair

    hb = len(hbond_frames.get(pair, set()))
    sb = len(salt_frames.get(pair, set()))
    hp = len(hydro_frames.get(pair, set()))
    total = len(frames)

    rows.append({
        "binder_resname": binder[0],
        "binder_resnum": binder[1],

        "receptor_resname": receptor[0],
        "receptor_resnum": receptor[1],

        "hbond_occ_percent": 100.0 * hb / NFRAMES,
        "saltbridge_occ_percent": 100.0 * sb / NFRAMES,
        "hydrophobic_occ_percent": 100.0 * hp / NFRAMES,

        # IMPORTANT:
        # union of frames, NOT sum of the three percentages
        "any_interaction_occ_percent": 100.0 * total / NFRAMES,

        "frames_with_any_interaction": total
    })


rows.sort(
    key=lambda x: x["any_interaction_occ_percent"],
    reverse=True
)


with open(outfile, "w", newline="") as f:

    fields = [
        "binder_resname",
        "binder_resnum",
        "receptor_resname",
        "receptor_resnum",
        "hbond_occ_percent",
        "saltbridge_occ_percent",
        "hydrophobic_occ_percent",
        "any_interaction_occ_percent",
        "frames_with_any_interaction"
    ]

    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)


print("\nDONE")
print(f"Output: {outfile}")
print(f"Total interacting residue pairs: {len(rows)}")

print("\nTop 30 persistent Binder1-TAS1R2 residue pairs:")
print(
    "Binder\tReceptor\tAny%\tHbond%\tSalt%\tHydrophobic%"
)

for r in rows[:30]:

    print(
        f"A:{r['binder_resname']}{r['binder_resnum']}\t"
        f"B:{r['receptor_resname']}{r['receptor_resnum']}\t"
        f"{r['any_interaction_occ_percent']:.1f}\t"
        f"{r['hbond_occ_percent']:.1f}\t"
        f"{r['saltbridge_occ_percent']:.1f}\t"
        f"{r['hydrophobic_occ_percent']:.1f}"
    )


# ==========================================================
# PER-RESIDUE IMPORTANCE SCORE
# ==========================================================

binder_residue_frames = defaultdict(set)
receptor_residue_frames = defaultdict(set)

for (binder, receptor), frames in all_frames.items():

    binder_residue_frames[binder].update(frames)
    receptor_residue_frames[receptor].update(frames)


with open(
    "binder1_binder_residue_interface_occupancy.csv",
    "w",
    newline=""
) as f:

    writer = csv.writer(f)

    writer.writerow([
        "binder_resname",
        "binder_resnum",
        "frames_interacting",
        "occupancy_percent"
    ])

    data = sorted(
        binder_residue_frames.items(),
        key=lambda x: len(x[1]),
        reverse=True
    )

    for residue, frames in data:

        writer.writerow([
            residue[0],
            residue[1],
            len(frames),
            100.0 * len(frames) / NFRAMES
        ])


with open(
    "binder1_receptor_residue_interface_occupancy.csv",
    "w",
    newline=""
) as f:

    writer = csv.writer(f)

    writer.writerow([
        "receptor_resname",
        "receptor_resnum",
        "frames_interacting",
        "occupancy_percent"
    ])

    data = sorted(
        receptor_residue_frames.items(),
        key=lambda x: len(x[1]),
        reverse=True
    )

    for residue, frames in data:

        writer.writerow([
            residue[0],
            residue[1],
            len(frames),
            100.0 * len(frames) / NFRAMES
        ])


print("\nAdditional outputs:")
print("binder1_binder_residue_interface_occupancy.csv")
print("binder1_receptor_residue_interface_occupancy.csv")
