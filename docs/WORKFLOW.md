# Computational workflow

The parameters below distinguish article-reported methods from settings directly visible in the supplied analysis scripts. This compact repository does not contain the original generation launch files, every intermediate model or the trajectories.

## 1. Receptor and backbone design — article-reported

Human TAS1R2: UniProt Q8TE23 (839 residues). Structural references were TAS1R2 chain B from the TAS1R2–TAS1R3 cryo-EM structure PDB 9NOR and the AlphaFold2 model of Q8TE23. The extracellular VFTD and CRD regions were the design target.

RFdiffusion generated 5,000 backbones across VFTD-directed and broader CRD-directed conditions, with one continuous 128-residue binder chain. The VFTD contact-bias positions were S40, K65, I67, L71, Y103, D142, N143, S144, E145, S165, A166, I167, S168, A187, H190, Y215, P277, D278, L279, Y282, E302, S303, W304, D307, V309, E382, R383 and V384. These are spatial design biases, not experimentally established energetic hotspots. Exact original seeds, checkpoint revisions, contig strings and CRD launch arguments are not reconstructed here.

## 2. Sequence optimization — article-reported

The `dl_binder_design` ProteinMPNN–FastRelax workflow masked the binder sequence, designed amino acids on its backbone, threaded them onto the complex, and used PyRosetta FastRelax before iterative redesign. Reported setting: `-seqs_per_struct=1`.

## 3. Complex prediction and initial selection

Boltz-1 with `--use_msa_server`; five predictions per binder–receptor pair. The reported filters were `confidence_score > 0.6` and `complex_pLDDT > 0.7`, followed by consideration of target-region engagement, nonredundant poses and energetics. These global confidence filters alone do not establish a reliable binding interface.

The included final comparison table contains Binder1–5 and brazzein, five predictions each. Source design indices 0–4 map to Binder1–5. An earlier table containing `Q8TE_A-amai` was not used for this brazzein comparison. Source JSON filesystem paths were removed from the published table; numeric values were preserved.

## 4. Structural preparation and static energies — article-reported

Schrödinger 2023.4 preparation included hydrogen-bond optimization and restrained minimization, with optional water removal. The top five complexes were scored using Uni-GBSA defaults. Static values in this repository were recovered from the original plotting notebook, not recalculated. ChimeraX 1.7 and PyMOL 3.0.0 were used for visualization.

## 5. MD — article-reported

| Parameter | Setting |
| --- | --- |
| Engine | Desmond, Schrödinger 2023-4 |
| Force field / water | OPLS4 / SPC |
| Periodic box | Cubic; minimum 12 Å solute-to-boundary distance |
| Ions | Neutralizing Na+ or Cl−, plus 0.15 M NaCl |
| Relaxation | Standard multistage Desmond relaxation |
| Production | NPT, 300 K, 1.01325 bar |
| Coordinate interval | 100 ps |
| Binder1–3 | 200 ns each |
| Advantame–TAS1R2 ECD, PDB 9OQ4 | 300 ns analyzed |
| Comparison windows | Common 0–200 ns; post-equilibration energetics targeted 20–200 ns |

The full simulated solvation systems and restart files are outside this lightweight release. Representative predicted PDBs are not claimed to be prepared production MD inputs.

## 6. Trajectory analysis — source-script definitions

- RMSD and RMSF describe coordinate deviations; Rg describes compactness. The binder curves describe the binder; advantame files separately describe receptor and ligand. Their atom selections are different and should not be treated as interchangeable measurements.
- The preserved Binder1 scripts select binder chain A and receptor chain B in the MD CMS. **This is an MD-system convention**; use `structures/chain_map.csv` for the separately supplied prediction PDBs.
- Hydrogen bonds, salt bridges and nonpolar contacts use Desmond `HydrogenBondFinder`, `SaltBridgeFinder` and `HydrophobicInterFinder`, respectively. The scripts use the installed library defaults rather than explicit distance/angle overrides; no numerical geometric cutoff has been invented in this companion.
- Contact counts are lengths of the atom-level event lists per frame. Residue-pair occupancy counts a pair at most once per frame and divides by the analyzed frame count. Combined occupancy uses the union of frames across contact types, not the sum of percentages. The preserved combination script hardcodes 2,001 frames for the 200-ns Binder1 trajectory.
- The contact scripts analyze the full input trajectory without a 20-ns discard. Consequently their exported occupancies must not silently be described as 20–200-ns occupancies.
- Interface area uses `SolventAccessibleSurfaceArea`, with the opposite partner excluded for the free-surface calculation. Total buried area is `(SASA_A_free − SASA_A_bound) + (SASA_B_free − SASA_B_bound)`; interface area is half this total. The source BSA script samples every tenth frame (~1 ns) and defines 20 ns for post-equilibration summaries.

## 7. Trajectory MM/GBSA — source launcher and result tables

The source launcher uses Schrödinger 2023-4 `thermal_mmgbsa.py`, `-start_frame 200 -end_frame 2000 -step_size 20 -NJOBS 8`, targeting 20–200 ns with 2-ns spacing. Binder ligand selection is `protein and chain.name "A"`; advantame selection is `res.pt A1CD`. The primary score column is `r_psp_MMGBSA_dG_Bind` (kcal/mol).

Each supplied result CSV has **90 rows**, although the launcher comment anticipates 91 snapshots. The result tables do not preserve a unique time for each row, so this release reports 90 observed results and the intended window; it does not fabricate an endpoint or a missing-frame explanation. Raw Prime exports contain repeated ligand-energy column names; they are preserved unchanged. The selected binding-score column is unique.

Summary plots show individual values and descriptive mean/sample SD, without new p-values or treating correlated trajectory frames as independent simulations. No new MD was run for this publication package.

## Methods references

- [Associated article](https://doi.org/10.1038/s41538-026-01202-x)
- [RFdiffusion](https://github.com/RosettaCommons/RFdiffusion)
- [ProteinMPNN–FastRelax binder design](https://github.com/nrbennet/dl_binder_design)
- [Boltz](https://github.com/jwohlwend/boltz)

The article and supplied research files are the basis for this historical workflow; current upstream defaults are not asserted to be the settings used in the original study.
