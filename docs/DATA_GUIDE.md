# Article and data guide

## Figure / analysis mapping

| Article topic | Included companion material | Scope |
| --- | --- | --- |
| Figure 1: design workflow | `docs/WORKFLOW.md` | Method summary, not the original campaign archive |
| Figure 2: confidence comparison | `data/boltz/metrics.csv`, `figures/boltz_summary.png` | Five predictions each for Binder1–5 and brazzein; new descriptive plot |
| Figure 3: structures / secondary structure | `structures/`, `data/secondary_structure/`, `figures/secondary_structure.svg` | Five representative complex PDBs, extracted binder sequences, SSE fractions |
| Figure 4: Binder1–3 MD | `data/md/binder*/`, `figures/binder_md_summary.png` | Rg, RMSD, RMSF and exported contacts / buried area |
| Figure 5: advantame reference | `data/md/advantame/`, `figures/advantame_*.svg` | Selected original 300-ns reference curves; not every article panel |
| Figure 6: energies | `data/energies/`, `figures/mmgbsa_summary.png` | Static plotting values and four trajectory energy tables are separate |
| Binder2 R37A/T55A decoy | Article only | No matching local decoy prediction table was identified; no substitute data created |

## Reading the numbers

| Field / quantity | Units and interpretation |
| --- | --- |
| `confidence_score`, `ptm`, `iptm` | Dimensionless prediction confidence; higher is generally better. ipTM concerns interface confidence, not measured binding affinity |
| `complex_plddt`, `complex_iplddt` | Local confidence in the source table's 0–1 scale; not 0–100 |
| `complex_pde`, `complex_ipde` | Predicted distance error in Å; lower is generally better |
| `ligand_iptm` | Source field retained; zeros for these protein complexes should not be interpreted as measured ligand affinity |
| Helix / Strand / Coil | Fractions from the original MDTraj/DSSP analysis, summing to one |
| RMSD / RMSF / Rg | Å; RMSD and RMSF depend on atom selection and alignment |
| Contact count | Number of detected atom-level events in a sampled frame, not the number of distinct residue pairs |
| Residue-pair occupancy | Fraction or percent of frames with at least one event for that pair; retain the source column's explicit unit |
| `total_BSA_A2` | Total buried solvent-accessible area, Å² |
| `interface_area_half_BSA_A2` | Half the total buried area, Å²; this is the article's interface-area convention |
| `r_psp_MMGBSA_dG_Bind` | Prime trajectory MM/GBSA score, kcal/mol; lower is more favorable within the same protocol |

See `data/summary/time_ranges.csv` for the actual exported time ranges and row counts. The advantame structural curves have 3,038 records ending at 303.652854 ns, beyond the article's nominal 300 ns; the raw exports are preserved. Its contact-count exports contain 1,800 records from 20.0004 to 199.903998 ns, whereas binder contact tables span the full 0–200 ns trajectory. Summary tables report both the requested window and the actual available start/end; a requested 0–200-ns window does not supply missing advantame contact records. Window selection uses a 0.01-ns tolerance for minor timestamp offsets. RMSF files are indexed by residue rather than time. Mean ± SD over snapshots describes observed variation, not uncertainty from independent simulations.

## Provenance and interpretation limits

`provenance/source_manifest.json` maps each copied or transformed source to its published path and SHA-256 hash. `structures/chain_map.csv` is obtained directly from CA records; `design_sequences.fasta` contains the unique 128-residue chain in each representative PDB. Each supplied PDB contains binder chain A (128 CA residues) and receptor chain B (839 CA residues): these are full-length receptor prediction complexes, whereas the article's MD systems used the extracellular domain. These PDBs are taken from the original figure preparation folder; their Boltz model indices and prediction seeds are not asserted.

Paper tables use several atom/residue interaction summaries. In particular, Tables 2 and 3 list different occupancies for some named contacts. The released CSVs retain their source definitions and are not relabeled as exact reconstructions of both paper tables. Full-trajectory occupancy and post-equilibration occupancy are distinct. The current summary script only computes descriptive numerical summaries and does not recompute contacts from coordinates.

Binder2 has the strongest mean ipTM among the five designs but the absolute interface confidence is modest. Favorable GBSA scores, structural persistence and interface contacts are complementary prioritization evidence, not proof of receptor activation. Larger protein interfaces can also receive more favorable endpoint energies than a small molecule, so comparisons with advantame require that context.

The collection intentionally omits full trajectories, CMS/checkpoint files, MSAs, repeated structures, oversized images and the publisher PDF. It supports inspection of the reported computational approach and selected numerical outputs, not a claim of complete independent reproduction.
