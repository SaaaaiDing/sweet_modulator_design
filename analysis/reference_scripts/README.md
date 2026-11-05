# Original Binder1 analysis scripts

These are unmodified copies from the research source directory, included to expose analysis definitions. They require the licensed Schrödinger Python environment and `binder1_200ns-out.cms` plus its trajectory. Inputs are intentionally not bundled. Scripts write output CSVs in the working directory.

Chain A is the MD binder; chain B is the receptor. Contact finders use library defaults. The combined-occupancy script requires the event CSVs produced by the three contact scripts and assumes 2,001 frames. The BSA script uses a stride of 10. See `../../docs/WORKFLOW.md` for window and interpretation details. These scripts are references, not a portable end-to-end runner.
