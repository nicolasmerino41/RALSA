# Degree model selection

Start with `outputs/01_model_selection_overview.png` and `outputs/02_reference_table.png`.
The findings, reference table and concise methods are in `outputs/RESULTS_AND_METHODS.md`.

Consumer degrees; four MLE candidates; 100 removals per treatment; 199 refitted bootstrap samples per KS test. No singleton-plus-geometric model.

Run: `julia --startup-file=no --project=PATH_TO_RALSA run.jl PATH_TO_RALSA`.
All outputs stay in this folder. Earlier analyses remain unchanged.
