# Claude Code entry point

Read `AGENTS.md` and follow its project conventions. For the current cluster retrieval task, then read `docs/PROJECT_STATUS.md`, `docs/cluster/CLAUDE_TASK.md` and `docs/cluster/WES_RECOVERY.md`.

Begin with existing WES run provenance and QC. The raw archive is intentionally not committed: absence from the clone does not mean a study file is missing. `docs/data/archived_input_inventory.tsv` lists inputs available locally before repository creation; `reports/audit_2026-09-05/wes_cluster_path_hints.csv` contains historical cluster clues. All new location findings need verification on the cluster.

Do not automatically run `scripts/run_all.sh`, restore the full R environment or submit analysis jobs merely to explore this checkout. `python3 scripts/check_checkout.py` needs only the Python standard library and checks the committed handoff plus processed release.
