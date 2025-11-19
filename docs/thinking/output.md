# Output artifacts

- GA run exports only the `main` transcript (intron+exon) as DNA and RNA multifasta files, which lets downstream tools compare both representations while keeping the fixed UTR context out of the inference. Each entry is annotated with its EFE and boundary probability summary.
- Per-generation statistics (mean/min of the α·sum_bpp, β·(−EFE), and α·sum_bpp + β·EFE terms) are plotted and saved alongside the multifasta. `config.yaml` now carries the default `outputs/` directory and `prefix`, and `run_ga.py` exposes CLI overrides so downstream callers can keep distinct runs organized.
