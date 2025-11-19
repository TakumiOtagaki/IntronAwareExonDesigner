# Output artifacts

- GA run exports only the `main` transcript (intron+exon) as a multifasta so downstream analyses can focus on the design space without the fixed UTR context. The file is annotated with EFE and boundary probabilities for each entry.
- Per-generation statistics (mean/min of the α·sum_bpp, β·(−EFE), and α·sum_bpp + β·EFE terms) are plotted and saved alongside the multifasta. `config.yaml` now carries the default `outputs/` directory and `prefix`, and `run_ga.py` exposes CLI overrides so downstream callers can keep distinct runs organized.
