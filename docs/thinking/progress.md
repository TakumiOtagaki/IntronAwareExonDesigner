# Project progress notes

- Evaluator/context setup: `intron_design.py` now parses RNA-aware multi-FASTA, tracks exon/intron segments, and rebuilds full DNA/RNA sequences; `SequenceEvaluator` decodes genotypes and scores them via ViennaRNA while keeping DNA output intact for downstream reporting.
- Step 1 verification: `tests/test_evaluator.py` uses the real context helper and a ViennaRNA stub so we can assert boundary scores/energetics early. Added `sys.path` insertion to keep imports resolvable in the test runner.
- GA driver (Step 2): `ga.py` implements a uniform crossover GA with tournament selection, mutation, best-result tracking, and a `cost_history` list; `tests/test_ga.py` provides a deterministic dummy evaluator to confirm the GA prefers low-cost genotypes.
- Configuration + runner: `config.yaml` captures window, GA, and output preferences; `run_ga.py` ties the real FASTA data to `SequenceEvaluator` and `GeneticAlgorithm`, prints cost history plus DNA/RNA snippets, and exposes CLI overrides.
- Operational notes: ViennaRNA (the `RNA` module) must be installed for `run_ga.py` to succeed; `uv` runs still face caching permissions. Tests (`uv run pytest`) currently pass after installing required dependencies locally in the uv environment.
