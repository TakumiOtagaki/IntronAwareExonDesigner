# IntronAwareExonDesigner

## Overview

IntronAwareExonDesigner is a small genetic-algorithm-driven pipeline that redesigns exon coding sequences within a larger transcript while keeping intron boundaries intact. Each candidate is evaluated with ViennaRNA to estimate both minimum free energy and pairing probabilities near intron splice sites, and the GA steers toward sequences that balance stable intron boundaries with low (more negative) folding energies.

## Requirements

Install the project dependencies listed in `pyproject.toml`, then ensure the ViennaRNA bindings are available in your Python environment. For example:

```bash
# brew install uv # for macbook user if uv is not installed.
uv sync
```


## Configuration

`config.yaml` defines the GA hyperparameters (`generations`, `population_size`, `mutation_rate`, `tournament_size`) and file paths for the `data/` FASTA inputs. Inspect that file to point to your transcript/UTR inputs or to change output directories/prefixes.

## Running a GA experiment

1. Prepare `data/chimera.fa` and `data/aminoacid.fa` following the FASTA conventions used by `IntronAwaredExonDesignerContext`.
2. Rewrite `config.yaml` as you like.
3. Run `run_ga.py`.

```bash
uv run python run_ga.py
```

or you can override the config.yaml settings with arguments, like;

```bash
uv run python run_ga.py --generations 10 --population 5
```


4. Results are written to `outputs/`:
   - `intron_design_dna_main_sequences.fa` / `intron_design_rna_main_sequences.fa` contain population members with `efe` and `sum_bpp` annotations.
   - `intron_design_stem_prob.tsv` summarizes per-base stem probabilities for later visualization.
      - You can utilize this tsv file to understand base pairing probability on forna visualization server (http://rna.tbi.univie.ac.at/forna/).
   - `intron_design_metrics.png` plots the GA statistics (α·sum_bpp, β·(-EFE), mix).

## example
## Example visualization

After running the GA with the default configuration, paste the best RNA sequence (`outputs/..._rna_main_sequences.fa`) into [Forna](http://rna.tbi.univie.ac.at/forna/) and upload the matching `intron_design_stem_prob.tsv` row via the “Colors” tab. This overlays stem probabilities onto the MFE structure—here is the result from a sample run:

![Forna example](img/example.png)

As you can see, the 3 leftmost and right most bases do not form base pair, and the probabilities of base pairing are low (white means low probability).

## Evaluating sequences

`evaluator.py` implements the ViennaRNA scoring pipeline used by `ga.py`. It extracts boundary indices within intron segments and sums pairing probabilities from `ViennaRNA.fold_compound.bpp()` together with intron-adjacent window energies to produce costs that drive the GA.



# citation
Forna
 - Kerpedjiev P, Hammer S, Hofacker IL (2015). Forna (force-directed RNA): Simple and effective online RNA secondary structure diagrams. Bioinformatics 31(20):3377-9.


# Future work
We plan to optimize the hairpin (or another loop-type) probability near intron splice sites, which might improve the splicing efficiency.
It can be utilized for the ensemble awared RNA inverse folding.
