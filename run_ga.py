"""Run a short GA epoch on the provided data set to inspect cost improvements."""

from __future__ import annotations

from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from random import Random
from typing import Any, Callable, Iterable, List, Sequence
import RNA
import yaml
from matplotlib import pyplot as plt

from Bio import SeqIO

from evaluator import EvaluationResult, EvaluatorConfig, SequenceEvaluator
from ga import GeneticAlgorithm, GeneticAlgorithmConfig
from intron_design import IntronAwaredExonDesignerContext


DEFAULT_CONFIG_PATH = Path("config.yaml")
DEFAULT_CHIIRMA_PATH = Path("data/chimera.fa")
DEFAULT_AMINO_PATH = Path("data/aminoacid.fa")


def _sequence_snippet(sequence: str, width: int = 80) -> str:
    """Return a truncated representation of a sequence with start/end context."""
    if len(sequence) <= width:
        return sequence
    segment = max(1, (width - 1) // 2)
    return f"{sequence[:segment]}…{sequence[-segment:]}"


@dataclass(frozen=True)
class _FullSequenceEntry:
    label: str
    result: EvaluationResult
    dna_sequence: str
    rna_sequence: str


def _multifasta_header(label: str, energy: float, boundary_score: float) -> str:
    return f">{label} | efe={energy:.4f} | sum_bpp={boundary_score:.4f}"


def _collect_full_sequence_entries(
    context: IntronAwaredExonDesignerContext,
    best_result: EvaluationResult,
    final_population: Iterable[EvaluationResult],
) -> List[_FullSequenceEntry]:
    entries: List[tuple[str, EvaluationResult]] = [
        ("best", best_result),
    ]
    entries.extend((f"final_population{idx}", result) for idx, result in enumerate(final_population))
    full_entries: List[_FullSequenceEntry] = []
    for label, result in entries:
        dna_seq = context.build_full_sequence(result.exon_sequence, uppercase=False)
        rna_seq = context.build_full_sequence(
            result.exon_sequence, uppercase=False, rna_output=True
        )
        full_entries.append(_FullSequenceEntry(label, result, dna_seq, rna_seq))
    return full_entries


def _write_full_multifasta(
    path: Path,
    entries: Iterable[_FullSequenceEntry],
    sequence_extractor: Callable[[_FullSequenceEntry], str],
) -> None:
    lines: List[str] = []
    for entry in entries:
        lines.append(
            _multifasta_header(entry.label, entry.result.energy, entry.result.boundary_pair_score)
        )
        lines.append(sequence_extractor(entry))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _sum_bpp_values(bpp_matrix, sequence_length: int) -> List[float]:
    totals = [0.0] * sequence_length
    if sequence_length == 0:
        return totals
    for i_one, j_one, prob in SequenceEvaluator._iter_bpp_entries(bpp_matrix):
        i = i_one - 1
        j = j_one - 1
        if 0 <= i < sequence_length:
            totals[i] += prob
        if 0 <= j < sequence_length:
            totals[j] += prob
    return totals


def _compute_sequence_stem_probabilities(task: tuple[str, str]) -> tuple[str, List[float]]:
    label, sequence = task
    fold_compound = RNA.fold_compound(sequence)
    _, mfe = fold_compound.mfe()
    fold_compound.exp_params_rescale(mfe)
    fold_compound.pf()
    bpp_matrix = fold_compound.bpp()
    values = _sum_bpp_values(bpp_matrix, len(sequence))
    return label, values


def _compute_stem_probabilities(
    entries: Sequence[_FullSequenceEntry],
    max_workers: int | None,
) -> dict[str, List[float]]:
    tasks = [(entry.label, entry.rna_sequence) for entry in entries]
    if not tasks:
        return {}
    if max_workers == 1:
        results = [_compute_sequence_stem_probabilities(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            results = list(executor.map(_compute_sequence_stem_probabilities, tasks))
    return {label: values for label, values in results}


def _write_combined_stem_probability_tsv(
    path: Path, entries: Sequence[_FullSequenceEntry], probabilities: dict[str, List[float]]
) -> None:
    sequence_length = max(len(entry.rna_sequence) for entry in entries)
    header = ["name"] + [str(idx) for idx in range(sequence_length)]
    lines = ["\t".join(header)]
    for entry in entries:
        values = probabilities.get(entry.label, [])
        rounded = [
            f"{values[idx]:.2f}" if idx < len(values) else ""
            for idx in range(sequence_length)
        ]
        lines.append("\t".join([entry.label] + rounded))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _plot_generation_metrics(metrics: Iterable[dict[str, float]], destination: Path) -> None:
    data = list(metrics)
    if not data:
        return
    generations = [entry["generation"] for entry in data]
    fields = [
        ("alpha_bpp", "Alpha × sum_bpp"),
        ("beta_minus_efe", "Beta × (-EFE)"),
        ("alpha_beta_mix", "Alpha × sum_bpp + Beta × (-EFE)"),
    ]
    fig, axes = plt.subplots(nrows=len(fields), ncols=1, figsize=(10, 8), sharex=True)
    if len(fields) == 1:
        axes = [axes]
    for ax, (prefix, label) in zip(axes, fields):
        mean_key = f"{prefix}_mean"
        min_key = f"{prefix}_min"
        ax.plot(generations, [entry[mean_key] for entry in data], label="mean", marker="o")
        ax.plot(generations, [entry[min_key] for entry in data], label="min", marker="s")
        ax.set_ylabel(label)
        ax.grid(True, linestyle="--", alpha=0.5)
        ax.legend()
    axes[-1].set_xlabel("Generation")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.suptitle("GA population metrics", fontsize=14)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination)
    plt.close(fig)


def _load_config(path: Path | None = None) -> dict[str, Any]:
    target = path or DEFAULT_CONFIG_PATH
    if not target.exists():
        return {}
    return yaml.safe_load(target.read_text(encoding="utf-8")) or {}


def _load_amino_sequence(path: Path | None = None) -> str:
    target = path or DEFAULT_AMINO_PATH
    record = SeqIO.read(target, "fasta")
    return str(record.seq)


def main() -> None:
    config = _load_config()
    parser = ArgumentParser(description="Run a short GA search to inspect cost improvements.")
    parser.add_argument("--generations", type=int, help="Override GA generations count.")
    parser.add_argument("--population", type=int, help="Override GA population size.")
    parser.add_argument("--mutation", type=float, help="Override GA mutation rate.")
    parser.add_argument("--tournament", type=int, help="Override GA tournament size.")
    parser.add_argument("--workers", type=int, help="Override GA parallel_workers setting.")
    parser.add_argument("--output-prefix", help="Override output file prefix.")
    parser.add_argument("--output-dir", help="Override output directory.")
    parser.add_argument("--seed", type=int, help="Deterministic RNG seed.")
    args = parser.parse_args()

    ga_section = config.get("ga", {})
    generations = args.generations if args.generations is not None else ga_section.get("generations", 5)
    population = args.population if args.population is not None else ga_section.get("population_size", 20)
    mutation_rate = args.mutation if args.mutation is not None else ga_section.get("mutation_rate", 0.05)
    tournament_size = args.tournament if args.tournament is not None else ga_section.get("tournament_size", 3)
    seed = args.seed
    workers_override = args.workers
    parallel_workers = (
        workers_override
        if workers_override is not None
        else ga_section.get("parallel_workers", config.get("parallel_workers"))
    )
    output_config = config.get("output", {})
    output_prefix = args.output_prefix or output_config.get("prefix", "ga_run")
    output_dir = Path(args.output_dir or output_config.get("directory", "outputs"))

    evaluator_config = EvaluatorConfig(
        window_upstream=config.get("window_upstream", 60),
        window_downstream=config.get("window_downstream", 30),
        alpha=config.get("weights", {}).get("alpha", 1.0),
        beta=config.get("weights", {}).get("beta", 1.0),
    )

    amino_seq = _load_amino_sequence()
    context = IntronAwaredExonDesignerContext(DEFAULT_CHIIRMA_PATH, amino_seq)
    evaluator = SequenceEvaluator(context, config=evaluator_config)

    initial_exon = context.get_exon_sequence()
    initial_dna_full = context.build_full_sequence(initial_exon)
    initial_rna_full = context.build_full_sequence(initial_exon, rna_output=True)
    baseline_boundary_score, baseline_energy = evaluator.evaluate_rna_sequence(initial_rna_full)
    initial_cost = evaluator.config.alpha * baseline_boundary_score - evaluator.config.beta * baseline_energy
    print("\nBaseline (initial sequence)")
    print(f"  Boundary score: {baseline_boundary_score:.4f}")
    print(f"  Energy: {baseline_energy:.4f}")
    print(f"  Cost: {initial_cost:.4f}")

    ga_config = GeneticAlgorithmConfig(
        population_size=population,
        generations=generations,
        mutation_rate=mutation_rate,
        tournament_size=tournament_size,
        parallel_workers=parallel_workers,
    )
    rng = Random(seed)
    ga = GeneticAlgorithm(evaluator, config=ga_config, rng=rng)


    best = ga.run(generations=generations)

    history = " → ".join(f"{cost:.4f}" for cost in ga.cost_history)
    snippet = _sequence_snippet(best.full_sequence)
    rna_full = context.build_full_sequence(best.exon_sequence, rna_output=True)

    print("GA run complete")
    print(f"Generations: {generations}; Population: {population}; Mutation: {mutation_rate:.3f}")
    print(f"Tournament size: {tournament_size}; RNG seed: {seed if seed is not None else 'random'}")
    print(f"number of workers: {parallel_workers or 'disabled'}")
    print(f"Cost history: {history}")
    print(f"Best cost: {best.cost:.4f}")
    print(f"Best genotype: {best.genotype}")
    print(f"Parallel workers: {parallel_workers or 'disabled'}")
    print(f"Best full (DNA) sequence snippet: {snippet}")
    print(f"Best full (RNA) sequence snippet: {_sequence_snippet(rna_full)}")
    print(f"Boundary score: {best.boundary_pair_score:.4f}; Energy: {best.energy:.4f}")

    output_dir.mkdir(parents=True, exist_ok=True)
    entries = _collect_full_sequence_entries(context, best, ga.final_population)
    dna_path = output_dir / f"{output_prefix}_dna_main_sequences.fa"
    rna_path = output_dir / f"{output_prefix}_rna_main_sequences.fa"
    _write_full_multifasta(dna_path, entries, lambda entry: entry.dna_sequence)
    _write_full_multifasta(rna_path, entries, lambda entry: entry.rna_sequence)

    stem_probabilities = _compute_stem_probabilities(entries, parallel_workers)
    stem_path = output_dir / f"{output_prefix}_stem_prob.tsv"
    _write_combined_stem_probability_tsv(stem_path, entries, stem_probabilities)

    metrics_path = output_dir / f"{output_prefix}_metrics.png"
    _plot_generation_metrics(ga.generation_metrics, metrics_path)
    print(f"DNA multifasta written to {dna_path}")
    print(f"RNA multifasta written to {rna_path}")
    print(f"Stem probability TSV written to {stem_path}")
    print(f"Generation metric plot written to {metrics_path}")


if __name__ == "__main__":
    main()
