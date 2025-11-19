"""Run a short GA epoch on the provided data set to inspect cost improvements."""

from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path
from random import Random
from typing import Any
import yaml


from Bio import SeqIO

from evaluator import EvaluatorConfig, SequenceEvaluator
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
    print(f"Cost history: {history}")
    print(f"Best cost: {best.cost:.4f}")
    print(f"Best genotype: {best.genotype}")
    print(f"Parallel workers: {parallel_workers or 'disabled'}")
    print(f"Best full (DNA) sequence snippet: {snippet}")
    print(f"Best full (RNA) sequence snippet: {_sequence_snippet(rna_full)}")
    print(f"Boundary score: {best.boundary_pair_score:.4f}; Energy: {best.energy:.4f}")


if __name__ == "__main__":
    main()
