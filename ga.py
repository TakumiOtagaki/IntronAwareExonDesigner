"""Basic genetic algorithm driver for exon design."""

from __future__ import annotations

from dataclasses import dataclass
from random import Random
from typing import Iterable, List, Protocol, Sequence, Tuple

from evaluator import EvaluationResult, GenotypeDecoder


class EvaluatorProtocol(Protocol):
    """Minimal evaluator protocol that exposes decoder and evaluation."""

    decoder: GenotypeDecoder

    def evaluate(self, genotype: Sequence[int]) -> EvaluationResult:
        ...


@dataclass(frozen=True)
class GeneticAlgorithmConfig:
    """Tunable GA parameters."""

    population_size: int = 20
    generations: int = 20
    mutation_rate: float = 0.05
    tournament_size: int = 3


class GeneticAlgorithm:
    """Standalone GA loop coordinating decoder, evaluator, and operators."""

    def __init__(
        self,
        evaluator: EvaluatorProtocol,
        config: GeneticAlgorithmConfig | None = None,
        rng: Random | None = None, # 乱数発生器
    ):
        self.evaluator = evaluator
        self.config = config or GeneticAlgorithmConfig()
        self.rng = rng or Random()
        self.population: List[List[int]] = []
        self.best_result: EvaluationResult | None = None
        self.cost_history: List[float] = []
        self._validate_config()

    @property
    def decoder(self) -> GenotypeDecoder:
        return self.evaluator.decoder

    def _validate_config(self) -> None:
        if self.config.population_size < 2:
            raise ValueError("population_size must be at least 2")
        if not 0.0 <= self.config.mutation_rate <= 1.0:
            raise ValueError("mutation_rate must be between 0 and 1")
        if self.config.tournament_size < 1:
            raise ValueError("tournament_size must be a positive integer")

    def run(self, generations: int | None = None) -> EvaluationResult:
        """Execute the GA for `generations` iterations and return the best solution."""
        self.cost_history = []
        self.population = self._initialize_population()
        evaluated = self._evaluate_population(self.population)
        self._record_best(evaluated)

        max_generations = generations if generations is not None else self.config.generations
        for _ in range(max_generations):
            self.population = self._breed_population(evaluated)
            evaluated = self._evaluate_population(self.population)
            self._record_best(evaluated)

        if self.best_result is None:
            raise RuntimeError("GA run did not produce any evaluation results.")
        return self.best_result

    def _initialize_population(self) -> List[List[int]]:
        return [self.decoder.random_genotype(rng=self.rng) for _ in range(self.config.population_size)]

    def _evaluate_population(
        self, population: Iterable[List[int]]
    ) -> List[Tuple[List[int], EvaluationResult]]:
        return [(genotype, self.evaluator.evaluate(genotype)) for genotype in population]

    def _record_best(self, evaluated: Iterable[Tuple[List[int], EvaluationResult]]) -> None:
        for _, result in evaluated:
            if self.best_result is None or result.cost < self.best_result.cost:
                self.best_result = result
        if self.best_result is not None:
            self.cost_history.append(self.best_result.cost)

    def _breed_population(
        self, evaluated: Sequence[Tuple[List[int], EvaluationResult]]
    ) -> List[List[int]]:
        children: List[List[int]] = []
        while len(children) < self.config.population_size:
            parent_a = self._select_parent(evaluated)
            parent_b = self._select_parent(evaluated)
            child = self._crossover(parent_a, parent_b)
            child = self._mutate(child)
            children.append(child)
        return children

    def _select_parent(
        self, evaluated: Sequence[Tuple[List[int], EvaluationResult]]
    ) -> List[int]:
        tournament_size = min(len(evaluated), self.config.tournament_size)
        candidates = self.rng.sample(list(evaluated), tournament_size)
        return min(candidates, key=lambda pair: pair[1].cost)[0]

    def _crossover(self, parent_a: Sequence[int], parent_b: Sequence[int]) -> List[int]:
        length = self.decoder.length
        child: List[int] = []
        for locus in range(length):
            if self.rng.random() < 0.5:
                child.append(parent_a[locus])
            else:
                child.append(parent_b[locus])
        return child

    def _mutate(self, genotype: Sequence[int]) -> List[int]:
        mutated = list(genotype)
        for locus in range(self.decoder.length):
            if self.rng.random() < self.config.mutation_rate:
                mutated[locus] = self.decoder.random_allele(locus, rng=self.rng)
        return mutated
