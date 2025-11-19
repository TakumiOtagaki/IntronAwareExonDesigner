"""Unit tests for the GA driver using deterministic dummy evaluator."""

from random import Random

from evaluator import EvaluationResult
from ga import GeneticAlgorithm, GeneticAlgorithmConfig


class DummyDecoder:
    """Simple decoder that exposes allele ranges without biological details."""

    def __init__(self, allele_counts):
        self._allele_counts = tuple(allele_counts)
        self.length = len(self._allele_counts)

    def random_genotype(self, rng=None):
        rng = rng or Random()
        return [rng.randrange(count) for count in self._allele_counts]

    def random_allele(self, locus, rng=None):
        rng = rng or Random()
        return rng.randrange(self._allele_counts[locus])


class DummyEvaluator:
    """Evaluator stub that returns cost equal to the sum of the genotype."""

    def __init__(self, decoder):
        self.decoder = decoder

    def evaluate(self, genotype):
        cost = sum(genotype)
        length = self.decoder.length * 3
        return EvaluationResult(
            genotype=tuple(genotype),
            exon_sequence="A" * length,
            full_sequence="A" * (length + 10),
            cost=float(cost),
            boundary_pair_score=0.0,
            energy=0.0,
        )


def test_genetic_algorithm_favors_low_cost_sequences():
    decoder = DummyDecoder([2, 2, 2])
    evaluator = DummyEvaluator(decoder)
    config = GeneticAlgorithmConfig(population_size=6, generations=5, mutation_rate=0.7, tournament_size=2)
    ga = GeneticAlgorithm(evaluator, config=config, rng=Random(0))

    best = ga.run()

    assert best.cost == 0.0
