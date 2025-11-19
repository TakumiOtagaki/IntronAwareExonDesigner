"""Evaluator components for the intron-aware exon design pipeline."""

from __future__ import annotations

import RNA

from collections import defaultdict
from dataclasses import dataclass, field
from random import Random
from typing import Dict, Iterator, List, Mapping, Sequence, Tuple

from Bio.Data import CodonTable

from intron_design import IntronAwaredExonDesignerContext


@dataclass(frozen=True)
class EvaluatorConfig:
    """Parameters for controlling the ViennaRNA evaluation and boundary emphases."""

    window_upstream: int = 60
    window_downstream: int = 30
    alpha: float = 1.0
    beta: float = 1.0
    flank: int = 3

    @property
    def window_size(self) -> int:
        return self.window_upstream + self.window_downstream


@dataclass
class CodonLookup:
    """Utility that exposes synonym codons for each amino-acid residue."""

    table_id: int = 1
    codons_by_aa: Dict[str, List[str]] = field(init=False)

    def __post_init__(self) -> None:
        table = CodonTable.unambiguous_dna_by_id[self.table_id]
        mapping: Dict[str, List[str]] = {}
        for codon, residue in table.forward_table.items():
            mapping.setdefault(residue, []).append(codon)
        self.codons_by_aa = {residue: sorted(codons) for residue, codons in mapping.items()}

    def get_codons(self, amino_acid: str) -> List[str]:
        """Return the ordered codon list for a given amino-acid symbol."""
        aa = amino_acid.upper()
        try:
            return self.codons_by_aa[aa]
        except KeyError as exc:
            raise ValueError(f"Unknown amino-acid '{amino_acid}' for codon lookup.") from exc


@dataclass
class GenotypeDecoder:
    """Translate genotype indices into the corresponding DNA codon sequence."""

    amino_acid_sequence: str
    codon_lookup: CodonLookup
    _codon_options: List[List[str]] = field(init=False)

    def __post_init__(self) -> None:
        residues = [
            char.upper()
            for char in self.amino_acid_sequence
            if char.isalpha() or char == "*"
        ]
        if not residues:
            raise ValueError("Amino-acid sequence must contain at least one designable residue.")
        self._codon_options = [self.codon_lookup.get_codons(residue) for residue in residues]

    @property
    def length(self) -> int:
        return len(self._codon_options)

    def decode(self, genotype: Sequence[int]) -> str:
        """Map each allele index into the DNA codon and concatenate."""
        if len(genotype) != self.length:
            raise ValueError(
                f"Genotype length {len(genotype)} does not match expected design length {self.length}."
            )
        codons: List[str] = []
        for locus, allele in enumerate(genotype):
            options = self._codon_options[locus]
            if not 0 <= allele < len(options):
                raise ValueError(
                    f"Allele {allele} at position {locus} is out of range for amino-acid "
                    f"with {len(options)} synonymous codons."
                )
            codons.append(options[allele])
        return "".join(codons)

    def random_genotype(self, rng: Random | None = None) -> List[int]:
        """Produce a random genotype respecting each amino-acid's codon count."""
        rng = rng or Random()
        return [rng.randrange(len(options)) for options in self._codon_options]


@dataclass(frozen=True)
class EvaluationResult:
    genotype: Tuple[int, ...]
    exon_sequence: str
    full_sequence: str
    cost: float
    boundary_pair_score: float
    energy: float


class SequenceEvaluator:
    """Wraps ViennaRNA evaluation logic so that genotypes can be compared."""

    def __init__(
        self,
        context: IntronAwaredExonDesignerContext,
        config: EvaluatorConfig | None = None,
        codon_lookup: CodonLookup | None = None,
    ):
        self.context = context
        self.config = config or EvaluatorConfig()
        self.codon_lookup = codon_lookup or CodonLookup()
        self.decoder = GenotypeDecoder(context.amino_acid_sequence, self.codon_lookup)
        self._sequence_cache: Dict[str, Tuple[float, float]] = {}

    @property
    def design_length(self) -> int:
        return self.decoder.length

    def evaluate(self, genotype: Sequence[int]) -> EvaluationResult:
        exon_design = self.decoder.decode(genotype)
        full_sequence = self.context.build_full_sequence(exon_design)
        boundary_score, energy = self._cached_evaluate(full_sequence)
        cost = self.config.alpha * boundary_score - self.config.beta * energy
        return EvaluationResult(
            genotype=tuple(genotype),
            exon_sequence=exon_design,
            full_sequence=full_sequence,
            cost=cost,
            boundary_pair_score=boundary_score,
            energy=energy,
        )

    def _cached_evaluate(self, sequence: str) -> Tuple[float, float]:
        if sequence in self._sequence_cache:
            return self._sequence_cache[sequence]
        metrics = self._evaluate_sequence(sequence)
        self._sequence_cache[sequence] = metrics
        return metrics

    def _evaluate_sequence(self, sequence: str) -> Tuple[float, float]:
        fold_compound = RNA.fold_compound(sequence)
        pf_result = fold_compound.pf()
        if isinstance(pf_result, tuple):
            _, energy = pf_result
        else:
            energy = float(pf_result)
        bpp_matrix = fold_compound.bpp()
        boundary_score = self._boundary_pair_probability(bpp_matrix)
        return boundary_score, energy

    def _boundary_pair_probability(self, bpp_matrix) -> float:
        weights: Dict[int, float] = defaultdict(float)
        for i, j, prob in self._iter_bpp_entries(bpp_matrix):
            weights[i - 1] += prob
            weights[j - 1] += prob

        score = 0.0
        for idx in self.context.get_boundary_indices(flank=self.config.flank):
            score += weights.get(idx, 0.0)
        return score

    @staticmethod
    def _iter_bpp_entries(bpp_matrix) -> Iterator[Tuple[int, int, float]]:
        if not bpp_matrix:
            return iter(())

        if isinstance(bpp_matrix, Mapping):
            return (
                (coords[0], coords[1], float(prob))
                for coords, prob in bpp_matrix.items()
                if isinstance(coords, tuple) and len(coords) == 2
            )

        def generator():
            for entry in bpp_matrix:
                if isinstance(entry, tuple) and len(entry) == 3:
                    yield int(entry[0]), int(entry[1]), float(entry[2])
                    continue

                i = getattr(entry, "i", None) or getattr(entry, "I", None)
                j = getattr(entry, "j", None) or getattr(entry, "J", None)
                prob = (
                    getattr(entry, "p", None)
                    or getattr(entry, "prob", None)
                    or getattr(entry, "probability", None)
                )
                if i is None or j is None or prob is None:
                    continue
                yield int(i), int(j), float(prob)

        return generator()
