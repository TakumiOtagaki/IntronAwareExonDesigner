"""Tests for the evaluator components using a lightweight ViennaRNA stub."""

import math

from evaluator import SequenceEvaluator


SAMPLE_FASTA = """\
>5utr
AAAA
>main
ATGaaaGGG
>3utr
CCCC
"""


def make_context(tmp_path):
    path = tmp_path / "sample.fa"
    path.write_text(SAMPLE_FASTA, encoding="utf-8")
    return IntronAwaredExonDesignerContext(str(path), "MG")


def test_decoder_returns_design_length(tmp_path):
    context = make_context(tmp_path)
    evaluator = SequenceEvaluator(context)
    genotype = evaluator.decoder.random_genotype()
    exon_design = evaluator.decoder.decode(genotype)
    assert len(exon_design) == context.design_length


def test_evaluator_cost_matches_stubbed_rna(tmp_path):
    context = make_context(tmp_path)
    evaluator = SequenceEvaluator(context)
    result = evaluator.evaluate([0, 0])

    assert math.isclose(result.boundary_pair_score, 0.6, rel_tol=1e-6)
    assert math.isclose(result.energy, -5.0, rel_tol=1e-6)
    assert math.isclose(result.cost, 5.6, rel_tol=1e-6)
    assert result.exon_sequence == "ATGGGA"
    assert result.full_sequence == "AAAAATGAAAGGACCCC"
