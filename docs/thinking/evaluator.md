# Evaluator design notes

- Step 1 focuses on translating amino-acid sequences into codons with Biopython's NCBI table to keep the decision space purely synonymous.
- ViennaRNA is imported lazily so that the module can still be inspected even if the bindings are unavailable; evaluation errors are surfaced with a clear runtime message.
- Boundary pair scoring aggregates probabilities for the intron termini that the context reports, while caching avoids repeated ViennaRNA runs for identical sequences.
- Base-pair indexes returned by ViennaRNA come either as 0-based arrays or 1-based tuple-like entries, so the evaluator now normalizes everything to 0-based positions before comparing against the boundary offsets to avoid missing edge probabilities.

## オブジェクト情報

- **Genotype（`Sequence[int]`）**: コドンごとに同義コドン候補のインデックスを保持した整数列。`GenotypeDecoder` が各遺伝子座を対応するコドンに展開し、`SequenceEvaluator` は入力された遺伝子型から塩基配列を復元して評価に回す。
- **codon_lookup（`CodonLookup`）**: Biopython の `CodonTable` から NCBI table 1 を参照し、1文字のアミノ酸記号からソート済みコドンリストを取り出す。`GenotypeDecoder` で allele 値がどのコドンに対応するか決めるルックアップテーブル。
- **exon_design（`str`）**: `GenotypeDecoder.decode` の結果として得られる、全エクソン領域（上流＋下流）の DNA 文字列。長さはコンテキストの `design_length`（エクソンの塩基数）と一致し、`rebuild_main_with_exons` で元のメイン配列に差し替えられる。
- **full_sequence（`str`）**: `SequenceEvaluator` で評価対象となる完全配列。`utr5 + main（exon_design を反映） + utr3` を連結し、必要であれば大文字化したもの。これを ViennaRNA の `fold_compound` に渡す。
- **EvaluationResult**: `evaluate` が吐き出す結果オブジェクト。`genotype`（入力値）、`exon_sequence`/`full_sequence`（再構築された配列）、`boundary_pair_score`（イントロン端でのペア確率）、`energy`（パーティション関数由来の自由エネルギー）、`cost`（目的関数スコア）をまとめて保持する。
