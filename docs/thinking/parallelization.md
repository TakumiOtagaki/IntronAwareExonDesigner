# GA evaluation parallelization notes

- **現状**: `GeneticAlgorithm._evaluate_population` が `[ (g, evaluator.evaluate(g)) ... ]` で連続実行されており、ViennaRNA 評価がボトルネックになっている（1個あたり 0.05 秒程度）。cache は `SequenceEvaluator` が保持しているが、プロセスを並列化するとプロセスごとに再構築される。
- **目的**: `ProcessPoolExecutor` を使って評価部分だけをネイティブに並列化し、GA のループ自体はそのまま一貫したまま結果を集約したい。

### 計画
1. `GeneticAlgorithm` に `WorkerEvaluator` 構造を導入し、`SequenceEvaluator` のコンテキスト（ファイルパス・アミノ酸配列・EvaluatorConfig）から各ワーカープロセスで独立にインスタンス化できるようにする。メインプロセスからは genotype と必要なメタ情報をまとめて渡す。
2. `_evaluate_population` を並列版 `parallel_evaluate_population` に差し替える。`ProcessPoolExecutor` の `map` あるいは `executor.submit` で各 genotype を評価させ、結果をリストに戻す。
3. ViennaRNA や Bio.Python がプロセス内でインポートされるよう、`initializer` で `SequenceEvaluator` を初期化するか、worker 関数内で遅延インポートする。キャッシュはプロセス境界で共有できない点をドキュメントに明記する。
4. 並列化はオプション（例: `GeneticAlgorithmConfig.parallel_workers`) として実装し、0 または 1 なら現在の連続実行ループへフォールバックする。
5. 必要なら `SequenceEvaluator` に `to_worker_state()` 的なメソッドを追加して pickle 可能な状態を渡す。スコア計算自体は純粋なので、worker 内で `evaluate` を呼べばよい。
6. 実装前にテストの仮装（`tests/test_ga.py`）を改修して ProcessPoolExecutor をモックまたは `ThreadPoolExecutor` に差し替えるだけで検証できるように整備しておく。

今後はこの doc を参考に段階的に実装していく。
