# スプライシング効率最適化システム 要件定義書

## 1. システム概要
イントロンおよびその周辺のエクソン領域の二次構造を制御し、スプライシング効率を向上させる（構造形成を抑制する）塩基配列を探索する。
**Hybrid GA-MCMC Strategy** を採用する。

## 2. データ構造と入力定義

### 2.1 入力データ (Config)
* **固定情報**:
    - 5utr, main(intron + exon), 3utr の multifasta file と アミノ酸配列を受け取る。intron は複数存在することも想定する.
    * `seq_upstream_aa`: 上流エクソンのアミノ酸配列（またはDNA配列 → アミノ酸に翻訳して保持）
    * `seq_intron`: イントロンの塩基配列（完全固定）
    * `seq_downstream_aa`: 下流エクソンのアミノ酸配列
    --> 例：`data/chimera.fa`
* **パラメータ**:
    * `window_size`: 構造計算を行う範囲（例: イントロンを中心に上流60nt + 下流30nt）
    * `weights`: 目的関数の重み係数 $\alpha$ (構造エネルギー), $\beta$ (末端対合確率)
    * `temperature`: MCMCの初期温度と減衰率
    --> `config.yaml` を作成する

### 2.2 遺伝子表現 (Genotype)

GAおよびMCMCで操作する対象。
* **形式**: 整数リスト `List[int]`
* **定義**: 各アミノ酸位置における「同義コドンリスト内のインデックス」。
    * 例: Met(ATG) は候補が1つなので常に `0`。Leu(TTA, TTG, CTT...) は候補が6つなので `0~5` の値をとる。
* **メリット**: 致死遺伝子（アミノ酸変異）を構造レベルで排除。

### 2.3 コドンテーブル (Lookup Table)
* **構造**: `Dict[str, List[str]]` (Amino Acid 1文字表記 -> コドン文字列のリスト)
    * 例: `{'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], ...}`
    -->. `Bio.Data.CodonTable.unambiguous_dna_by_id[1]` が NCBI table 1 であり、これを jsonl に変換しておく.
---

## 3. コアロジック: 配列構築と評価関数

### 3.1 配列構築 (Decoder)
Genotype（インデックス列）から表現型（塩基配列文字列）へ変換する処理。
* `decode(genotype) -> str`:
    1.  上流コドンインデックス列を塩基配列に変換。
    2.  イントロン配列（固定）を結合。
    3.  下流コドンインデックス列を塩基配列に変換。
    4.  結合して完全なターゲット配列を返す。

### 3.2 目的関数 (Cost Function)
最小化問題として定義する。スコアが低いほど良い。

$$Cost = \alpha \cdot (\sum_{pos \in Ends} P_{pair}(pos)) - \beta \cdot (Energy_{ensemble})$$

* **計算プロセス**:
    1.  **ViennaRNA Partition Function**:
        * `fc = RNA.fold_compound(sequence)`
        * `(propensity, energy) = fc.pf()`
        * `bpp_matrix = fc.bpp()`
    2.  **項1: イントロン端の対合確率**:
        * イントロン5'端 3塩基、3'端 3塩基のインデックスを特定。
        * `bpp_matrix` から該当箇所の対合確率（行または列の和）を取得し合計する。
    3.  **項2: 構造安定性 (Unstability Reward)**:
        * `energy` ($\Delta G$) は通常、安定なほど負の大きな値をとる（例: -20.5 kcal/mol）。
        * 構造を壊したい = $\Delta G$ を大きく（0に近く）したい。
        * 最小化問題にするため、マイナスをかけて符号を反転させる（$- \Delta G$）。
        * ※ $\beta$ は正の数。

---

## 4. 最適化エンジン仕様

### 4.1 Genetic Algorithm (GA) フェーズ
大域的な探索を行う。

* **個体群管理**:
    * `Population`: `List[Genotype]` (サイズ例: 100)
* **並列化 (必須)**:
    * 評価関数（ViennaRNA部分）がボトルネックとなるため、`ProcessPoolExecutor` を用いて個体群の評価を並列化する。
* **選択 (Selection)**:
    * トーナメント選択 (サイズ=3〜5)。
* **交叉 (Crossover)**:
    * **Uniform Crossover**: 各コドン位置について、50%の確率で親A、50%で親Bの遺伝子を選択。
* **突然変異 (Mutation)**:
    1.  **Random Mutation**: 確率 $P_m$ でランダムな同義コドンに変更。
    2.  **Smart Mutation (Destabilizer)**:
        * 確率 $P_s$ で発動。
        * 評価時に `bpp` が高かった（強くペアを組んでいる）エクソン領域の塩基を特定。
        * その塩基を含むコドンを、「GC含量が低い」または「元の塩基と異なる」同義コドンへ強制変更。
    - まずは random mutation だけの実装をしましょう.

### 4.2 MCMC / Simulated Annealing フェーズ
GAのベスト解を初期値として局所最適化を行う。

* **状態遷移**:
    * ランダムに1つのコドン座を選び、別の同義コドンへ変更。
* **受容確率 (Metropolis法)**:
    * $\Delta Cost = Cost_{new} - Cost_{old}$
    * $\Delta Cost < 0$ (改善) なら即受容。
    * $\Delta Cost \ge 0$ (改悪) なら確率 $\exp(-\Delta Cost / T)$ で受容。
* **クーリング**:
    * ステップごとに $T \leftarrow T \times decay\_rate$ (例: 0.995)

---

## 5. 実装上の考慮事項（Non-functional Requirements）

1.  **キャッシュ (Memoization)**:
    * 同じ配列（またはGenotype）が何度も評価されるのを防ぐため、LRUキャッシュまたは辞書を用いて `Sequence -> (Cost, BPP)` の結果を保持する。
2.  **インデックス管理**:
    * `Upstream` + `Intron` + `Downstream` を結合した際、イントロンの開始・終了位置が全体の何番目の塩基になるかを正確に追跡するオフセット管理クラスを作る。
3.  **ViennaRNAスレッドセーフ性**:
    * Pythonの `multiprocessing` を使う場合、各プロセス内で `import RNA` してインスタンス化することを推奨（Global Interpreter Lock回避とCライブラリの競合回避のため）。

---

### 開発順序

1.  **Step 1: Evaluatorの実装**
    * アミノ酸配列とイントロンを与えて、ランダムなコドンでDNA配列を作り、ViennaRNAでスコアを返す関数だけ作る。
    * テスト：わざと `GGGG...` などを入れてスコアが悪化（安定化）することを確認。
2.  **Step 2: GAクラスの実装**
    * Evaluatorを組み込み、並列処理なしで動作確認。
    * 続いて並列化を実施
3.  **Step 3: 並列化とMCMCの追加**
    * `concurrent.futures` で高速化し、最後にMCMCループを足す。
