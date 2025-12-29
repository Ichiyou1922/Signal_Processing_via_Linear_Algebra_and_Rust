# Signal Processing via Linear Algebra & Rust
> **Reconstruction of Signal Processing Theory based on Vector Space Structures**

## 1. Project Overview
本プロジェクトの目的は，大学の「信号処理論」における主要トピック（DFT, 窓関数, FFT, フィルタ）を，既存ライブラリに頼らず**線形代数（内積空間・基底変換・作用素）の視点から再構築・実装**することである．

単なる数値計算ではなく，信号を**「複素ベクトル空間 $\mathbb{C}^N$ 上の点」**として捉え直し，その幾何学的構造（直交性，回転，射影）をRustによるスクラッチ実装を通じて検証する．

## 2. Report Structure (The Logic)
レポートは以下の論理構成で記述し，線形代数の概念がいかに信号処理の実体を支えているかを論証する．

### Chapter 1: 信号空間と内積 (The Vector Space of Signals)
* **Concept**: 離散信号 $x[n]$ の $\mathbb{C}^N$ ベクトルとしての再定義．
* **Core Logic**: 「似ている」ことの数学的定義としての**エルミート内積**とノルム．
* **Implementation**: 複素数構造体と内積演算の定義．

### Chapter 2: DFTの本質と基底変換 (DFT as Change of Basis)
* **Concept**: DFTを「時間領域基底」から「周波数領域基底（直交基底）」への**座標変換**として定式化する．
* **Core Logic**: 離散指数関数系 $\{e^{j\omega n}\}$ の**直交性**の証明．DFT行列（ヴァンデルモンド行列）のユニタリ性．
* **Implementation**: 直交基底生成とDFT（行列ベクトル積）の実装．

### Chapter 3: 窓関数とスペクトル漏れ (Window Functions & Leakage)
* **Concept**: 有限区間切り出しによる「基底の直交性の破壊」としてのスペクトル漏れ．
* **Core Logic**: 窓関数を**対角行列 $W$ による線形変換（作用素）**として記述する．
* **Implementation**: ハニング窓等の行列生成と信号への作用．漏れの可視化．

### Chapter 4: 高速フーリエ変換 (FFT: Matrix Factorization)
* **Concept**: 計算量 $O(N^2)$ から $O(N \log N)$ への短縮．
* **Core Logic**: DFT行列の**因数分解（疎行列の積への分解）**と再帰構造．
* **Implementation**: Cooley-Tukeyアルゴリズム（再帰）の実装．

### Chapter 5: ディジタルフィルタ (Digital Filters as Operators)
* **Concept**: 周波数領域における操作と時間領域への書き戻し．
* **Core Logic**: 線形時不変システム（LTI）としてのフィルタ．畳み込みと積の双対性．
* **Implementation**: 理想低域通過フィルタ（LPF）の設計とIDFTの実装．

---

## 3. Implementation & Learning Roadmap (Action Plan)

各フェーズにおいて，**「線形代数の理論学習（Note）」**と**「Rustによる検証（Code）」**をセットで行う．

### Phase 1: 空間の定義 (The Foundation)
- [ ] **Linear Algebra**: `Inner_Product.md` の完了（複素内積，共役対称性，ノルム）．
- [ ] **Rust**:
    - `struct Complex` (複素数型) の実装．
    - `struct Signal` (ベクトル型) の実装．
    - `fn inner_product` (内積 $<x, y> = x^H y$) の実装．

### Phase 2: 変換の実装 (The Transform)
- [ ] **Linear Algebra**: DFT基底の直交性を手計算で証明する．
- [ ] **Rust**:
    - 基底生成メソッド `fn basis(N, k)` の実装．
    - 素朴なDFT `fn dft(signal)` の実装（行列演算として）．
- [ ] **Experiment**: 正弦波の分解実験（スペクトルのピーク確認）．

### Phase 3: 現実との対峙 (The Artifacts)
- [ ] **Linear Algebra**: 射影と部分空間，作用素としての対角行列．
- [ ] **Rust**:
    - 窓関数生成 `fn hanning(N)` の実装．
    - 窓適用メソッドの実装．
- [ ] **Experiment**: 矩形窓 vs ハニング窓のスペクトル比較（漏れの観測）．

### Phase 4: 高速化と応用 (Optimization & Application)
- [ ] **Linear Algebra**: 行列のブロック分解と回転因子．
- [ ] **Rust**:
    - 再帰的FFT `fn fft(signal)` の実装．
    - IDFT（逆変換）の実装．
    - フィルタリング処理の実装．

## 4. Environment
* Language: Rust (Edition 2021)
* Dependencies: `rand` (for signal generation only), `num-traits` (optional)
    * *Note: `rustfft`等のFFTライブラリは使用禁止（Scracth Build）．*Rust
