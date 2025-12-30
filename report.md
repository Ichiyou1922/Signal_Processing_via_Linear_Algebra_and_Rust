---
title: "フーリエ変換・窓関数・高速フーリエ変換"
author: "24c1121 茂木一葉"
date: "2025-12-30"
documentclass: ltjsarticle
header-includes:
  - \usepackage{amsmath}
  - \usepackage{amssymb}
  - \usepackage{bm}
  - \usepackage{graphicx}
---

# 1. 目的 (Objective)
本レポートの目的は2つ．

- 信号処理論において学んだ，フーリエ変換・窓関数・高速フーリエ変換を線形代数の視点から捉え直すこと(理論)．
- 離散フーリエ変換・高速フーリエ変換の計算量の差について記述する(アルゴリズム)．

# 2. 理論 (Theory)

本実験で使用する基本原理について，線形代数的な観点から記述する．

## 2.1 フーリエ変換における基底変換

長さ $N$ の任意の離散信号 $\mathbf{x} \in \mathbb{C}^N$ は，「時刻 $t=k$ のとき値が $1$」であるインパルス信号（標準基底） $\mathbf{e}_k$ を用いた線形結合で表せる．

$$
\mathbf{x} = x[0]\mathbf{e}_0 + x[1]\mathbf{e}_1 + \dots + x[N-1]\mathbf{e}_{N-1}
$$

フーリエ変換の本質は，この「時間ごとのインパルス」という基底を，「周波数ごとの正弦波ベクトル（複素正弦波）」という別の直交基底に取り替える座標変換である．

信号 $\mathbf{x}$ の中に特定の周波数成分（基底ベクトル $\mathbf{u}_k$）がどれだけ含まれているかを知るためには，信号と基底の内積を取ればよい．複素ベクトル空間におけるエルミート内積を以下で定義する．

$$
\langle \mathbf{a}, \mathbf{b} \rangle = \sum_{n=0}^{N-1} a[n] \cdot \overline{b[n]}
$$

もし基底ベクトル $\mathbf{u}_k$ が正規化（長さ $1$）されていれば，信号 $\mathbf{x}$ との内積は，$\mathbf{x}$ の $\mathbf{u}_k$ 方向への射影成分（影の長さ）を与える．

$$
c_k = \langle \mathbf{x}, \mathbf{u}_k \rangle
$$

連続時間におけるフーリエ係数 $\int x(t) e^{-i \omega t} dt$ は，無限次元ベクトル $x(t)$ と基底ベクトル $e^{i\omega t}$ の内積操作そのものである．

## 2.2 リーマン和による連続への移行

離散和 $\sum$ から積分 $\int$ への移行において，微小時間 $dt$ は内積の定義において不可欠な役割を担う．無限次元の関数空間において，単に $f(t) \cdot \overline{g(t)}$ を無限個足し合わせれば値は発散しうるためである．

区間 $T$ を $N$ 等分した幅を $\Delta t = T/N$ とする．連続信号の内積は，この微小幅 $\Delta t$ を重みとしたリーマン和の極限として定義される．

$$
\langle f, g \rangle = \lim_{N \to \infty} \sum_{n=0}^{N-1} f(t_n) \cdot \overline{g(t_n)} \cdot \Delta t = \int_{0}^{T} f(t) \cdot \overline{g(t)} dt
$$

## 2.3 三角関数の直交性と正規化

インパルス基底を三角関数系 $\{1, \cos t, \sin t, \cos 2t, \dots\}$ に変換する場合，基底の「直交性」と「ノルム（長さ）」を確認する必要がある．
三角関数系は区間 $[0, 2\pi]$ において互いに直交する．

2.1節では基底ベクトルが正規化済みであると仮定したが，一般に基底ベクトル $\mathbf{v}$ の長さは $1$ とは限らない．その場合，射影成分（フーリエ係数）を求めるには，基底自身の長さの二乗（ノルムの二乗）で除算する必要がある．

$$
c_k = \frac{\langle \mathbf{x}, \mathbf{v} \rangle}{\langle \mathbf{v}, \mathbf{v} \rangle} = \frac{\langle \mathbf{x}, \mathbf{v} \rangle}{\|\mathbf{v}\|^2}
$$

交流成分 $\cos(nt), \sin(nt)$ の区間 $[0, 2\pi]$ におけるノルムの二乗は以下のようになる．

$$
\|\cos(nt)\|^2 = \int_{0}^{2\pi} \cos^2(nt) dt = \pi \quad (\sin(nt)\text{も同様})
$$

よって長さは $\sqrt{\pi}$ である．一方，直流成分（$1$）のノルムの二乗は，

$$
\|1\|^2 = \int_{0}^{2\pi} 1^2 dt = 2\pi
$$

となり，長さは $\sqrt{2\pi}$ である．これらを射影公式に代入することで，実フーリエ級数の係数公式が導かれる．

$$
a_n = \frac{1}{\pi}\int_{-\pi}^{\pi} f(t)\cos(nt) dt \quad (n \ge 1)
$$
$$
a_0 = \frac{1}{2\pi}\int_{-\pi}^{\pi} f(t) \cdot 1 dt
$$

## 2.4 オイラーの公式と複素指数関数基底

実フーリエ級数では $\cos, \sin$ という2種類の基底を管理する必要があるが，オイラーの公式 $e^{i\theta} = \cos\theta + i\sin\theta$ を用いることで，これらを「複素平面上を回転するベクトル」として統合できる．

基底を $u_n(t) = e^{int}$ とすれば，区間 $[0, 2\pi]$ におけるノルムの二乗は，

$$
\| e^{int} \|^2 = \int_{0}^{2\pi} e^{int} \overline{e^{int}} dt = \int_{0}^{2\pi} 1 dt = 2\pi
$$

となり，全ての整数 $n$ において一定値 $2\pi$ をとる．これにより，複素フーリエ係数 $c_n$ の公式は統一的に記述できる．

$$
c_n = \frac{1}{2\pi}\int_{0}^{2\pi} f(t)e^{-int} dt
$$

負の周波数 $n$ まで拡張することで，$\cos(nt) = \frac{e^{int}+e^{-int}}{2}$ のように，正負の回転の合成としてすべての実信号を表現可能となる．

## 2.5 離散化と直交性

連続信号 $x(t)$ を $N$ 点にサンプリングしたベクトル $\mathbf{x}=[x[0], \dots, x[N-1]]$ を考える．同様に基底関数 $e^{i\omega t}$ も離散化する．
周波数 $k$ の基底ベクトル $\mathbf{w}_k$ の第 $n$ 成分は以下のように書ける．

$$
w_k[n] = e^{i\frac{2\pi}{N}kn} \quad (n=0, 1, \dots, N-1)
$$

異なる周波数 $k$ と $l$ ($k \neq l$) を持つ2つの基底ベクトルの内積を計算する．

$$
\langle \mathbf{w}_k , \mathbf{w}_l \rangle = \sum_{n=0}^{N-1} e^{i\frac{2\pi}{N}kn} \cdot e^{-i\frac{2\pi}{N}ln} = \sum_{n=0}^{N-1} e^{i\frac{2\pi}{N}(k-l)n}
$$

ここで $r = e^{i\frac{2\pi}{N}(k-l)}$ と置くと，この式は初項 $1$，公比 $r$ の等比級数の和となる．$k \neq l$ のとき $r \neq 1$ であるが，

$$
r^N = \left(e^{i\frac{2\pi}{N}(k-l)}\right)^N = e^{i2\pi (k-l)} = 1
$$

となるため，等比級数の和の公式より，

$$
\text{Sum} = \frac{1-r^N}{1-r} = \frac{1-1}{1-r} = 0
$$

よって，離散化しても異なる周波数の基底ベクトルは互いに直交することが示された．

## 2.6 DFT行列

信号ベクトル $\mathbf{x}$ を直交基底 $\mathbf{w}_k$ の線形結合で表した際の係数（スペクトル）を $\mathbf{X}[k]$ とする．

$$
\mathbf{X}[k] = \langle \mathbf{x}, \mathbf{w}_k \rangle = \sum_{n=0}^{N-1} x[n] e^{-i\frac{2\pi}{N}kn}
$$

これを行列形式で記述すると，DFT行列 $\mathbf{F}$ が導かれる．回転因子 $W_N = e^{-i\frac{2\pi}{N}}$ を用いると，

$$
\begin{bmatrix} X[0] \\ X[1] \\ \vdots \\ X[N-1] \end{bmatrix} 
= 
\begin{bmatrix} 
1 & 1 & 1 & \dots & 1 \\ 
1 & W_N & W_N^2 & \dots & W_N^{N-1} \\ 
1 & W_N^2 & W_N^4 & \dots & W_N^{2(N-1)} \\ 
\vdots & \vdots & \vdots & \ddots & \vdots \\
1 & W_N^{N-1} & W_N^{2(N-1)} & \dots & W_N^{(N-1)(N-1)}
\end{bmatrix} 
\begin{bmatrix} x[0] \\ x[1] \\ \vdots \\ x[N-1] \end{bmatrix}
$$

この行列 $\mathbf{F}$ は，信号空間を時間軸から周波数軸へ回転させる座標変換行列である．$\mathbf{F}$ はヴァンデルモンド行列の一種であり，かつ定数倍の違いを除いてユニタリ行列の性質を持つ．

## 2.7 スペクトル漏れと窓関数

DFTは信号が周期的である（$x[0]$ と $x[N-1]$ が滑らかに接続される）ことを暗に仮定している．
しかし，実際の信号を有限区間で切り出すと，始端と終端に不連続なズレが生じることが多い．DFTはこの不連続性を表現するために，本来存在しない高周波成分まで動員して波形を合成しようとする．これが目的の周波数以外にエネルギーが拡散する「スペクトル漏れ」である．

これを抑制するために「窓関数」を用いる．ハニング窓（Hanning Window）は，信号の両端を滑らかにゼロに減衰させることで不連続性を解消する．

$$
w[n] = 0.5 - 0.5\cos\left(\frac{2\pi n}{N-1}\right)
$$

これにより，振幅分解能は多少犠牲になるが，スペクトル漏れによる信号の汚染を防ぐことができる．

## 2.8 逆変換（IDFT）とユニタリ性

周波数領域 $\mathbf{X}$ から時間領域 $\mathbf{x}$ を復元する逆変換は，行列 $\mathbf{F}$ の逆行列を用いて記述される．

$$
\mathbf{x} = \mathbf{F}^{-1}\mathbf{X}
$$

通常，逆行列の計算は計算コストが高いが，DFT行列の持つユニタリ性（に近い性質）により，以下の関係が成り立つ．
行列 $\mathbf{F}$ とそのエルミート共役（随伴行列）$\mathbf{F}^{\dagger}$ の積を考える．

$$
(\mathbf{F}^{\dagger}\mathbf{F})_{mn} = \sum_{k=0}^{N-1} e^{i\frac{2\pi}{N}k(m-n)} = N \delta_{mn}
$$

ここで $\delta_{mn}$ はクロネッカーのデルタである．これより $\mathbf{F}^{\dagger}\mathbf{F} = N \mathbf{I}$ （$\mathbf{I}$ は単位行列）となるため，

$$
\mathbf{F}^{-1} = \frac{1}{N}\mathbf{F}^{\dagger}
$$

すなわち，逆DFT行列は，元の行列の複素共役を取り，係数 $1/N$ を掛けるだけで求まる．これがIDFTの計算原理である．

# 3. 方法 (Methods)
## 3.1 実験環境・使用機器 (Environment & Equipment)
* OS: 
* Compiler: 
* Hardware: 

## 3.2 手順 (Procedure)
以下の手順で実験を行った．

# 4. 結果 (Results)
得られたデータを以下に示す．

# 5. 考察 (Discussion)
結果より，以下の知見が得られた．

# 6. 結論 (Conclusion)
以上の結果より，本実験の目的は達成されたと結論づけられる．

# 参考文献 (References)
[1] 
