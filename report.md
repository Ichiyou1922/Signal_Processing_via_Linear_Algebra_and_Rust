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
本実験で使用する基本原理について述べる．
## 2.1 フーリエ変換における基底変換
任意の信号 $\mathbf{x}$ は「時刻t=kのとき値が1」であるインパルス信号 $\mathbf{e_k}$ を用いて線型結合で表せる．
```math
\mathbf{x}=x[0]\mathbf{e}_0+x[1]\mathbf{e}_1+\dots +x[N-1]\mathbf{e}_{N-1}\in \mathbb{C}^N
```

フーリエ変換では，「時間ごとのインパルス」という基底を「周波数ごとの正弦波ベクトル」に取り替える．

信号 $\mathbf{x}$ の中に特定の成分(例: 基本波 $\mathbf{v}$ )がどれだけ含まれているのか知るためにはこの2つのベクトルの内積を取れば良い．

エルミート内積を以下で定義する．
```math
\langle \mathbf{a}, \mathbf{b} \rangle = \Sum_{n=0}^{N-1}a[n]\cdot \bar b[n]
```

もし基底ベクトル $\mathbf{u}_k$ が長さ1(正規化済み)であれば，信号 $\mathbf{x}$ との内積は， $\mathbf{x}$ の $\mathbf{u}_k$ 方向への影の長さ(成分の大きさ)を与える．
```math
成分c_k=\langle \mathbf{x}, \mathbf{u}_k \rangle
```

フーリエ係数 $\int x(t) e^{-i \omega t} dt$ とは無限次元ベクトル $x(t)$ と基底ベクトル $e^{-i\omega t}$ の内積を表していると言える．

## 2.2 リーマン和
離散和 $\Sum$ から積分 $\int$ への移行において， $dt$ は不可欠な役割を担っている．

無限次元の関数空間において，単に $f(t)\cdot g(t)$ を無限個足し合わせれば，値は発散しうる．

これを防ぐのが「リーマン和」の思想である．

区間TをN等分した幅を $\Delta t = T/N$ とする．

連続信号の内積は，この微小幅 $\Delta t$ をかけた長方形の面積の総和の極限として定義される．
```math
\langle f, g \rangle = \Sum_{n=0}^{N-1}f(t_n)\cdot \bar g(t_n)\cdot \Delta t = \int_{0}^{T}f(t)\cdot \bar g(t) dt
```

## 2.3 三角関数の直交性と長さからフーリエ係数の公式へ
インパルス基底を三角関数系{1, cost, sint, cos2t, \dots}に変換したい．ここで，三角関数系は直交しているのか，長さはどうなっているのかという点が重要である．

三角関数系は区間 $[0, 2\pi]$ において直交していることが知られているため，この点は問題ない．また後述するオイラーの公式による三角関数の指数変換により簡単に示すことができる．

2.1においては影の長さを求める際に基底ベクトルが正規化済みだと仮定した．しかし多くの場合，基底ベクトルの長さは1ではない．よって，より一般的に影の長さを求めるためには以下のように記述する必要がある．
```math
成分c_k = \frac{\langle \mathbf{x}, \mathbf{v}}{\langle \mathbf{v}, \mathbf{v} \rangle}=\frac{\langle \mathbf{x}, \mathbf{v}}{\|v \| ^2}
```

交流成分( $cos(nt), sin(nt)$ )の長さについて調べる．
```math
\|cosnt \| ^2 = \int_{-\infty}^{\infty}cos^2(nt) dt=\pi \quad (\text{sin(nt)も同様})
```

よって長さは $\sqrt{\pi}$ ．

直流成分(1)の長さについて調べる．
```math
\|1 \| ^2 = \int_{-\infty}^{\infty}1^2 dt = 2\pi
```

よって長さは $\sqrt{2\pi$ ．

前述した射影公式に代入することでフーリエ係数の公式は導かれる．
```math
a_n = \frac{1}{\pi}\int_{-\pi}^{\pi}f(t)cos(nt)dt \quad (n \ge 1)
```

```math
a_0 = \frac{1}{2\pi}\int_{-\pi}^{\pi}f(t)\cdot 1 dt
```

## 2.4 オイラーの公式 
実フーリエ係数では， $cos(nt), sin(nt)$ という2つの基底を管理しなければならなかった．しかし，オイラーの公式( $e^{i\theta}=cos\theta + isin\theta$ )はこれらを「回転するベクトル」として統合する．

基底を $u_n (t)=e^{int}$ とすれば，区間 $[0, 2\pi]$ における長さは，
```math
\| e^{int} \| = \int_{0}^{2\pi}e^{int} \bar e^{int} dt=\int_{0}^{2\pi}a dt=2\pi
```

驚くべきことに，複素指数関数基底はすべてのnで長さが $\2pi$ で一定である．故に係数 $c_n$ の公式は，
```math
c_n=\frac{1}{2\pi}\int_{0}^{2\pi}f(t)e^{-int} dt 
```

として統一される．更にこれを利用して複素フーリエ級数は，
```math
f(t)=\Sum_{n=-\infty}^{\infty}c_n e^{int}
```

として表せる．

負の整数nまで考えているが，これにより三角関数系に含まれるすべての成分をまとめて記述できる．なぜならすべての三角関数は $cos(nt)=\frac{e^{int}+e^{-int}}{2}, sin(nt)=\frac{e^{int}-e^{-int}}{2i}$ と記述でき，正の回転，負の回転の組み合わせにより一意に決定されるためである．

## 2.5 離散化
連続信号 $x(t)$ を区間 $[0, 2\pi]$ で考える変わりに，これをN個の点にサンプリングしたベクトル $\mathbf{x}=[x[0], \dots, x[N-1]]$ を考える．

基底関数 $e^{i\omega t}$ も同様に離散化する．

周波数kの基底ベクトル $\mathbf{w}_k$ の第n成分は以下のように書ける．
```math
w_k[n]=e^{i\frac{2\pi}{N}kn}\quad (n=0, 1, \dots, N-1)
```

異なる周波数kとl( $k\neq l$ )を持つ2つの基底ベクトルの内積を取ると，
```math
\langle \mathbf{w}_k , \mathbf{w}_l \rangle = \Sum_{n=0}^{N-1}e^{i\frac{2\pi}{N}kn}\cdot \bar e^{i\frac{2\pi}{N}ln}=\Sum_{n=0}^{N-1}e^{i\frac{2\pi}{N}(k-l)n}
```

ここで $r=e^{i\frac{2\pi}{N}(l-k)}$ と置くと，この式は，初項1, 項比rの等比級数の和になる．
```math
\text{Sum}=1+e+r^2+\dots + r^{N-1}=\frac{1-r^N}{1-r}
```

ここで
```math
r^N=(e^{i\frac{2\pi}{N}(k-l)})^N = e^{i2\pi (k-l)}=1
```

よって，内積はゼロであり，離散化しても異なる周波数の波は直交することがわかる．

## 2.6 DFT行列
任意の信号ベクトル $\mathbf{x}$ を，これらの直交基底 $\mathbf{w}_k$ の線型結合で表したい．係数（スペクトル）を $\mathbf{X}[k]$ とすると，
```math
\mathbf{X}[k]=\langle \mathbf{x}, mathbf{w}_k \rangle = \Sum_{n=0}^{N-1}x[n]e^{-i\frac{2\pi}{N}kn}
```

これを行列形式で書くと，巨大な行列 $\mathbf{F}$ が現れる．
```math
\begin{bmatrix} X[0] \\ \vdots \\ X[N-1] \end{bmatrix} = \begin{bmatrix} 1 & 1 & 1 & \dots \\ 1 & W & W^2 & \dots \\ 1 & W^2 & W^4 & \dots \\ \vdots & \vdots & \vdots & \ddots \end{bmatrix} \begin{bmatrix} x[0] \\ \vdots \\ x[N-1] \end{bmatrix}
```

この行列 $\mathbf{F}$ こそが，信号空間を時間軸から，周波数軸へ回転させる座標変換行列（DFT行列）である．また，この行列はヴァンデルモンド行列であり，かつユニタリ行列でもある．この性質はIDFT(逆離散フーリエ変換)を考える際に効力を発揮する．

## 2.7 スペクトル漏れと窓関数
DFT行列の周期構造は，信号ベクトル $\mathbf{x}$ が $x_0$, $x_{N-1}$ において，なめらかに接続されていることを要請する．

しかし，多くの信号はこの要請を満たせない．このとき $x_o$ と $x_{N-1}$ にはズレが発生し，DFTはこのズレの造形をすべての周波数の基底ベクトルを合成することで再現する．これが目的の周波数以外にエネルギーが飛び散るスペクトル漏れの実態である．

このズレを抹消するためにハニング窓を使用する．ハニング窓への要求は信号の両端を強制的にゼロに押しつぶすことである．これにより，元の信号の振幅情報は一部失われるが，スペクトル漏れによる汚染からは守られる．この要求を満たす数式として，ハニング窓 
```math
\mathbf{w}(n)=0.5-0.5cos(\frac{2\pi n}{N-1})
```

が定義される．

## 2.8 逆変換（IDFT）
$\mathbf{X}$ から $\mathbf{x}$ を求めたい．数式的には，
```math
\mathbf{x}=\mathbf{F}^{-1}\mathbf{X}
```

により求められるはずだ．

通常巨大な行列の逆行列を求めることは困難であるが，DFT行列の持つユニタリ性により逆行列に近い行列を簡単に求めることができる．

$\mathbf{F}$ とそのエルミート共役 $\mathbf{F}^{\dagger}$ を掛けてみる．
```math
(F^{\dagger}F)_{mn}=\Sum_{k=0}^{N-1}e^{i\frac{2\pi}{N}k(m-n)}=\begin{bmatrix} N & 0 & \dots \\ 0 & N & \dots \\ \vdots & \vdots & \ddots \end{bmatrix} = N I
```

両辺をNで割れば，
```math
\frac{1}{N}\mathbf{F}^\dagger \mathbf{F} = I \Rightarrow \mathbf{F}^{-1}=\frac{1}{N}\mathbf{F}^\dagger
```

すなわち，逆DFT行列は，元の行列の共役を取り， $\frac{1}{N}$ を掛けたものに等しい．
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
