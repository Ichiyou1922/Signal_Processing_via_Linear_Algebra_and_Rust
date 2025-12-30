use std::ops::{Add, Mul, Sub}; // Subを追加
use std::f64::consts::PI;

fn main() {
    // 1. 信号生成 (Original Signal)
    let n = 8; // 検証用なので少なめでOK
    let mut sig_data = Vec::new();
    for i in 0..n {
        let t = i as f64;
        // 複雑な信号を作る: 2.5Hz + 1.0Hz の合成波
        let val = (2.0 * PI * 2.5 * t / n as f64).sin() + 0.5 * (2.0 * PI * 1.0 * t / n as f64).cos();
        sig_data.push(Complex::new(val, 0.0));
    }
    let original = Signal::new(sig_data);

    // 2. DFT (Analysis)
    // 時間領域 -> 周波数領域
    let spectrum = original.dft();

    // 3. IDFT (Synthesis)
    // 周波数領域 -> 時間領域
    // スペクトルから元の波形を復元する
    let reconstructed = spectrum.idft();

    // 4. 検証 (Verification)
    // 元の信号と復元信号の差（誤差）を確認する
    println!("Index | Original (Re) | Recon (Re) | Error");
    println!("------+---------------+------------+-------------");
    
    let mut max_error = 0.0;
    for i in 0..n {
        let orig_re = original.data[i].re;
        let recon_re = reconstructed.data[i].re;
        let error = (orig_re - recon_re).abs();
        
        if error > max_error {
            max_error = error;
        }

        println!("{:5} | {:13.8} | {:10.8} | {:.4e}", i, orig_re, recon_re, error);
    }

    println!("-------------------------------------------------");
    println!("Max Reconstruction Error: {:.4e}", max_error);
    
    if max_error < 1e-10 {
        println!("Result: SUCCESS. F^(-1) * F = I is proved.");
    } else {
        println!("Result: FAILED. Check the implementation.");
    }
}

// --- 構造体定義 ---

#[derive(Debug, Clone, Copy)]
struct Complex {
    re: f64,
    im: f64,
}

impl Complex {
    fn new(re: f64, im: f64) -> Complex {
        Complex { re, im }
    }
    fn conjugate(&self) -> Complex {
        Complex { re: self.re, im: -1.0 * self.im }
    }
    fn norm(&self) -> f64 {
        (self.re * self.re + self.im * self.im).sqrt()
    }
}

impl Mul for Complex {
    type Output = Complex;
    fn mul(self, other: Complex) -> Complex {
        Complex {
            re: self.re * other.re - self.im * other.im,
            im: self.re * other.im + self.im * other.re,
        }
    }
}

impl Add for Complex {
    type Output = Complex;
    fn add(self, other: Complex) -> Complex {
        Complex {
            re: self.re + other.re,
            im: self.im + other.im,
        }
    }
}

// Subトレイトの実装（誤差計算用）
impl Sub for Complex {
    type Output = Complex;
    fn sub(self, other: Complex) -> Complex {
        Complex {
            re: self.re - other.re,
            im: self.im - other.im,
        }
    }
}

#[derive(Debug, Clone)]
struct Signal {
    data: Vec<Complex>,
}

impl Signal {
    fn new(data: Vec<Complex>) -> Signal {
        Signal { data }
    }

    fn inner_product(&self, other: &Signal) -> Result<Complex, String> {
        let mut sum = Complex::new(0.0, 0.0);
        for i in 0..self.data.len() {
            sum = sum + self.data[i] * other.data[i].conjugate();
        }
        Ok(sum)
    }

    fn basis(n: usize, k: usize) -> Signal {
        let mut data = Vec::with_capacity(n);
        for i in 0..n {
            let theta = 2.0 * PI * (k as f64) * (i as f64) / (n as f64);
            data.push(Complex::new(theta.cos(), theta.sin()));
        }
        Signal::new(data)
    }

    fn dft(&self) -> Signal {
        let n = self.data.len();
        let mut spectrum_data = Vec::with_capacity(n);
        for k in 0..n {
            let e_k = Signal::basis(n, k);
            let x_k = self.inner_product(&e_k).unwrap(); 
            spectrum_data.push(x_k);
        }
        Signal::new(spectrum_data)
    }
    
    // IDFT: x[n] = (1/N) * sum( X[k] * e_k )
    // スペクトル成分を重みとして，基底ベクトルを足し合わせる
    fn idft(&self) -> Signal {
        let n = self.data.len();
        // ゼロで初期化した信号ベクトルを作成
        let mut reconstructed_data = vec![Complex::new(0.0, 0.0); n];

        for k in 0..n {
            let X_k = self.data[k];       // スペクトル係数
            let e_k = Signal::basis(n, k); // 基底ベクトル (波)

            // 線形結合: vec += X_k * e_k
            for i in 0..n {
                // X_k * e_k[i]
                let term = X_k * e_k.data[i];
                reconstructed_data[i] = reconstructed_data[i] + term;
            }
        }

        // 最後に 1/N で正規化
        for i in 0..n {
            reconstructed_data[i].re /= n as f64;
            reconstructed_data[i].im /= n as f64;
        }

        Signal::new(reconstructed_data)
    }
}

// 既存のMul for Signal等は今回は使わないが，残しておいても問題ない
impl Mul for Signal {
    type Output = Signal;
    fn mul(self, other: Signal) -> Signal {
         let mut new_data = Vec::with_capacity(self.data.len());
        for i in 0..self.data.len() {
            new_data.push(self.data[i] * other.data[i]);
        }
        Signal::new(new_data)
    }
}
