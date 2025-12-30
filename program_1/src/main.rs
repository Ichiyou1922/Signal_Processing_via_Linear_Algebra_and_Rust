use std::ops::{Add, Mul};
use std::f64::consts::{PI};

fn main() {
    println!("Hello, world!");
    
    let mut sig_data = Vec::new();
    let n = 8;
    for i in 0..n {
        let t = i as f64;
        sig_data.push(Complex::new((2.0 * PI * 2.5 * t / 8.0).sin(), 0.0));
    }
    let raw_signal = Signal::new(sig_data);

    // ハニング窓を適用
    let hanning_window = Signal::create_hanning(n);
    let windowed_signal = raw_signal.clone() * hanning_window;

    // DFT Transform
    let spectrum = windowed_signal.dft();

    println!("--- Spectrum (With Hanning Window) ---");
    for k in 0..n {
        let amp = spectrum.data[k].norm(); // 振幅表示
        println!("k={}: Re={:.4}, Im={:.4}, Mag={:.4}",
            k, spectrum.data[k].re, spectrum.data[k].im, amp);
    }
}

#[derive(Debug, Clone, Copy)]
struct Complex {
    re: f64,
    im: f64,
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

impl Complex {
    fn new(re: f64, im: f64) -> Complex {
        Complex {
            re,
            im,
        }
    }

    fn conjugate(&self) -> Complex {
        Complex {
            re: self.re,
            im: -1.0 * self.im,
        }
    }

    fn norm(&self) -> f64 {
        (self.re * self.re + self.im * self.im).sqrt()
    }
}

#[derive(Debug, Clone)]
struct Signal {
    data: Vec<Complex>,
}

// アダマール積
impl Mul for Signal {
    type Output = Signal;

    fn mul(self, other: Signal) -> Signal {
        if self.data.len() != other.data.len() {
            panic!("Vector dimension mismatch in Hadamard product");
        }

        let mut new_data = Vec::with_capacity(self.data.len());
        for i in 0..self.data.len() {
            new_data.push(self.data[i] * other.data[i]);
        }
        Signal::new(new_data)
    }
}

impl Signal {
    fn new(data: Vec<Complex>) -> Signal {
        Signal {
            data,
        }
    }

    fn inner_product(&self, other: &Signal) -> Result<Complex, String> {
        if self.data.len() != other.data.len() {
           return Err("Vector sizes differ".to_string());
        }

        let mut sum = Complex::new(0.0, 0.0);

        for i in 0..self.data.len() {
            let term = self.data[i] * other.data[i].conjugate();
            sum = sum + term;
        }
        Ok(sum)
    }

    fn basis(n: usize, k: usize) -> Signal {
        // 周波数 k の基底ベクトル e_k を生成
        let mut data = Vec::with_capacity(n);
        for i in 0..n {
            let theta = 2.0 * PI * (k as f64) * (i as f64) / (n as f64);
            data.push(Complex::new(theta.cos(), theta.sin()));
        }
        Signal::new(data)
    }

    // DFT: X[k] = <x, e_k>
    // 信号xを基底e_kに射影する
    fn dft(&self) -> Signal {
        let n = self.data.len();
        let mut spectrum_data = Vec::with_capacity(n);

        for k in 0..n {
            let e_k = Signal::basis(n, k);
            let x_k = self.inner_product(&e_k).unwrap(); // <x, e_k>
            spectrum_data.push(x_k);
        }

        Signal::new(spectrum_data)
    }

    // ハニング窓の生成
    // w[n] = 0.5 - 0.5 * cos(2πn / (N-1))
    fn create_hanning(n: usize) -> Signal {
        let mut hanning = Vec::with_capacity(n);
        for i in 0..n {
            hanning.push(Complex::new(0.5 - 0.5 * (2.0 * PI * i as f64 / (n as f64 - 1.0)).cos(), 0.0));
        }
        Signal::new(hanning)
    }
}
