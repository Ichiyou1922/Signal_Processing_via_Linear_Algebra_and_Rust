use std::ops::{Add, Mul};
use std::f64::consts::{PI};

fn main() {
    println!("Hello, world!");
    let mut sig = Vec::new();
    
    for i in 0..8 {
        sig.push(Complex::new((2.0 * PI * 2.5 * i as f64 / 8.0).sin(), 0.0));
    }
    let mut spectrum = Signal::new(sig);
    spectrum = spectrum * Signal::create_hanning(sig.data.len());
    spectrum = Signal::dft(&spectrum);
    for i in 0..8 {
        println!("data{}: ({})+i({})", i, spectrum.data[i].re, spectrum.data[i].im);
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
}

#[derive(Debug, Clone, Copy)]
struct Signal {
    data: Vec<Complex>,
}

impl Mul for Signal {
    type Output: Signal;
    fn mul (self, other: Signal) -> Signal {
        for i in 0..self.data.len() {
            self.data[i] = self.data[i] * other.data[i];
        }
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
            data.push(Complex::new((2.0 * PI * k as f64 * i as f64 / n as f64).cos(), (2.0 * PI * k as f64 * i as f64 / n as f64).sin()));
        }

        Signal::new(data)
    }

    // DFT: X[k] = <x, e_k>
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

    fn create_hanning(&self, n: usize) -> Signal {
        let mut hanning = Vec::with_capacity(n);
        for i in 0..n {
            hanning.push(Complex::new(0.5 - 0.5 * 2.0 * PI * i as f64 / (n as f64 - 1.0), 0.0));
        }

        Signal::new(hanning)
    }
}
