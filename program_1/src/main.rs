use std::process::exit;
use std::ops::{Add, AddAssign, Mul};

fn main() {
    println!("Hello, world!");
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

struct Signal {
    data: Vec<Complex>,
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
}
