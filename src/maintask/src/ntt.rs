#![allow(dead_code, unused_imports)]

use std::{env, ops};
use std::cmp::min;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use rand::{Rng, RngCore};
use rand_chacha::ChaChaRng;
use rand_chacha::rand_core::SeedableRng;

const MOD: u128 = 0x1000000000000e001;
const R2: u128 = 0x9612ae0498038001;
const M2: u128 = 0x7f36f589911a834e69f0ab7f3c00dfff;
const PRIM_ROOT: u128 = 13;
const NTH_ROOT: [Mint; 13] = {
    let mut nth_root = [Mint::new(0); 13];
    nth_root[12] = Mint::new(0x1f924307b4d8a243);
    let mut i = 12;
    while i > 0 {
        nth_root[i - 1] = Mint(Mint::modulo(Mint::mul_reduce(nth_root[i].0, nth_root[i].0)));
        i -= 1;
    }
    nth_root
};

#[derive(Clone, Copy)]
struct Mint(u128);
impl Mint {
    const fn new(x: u128) -> Self {
        Mint(Self::mul_reduce(x % MOD, R2))
    }

    fn new_i64(x: i64) -> Self {
        let v = if x >= 0 { x as u128 } else { (x as i128 + MOD as i128) as u128 };
        Mint(Self::mul_reduce(v, R2))
    }

    const fn high(x: u128, y: u128) -> u128 {
        let xh = x >> 64;
        let yh = y >> 64;
        let xl = x & 0xFFFFFFFFFFFFFFFF;
        let yl = y & 0xFFFFFFFFFFFFFFFF;
        let a = xh * yl + (xl * yl >> 64);
        let b = xl * yh + (a & 0xFFFFFFFFFFFFFFFF);
        xh * yh + (a >> 64) + (b >> 64)
    }

    const fn reduce(x: u128) -> u128 {
        Self::high(M2 * x, MOD) + (x != 0) as u128
    }

    const fn mul_reduce(x: u128, y: u128) -> u128 {
        Self::high(x, y) + Self::high(M2 * x * y, MOD) + (x * y != 0) as u128
    }

    const fn modulo(x: u128) -> u128 {
        if x < x - MOD { x } else { x - MOD }
    }

    fn pow(self, n: u128) -> Mint {
        let mut z = self;
        let mut r = Mint::new(1);
        let mut n = n;
        while n > 0 {
            if (n & 1) == 1 {
                r = r * z;
            }
            z = z * z;
            n >>= 1;
        }
        r
    }

    fn inv(self) -> Mint {
        self.pow(MOD - 2)
    }
}
impl ops::Add<Mint> for Mint {
    type Output = Mint;

    fn add(self, rhs: Mint) -> Mint {
        Mint(Self::modulo(self.0 + rhs.0))
    }
}
impl ops::Sub<Mint> for Mint {
    type Output = Mint;

    fn sub(self, rhs: Mint) -> Mint {
        Mint(Self::modulo(self.0 - rhs.0 + MOD))
    }
}
impl ops::Mul<Mint> for Mint {
    type Output = Mint;

    fn mul(self, rhs: Mint) -> Mint {
        Mint(Self::modulo(Self::mul_reduce(self.0, rhs.0)))
    }
}
impl From<Mint> for u128 {
    fn from(value: Mint) -> Self {
        Mint::reduce(value.0)
    }
}
impl From<Mint> for i64 {
    fn from(value: Mint) -> Self {
        let x = Mint::reduce(value.0);
        if x < MOD / 2 { x as i64 } else { (x as i128 - MOD as i128) as i64 }
    }
}

fn nth_root(n: usize) -> Mint {
    NTH_ROOT[n.ilog2() as usize]
}

fn ntt<const N: usize>(a: &mut Vec<Mint>) {
    let mut width = 2 * N;
    let mut offset = N;
    while width > 1 {
        let w = nth_root(width);
        for top in (0..2*N).step_by(width) {
            let mut root = Mint::new(1);
            for i in top..top+offset {
                let (c0, c1) = (a[i], a[i + offset]);
                a[i] = c0 + c1;
                a[i + offset] = (c0 - c1) * root;
                root = root * w;
            }
        }
        width >>= 1;
        offset >>= 1;
    }
}

fn intt<const N: usize>(a: &mut Vec<Mint>) {
    let mut width = 2;
    let mut offset = 1;
    while width <= 2 * N {
        let w = nth_root(width).inv();
        for top in (0..2*N).step_by(width) {
            let mut root = Mint::new(1);
            for i in top..top+offset {
                let (c0, c1) = (a[i], a[i + offset]);
                a[i] = c0 + c1 * root;
                a[i + offset] = c0 - c1 * root;
                root = root * w;
            }
        }
        width <<= 1;
        offset <<= 1;
    }
    let b = Mint::new(2 * N as u128).inv();
    for i in 0..2*N {
        a[i] = a[i] * b;
    }
}

pub fn polymul<const N: usize>(a: [i32; N], b: [i32; N]) -> [i64; N] {
    let mut a: Vec<Mint> = a.iter().map(|&v| Mint::new_i64(v as i64)).collect();
    let mut b: Vec<Mint> = b.iter().map(|&v| Mint::new_i64(v as i64)).collect();
    a.resize(2 * N, Mint::new(0));
    b.resize(2 * N, Mint::new(0));
    ntt::<N>(&mut a);
    ntt::<N>(&mut b);
    for i in 0..2*N {
        a[i] = a[i] * b[i];
    }
    intt::<N>(&mut a);
    let mut c = [Mint::new(0); N];
    for i in 0..2*N {
        if i < N {
            c[i] = c[i] + a[i];
        } else {
            c[i - N] = c[i - N] - a[i];
        }
    }
    c.map(|v| v.into())
}
