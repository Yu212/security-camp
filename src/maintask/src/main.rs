use std::ops::{Add, Sub};
use rand_chacha::ChaChaRng;

const MU: Torus32 = 1 << 29;
const k: usize = 2;
const N: usize = 512;
const kN: usize = k * N;
const Bgbit: usize = 8;
const l: usize = 2;
const t: usize = 5;
const basebit: usize = 2;

fn main() {
}

fn sample_extract_index(trlwe: &TRLWE) -> TLWELv1 {
    todo!()
}

fn decomposition<const n: usize>(a: Torus32) -> [i8; n] {
    todo!()
}

fn external_product(c: &TRGSW, trlwe: &TRLWE) -> TRLWE {
    todo!()
}

fn cmux(trgsw: &TRGSW, trlwe0: &TRLWE, trlwe1: &TRLWE) -> TRLWE {
    todo!()
}

fn blind_rotate(tlwe: &TLWELv0, bk: &Vec<TRGSW>, test_vec: &TRLWE) -> TRLWE {
    todo!()
}

fn gate_bootstrapping_tlwe_to_tlwe(tlwe: &TLWELv0, bk: &Vec<TRGSW>) -> TLWELv1 {
    todo!()
}

fn identity_key_switching(tlwe: &TLWELv1, ks: &TLWELv0) -> TLWELv0 {
    todo!()
}

fn hom_nand(x: &TLWELv1, y: &TLWELv1) -> TLWELv1 {
    todo!()
}

fn gen_s(rng: &mut ChaChaRng) -> Vec<bool> {
    todo!()
}

fn encrypt(m: bool, s: &Vec<bool>) -> TLWELv1 {
    todo!()
}

fn decrypt(tlwe: &TLWELv1, s: &Vec<bool>) -> bool {
    todo!()
}

type Torus16 = u16;
type Torus32 = u32;
type TorusPoly = [Torus32; N];

type TLWELv0 = TLWE<636>;
type TLWELv1 = TLWE<kN>;
struct TLWE<const M: usize> {
    a: [Torus16; M],
    b: Torus16,
}
impl<const M: usize> TLWE<M> {
    pub fn new(a: [Torus16; M], b: Torus16) -> Self {
        Self {
            a,
            b,
        }
    }
}
impl<const M: usize> Add<TLWE<M>> for TLWE<M> {
    type Output = TLWE<M>;

    fn add(self, rhs: Self) -> Self::Output {
        let mut a = [0; M];
        for i in 0..M {
            a[i] = self.a[i] + rhs.a[i];
        }
        Self::new(a, self.b + rhs.b)
    }
}
impl<const M: usize> Sub<TLWE<M>> for TLWE<M> {
    type Output = TLWE<M>;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut a = [0; M];
        for i in 0..M {
            a[i] = self.a[i] - rhs.a[i];
        }
        Self::new(a, self.b - rhs.b)
    }
}

struct TRLWE {
    a0: TorusPoly,
    a1: TorusPoly,
    b: TorusPoly,
}
impl TRLWE {
    pub fn new(a0: TorusPoly, a1: TorusPoly, b: TorusPoly) -> Self {
        return Self {
            a0,
            a1,
            b,
        }
    }
    pub fn new_zero() -> Self {
        todo!()
    }
}
impl Add<TRLWE> for TRLWE {
    type Output = TRLWE;

    fn add(self, rhs: Self) -> Self::Output {
        let mut a0 = [0; N];
        let mut a1 = [0; N];
        let mut b = [0; N];
        for i in 0..N {
            a0[i] = self.a0[i] + rhs.a0[i];
            a1[i] = self.a1[i] + rhs.a1[i];
            b[i] = self.b[i] + rhs.b[i];
        }
        Self::new(a0, a1, b)
    }
}
impl Sub<TRLWE> for TRLWE {
    type Output = TRLWE;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut a0 = [0; N];
        let mut a1 = [0; N];
        let mut b = [0; N];
        for i in 0..N {
            a0[i] = self.a0[i] - rhs.a0[i];
            a1[i] = self.a1[i] - rhs.a1[i];
            b[i] = self.b[i] - rhs.b[i];
        }
        Self::new(a0, a1, b)
    }
}

struct TRGSW {
    mat: [[TorusPoly; 3]; l * 3],
}
impl TRGSW {
    pub fn new(mu: Torus32) -> Self {
        todo!()
    }
}
