use std::ops;
use rand_chacha::ChaChaRng;

const MU: u32 = 1 << 29;
const k: usize = 2;
const N: usize = 512;
const kN: usize = k * N;
const Bgbit: usize = 8;
const l: usize = 2;
const t: usize = 5;
const basebit: usize = 2;

fn main() {
}

fn sample_extract_index(trlwe: TRLWE) -> TLWELv1 {
    todo!()
}

fn decomposition(a: usize) -> [i8; l] {
    todo!()
}

fn external_product(c: TRGSW, trlwe: TRLWE) -> TRLWE {
    todo!()
}

fn cmux(trgsw: TRGSW, trlwe0: TRLWE, trlwe1: TRLWE) -> TRLWE {
    todo!()
}

fn blind_rotate(tlwe: TLWELv0, bk: Vec<TRGSW>, test_vec: TRLWE) -> TRLWE {
    todo!()
}

fn gate_bootstrapping_tlwe_to_tlwe(tlwe: TLWELv0, bk: Vec<TRGSW>) -> TLWELv1 {
    todo!()
}

fn identity_key_switching(tlwe: TLWELv1, ks: TLWELv0) -> TLWELv0 {
    todo!()
}

fn hom_nand(x: TLWELv1, y: TLWELv1) -> TLWELv1 {
    todo!()
}

fn encrypt(m: bool, s: &Vec<bool>) -> TLWELv1 {
    todo!()
}

fn decrypt(tlwe: TLWELv1, s: &Vec<bool>) -> bool {
    todo!()
}

type TLWELv0 = TLWE<636>;
type TLWELv1 = TLWE<kN>;
struct TLWE<const M: usize> {
    a: [u16; M],
    b: u32,
}
impl<const M: usize> ops::Add<TLWE<M>> for TLWE<M> {
    type Output = TLWE<M>;

    fn add(self, rhs: Self) -> Self::Output {
        todo!()
    }
}
impl<const M: usize> ops::Sub<TLWE<M>> for TLWE<M> {
    type Output = TLWE<M>;

    fn sub(self, rhs: Self) -> Self::Output {
        todo!()
    }
}

struct TRLWE {
    a0: [u32; N],
    a1: [u32; N],
    b: [u32; N],
}
impl ops::Add<TRLWE> for TRLWE {
    type Output = TRLWE;

    fn add(self, rhs: Self) -> Self::Output {
        todo!()
    }
}
impl ops::Sub<TRLWE> for TRLWE {
    type Output = TRLWE;

    fn sub(self, rhs: Self) -> Self::Output {
        todo!()
    }
}

struct TRGSW {
}
