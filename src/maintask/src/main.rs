#![allow(non_upper_case_globals)]

mod ntt;

use std::array;
use std::ops::{Add, Sub};
use arrayref::array_ref;
use num_traits::Unsigned;
use rand::{Rng, RngCore};
use rand::distributions::Standard;
use rand::prelude::Distribution;
use rand_chacha::ChaChaRng;
use rand_chacha::rand_core::SeedableRng;

const MU16: Torus16 = 1 << 13;
const MU32: Torus32 = 1 << 29;
const ETA: usize = 40;
const k: usize = 2;
const N: usize = 512;
const kN: usize = k * N;
const l: usize = 2;
const t: usize = 5;
const base: usize = 4;
const basebit: usize = 2;
const Nbit: usize = 9;
const n: usize = 636;

fn  main() {
    for itr in 0..50 {
        let timer = std::time::Instant::now();
        let v0 = (itr & 1) == 0;
        let v1 = (itr & 2) == 0;
        let mut rng = ChaChaRng::from_entropy();
        let s0 = gen_s::<n>(&mut rng);
        let s1 = gen_s::<kN>(&mut rng);
        let ct1 = encrypt(v0, &s1);
        let ct2 = encrypt(v1, &s1);
        let ct3 = hom_nand(&ct1, &ct2, &gen_bk(s0, s1), gen_ks(s0, s1));
        let pt0 = decryptlv1(&ct3, &s1);
        let pt = decrypt(&ct3, &s1);
        assert_eq!(!(v0 & v1), pt);
        println!("{:?}: {:?} nand {:?} = {:?} ({:016b})", itr, v0, v1, pt, pt0);
        println!("{:?}", timer.elapsed());
    }
}

fn torus32_to_16(a: Torus32) -> Torus16 {
    (a >> 16) as u16
}

fn sample_extract_index(trlwe: &TRLWE, x: usize) -> TLWELv1 {
    let mut a = [0; kN];
    for j in 0..=x {
        a[j] = trlwe.a0[x - j];
        a[N + j] = trlwe.a1[x - j];
    }
    for j in x+1..N {
        a[j] = trlwe.a0[N + x - j].wrapping_neg();
        a[N + j] = trlwe.a1[N + x - j].wrapping_neg();
    }
    TLWELv1::new(a, trlwe.b[x])
}

#[test]
fn test_sample_extract_index() {
    let mut rng = ChaChaRng::from_entropy();
    let s0 = gen_s::<n>(&mut rng);
    let s1 = gen_s::<kN>(&mut rng);
    let test_vec = TRLWE::new_test_vec();
    let bk = gen_bk(s0, s1);
    let lv0 = encrypt_lv0(MU16.wrapping_neg(), &s0);
    let trlwe = blind_rotate(&lv0, &bk, &test_vec);
    let dec0 = trlwe.decrypt(&s1);
    for i in 0..N {
        let ext = sample_extract_index(&trlwe, i);
        let dec1 = decryptlv1(&ext, &s1);
        assert_torus32_eq(dec0[i], dec1, 27);
    }
}

fn decomposition<const m: usize>(a: u32) -> [i8; m] {
    let a2 = a as i64 + 0x80808080;
    let mut r = [0; m];
    for i in 0..m {
        r[i] = ((a2 >> (24 - 8 * i) & 255) - 128) as i8;
    }
    r
}

fn decomposition_poly<const m: usize>(a: &TorusPoly) -> [[i8; N]; m] {
    let mut r = [[0; N]; m];
    for i in 0..N {
        let decomp = decomposition::<m>(a[i]);
        for j in 0..m {
            r[j][i] = decomp[j];
        }
    }
    r
}

fn polymul_i8(a: &[i8; N], b: &TorusPoly) -> TorusPoly {
    let a: [i32; N] = a.map(|v| v as i32);
    let b: [i32; N] = b.map(|v| v as i32);
    let c = ntt::polymul(a, b);
    Box::new(c.map(|v| v as Torus32))
}

fn polymul_bool(a: &[bool; N], b: &TorusPoly) -> TorusPoly {
    let a: [i32; N] = a.map(|v| v as i32);
    let b: [i32; N] = b.map(|v| v as i32);
    let c = ntt::polymul(a, b);
    Box::new(c.map(|v| v as Torus32))
}

fn external_product(c: &TRGSW, trlwe: &TRLWE) -> TRLWE {
    let a0_decomp = decomposition_poly::<l>(&trlwe.a0);
    let a1_decomp = decomposition_poly::<l>(&trlwe.a1);
    let b_decomp = decomposition_poly::<l>(&trlwe.b);
    let mut r = [[0; N]; 3];
    for i in 0..3 {
        for j in 0..l {
            let mul0 = polymul_i8(&a0_decomp[j], &c.mat[j][i]);
            let mul1 = polymul_i8(&a1_decomp[j], &c.mat[j+l][i]);
            let mul2 = polymul_i8(&b_decomp[j], &c.mat[j+l+l][i]);
            for x in 0..N {
                r[i][x] += mul0[x] + mul1[x] + mul2[x];
            }
        }
    }
    TRLWE {
        a0: Box::new(r[0]),
        a1: Box::new(r[1]),
        b: Box::new(r[2]),
    }
}

#[test]
fn test_external_product() {
    for val in [0, 1] {
        let mut rng = ChaChaRng::from_entropy();
        let s = gen_s::<kN>(&mut rng);
        let mut pt = [0; N];
        for i in 0..N {
            pt[i] = i as Torus32 * MU32;
        }
        let trlwe = TRLWE::encrypt(&Box::new(pt), s);
        let c = TRGSW::new(val, s);
        let res = external_product(&c, &trlwe);
        let pt1 = trlwe.decrypt(&s);
        let pt2 = res.decrypt(&s);
        for i in 0..N {
            assert_torus32_eq(pt1[i] * val, pt2[i], 27);
        }
    }
}

fn cmux(trgsw: &TRGSW, trlwe1: &TRLWE, trlwe0: &TRLWE) -> TRLWE {
    &external_product(trgsw, &(trlwe1 - trlwe0)) + trlwe0
}

fn blind_rotate(tlwe: &TLWELv0, bk: &[TRGSW; n], test_vec: &TRLWE) -> TRLWE {
    let sh = 16 - Nbit - 1;
    let b = (tlwe.b >> sh) as isize;
    let mut acc = test_vec.mul_x(-b);
    for i in 0..n {
        let a = (tlwe.a[i] + (1 << sh - 1) >> sh) as isize;
        acc = cmux(&bk[i], &acc.mul_x(a), &acc);
    }
    acc
}

#[test]
fn test_blind_rotate() {
    for (val, expect) in [(MU16, MU32), (3 * MU16, MU32), (5 * MU16, 7 * MU32), (7 * MU16, 7 * MU32)] {
        let mut rng = ChaChaRng::from_entropy();
        let test_vec = TRLWE::new_test_vec();
        let s0 = gen_s::<n>(&mut rng);
        let s1 = gen_s::<kN>(&mut rng);
        let bk = gen_bk(s0, s1);
        let lv0 = encrypt_lv0(val, &s0);
        let trlwe = blind_rotate(&lv0, &bk, &test_vec);
        let dec = trlwe.decrypt(&s1);
        assert_torus32_eq(dec[0], expect, 27);
    }
}

fn gate_bootstrapping_tlwe_to_tlwe(tlwe: &TLWELv0, bk: &[TRGSW; n]) -> TLWELv1 {
    let test_vec = TRLWE::new_test_vec();
    let trlwe = blind_rotate(tlwe, bk, &test_vec);
    sample_extract_index(&trlwe, 0)
}

#[test]
fn test_gate_bootstrapping_tlwe_to_tlwe() {
    for (val, expected) in [(MU16, MU16), (3 * MU16, MU16), (5 * MU16, 7 * MU16), (7 * MU16, 7 * MU16)] {
        let mut rng = ChaChaRng::from_entropy();
        let s0 = gen_s::<n>(&mut rng);
        let s1 = gen_s::<kN>(&mut rng);
        let bk = gen_bk(s0, s1);
        let lv0 = encrypt_lv0(val, &s0);
        let lv1 = gate_bootstrapping_tlwe_to_tlwe(&lv0, &bk);
        let dec = decryptlv1(&lv1, &s1);
        assert_torus16_eq(torus32_to_16(dec), expected, 11);
    }
}

fn identity_key_switching(tlwe: &TLWELv1, ks: [[[Box<TLWELv0>; base-1]; t]; kN]) -> TLWELv0 {
    let mut r = TLWELv0::new([0; n], torus32_to_16(tlwe.b));
    let mask = (base - 1) as u32;
    for i in 0..kN {
        for m in 0..t {
            let sh = 32 - (m + 1) * basebit;
            let o = ((tlwe.a[i] + (1 << (31 - t * basebit))) >> sh & mask) as usize;
            if o > 0 {
                r = &r - &ks[i][m][o-1];
            }
        }
    }
    r
}

#[test]
fn test_identity_key_switching() {
    for val in [MU32, 3 * MU32, 5 * MU32, 7 * MU32] {
        let mut rng = ChaChaRng::from_entropy();
        let s0 = gen_s::<n>(&mut rng);
        let s1 = gen_s::<kN>(&mut rng);
        let ks = gen_ks(s0, s1);
        let lv1 = encrypt_lv1(val, &s1);
        let lv0 = identity_key_switching(&lv1, ks);
        let dec0 = decryptlv0(&lv0, &s0);
        assert_torus16_eq(dec0, torus32_to_16(val), 11);
    }
}

fn hom_nand(x: &TLWELv1, y: &TLWELv1, bk: &[TRGSW; n], ks: [[[Box<TLWELv0>; base-1]; t]; kN]) -> TLWELv1 {
    let lv1 = &(&TLWELv1::new_mu() - x) - y;
    let lv0 = identity_key_switching(&lv1, ks);
    gate_bootstrapping_tlwe_to_tlwe(&lv0, bk)
}

#[test]
fn test_hom_nand() {
    for (v0, v1) in [(false, false), (false, true), (true, false), (true, true)] {
        let mut rng = ChaChaRng::from_entropy();
        let s0 = gen_s::<n>(&mut rng);
        let s1 = gen_s::<kN>(&mut rng);
        let ct1 = encrypt(v0, &s1);
        let ct2 = encrypt(v1, &s1);
        let ct3 = hom_nand(&ct1, &ct2, &gen_bk(s0, s1), gen_ks(s0, s1));
        let pt = decrypt(&ct3, &s1);
        assert_eq!(!(v0 & v1), pt);
    }
}

fn gen_bk(s0: [bool; n], s1: [bool; kN]) -> [TRGSW; n] {
    s0.map(|b| TRGSW::new(b as u32, s1))
}

fn gen_ks(s0: [bool; n], s1: [bool; kN]) -> [[[Box<TLWELv0>; base-1]; t]; kN] {
    array::from_fn(|i| {
        array::from_fn(|m| {
            let sh = 16 - (m + 1) * basebit;
            array::from_fn(|o| {
                let pt = (s1[i] as Torus16 * (o + 1) as Torus16) << sh;
                Box::new(encrypt_lv0(pt, &s0))
            })
        })
    })
}

fn gen_nonce<T, const m: usize>(rng: &mut ChaChaRng) -> [T; m] where Standard: Distribution<T> {
    [0; m].map(|_| rng.gen())
}

fn gen_cbd(rng: &mut ChaChaRng) -> Torus32 {
    let mask = (1 << ETA) - 1;
    let a = (rng.next_u64() & mask).count_ones();
    let b = (rng.next_u64() & mask).count_ones();
    (a - b) as Torus32
}

fn gen_s<const m: usize>(rng: &mut ChaChaRng) -> [bool; m] {
    [0; m].map(|_| (rng.next_u32() & 1) == 1)
}

fn encrypt(m: bool, s: &[bool; kN]) -> TLWELv1 {
    let mut rng = ChaChaRng::from_entropy();
    let a = gen_nonce::<Torus32, kN>(&mut rng);
    let e = gen_cbd(&mut rng);
    let prod = (0..kN).filter(|&i| s[i]).map(|i| a[i]).sum::<Torus32>();
    if m {
        TLWELv1::new(a, prod + MU32 + e)
    } else {
        TLWELv1::new(a, prod - MU32 + e)
    }
}
fn encrypt_lv0(m: Torus16, s: &[bool; n]) -> TLWELv0 {
    let mut rng = ChaChaRng::from_entropy();
    let a = gen_nonce::<Torus16, n>(&mut rng);
    let e = gen_cbd(&mut rng) as Torus16;
    let prod = (0..n).filter(|&i| s[i]).map(|i| a[i]).sum::<Torus16>();
    TLWELv0::new(a, prod + m + e)
}
fn encrypt_lv1(m: Torus32, s: &[bool; kN]) -> TLWELv1 {
    let mut rng = ChaChaRng::from_entropy();
    let a = gen_nonce::<Torus32, kN>(&mut rng);
    let e = gen_cbd(&mut rng);
    let prod = (0..kN).filter(|&i| s[i]).map(|i| a[i]).sum::<Torus32>();
    TLWELv1::new(a, prod + m + e)
}

fn decrypt(tlwe: &TLWELv1, s: &[bool; kN]) -> bool {
    let prod = (0..kN).filter(|&i| s[i]).map(|i| tlwe.a[i]).sum::<Torus32>();
    ((tlwe.b - prod) >> 31 & 1) == 0
}
fn decryptlv0(tlwe: &TLWELv0, s: &[bool; n]) -> Torus16 {
    let prod = (0..n).filter(|&i| s[i]).map(|i| tlwe.a[i]).sum::<Torus16>();
    tlwe.b - prod
}
fn decryptlv1(tlwe: &TLWELv1, s: &[bool; kN]) -> Torus32 {
    let prod = (0..kN).filter(|&i| s[i]).map(|i| tlwe.a[i]).sum::<Torus32>();
    tlwe.b - prod
}

type Torus16 = u16;
type Torus32 = u32;
type TorusPoly = Box<[Torus32; N]>;

type TLWELv0 = TLWE<Torus16, n>;
type TLWELv1 = TLWE<Torus32, kN>;
#[derive(Debug)]
struct TLWE<T, const M: usize> {
    a: [T; M],
    b: T,
}
impl<T, const M: usize> TLWE<T, M> {
    pub fn new(a: [T; M], b: T) -> Self {
        Self {
            a,
            b,
        }
    }
}
impl<const M: usize> TLWE<Torus16, M> {
    pub fn new_mu() -> Self {
        Self::new([0; M], MU16)
    }
}
impl<const M: usize> TLWE<Torus32, M> {
    pub fn new_mu() -> Self {
        Self::new([0; M], MU32)
    }
}
impl<'a, 'b, T, const M: usize> Add<&'b TLWE<T, M>> for &'a TLWE<T, M> where T: Unsigned + Copy + Default {
    type Output = TLWE<T, M>;

    fn add(self, rhs: &'b TLWE<T, M>) -> Self::Output {
        let mut a = [T::default(); M];
        for i in 0..M {
            a[i] = self.a[i] + rhs.a[i];
        }
        TLWE::new(a, self.b + rhs.b)
    }
}
impl<'a, 'b, T, const M: usize> Sub<&'b TLWE<T, M>> for &'a TLWE<T, M> where T: Unsigned + Copy + Default {
    type Output = TLWE<T, M>;

    fn sub(self, rhs: &'b TLWE<T, M>) -> Self::Output {
        let mut a = [T::default(); M];
        for i in 0..M {
            a[i] = self.a[i] - rhs.a[i];
        }
        TLWE::new(a, self.b - rhs.b)
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
    pub fn mul_x(&self, d: isize) -> Self {
        let mut a0 = [0; N];
        let mut a1 = [0; N];
        let mut b = [0; N];
        for i in 0..N {
            let j = (i as isize - d).rem_euclid(2 * N as isize) as usize;
            if j < N {
                a0[i] = self.a0[j];
                a1[i] = self.a1[j];
                b[i] = self.b[j];
            } else {
                a0[i] = self.a0[j-N].wrapping_neg();
                a1[i] = self.a1[j-N].wrapping_neg();
                b[i] = self.b[j-N].wrapping_neg();
            }
        }
        TRLWE::new(Box::new(a0), Box::new(a1), Box::new(b))
    }
    pub fn new_zero(s: [bool; kN]) -> Self {
        let mut rng = ChaChaRng::from_entropy();
        let a0 = Box::new(gen_nonce::<Torus32, N>(&mut rng));
        let a1 = Box::new(gen_nonce::<Torus32, N>(&mut rng));
        let e: TorusPoly = Box::new(array::from_fn(|_| gen_cbd(&mut rng)));
        let s0 = array_ref!(s, 0, N);
        let s1 = array_ref!(s, N, N);
        let as0 = polymul_bool(s0, &a0);
        let as1 = polymul_bool(s1, &a1);
        let mut b = [0; N];
        for i in 0..N {
            b[i] = as0[i] + as1[i] + e[i];
        }
        TRLWE::new(a0, a1, Box::new(b))
    }
    pub fn new_test_vec() -> Self {
        let a0 = [0; N];
        let a1 = [0; N];
        let b = [MU32; N];
        TRLWE::new(Box::new(a0), Box::new(a1), Box::new(b))
    }
    pub fn encrypt(poly: &TorusPoly, s: [bool; kN]) -> Self {
        let mut rng = ChaChaRng::from_entropy();
        let a0 = Box::new(gen_nonce::<Torus32, N>(&mut rng));
        let a1 = Box::new(gen_nonce::<Torus32, N>(&mut rng));
        let e: TorusPoly = Box::new(array::from_fn(|_| gen_cbd(&mut rng)));
        let s0 = array_ref!(s, 0, N);
        let s1 = array_ref!(s, N, N);
        let as0 = polymul_bool(s0, &a0);
        let as1 = polymul_bool(s1, &a1);
        let mut b = [0; N];
        for i in 0..N {
            b[i] = as0[i] + as1[i] + e[i] + poly[i];
        }
        TRLWE::new(a0, a1, Box::new(b))
    }
    fn decrypt(&self, s: &[bool; kN]) -> TorusPoly {
        let s0 = array_ref!(s, 0, N);
        let s1 = array_ref!(s, N, N);
        let as0 = polymul_bool(s0, &self.a0);
        let as1 = polymul_bool(s1, &self.a1);
        let mut pt = [0; N];
        for i in 0..N {
            pt[i] = self.b[i] - as0[i] - as1[i];
        }
        Box::new(pt)
    }
}
impl <'a, 'b> Add<&'b TRLWE> for &'a TRLWE {
    type Output = TRLWE;

    fn add(self, rhs: &'b TRLWE) -> Self::Output {
        let mut a0 = [0; N];
        let mut a1 = [0; N];
        let mut b = [0; N];
        for i in 0..N {
            a0[i] = self.a0[i] + rhs.a0[i];
            a1[i] = self.a1[i] + rhs.a1[i];
            b[i] = self.b[i] + rhs.b[i];
        }
        TRLWE::new(Box::new(a0), Box::new(a1), Box::new(b))
    }
}
impl <'a, 'b> Sub<&'b TRLWE> for &'a TRLWE {
    type Output = TRLWE;

    fn sub(self, rhs: &'b TRLWE) -> Self::Output {
        let mut a0 = [0; N];
        let mut a1 = [0; N];
        let mut b = [0; N];
        for i in 0..N {
            a0[i] = self.a0[i] - rhs.a0[i];
            a1[i] = self.a1[i] - rhs.a1[i];
            b[i] = self.b[i] - rhs.b[i];
        }
        TRLWE::new(Box::new(a0), Box::new(a1), Box::new(b))
    }
}

struct TRGSW {
    mat: [[TorusPoly; 3]; l * 3],
}
impl TRGSW {
    pub fn new(mu: u32, s: [bool; kN]) -> Self {
        let mut mat = array::from_fn(|_| {
            let zero_enc = TRLWE::new_zero(s);
            [zero_enc.a0, zero_enc.a1, zero_enc.b]
        });
        for i in 0..l {
            mat[i][0][0] += mu << (24 - i * 8);
            mat[i+l][1][0] += mu << (24 - i * 8);
            mat[i+l+l][2][0] += mu << (24 - i * 8);
        }
        TRGSW {
            mat,
        }
    }
}

fn assert_torus16_eq(a: Torus16, b: Torus16, err: usize) {
    assert!(((a - b) as i16).abs() < 1 << err, "{a:016b} != {b:016b}");
}
fn assert_torus32_eq(a: Torus32, b: Torus32, err: usize) {
    assert!(((a - b) as i32).abs() < 1 << err, "{a:032b} != {b:032b}");
}
