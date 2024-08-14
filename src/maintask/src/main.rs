#![allow(non_upper_case_globals)]

mod ntt;

use std::array;
use std::f64::consts::PI;
use std::ops::{Add, DerefMut, Mul, Sub};
use std::sync::{LazyLock, Mutex};
use std::time::Duration;
use arrayref::array_ref;
use concrete_fft::c64;
use concrete_fft::ordered::{Method, Plan};
use dyn_stack::{GlobalPodBuffer, PodStack, ReborrowMut};
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
const kN: usize = k*N;
const l: usize = 2;
const t: usize = 5;
const base: usize = 4;
const basebit: usize = 2;
const Nbit: usize = 9;
const n: usize = 636;

fn  main() {
    let mut rng = ChaChaRng::from_entropy();
    let s0 = gen_s::<n>(&mut rng);
    let s1 = gen_s::<kN>(&mut rng);
    let bk = gen_bk(s0, s1);
    let ks = gen_ks(s0, s1);
    let timer_all = std::time::Instant::now();
    for itr in 0..50 {
        let timer = std::time::Instant::now();
        let v1 = (itr & 1) == 0;
        let v2 = (itr & 2) == 0;
        let v3 = (itr & 4) == 0;
        let ct1 = encrypt(v1, &s1);
        let ct2 = encrypt(v2, &s1);
        let ct3 = encrypt(v3, &s1);
        let (ct_s, ct_c) = hom_full_adder(&ct1, &ct2, &ct3, &bk, &ks);
        let pt_s = decrypt(&ct_s, &s1);
        let pt_c = decrypt(&ct_c, &s1);
        assert_eq!(v1 ^ v2 ^ v3, pt_s);
        assert_eq!((v1 & v2) | (v1 & v3) | (v2 & v3), pt_c);
        println!("{:?} | {}: {} + {} + {} = {} {}", timer.elapsed(), itr, v1, v2, v3, pt_c, pt_s);
    }
    println!("{:?}", timer_all.elapsed());
}

fn torus32_to_16(a: Torus32) -> Torus16 {
    (a + 0x8000 >> 16) as u16
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
    TLWELv1::new(Box::new(a), trlwe.b[x])
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
    let dec1 = trlwe.decrypt(&s1);
    for i in 0..N {
        let ext = sample_extract_index(&trlwe, i);
        let dec2 = decryptlv1(&ext, &s1);
        assert_torus32_eq(dec1[i], dec2, 27);
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

static FFT: LazyLock<Mutex<FFT>> = LazyLock::new(|| {
    let plan = Plan::new(N / 2, Method::Measure(Duration::from_millis(10)));
    let scratch_memory = GlobalPodBuffer::new(plan.fft_scratch().unwrap());
    Mutex::new(FFT {
        plan,
        scratch_memory,
    })
});
struct FFT {
    plan: Plan,
    scratch_memory: GlobalPodBuffer,
}
type TorusPolyFFT = Box<[c64; N/2]>;
impl FFT {
    fn torus32_to_f64(a: Torus32) -> f64 {
        a as f64 / 4294967296.
    }
    fn f64_to_torus32(a: f64) -> Torus32 {
        (a * 4294967296.) as i64 as Torus32
    }
    fn fft(&mut self, a: &TorusPoly) -> TorusPolyFFT {
        let mut stack = PodStack::new(&mut self.scratch_memory);
        let mut ac: [c64; N/2] = [Default::default(); N/2];
        for i in 0..N/2 {
            let w = c64::cis(PI * i as f64 / N as f64);
            ac[i] = c64::new(Self::torus32_to_f64(a[i]), Self::torus32_to_f64(a[i + N/2])) * w;
        }
        self.plan.fwd(&mut ac, stack.rb_mut());
        Box::new(ac)
    }
    fn ifft(&mut self, ac: &mut TorusPolyFFT) -> TorusPoly {
        let mut stack = PodStack::new(&mut self.scratch_memory);
        self.plan.inv(ac.deref_mut(), stack.rb_mut());
        let mut c = [0; N];
        for i in 0..N/2 {
            let w = c64::cis(-PI * i as f64 / N as f64);
            ac[i] *= w;
            c[i] = Self::f64_to_torus32(ac[i].re / (N/2) as f64);
            c[i + N/2] = Self::f64_to_torus32(ac[i].im / (N/2) as f64);
        }
        Box::new(c)
    }
    fn polymul_i8(&mut self, a: &[i8; N], bc: &TorusPolyFFT) -> TorusPolyFFT {
        let mut stack = PodStack::new(&mut self.scratch_memory);
        let mut ac: [c64; N/2] = array::from_fn(|i| {
            let w = c64::cis(PI * i as f64 / N as f64);
            c64::new(a[i] as f64, a[i + N/2] as f64) * w
        });
        self.plan.fwd(&mut ac, stack.rb_mut());
        for i in 0..N/2 {
            ac[i] = ac[i] * bc[i];
        }
        Box::new(ac)
    }
    fn polymul_bool(&mut self, a: &[bool; N], b: &TorusPoly) -> TorusPoly {
        let bc = self.fft(b);
        let mut ac = self.polymul_i8(&a.map(|v| v as i8), &bc);
        self.ifft(&mut ac)
    }
}

#[allow(unused)]
fn polymul_i8(a: &[i8; N], b: &TorusPoly) -> TorusPoly {
    let a: [i32; N] = a.map(|v| v as i32);
    let b: [i32; N] = b.map(|v| v as i32);
    let c = ntt::polymul(a, b);
    Box::new(c.map(|v| v as Torus32))
}

#[allow(unused)]
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
    let mut r: [[c64; N/2]; 3] = [[Default::default(); N/2]; 3];
    let mut fft = FFT.lock().unwrap();
    for i in 0..3 {
        for j in 0..l {
            let mul0 = fft.polymul_i8(&a0_decomp[j], &c.mat[j][i]);
            let mul1 = fft.polymul_i8(&a1_decomp[j], &c.mat[j+l][i]);
            let mul2 = fft.polymul_i8(&b_decomp[j], &c.mat[j+l+l][i]);
            for x in 0..N/2 {
                r[i][x] += mul0[x] + mul1[x] + mul2[x];
            }
        }
    }
    TRLWE {
        a0: fft.ifft(&mut Box::new(r[0])),
        a1: fft.ifft(&mut Box::new(r[1])),
        b: fft.ifft(&mut Box::new(r[2])),
    }
}

#[test]
fn test_external_product() {
    let mut rng = ChaChaRng::from_entropy();
    for val in [0, 1] {
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

fn blind_rotate(tlwe: &TLWELv0, bk: &BKey, test_vec: &TRLWE) -> TRLWE {
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
    let mut rng = ChaChaRng::from_entropy();
    for (val, expect) in [(MU16, MU32), (3 * MU16, MU32), (5 * MU16, 7 * MU32), (7 * MU16, 7 * MU32)] {
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

fn gate_bootstrapping_tlwe_to_tlwe(tlwe: &TLWELv0, bk: &BKey) -> TLWELv1 {
    let test_vec = TRLWE::new_test_vec();
    let trlwe = blind_rotate(tlwe, bk, &test_vec);
    sample_extract_index(&trlwe, 0)
}

#[test]
fn test_gate_bootstrapping_tlwe_to_tlwe() {
    let mut rng = ChaChaRng::from_entropy();
    for (val, expected) in [(MU16, MU16), (3 * MU16, MU16), (5 * MU16, 7 * MU16), (7 * MU16, 7 * MU16)] {
        let s0 = gen_s::<n>(&mut rng);
        let s1 = gen_s::<kN>(&mut rng);
        let bk = gen_bk(s0, s1);
        let lv0 = encrypt_lv0(val, &s0);
        let lv1 = gate_bootstrapping_tlwe_to_tlwe(&lv0, &bk);
        let dec = decryptlv1(&lv1, &s1);
        assert_torus16_eq(torus32_to_16(dec), expected, 11);
    }
}

fn identity_key_switching(tlwe: &TLWELv1, ks: &KSKey) -> TLWELv0 {
    let mut r = TLWELv0::new(Box::new([0; n]), torus32_to_16(tlwe.b));
    for i in 0..kN {
        let a2 = tlwe.a[i] as i64 + (1 << (31 - t * basebit)) + 0xaa800000;
        for m in 0..t {
            let o = ((a2 >> (32 - (m + 1) * basebit) & 3) - 2) as i8;
            if o > 0 {
                r = &r - &ks[i][m][(o-1) as usize];
            } else if o < 0 {
                r = &r + &ks[i][m][(-o-1) as usize];
            }
        }
    }
    r
}

#[test]
fn test_identity_key_switching() {
    let mut rng = ChaChaRng::from_entropy();
    for val in [MU32, 3 * MU32, 5 * MU32, 7 * MU32] {
        let s0 = gen_s::<n>(&mut rng);
        let s1 = gen_s::<kN>(&mut rng);
        let ks = gen_ks(s0, s1);
        let lv1 = encrypt_lv1(val, &s1);
        let lv0 = identity_key_switching(&lv1, &ks);
        let dec0 = decryptlv0(&lv0, &s0);
        assert_torus16_eq(dec0, torus32_to_16(val), 11);
    }
}

fn hom_nand(x: &TLWELv1, y: &TLWELv1, bk: &BKey, ks: &KSKey) -> TLWELv1 {
    let lv1 = &(&TLWELv1::new_const(MU32) - x) - y;
    let lv0 = identity_key_switching(&lv1, ks);
    gate_bootstrapping_tlwe_to_tlwe(&lv0, bk)
}
fn hom_and(x: &TLWELv1, y: &TLWELv1, bk: &BKey, ks: &KSKey) -> TLWELv1 {
    let lv1 = &(x + y) - &TLWELv1::new_const(MU32);
    let lv0 = identity_key_switching(&lv1, ks);
    gate_bootstrapping_tlwe_to_tlwe(&lv0, bk)
}
fn hom_or(x: &TLWELv1, y: &TLWELv1, bk: &BKey, ks: &KSKey) -> TLWELv1 {
    let lv1 = &(x + y) + &TLWELv1::new_const(MU32);
    let lv0 = identity_key_switching(&lv1, ks);
    gate_bootstrapping_tlwe_to_tlwe(&lv0, bk)
}
fn hom_xor(x: &TLWELv1, y: &TLWELv1, bk: &BKey, ks: &KSKey) -> TLWELv1 {
    let lv1 = &(&(x + y) + &TLWELv1::new_const(MU32)) * 2;
    let lv0 = identity_key_switching(&lv1, ks);
    gate_bootstrapping_tlwe_to_tlwe(&lv0, bk)
}
fn hom_full_adder(x: &TLWELv1, y: &TLWELv1, c0: &TLWELv1, bk: &BKey, ks: &KSKey) -> (TLWELv1, TLWELv1) {
    let xor = hom_xor(x, y, bk, ks);
    let s = hom_xor(&xor, c0, bk, ks);
    let t1 = hom_and(&xor, c0, bk, ks);
    let t2 = hom_and(x, y, bk, ks);
    let c1 = hom_or(&t1, &t2, bk, ks);
    (s, c1)
}

#[test]
fn test_hom_op() {
    let mut rng = ChaChaRng::from_entropy();
    let s0 = gen_s::<n>(&mut rng);
    let s1 = gen_s::<kN>(&mut rng);
    let bk = gen_bk(s0, s1);
    let ks = gen_ks(s0, s1);
    let ops: &[(fn (bool, bool) -> bool, fn (&TLWELv1, &TLWELv1, &BKey, &KSKey) -> TLWELv1)] = &[
        (|x, y| !(x & y), hom_nand),
        (|x, y| x & y, hom_and),
        (|x, y| x | y, hom_or),
        (|x, y| x ^ y, hom_xor)
    ];
    for (op, func) in ops {
        for (v1, v2) in [(false, false), (false, true), (true, false), (true, true)] {
            let ct1 = encrypt(v1, &s1);
            let ct2 = encrypt(v2, &s1);
            let ct3 = func(&ct1, &ct2, &bk, &ks);
            let pt = decrypt(&ct3, &s1);
            assert_eq!(op(v1, v2), pt);
        }
    }
}

#[test]
fn test_full_adder() {
    for (v1, v2, v3) in [(false, false, false), (false, true, false), (true, false, false), (true, true, false), (false, false, true), (false, true, true), (true, false, true), (true, true, true)] {
        let mut rng = ChaChaRng::from_entropy();
        let s0 = gen_s::<n>(&mut rng);
        let s1 = gen_s::<kN>(&mut rng);
        let ct1 = encrypt(v1, &s1);
        let ct2 = encrypt(v2, &s1);
        let ct3 = encrypt(v3, &s1);
        let bk = gen_bk(s0, s1);
        let ks = gen_ks(s0, s1);
        let (ct_s, ct_c) = hom_full_adder(&ct1, &ct2, &ct3, &bk, &ks);
        let pt_s = decrypt(&ct_s, &s1);
        let pt_c = decrypt(&ct_c, &s1);
        assert_eq!(v1 ^ v2 ^ v3, pt_s);
        assert_eq!((v1 & v2) | (v1 & v3) | (v2 & v3), pt_c);
    }
}

fn gen_bk(s0: [bool; n], s1: [bool; kN]) -> BKey {
    s0.map(|b| TRGSW::new(b as u32, s1))
}

fn gen_ks(s0: [bool; n], s1: [bool; kN]) -> KSKey {
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

fn gen_nonce<T, const m: usize>(rng: &mut ChaChaRng) -> Box<[T; m]> where Standard: Distribution<T> {
    Box::new([0; m].map(|_| rng.gen()))
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
type BKey = [TRGSW; n];
type KSKey = [[[Box<TLWELv0>; (1<<basebit)/2]; t]; kN];

type TLWELv0 = TLWE<Torus16, n>;
type TLWELv1 = TLWE<Torus32, kN>;
#[derive(Debug)]
struct TLWE<T, const M: usize> {
    a: Box<[T; M]>,
    b: T,
}
impl<T, const M: usize> TLWE<T, M> {
    pub fn new(a: Box<[T; M]>, b: T) -> Self {
        Self {
            a,
            b,
        }
    }
}
impl<const M: usize> TLWE<Torus32, M> {
    pub fn new_const(x: Torus32) -> Self {
        Self::new(Box::new([0; M]), x)
    }
}
impl<'a, 'b, T, const M: usize> Add<&'b TLWE<T, M>> for &'a TLWE<T, M> where T: Unsigned + Copy + Default {
    type Output = TLWE<T, M>;

    fn add(self, rhs: &'b TLWE<T, M>) -> Self::Output {
        let mut a = [T::default(); M];
        for i in 0..M {
            a[i] = self.a[i] + rhs.a[i];
        }
        TLWE::new(Box::new(a), self.b + rhs.b)
    }
}
impl<'a, 'b, T, const M: usize> Sub<&'b TLWE<T, M>> for &'a TLWE<T, M> where T: Unsigned + Copy + Default {
    type Output = TLWE<T, M>;

    fn sub(self, rhs: &'b TLWE<T, M>) -> Self::Output {
        let mut a = [T::default(); M];
        for i in 0..M {
            a[i] = self.a[i] - rhs.a[i];
        }
        TLWE::new(Box::new(a), self.b - rhs.b)
    }
}
impl<'a, T, const M: usize> Mul<T> for &'a TLWE<T, M> where T: Unsigned + Copy + Default {
    type Output = TLWE<T, M>;

    fn mul(self, rhs: T) -> Self::Output {
        let mut a = [T::default(); M];
        for i in 0..M {
            a[i] = self.a[i] * rhs;
        }
        TLWE::new(Box::new(a), self.b * rhs)
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
        let a0 = gen_nonce::<Torus32, N>(&mut rng);
        let a1 = gen_nonce::<Torus32, N>(&mut rng);
        let e: TorusPoly = Box::new(array::from_fn(|_| gen_cbd(&mut rng)));
        let s0 = array_ref!(s, 0, N);
        let s1 = array_ref!(s, N, N);
        let mut fft = FFT.lock().unwrap();
        let as0 = fft.polymul_bool(s0, &a0);
        let as1 = fft.polymul_bool(s1, &a1);
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
        let a0 = gen_nonce::<Torus32, N>(&mut rng);
        let a1 = gen_nonce::<Torus32, N>(&mut rng);
        let e: TorusPoly = Box::new(array::from_fn(|_| gen_cbd(&mut rng)));
        let s0 = array_ref!(s, 0, N);
        let s1 = array_ref!(s, N, N);
        let mut fft = FFT.lock().unwrap();
        let as0 = fft.polymul_bool(s0, &a0);
        let as1 = fft.polymul_bool(s1, &a1);
        let mut b = [0; N];
        for i in 0..N {
            b[i] = as0[i] + as1[i] + e[i] + poly[i];
        }
        TRLWE::new(a0, a1, Box::new(b))
    }
    fn decrypt(&self, s: &[bool; kN]) -> TorusPoly {
        let s0 = array_ref!(s, 0, N);
        let s1 = array_ref!(s, N, N);
        let mut fft = FFT.lock().unwrap();
        let as0 = fft.polymul_bool(s0, &self.a0);
        let as1 = fft.polymul_bool(s1, &self.a1);
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
    mat: [[TorusPolyFFT; 3]; l * 3],
}
impl TRGSW {
    pub fn new(mu: u32, s: [bool; kN]) -> Self {
        let mat = array::from_fn(|i| {
            let mut zero_enc = TRLWE::new_zero(s);
            if i < l {
                zero_enc.a0[0] += mu << (24 - i * 8);
            } else if i < l * 2 {
                zero_enc.a1[0] += mu << (24 - (i-l) * 8);
            } else {
                zero_enc.b[0] += mu << (24 - (i-l-l) * 8);
            }
            let mut fft = FFT.lock().unwrap();
            [fft.fft(&zero_enc.a0), fft.fft(&zero_enc.a1), fft.fft(&zero_enc.b)]
        });
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
