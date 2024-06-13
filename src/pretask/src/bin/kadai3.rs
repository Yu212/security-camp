#![allow(dead_code, unused_imports)]

use rand::{Rng, RngCore};
use rand_chacha::ChaChaRng;
use rand_chacha::rand_core::SeedableRng;

const N: usize = 512;
const Q: usize = 0x10000;
const ETA: usize = 40;

fn gen_nonce(rng: &mut ChaChaRng) -> Vec<u16> {
    (0..N).map(|_| rng.gen_range(0..Q) as u16).collect()
}

fn gen_cbd(rng: &mut ChaChaRng) -> u16 {
    let mask = !(!0 << ETA);
    // 40ビットの整数値の立っているビットの数から二項分布を生成
    let a = (rng.next_u64() & mask).count_ones();
    let b = (rng.next_u64() & mask).count_ones();
    (a - b) as u16
}

fn gen_s(rng: &mut ChaChaRng) -> Vec<bool> {
    (0..N).map(|_| (rng.next_u32() & 1) == 1).collect::<Vec<bool>>()
}

// 暗号化
fn encrypt(m: bool, s: &Vec<bool>) -> (Vec<u16>, u16) {
    let mut rng = ChaChaRng::from_entropy();
    let a = gen_nonce(&mut rng);
    let e = gen_cbd(&mut rng);
    let prod = (0..N).filter(|&i| s[i]).map(|i| a[i]).sum::<u16>();
    (a, prod + e + m as u16 * (Q / 2) as u16)
}

// 復号
fn decrypt(a: Vec<u16>, b: u16, s: &Vec<bool>) -> bool {
    let prod = (0..N).filter(|&i| s[i]).map(|i| a[i]).sum::<u16>();
    // e が負の時最上位ビットが反転してしまうため1/4だけずらす
    let v = b - prod + (Q / 4) as u16;
    // 最上位ビットを取得
    (v >> 15 & 1) == 1
}

// 暗号文二つの加算
fn add_ciphertext(a0: Vec<u16>, b0: u16, a1: Vec<u16>, b1: u16) -> (Vec<u16>, u16) {
    (a0.iter().zip(a1).map(|t| t.0 + t.1).collect(), b0 + b1)
}

fn main() {
}

#[test]
fn test_lwe_0() {
    // テストのため乱数を固定
    let mut rng = ChaChaRng::seed_from_u64(0);
    let s = gen_s(&mut rng);
    let (a, b) = encrypt(false, &s);
    let pt = decrypt(a, b, &s);
    assert_eq!(pt, false);
}

#[test]
fn test_lwe_1() {
    let mut rng = ChaChaRng::seed_from_u64(0);
    let s = gen_s(&mut rng);
    let (a, b) = encrypt(true, &s);
    let pt = decrypt(a, b, &s);
    assert_eq!(pt, true);
}

// すべての組み合わせで加法準同型性が成り立つこと確認
#[test]
fn test_lwe_additive() {
    let mut rng = ChaChaRng::seed_from_u64(0);
    for (x, y) in [(false, false), (false, true), (true, false), (true, true)] {
        let s = gen_s(&mut rng);
        let (a0, b0) = encrypt(x, &s);
        let (a1, b1) = encrypt(y, &s);
        let (a2, b2) = add_ciphertext(a0, b0, a1, b1);
        let pt = decrypt(a2, b2, &s);
        assert_eq!(pt, x ^ y);
    }
}
