#![allow(dead_code, unused_imports)]

use std::{env, ops};
use std::cmp::min;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use rand::{Rng, RngCore};
use rand_chacha::ChaChaRng;
use rand_chacha::rand_core::SeedableRng;

// NNT を用いる
// 複数のmodで答えを求めて garner で復元する方法が一般的だが、ここでは法として2^64以上の素数を用いて一度のNTTで答えを求める
// オーバーフローしない場合正しい値が得られることを保証できる
// この際二数の積がu128に収まらないことがあるため、モンゴメリ乗算を用いて剰余をとる
// 参考(自分の記事): https://yu212.hatenablog.com/entry/2023/12/14/203400

// この課題に限らず、Cargo.tomlにて overflow-checks = false としてオーバーフローがエラーになることを抑制している。

const N: usize = 1024;
const MOD: u128 = 0x1000000000000e001; // MOD > 2^64, prime
const R2: u128 = 0x9612ae0498038001; // R2 = (2^128)^2 = 2^256 % MOD
const M2: u128 = 0x7f36f589911a834e69f0ab7f3c00dfff; // MOD * M2 = -1
const PRIM_ROOT: u128 = 13; // 1の原始 MOD 乗根

// 内部では x の代わりに x * 2^128 % MOD を持つ
#[derive(Clone, Copy)]
struct Mint(u128);
impl Mint {
    fn new(x: u128) -> Self {
        // x * 2^-128 * 2^256 % MOD = x * 2^128 % MOD
        Mint(Self::mul_reduce(x % MOD, R2))
    }

    fn new_i64(x: i64) -> Self {
        let v = if x >= 0 { x as u128 } else { (x as i128 + MOD as i128) as u128 };
        Mint(Self::mul_reduce(v, R2))
    }

    // x * y >> 128 を返す
    fn high(x: u128, y: u128) -> u128 {
        let xh = x >> 64;
        let yh = y >> 64;
        let xl = x & 0xFFFFFFFFFFFFFFFF;
        let yl = y & 0xFFFFFFFFFFFFFFFF;
        let a = xh * yl + (xl * yl >> 64);
        let b = xl * yh + (a & 0xFFFFFFFFFFFFFFFF);
        xh * yh + (a >> 64) + (b >> 64)
    }

    // x のモンゴメリ表現 x * 2^-128 % MOD を返す
    fn reduce(x: u128) -> u128 {
        Self::high(M2 * x, MOD) + (x != 0) as u128
    }

    // x * y のモンゴメリ表現 x * y * 2^-128 % MOD を返す
    fn mul_reduce(x: u128, y: u128) -> u128 {
        Self::high(x, y) + Self::high(M2 * x * y, MOD) + (x * y != 0) as u128
    }

    // x in [0, MOD * 2) を MOD で割った余りを返す
    // 符号なし整数なので x < MOD だとオーバーフローして x < x - MOD になることを利用する
    fn modulo(x: u128) -> u128 {
        min(x, x - MOD)
    }

    // self ^ n を返す
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

    // self ^ -1 を返す
    // フェルマーの小定理より self ^ (MOD - 2) % MOD = self ^ -1 % MOD
    fn inv(self) -> Mint {
        self.pow(MOD - 2)
    }
}
// 各演算子のオーバーロード
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
        // (a * 2^128 * b * 2^128) * 2^-128 % MOD = a * b * 2^128 % MOD
        Mint(Self::modulo(Self::mul_reduce(self.0, rhs.0)))
    }
}
// u128, i64への変換
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

// MOD の原始 n 乗根
// メモ化することもできるがここでは毎回求めている
fn nth_root(n: usize) -> Mint {
    Mint::new(PRIM_ROOT).pow((MOD - 1) / n as u128)
}

// NTT (周波数間引き)
// ref: https://tayu0110.hatenablog.com/entry/2023/05/06/023244
fn ntt(a: &mut Vec<Mint>) {
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

// inverse NTT (時間間引き)
fn intt(a: &mut Vec<Mint>) {
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
    for i in 0..2*N {
        a[i] = a[i] * Mint::new(2 * N as u128).inv();
    }
}

// NTT O(N log N)
fn polymul(a: [i32; N], b: [i32; N]) -> [i64; N] {
    // 扱いやすい型に変換
    let mut a: Vec<Mint> = a.iter().map(|&v| Mint::new_i64(v as i64)).collect();
    let mut b: Vec<Mint> = b.iter().map(|&v| Mint::new_i64(v as i64)).collect();
    a.resize(2 * N, Mint::new(0));
    b.resize(2 * N, Mint::new(0));
    ntt(&mut a);
    ntt(&mut b);
    // フーリエ変換したものの各点積
    for i in 0..2*N {
        a[i] = a[i] * b[i];
    }
    intt(&mut a);
    let mut c = [Mint::new(0); N];
    for i in 0..2*N {
        // mod X^n + 1 において X^n = -1 のため、
        // i+j >= N のとき、X^{i+j} = X^N * X^{i+j-N} = -X^{i+j-N}
        if i < N {
            c[i] = c[i] + a[i];
        } else {
            c[i - N] = c[i - N] - a[i];
        }
    }
    c.map(|v| v.into())
}

// naive O(N^2)
fn polymul_naive(a: [i32; N], b: [i32; N]) -> [i64; N] {
    let mut c = [0; N];
    for i in 0..N {
        for j in 0..N {
            if i + j < N {
                c[i + j] += a[i] as i64 * b[j] as i64;
            } else {
                c[i + j - N] -= a[i] as i64 * b[j] as i64;
            }
        }
    }
    c
}

// naive with mod
fn polymul_naive_mod(a: [i32; N], b: [i32; N]) -> [i64; N] {
    let a: Vec<Mint> = a.iter().map(|&v| Mint::new_i64(v as i64)).collect();
    let b: Vec<Mint> = b.iter().map(|&v| Mint::new_i64(v as i64)).collect();
    let mut c = [Mint::new(0); N];
    for i in 0..N {
        for j in 0..N {
            if i + j < N {
                c[i + j] = c[i + j] + a[i] * b[j];
            } else {
                c[i + j - N] = c[i + j - N] - a[i] * b[j];
            }
        }
    }
    c.map(|v| v.into())
}

fn main() {
    let args: Vec<String> = env::args().collect();

    let mut a : [i32; N] = [0; N];
    let mut b : [i32; N] = [0; N];

    {
        let path =  Path::new(&args[1]);
        let display = path.display();
        //パスを指定してファイルを開く
        let f = match File::open(&path){
            Err(why) => panic!("couldn't open {}: {}", display,
                               why.to_string()),
            Ok(file) => file,
        };
        let mut reader = BufReader::new(f);
        for i in a.iter_mut() {
            let mut buffer = String::new();
            reader.read_line(&mut buffer).unwrap();
            *i = buffer.trim().parse::<i32>().unwrap();
        }
    }

    {
        let path =  Path::new(&args[2]);
        let display = path.display();
        //パスを指定してファイルを開く
        let f = match File::open(&path){
            Err(why) => panic!("couldn't open {}: {}", display,
                               why.to_string()),
            Ok(file) => file,
        };
        let mut reader = BufReader::new(f);
        for i in b.iter_mut() {
            let mut buffer = String::new();
            reader.read_line(&mut buffer).unwrap();
            *i = buffer.trim().parse::<i32>().unwrap();
        }
    }

    let res: [i64; N] = polymul(a,b);
    for i in &res {
        println!("{}", i);
    }
}

// テスト
#[test]
fn test_polymul_naive() {
    let mut a = [0; N];
    let mut b = [0; N];
    let mut c = [0; N];
    a[0] = -1; a[1] = 3;
    b[0] = 1; b[1] = -1;
    c[0] = -1; c[1] = 4; c[2] = -3;
    // (3X-1)⋅(-X+1) = -3X^2+4X-1
    assert_eq!(polymul_naive(a, b), c);
}

#[test]
fn test_polymul_ntt() {
    let mut a = [0; N];
    let mut b = [0; N];
    let mut c = [0; N];
    a[0] = -1; a[1] = 3;
    b[0] = 1; b[1] = -1;
    c[0] = -1; c[1] = 4; c[2] = -3;
    assert_eq!(polymul(a, b), c);
}

#[test]
fn test_polymul() {
    let mut rng = ChaChaRng::seed_from_u64(0);
    for _ in 0..100 {
        let mut a = [0; N];
        let mut b = [0; N];
        for i in 0..N {
            a[i] = rng.gen_range(-1000000..1000000);
            b[i] = rng.gen_range(-1000000..1000000);
        }
        assert_eq!(polymul(a, b), polymul_naive(a, b));
    }
}

#[test]
fn test_polymul_mod() {
    let mut rng = ChaChaRng::seed_from_u64(0);
    for _ in 0..100 {
        let mut a = [0; N];
        let mut b = [0; N];
        for i in 0..N {
            a[i] = rng.next_u32() as i32;
            b[i] = rng.next_u32() as i32;
        }
        assert_eq!(polymul(a, b), polymul_naive_mod(a, b));
    }
}

// Mod Intのテスト
#[test]
fn test_mint_add() {
    let mut rng = ChaChaRng::seed_from_u64(0);
    for _ in 0..1000000 {
        let a = MOD - rng.next_u64() as u128;
        let b = rng.next_u64() as u128;
        let c: u128 = (Mint::new(a) + Mint::new(b)).into();
        assert_eq!(c, (a + b) % MOD);
    }
}

#[test]
fn test_mint_mul() {
    let mut rng = ChaChaRng::seed_from_u64(0);
    for _ in 0..1000000 {
        let a = MOD - rng.next_u64() as u128;
        let b = rng.next_u64() as u128;
        let c: u128 = (Mint::new(a) * Mint::new(b)).into();
        assert_eq!(c, a * b % MOD);
    }
}

#[test]
fn test_mint_inv() {
    let mut rng = ChaChaRng::seed_from_u64(0);
    for _ in 0..1000000 {
        let a = rng.next_u64() as u128;
        let c: u128 = (Mint::new(a) * Mint::new(a).inv()).into();
        assert_eq!(c, 1);
    }
}

#[test]
fn test_nth_root() {
    for i in 0..10 {
        let a = nth_root(1 << i);
        let b: u128 = a.pow(1 << i).into();
        assert_eq!(b, 1);
    }
}
