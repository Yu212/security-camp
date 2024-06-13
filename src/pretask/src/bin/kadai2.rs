#![allow(dead_code, unused_imports)]

use std::io::BufRead;

const L: usize = 4;

fn decomposition_naive(n: u16) -> [i8; L] {
    let a4 = n as i16;
    // 残りの3桁で表せる範囲 [-0x888, 0x777] に収まるように最上位桁を選ぶ これは一意に決まる
    let a3 = (-8..8).map(|a| (a, a4 - a * 0x1000)).find(|&(_, b)| -0x888 <= b && b <= 0x777).unwrap();
    // 次の桁も同様
    let a2 = (-8..8).map(|a| (a, a3.1 - a * 0x100)).find(|&(_, b)| -0x88 <= b && b <= 0x77).unwrap();
    let a1 = (-8..8).map(|a| (a, a2.1 - a * 0x10)).find(|&(_, b)| -0x8 <= b && b <= 0x7).unwrap();
    [a1.1 as i8, a1.0 as i8, a2.0 as i8, a3.0 as i8]
}

fn decomposition(n: u16) -> [i8; L] {
    // あらかじめ 0x8888 を足すと16進数各桁が変換後の各桁に対応する
    let n2 = n as i32 + 0x8888;
    // [0, 16) -> [-8, 8) の変換
    [n2, n2 >> 4, n2 >> 8, n2 >> 12].map(|x| ((x & 15) - 8) as i8)
}

// decomposition の逆操作
fn decomposition_inv(a: [i8; L]) -> u16 {
    let a32: [i16; L] = a.map(i16::from);
    (a32[0] + a32[1] * 0x10 + a32[2] * 0x100 + a32[3] * 0x1000) as u16
}

fn main() {
    // https://programming-idioms.org/idiom/120/read-integer-from-stdin/1906/rust
    let a: u16 = std::io::stdin()
        .lock()
        .lines()
        .next()
        .expect("stdin should be available")
        .expect("couldn't read from stdin")
        .trim()
        .parse::<u16>()
        .expect("input was not an integer");
    let res: [i8; L] = decomposition(a);
    for i in &res {
        println!("{}", i);
    }
}

// すべての入力値で正しく変換し元に戻ることを確かめる
#[test]
fn test_decomposition_naive() {
    for i in 0..u16::MAX {
        assert_eq!(i, decomposition_inv(decomposition_naive(i)));
    }
}

#[test]
fn test_decomposition() {
    for i in 0..u16::MAX {
        assert_eq!(i, decomposition_inv(decomposition(i)));
    }
}
