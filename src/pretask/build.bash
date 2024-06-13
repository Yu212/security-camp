# 同一ディレクトリに配布されたフォルダを置いて実行

# テストの実行
cargo test

# コンパイルして配布されたフォルダ内に配置
cargo build --package security-camp-pretask --bin kadai1 --out-dir kadai2/kadai2-1 -Z unstable-options
cargo build --package security-camp-pretask --bin kadai2 --out-dir kadai2/kadai2-2 -Z unstable-options
cargo build --package security-camp-pretask --bin kadai3 --out-dir kadai2/ -Z unstable-options

# 配布された testRust.bash の実行
cd ./kadai2/kadai2-1/
mv kadai1 kadai2-1
./testRust.bash
cd ../kadai2-2/
mv kadai2 kadai2-2
./testRust.bash
