# fuller_md_fortran — フラーレン結晶 NPT 分子動力学 (Fortran 95版)

C60フラーレン結晶の NPT 分子動力学シミュレーションコード（Fortran 95版）。
[C++版 fuller_md](https://github.com/focusnishikawa/fuller_md) からの移植。

## ディレクトリ構造

```
fuller_md_fortran/
├── README.md              ← このファイル
├── FullereneLib/          ← フラーレン分子座標データ (.cc1)
│   ├── C60-76/            ←   C60(Ih), C70(D5h), C72(D6d), C74(D3h), C76(D2,Td)
│   └── C84/               ←   C84 No.01〜No.24 (24異性体)
├── src/                   ← ソースコード・スクリプト
│   ├── Build_fuller.sh    ←   ビルドスクリプト
│   └── fuller_LJ_npt_md_core_serial.f90   [1] LJ剛体コア版 (Serial)
└── bin/                   ← 実行モジュール出力先 (ビルド時に自動作成)
```

## ソースファイル一覧

### コア版 [1] — 学習・ベンチマーク用

パラメータはソースコード内に固定（T=300K, P=0GPa, dt=1fs, 1000ステップ）。
引数は `nc`（セルサイズ）のみ。

| # | ファイル | 説明 |
|---|---------|------|
| 1 | `fuller_LJ_npt_md_core_serial.f90` | LJ剛体 Serial版。C++版からの完全移植 |

## ビルド方法

### 前提条件

- **gfortran** (GCC Fortran コンパイラ)
- macOSの場合: `brew install gcc`

### ビルド

```bash
cd fuller_md_fortran
src/Build_fuller.sh
```

環境変数 `FC` でコンパイラを明示指定可能:

```bash
FC=/opt/gcc/bin/gfortran src/Build_fuller.sh
```

### 生成される実行モジュール

| 実行モジュール | ソース | モード |
|---------------|--------|--------|
| `fuller_LJ_core_serial` | [1] | Serial |

## 実行方法

```bash
cd fuller_md_fortran

# デフォルト (3x3x3, N=108分子, 1000ステップ)
bin/fuller_LJ_core_serial

# セルサイズ指定 (5x5x5, N=500分子)
bin/fuller_LJ_core_serial 5
bin/fuller_LJ_core_serial --cell=5
```

## 物理モデル

- **アンサンブル**: NPT (等温等圧)
- **熱浴**: Nose-Hoover chain
- **圧力制御**: Parrinello-Rahman
- **時間積分**: Velocity-Verlet
- **周期境界条件**: 3次元 (triclinic cell)
- **剛体回転**: 四元数表現
- **分子間力**: Lennard-Jones (C-C, sigma=3.4A, epsilon=2.63meV)
- **近接リスト**: 対称フルリスト

## 単位系

| 物理量 | 単位 |
|--------|------|
| 距離 | A (オングストローム) |
| 質量 | amu (原子質量単位) |
| エネルギー | eV (電子ボルト) |
| 時間 | fs (フェムト秒) |
| 温度 | K (ケルビン) |
| 圧力 | GPa (ギガパスカル) |

## C++版からの移植ノート

- 0始まり配列 → 1始まり配列に全て変換
- `std::mt19937` → Fortran `RANDOM_NUMBER` + Box-Muller法
- `std::chrono` → `CPU_TIME`
- `std::clamp` → `clamp_val` 関数
- flat 1D配列レイアウトを維持（将来のOpenACC対応に備えて）
- MODULE構造で全サブルーチンを格納

## 動作環境

- macOS (Homebrew GCC / gfortran)
- Linux (gfortran)

## ライセンス

本プロジェクトは [BSD 3-Clause License](LICENSE) の下で公開されています。

## 他言語・他バージョン

| リポジトリ | 言語 | 説明 |
|-----------|------|------|
| [fuller_md](https://github.com/focusnishikawa/fuller_md) | C++ (日本語) | オリジナル |
| [fuller_md_en](https://github.com/focusnishikawa/fuller_md_en) | C++ (English) | C++版 英語 |
| [fuller_md_Julia](https://github.com/focusnishikawa/fuller_md_Julia) | Julia (English) | Julia版 英語 |
| [fuller_md_Julia_ja](https://github.com/focusnishikawa/fuller_md_Julia_ja) | Julia (日本語) | Julia版 日本語 |
| [fuller_md_fortran](https://github.com/focusnishikawa/fuller_md_fortran) | Fortran 95 (日本語) | このリポジトリ |
| [fuller_md_fortran_en](https://github.com/focusnishikawa/fuller_md_fortran_en) | Fortran 95 (English) | Fortran版 英語 |
