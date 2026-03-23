# fuller_md C++ → Fortran 95 移植 エラー報告書

**プロジェクト**: fuller_md (フラーレン結晶 NPT 分子動力学シミュレーション)
**報告日**: 2026-03-23
**対象**: C++版 → Fortran 95版 移植における全エラーの記録と対処

---

## 1. エラー一覧

| # | 分類 | ファイル | 重要度 | 修正コミット |
|---|------|---------|--------|------------|
| 1 | コンパイルエラー | MMMD版 (f90) | 高 | `9720641` |
| 2 | 移植バグ | AIREBO版 (f90) | 高 | `123262d` |
| 3 | 移植バグ | MMMD版 (f90) | 高 | `ba29af0` |
| 4 | 移植漏れ | LJフル版 (f90) | 中 | `ba29af0` |
| 5 | 名前衝突 | AIREBO版 (f90) | 中 | `123262d` |
| 6 | C++版は正常 | MMMD版 (cpp) | — | 修正不要 |

---

## 2. 各エラーの詳細

### エラー#1: OpenMP ATOMIC 構文違反（MMMD版）

**ファイル**: `fuller_LJ_npt_mmmd_serial_omp_acc.f90` 655行目

**症状**: gfortran-15 -fopenmp でコンパイルエラー
```
Error: !$OMP ATOMIC var = var op expr not mathematically equivalent
to var = var op (expr) at (1)
```

**原因**: Fortranの`!$OMP ATOMIC`は`var = var op expr`の形式で**1つの演算子**のみ許可。
```fortran
! NG: 2つの減算演算子
!$OMP ATOMIC
Fv((j-1)*3+c) = Fv((j-1)*3+c) - fi_c - fk_c
```

**修正**:
```fortran
! OK: 括弧で1つの演算に統合
!$OMP ATOMIC
Fv((j-1)*3+c) = Fv((j-1)*3+c) - (fi_c + fk_c)
```

**C++版との比較**: C++版は`F[j*3+c]-=fi+fk;`で、`-=`は1つの複合代入演算子のため問題なし。Fortranには`-=`演算子がないため展開時に2演算になった移植固有のエラー。

**影響範囲**: OpenMPビルド時のみ発生。Serialビルドでは`!$OMP`行はコメント扱いのためエラーにならない。

---

### エラー#2: AIREBO版 3体力項の欠落

**ファイル**: `fuller_airebo_npt_md_serial_omp_acc.f90`（初版）

**症状**: 計算結果は出るが、REBO-IIの3体力（db_ij/dr_k, db_ji/dr_l）が未実装のため物理的に不正確な結果。

**原因**: 初回移植時（`ba29af0`）のAIREBO版(484行)は骨格のみの簡略版で、結合次数の角度依存3体力項が省略されていた。

**修正**: 完全版(1258行)に置き換え（`123262d`）。
- 3体力項`db_ij/dr_k`と`db_ji/dr_l`を完全実装
- ビリアルテンソル9成分全て計算（非立方晶セル対応）
- `!$OMP ATOMIC`による並列安全な力蓄積

**検証**: 修正後のREBO-IIエネルギーが物理的に妥当な値（C60で約-450 eV/mol）を出力することを確認。

---

### エラー#3: MMMD版 トポロジー構築の不完全さ

**ファイル**: `fuller_LJ_npt_mmmd_serial_omp_acc.f90`（初版）

**症状**: 初版(462行)ではBond力のみ実装で、Angle/Dihedral/Improper力場が欠落。

**原因**: 初回移植時の簡略化。C++版(1291行)の全分子力学力場を移植するには、トポロジー構築（隣接リストからAngle/Dihedral/Improperを自動生成）の完全な移植が必要だった。

**修正**: 完全版(1378行)に置き換え。
- Bond stretching (調和ポテンシャル)
- Angle bending (調和ポテンシャル)
- Dihedral torsion (Bekker解析的力)
- Improper dihedral
- LJ分子間力（分子内排除）
- .cc1ファイルからの結合トポロジー読み込み

---

### エラー#4: LJフル版 CLI引数解析・リスタートの不足

**ファイル**: `fuller_LJ_npt_md_serial_omp_acc.f90`（初版）

**症状**: 初版(617行)では基本的なCLI引数とシミュレーションは動作するが、一部のランタイムオプション（`--warmup_mon`, `--ofile`等）やリスタート復元の一部が未実装。

**修正**: 完全版(1822行)に置き換え。
- 全ランタイムオプション対応
- .cc1ファイル読み込み（C60〜C84全異性体）
- OVITO XYZ出力
- リスタート保存/復元
- Cold start + Warmup（速度リスケーリング）
- abort.md/stop.md停止制御
- C76/C84のフラーレン名解決

---

### エラー#5: Fortran大文字小文字非区別による変数名衝突（AIREBO版）

**ファイル**: `fuller_airebo_npt_md_serial_omp_acc.f90`

**症状**: gfortranコンパイルエラー
```
Error: Symbol 'st' at (1) already has basic type of REAL
```

**原因**: Fortranは大文字小文字を区別しない。C++版から移植した際に:
- `sT` (double, sum of Temperature) と
- `st` (string, crystal structure)
が同一変数として衝突。

**修正**:
```fortran
! 修正前: sT と st が衝突
double precision :: sT, sP, sa
character(len=256) :: st

! 修正後: 明確に分離
double precision :: sum_T, sum_P, sum_a
character(len=256) :: cryst_st
```

**教訓**: Fortranへの移植時は、大文字小文字のみで区別されるC++変数名を全て検出し、一意な名前に変更する必要がある。

---

### エラー#6: C++版の確認結果（問題なし）

**確認対象**: C++版全5ファイルの`#pragma omp atomic`直後の演算

**結果**: C++版には同様のエラーなし。

C++の`-=`は1つの複合代入演算子であり、`F[j]-=fi+fk`は`F[j] = F[j] - (fi+fk)`と等価。OpenMP仕様上、`x op= expr`形式は右辺の`expr`が任意の式でも合法。

| 言語 | コード | OpenMP判定 |
|------|--------|-----------|
| C++ | `F[j] -= fi + fk;` | OK（`-=`は1つの複合代入） |
| Fortran | `Fv(j) = Fv(j) - fi_c - fk_c` | NG（`-`が2回 = 2演算） |
| Fortran (修正) | `Fv(j) = Fv(j) - (fi_c + fk_c)` | OK（`-`が1回 = 1演算） |

---

## 3. 移植時の共通課題と教訓

### 3.1 配列インデックス変換 (0始まり → 1始まり)

C++の全ての配列アクセス`arr[i*3+a]`をFortranの`arr((i-1)*3+a)`に変換する必要がある。約4800行のC++コードに数百箇所のインデックス変換が存在し、1箇所でも誤ると実行時に不正な結果やセグメンテーション違反が発生する。

### 3.2 Fortranの大文字小文字非区別

C++では`sT`(温度合計)と`st`(結晶構造文字列)は異なる変数だが、Fortranでは同一。移植時に全変数名を走査して衝突を検出・解消する必要がある。

### 3.3 OpenMP ATOMIC の言語差異

C++の`-=`, `+=`は1つの複合代入演算子だが、Fortranには存在しない。`var = var - a - b`と展開すると2演算になり、`!$OMP ATOMIC`違反になる。括弧`var = var - (a + b)`で解決。

### 3.4 乱数生成器の差異

C++の`std::mt19937` + `std::normal_distribution`はFortranに直接対応がない。`RANDOM_NUMBER`(一様分布) + Box-Muller法(正規分布変換)で代替。乱数シード初期化方法も異なるため、C++版と完全に同一の乱数列は再現できない。

### 3.5 文字列処理の差異

C++の`std::string`, `std::map`, `std::istringstream`等の文字列操作はFortranでは`CHARACTER`型と`READ`文で代替する必要があり、CLI引数解析や.cc1ファイルパースの実装が大幅に異なる。

---

## 4. 最終状態

| # | ファイル | 行数 | Serial | OpenMP | 状態 |
|---|---------|------|--------|--------|------|
| 1 | `fuller_LJ_npt_md_core_serial.f90` | 826 | OK | — | 正常 |
| 2 | `fuller_LJ_npt_md_core_serial_omp_acc.f90` | 939 | OK | OK | 正常 |
| 3 | `fuller_LJ_npt_md_serial_omp_acc.f90` | 1822 | OK | OK | 正常 |
| 4 | `fuller_LJ_npt_mmmd_serial_omp_acc.f90` | 1378 | OK | OK | 正常 |
| 5 | `fuller_airebo_npt_md_serial_omp_acc.f90` | 1258 | OK | OK | 正常 |

**合計: 6,223行 (C++ 4,773行からの移植)**

全ファイルがgfortran-15でSerial/OpenMP両モードでコンパイル・動作確認済み。
