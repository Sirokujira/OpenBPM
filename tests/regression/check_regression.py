#!/usr/bin/env python3
"""OpenBPM 数値回帰チェック (numerical regression check).

OpenBPM の HDF5 結果ファイル (``time_series_data.h5``) を読み込み、
最終電界から安定したスカラー指標を計算して、コミット済みの golden 値
(``reference.json``) と許容誤差の範囲内で比較する。

一部のサンプルでは、ReadMe.md に記載された解析理論値との物理的な
整合性も追加で検証する。

使い方
------
golden 値との比較 (CI での用途)::

    check_regression.py --name freespace --h5 time_series_data.h5

サンプルの golden 値を登録・更新する (known-good ビルドでローカル実行後に commit)::

    check_regression.py --name freespace --h5 time_series_data.h5 --update

数値の再現性に関する注意:
    OpenMP の並列化により浮動小数点演算の加算順序がスレッド数に依存して変わるため、
    golden 値の生成・検証ともに ``OMP_NUM_THREADS=1`` で solver を実行すること。

比較指標:
    電力 (power), ピーク強度 (peak), 重心 (cx, cy), ビーム幅 (wx, wy)
    各指標の定義は ``metrics()`` 関数を参照。

物理理論チェック (curated subset):
    - freespace     : ガウシアンビーム回折 w(z) = w0·√(1+(z/zR)²)
    - tilt_wideangle: 広角 BPM で重心変位が z·tan(20°) に収束すること
                      (近軸 z·sin(20°) との区別)
"""
import argparse
import json
import math
import os
import sys

import h5py
import numpy as np

# reference.json のデフォルトパス (このスクリプトと同じディレクトリ)
REF_DEFAULT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "reference.json")

# 比較する指標の名前リスト (reference.json の "metrics" キーと対応)
METRIC_FIELDS = ("power", "peak", "cx", "cy", "wx", "wy")

# 許容誤差のデフォルト値
# rtol  : 相対許容誤差 (1%)
# atol_m: メートル単位の指標 (重心・ビーム幅) に対する絶対下限 (0.01 um)
#         これにより、ゼロ近傍の重心に対して相対誤差が発散するのを防ぐ
DEFAULT_RTOL = 0.01
DEFAULT_ATOL_M = 1.0e-8

# 絶対誤差下限を適用するメートル単位の指標
METRES_METRICS = ("cx", "cy", "wx", "wy")


def load(h5path):
    """HDF5 ファイルを読み込んで強度・メタデータ・Ixz を返す。

    Returns
    -------
    intensity : ndarray, shape (Ny, Nx)
        最終電界の強度 |E|² = Er² + Ei²。Efinal は [iy][ix] 順で格納されている。
    md : dict[str, ndarray]
        metadata グループの全キーを 1 次元配列に変換した辞書。
        空間座標は Xc (長さ Nx), Yc (長さ Ny), Zn (長さ Nz+1)。
    ixz : ndarray or None
        伝搬マップ |E(x, y=Ny/2, z)|² (shape: Nz × Nx)。存在しない場合は None。
    """
    with h5py.File(h5path, "r") as f:
        er = f["field/Efinal_r"][()].astype(np.float64)
        ei = f["field/Efinal_i"][()].astype(np.float64)
        intensity = er * er + ei * ei
        md = {k: np.array(f["metadata"][k]).ravel() for k in f["metadata"]}
        ixz = f["field/Ixz"][()].astype(np.float64) if "field/Ixz" in f else None
    return intensity, md, ixz


def metrics(intensity, md):
    """最終電界の強度から安定したスカラー指標を計算する。

    指標の定義:
        power : 総出力パワー = Σ I(x,y)
        peak  : ピーク強度 = max I(x,y)
        cx, cy: 強度重心 [m] = Σ(x·I) / Σ I  (ビームの横ずれを検出)
        wx, wy: 1/e² ビーム径 [m] = 2·√(2次モーメント)
                  ガウシアンビームの場合は 1/e² 半径と一致する

    Parameters
    ----------
    intensity : ndarray, shape (Ny, Nx)
        load() が返す強度配列。
    md : dict
        load() が返すメタデータ辞書。Xc (len Nx), Yc (len Ny) を使用。

    Notes
    -----
    Efinal は [iy][ix] 順で格納されるため、Xc を axis=1 (列方向) に、
    Yc を axis=0 (行方向) にブロードキャストして計算する。
    """
    xc = md["Xc"]  # shape (Nx,) : x 方向セル中心座標 [m]
    yc = md["Yc"]  # shape (Ny,) : y 方向セル中心座標 [m]
    total = intensity.sum()
    cx = float((intensity * xc[None, :]).sum() / total)
    cy = float((intensity * yc[:, None]).sum() / total)
    # 2次モーメントの平方根 → 1/e² 半径 (= 2·σ)
    wx = 2.0 * math.sqrt(float((intensity * (xc[None, :] - cx) ** 2).sum() / total))
    wy = 2.0 * math.sqrt(float((intensity * (yc[:, None] - cy) ** 2).sum() / total))
    return {
        "power": float(total),
        "peak": float(intensity.max()),
        "cx": cx,
        "cy": cy,
        "wx": wx,
        "wy": wy,
    }


# ============================================================
# 物理理論チェック (curated subset)
#
# 各関数は (label, measured, theory, rtol) のタプルリストを返す。
# measured と theory は同じ単位の浮動小数点数。
# rtol は許容相対誤差 (0.0 ~ 1.0)。
#
# これらは golden 値との一致ではなく、ReadMe.md に記載された
# 解析解と solver の結果が整合することを検証する。
# ============================================================

def _physics_freespace(intensity, md, m):
    """ガウシアンビーム回折の解析解との比較 (freespace サンプル用).

    理論式: w(z) = w0 · √(1 + (z/zR)²)
            zR = π·w0²·n0 / λ  (レイリー長)

    ReadMe.md の記述:
        freespace.ofd : 均一媒質中のガウシアンビーム回折。
        解析解 w(z)=w0·sqrt(1+(z/zR)^2) と比較可能。
    """
    w0 = md["beam_w0"][0]   # 入射ビームウェスト [m]
    lam = md["lambda"][0]   # 波長 [m]
    n0 = md["n_0"][0]       # 参照屈折率
    length = md["Zn"][-1] - md["Zn"][0]  # 伝搬長 z=L [m]
    z_r = math.pi * w0 ** 2 * n0 / lam  # レイリー長 [m]
    w_theory = w0 * math.sqrt(1.0 + (length / z_r) ** 2)
    return [("output beam width w(L) [um]", m["wx"] * 1e6, w_theory * 1e6, 0.03)]


def _physics_tilt(intensity, md, m):
    """広角 BPM のビームチルト変位の解析値との比較 (tilt_wideangle サンプル用).

    広角 BPM (Pade(1,1)) では傾き θ のビームが横方向に z·tan(θ) だけ変位する。
    近軸近似では z·sin(θ) に留まるため、20° では約 6% の誤差が生じる。

    ReadMe.md の記述:
        tilt_wideangle.ofd : 広角 BPM (Pade(1,1))。
        20 度チルトビームの横変位が厳密値 z·tan(20)=21.8um へ収束
        (近軸は z·sin(20)=20.5um で頭打ち)。
    """
    tilt_deg = 20.0
    length = md["Zn"][-1] - md["Zn"][0]  # 伝搬長 [m]
    x0 = md["beam_x0"][0]                # 入射ビーム中心 x 座標 [m]
    disp_theory = length * math.tan(math.radians(tilt_deg))  # 厳密値 [m]
    disp_meas = m["cx"] - x0  # 重心変位 = 出力重心 - 入射中心 [m]
    return [("centroid shift = z*tan(20) [um]", disp_meas * 1e6, disp_theory * 1e6, 0.05)]


# 物理チェック関数の登録テーブル: サンプル名 → チェック関数
PHYSICS_CHECKS = {
    "freespace": _physics_freespace,
    "tilt_wideangle": _physics_tilt,
}


def compare(name, m, ref):
    """指標 m を reference.json の golden 値と比較する。

    比較条件: |got - want| ≤ atol + rtol·|want|
        atol はメートル単位の指標に対してのみ適用する絶対下限。

    Returns
    -------
    ok : bool
        全指標が許容誤差内であれば True。
    lines : list[str]
        各指標の OK/FAIL と数値を記したログ行。
    """
    tol = ref.get("_meta", {}).get("tol", {})
    rtol = tol.get("rtol", DEFAULT_RTOL)
    atol_m = tol.get("atol_m", DEFAULT_ATOL_M)
    golden = ref["samples"].get(name, {}).get("metrics")
    lines = []
    ok = True
    if golden is None:
        return False, [f"  no golden reference for sample '{name}'"]
    for key in METRIC_FIELDS:
        got = m[key]
        want = golden[key]
        # メートル単位の指標には絶対誤差下限を加算して 0 近傍での誤判定を防ぐ
        atol = atol_m if key in METRES_METRICS else 0.0
        diff = abs(got - want)
        limit = atol + rtol * abs(want)
        passed = diff <= limit
        ok = ok and passed
        lines.append(
            f"  [{'OK' if passed else 'FAIL'}] {key:6s} got={got:+.6e} ref={want:+.6e} "
            f"diff={diff:.2e} limit={limit:.2e}"
        )
    return ok, lines


def run_physics(name, intensity, md, m):
    """解析理論値との物理チェックを実行する (登録済みサンプルのみ)。

    Returns
    -------
    ok : bool
        全チェックが許容誤差内であれば True。登録なしサンプルは常に True。
    lines : list[str]
        チェック結果のログ行。
    """
    check = PHYSICS_CHECKS.get(name)
    if check is None:
        return True, []  # 対象外サンプルはスキップ
    lines = ["  physics (vs analytic theory):"]
    ok = True
    for label, meas, theory, rtol in check(intensity, md, m):
        rel = abs(meas - theory) / abs(theory) if theory else 0.0
        passed = rel <= rtol
        ok = ok and passed
        lines.append(
            f"    [{'OK' if passed else 'FAIL'}] {label}: meas={meas:.4f} "
            f"theory={theory:.4f} rel={rel:.2%} tol={rtol:.0%}"
        )
    return ok, lines


def main():
    ap = argparse.ArgumentParser(
        description="OpenBPM 数値回帰チェック",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "例:\n"
            "  OMP_NUM_THREADS=1 ./bin/obpm data/freespace.ofd\n"
            "  python check_regression.py --name freespace --h5 time_series_data.h5\n"
            "\n"
            "golden 値を更新する場合:\n"
            "  python check_regression.py --name freespace --h5 time_series_data.h5 --update\n"
        ),
    )
    ap.add_argument("--name", required=True, help="サンプル名 (.ofd の拡張子なし basename)")
    ap.add_argument("--h5",   required=True, help="time_series_data.h5 のパス")
    ap.add_argument("--ref",  default=REF_DEFAULT, help="reference.json のパス (省略時は同ディレクトリ)")
    ap.add_argument("--update", action="store_true",
                    help="現在の指標を reference.json に書き込む (golden 値の更新)")
    args = ap.parse_args()

    intensity, md, _ = load(args.h5)
    m = metrics(intensity, md)

    if args.update:
        # --- golden 値の更新モード ---
        if os.path.exists(args.ref):
            with open(args.ref) as fh:
                ref = json.load(fh)
        else:
            # reference.json が存在しない場合は新規作成する
            ref = {
                "_meta": {
                    "omp_threads": 1,  # OMP_NUM_THREADS=1 で生成した値であることを記録
                    "tol": {"rtol": DEFAULT_RTOL, "atol_m": DEFAULT_ATOL_M},
                },
                "samples": {},
            }
        ref.setdefault("samples", {})[args.name] = {"metrics": m}
        with open(args.ref, "w") as fh:
            json.dump(ref, fh, indent=2, sort_keys=True)
            fh.write("\n")
        print(f"updated reference for '{args.name}' in {args.ref}")
        return 0

    # --- 照合モード ---
    with open(args.ref) as fh:
        ref = json.load(fh)

    print(f"== {args.name} ==")
    ok_metrics, lines = compare(args.name, m, ref)
    for line in lines:
        print(line)
    ok_physics, plines = run_physics(args.name, intensity, md, m)
    for line in plines:
        print(line)

    if ok_metrics and ok_physics:
        print("  PASS")
        return 0

    print(f"  REGRESSION DETECTED for '{args.name}'", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
