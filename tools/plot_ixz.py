#!/usr/bin/env python3
"""OpenBPM の出力 (time_series_data.h5) を可視化する後処理スクリプト。

生成物 (入力 HDF5 と同じディレクトリに出力):
  - <prefix>_ixz.png    : 伝搬マップ |E(x, y=Ny/2, z)|^2 (/field/Ixz)
  - <prefix>_final.png  : 最終電界 |E(x,y)|^2 と屈折率分布 (/field/Efinal_*, n_out_r)
  - <prefix>_trace.png  : z ごとのスカラー推移 (/trace : パワー・ピーク強度・重心・幅)
  - <prefix>_overlap.png: 各導波モードへのパワー占有率の z 推移 (/trace/overlap。
                          `modes = <nModes>` 指定時のみ)
  - <prefix>_modes.png  : 導波モード形状と neff (/modes がある場合のみ。
                          ソルバ入力で `modes = <nModes>` を指定して解析する)
  - <prefix>_prop.gif   : |E(x,y)|^2 の伝搬アニメーション (/field/frames がある場合のみ。
                          ソルバ入力で `frames = <interval>` を指定して記録する)
  - <prefix>_phase.gif  : 位相 arg(E) の伝搬アニメーション (/field/frames_r/_i がある
                          場合のみ。`frames = <interval> complex` で記録する)

使い方:
  python3 tools/plot_ixz.py time_series_data.h5 [--prefix out] [--db] [--fps 15]

依存: h5py, numpy, matplotlib (GIF 出力は matplotlib.animation の PillowWriter)
"""
import argparse
import os
import sys

import h5py
import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    print("error: matplotlib が必要です (pip install matplotlib)", file=sys.stderr)
    sys.exit(1)


def read_scalar(f, name, default=None):
    if name in f:
        return float(np.asarray(f[name]).reshape(-1)[0])
    return default


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("h5file", help="time_series_data.h5")
    ap.add_argument("--prefix", default=None, help="出力ファイル名の接頭辞 (既定: 入力名)")
    ap.add_argument("--db", action="store_true", help="強度を dB スケールで表示")
    ap.add_argument("--fps", type=int, default=15, help="GIF のフレームレート")
    args = ap.parse_args()

    prefix = args.prefix or os.path.splitext(args.h5file)[0]

    with h5py.File(args.h5file, "r") as f:
        ixz = f["/field/Ixz"][:] if "/field/Ixz" in f else None
        er = f["/field/Efinal_r"][:]
        ei = f["/field/Efinal_i"][:]
        n_out = f["/field/n_out_r"][:] if "/field/n_out_r" in f else None
        frames = f["/field/frames"][:] if "/field/frames" in f else None
        # 各フレームの z 位置 (ソルバーが座標として書き出す)。
        # 無い場合は frame_interval と grid_dz から従来どおり計算する。
        frames_z = f["/field/frames_z"][:] if "/field/frames_z" in f else None
        frames_c = None
        if "/field/frames_r" in f and "/field/frames_i" in f:
            frames_c = f["/field/frames_r"][:] + 1j * f["/field/frames_i"][:]
        dx = read_scalar(f, "/metadata/grid_dx")
        dz = read_scalar(f, "/metadata/grid_dz")
        dy = read_scalar(f, "/metadata/grid_dy")
        frame_interval = read_scalar(f, "/metadata/frame_interval", 0)
        trace = {}
        if "/trace" in f:
            for k in ("z", "power", "peak", "centroid_x", "centroid_y", "width_x", "width_y"):
                if f"/trace/{k}" in f:
                    trace[k] = f[f"/trace/{k}"][:]
        # モード結合率 (nModes x ntr)。`modes = <n>` 指定時のみ存在する
        overlap = f["/trace/overlap"][:] if "/trace/overlap" in f else None
        mode_fields = []
        mode_neff = []
        if "/modes" in f:
            mode_neff = list(np.asarray(f["/modes/neff"]))
            for m in range(len(mode_neff)):
                mode_fields.append(f[f"/modes/mode{m + 1}"][:])

    Ifinal = er ** 2 + ei ** 2
    Ny, Nx = Ifinal.shape

    def scale(I):
        if args.db:
            ref = I.max() if I.max() > 0 else 1.0
            return 10 * np.log10(np.maximum(I / ref, 1e-8))
        return I

    def extent_xy():
        if dx and dy:
            return [-Nx / 2 * dx * 1e6, Nx / 2 * dx * 1e6,
                    -Ny / 2 * dy * 1e6, Ny / 2 * dy * 1e6]
        return None

    unit = "um" if dx else "cell"

    # --- 伝搬マップ ---
    if ixz is not None:
        Nz = ixz.shape[0]
        fig, ax = plt.subplots(figsize=(8, 5))
        ext = ([0, Nz * dz * 1e6, -Nx / 2 * dx * 1e6, Nx / 2 * dx * 1e6]
               if (dx and dz) else None)
        im = ax.imshow(scale(ixz).T, origin="lower", aspect="auto", extent=ext,
                       cmap="inferno")
        ax.set_xlabel(f"z [{unit}]")
        ax.set_ylabel(f"x [{unit}]")
        ax.set_title("|E(x, y=Ny/2, z)|^2" + (" [dB]" if args.db else ""))
        fig.colorbar(im, ax=ax)
        fig.tight_layout()
        fig.savefig(f"{prefix}_ixz.png", dpi=150)
        plt.close(fig)
        print(f"wrote {prefix}_ixz.png")

    # --- 最終電界と屈折率分布 ---
    ncols = 2 if n_out is not None else 1
    fig, axes = plt.subplots(1, ncols, figsize=(6 * ncols, 5))
    axes = np.atleast_1d(axes)
    im = axes[0].imshow(scale(Ifinal), origin="lower", extent=extent_xy(),
                        cmap="inferno")
    axes[0].set_title("|E_final(x, y)|^2" + (" [dB]" if args.db else ""))
    axes[0].set_xlabel(f"x [{unit}]")
    axes[0].set_ylabel(f"y [{unit}]")
    fig.colorbar(im, ax=axes[0])
    if n_out is not None:
        im = axes[1].imshow(n_out, origin="lower", extent=extent_xy(), cmap="viridis")
        axes[1].set_title("n(x, y)")
        axes[1].set_xlabel(f"x [{unit}]")
        fig.colorbar(im, ax=axes[1])
    fig.tight_layout()
    fig.savefig(f"{prefix}_final.png", dpi=150)
    plt.close(fig)
    print(f"wrote {prefix}_final.png")

    # --- z ごとのスカラー推移 (/trace) ---
    if trace.get("z") is not None and len(trace["z"]) > 1:
        zum = trace["z"] * 1e6
        panels = [
            # ラベルは英語 (matplotlib の既定フォントは日本語グリフを持たない)
            ("power",      "cross-section power",  1.0,  None),
            ("peak",       "peak |E|^2",           1.0,  None),
            ("centroid_x", "centroid [um]",        1e6,  "centroid_y"),
            ("width_x",    "beam width 2sigma [um]", 1e6, "width_y"),
        ]
        avail = [p for p in panels if p[0] in trace]
        fig, axes = plt.subplots(len(avail), 1, figsize=(8, 2.6 * len(avail)), sharex=True)
        axes = np.atleast_1d(axes)
        for ax, (key, label, sc, pair) in zip(axes, avail):   # sc : scale() 関数を隠さない名前
            ax.plot(zum, trace[key] * sc, label=key)
            if pair and pair in trace:
                ax.plot(zum, trace[pair] * sc, label=pair)
                ax.legend(fontsize=8)
            ax.set_ylabel(label)
            ax.grid(alpha=0.3)
        axes[-1].set_xlabel("z [um]")
        axes[0].set_title("BPM propagation traces (/trace)")
        fig.tight_layout()
        fig.savefig(f"{prefix}_trace.png", dpi=150)
        plt.close(fig)
        print(f"wrote {prefix}_trace.png")

    # --- モード結合率の z 推移 (/trace/overlap) ---
    if overlap is not None and trace.get("z") is not None and overlap.shape[1] > 1:
        zum = trace["z"] * 1e6
        fig, ax = plt.subplots(figsize=(8, 4))
        for m in range(overlap.shape[0]):
            ne = f"  (neff = {mode_neff[m]:.6f})" if m < len(mode_neff) else ""
            ax.plot(zum, overlap[m], label=f"mode {m + 1}{ne}")
        if overlap.shape[0] > 1:
            ax.plot(zum, overlap.sum(axis=0), "k--", lw=1, label="sum")
        ax.set_xlabel("z [um]")
        ax.set_ylabel("power fraction")
        ax.set_title("Mode coupling ratio (/trace/overlap)")
        ax.set_ylim(0, 1.05)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(f"{prefix}_overlap.png", dpi=150)
        plt.close(fig)
        print(f"wrote {prefix}_overlap.png ({overlap.shape[0]} modes)")

    # --- 導波モード形状 (/modes) ---
    if mode_fields:
        nm = len(mode_fields)
        ncols_m = min(nm, 3)
        nrows_m = (nm + ncols_m - 1) // ncols_m
        fig, axes = plt.subplots(nrows_m, ncols_m,
                                 figsize=(5 * ncols_m, 4.5 * nrows_m))
        axes = np.atleast_1d(axes).ravel()
        for m, (fld, ne) in enumerate(zip(mode_fields, mode_neff)):
            # モードは符号付きの実数場なので発散配色 (0 中心) で描く
            # (べき乗法の符号は任意のため、最大振幅点が正になるよう正規化)
            if fld.flat[np.abs(fld).argmax()] < 0:
                fld = -fld
            vmax = np.abs(fld).max()
            im = axes[m].imshow(fld, origin="lower", extent=extent_xy(),
                                cmap="RdBu_r", vmin=-vmax, vmax=vmax)
            axes[m].set_title(f"mode {m + 1}  (neff = {ne:.6f})")
            axes[m].set_xlabel(f"x [{unit}]")
            axes[m].set_ylabel(f"y [{unit}]")
            fig.colorbar(im, ax=axes[m])
        for m in range(nm, len(axes)):
            axes[m].set_visible(False)
        fig.tight_layout()
        fig.savefig(f"{prefix}_modes.png", dpi=150)
        plt.close(fig)
        print(f"wrote {prefix}_modes.png ({nm} modes)")

    # --- 伝搬アニメーション ---
    if frames is not None:
        from matplotlib.animation import FuncAnimation, PillowWriter
        nfr = frames.shape[0]
        vmax = frames.max()
        fig, ax = plt.subplots(figsize=(6, 5))
        im = ax.imshow(scale(frames[0]), origin="lower", extent=extent_xy(),
                       cmap="inferno",
                       vmin=(-80 if args.db else 0),
                       vmax=(0 if args.db else vmax))
        ax.set_xlabel(f"x [{unit}]")
        ax.set_ylabel(f"y [{unit}]")
        fig.colorbar(im, ax=ax)

        def frame_z_um(i):
            if frames_z is not None:
                return frames_z[i] * 1e6, unit
            if dz and frame_interval:
                return i * frame_interval * dz * 1e6, unit
            return i, "frame"

        def update(i):
            im.set_data(scale(frames[i]))
            z, zu = frame_z_um(i)
            ax.set_title(f"|E(x, y)|^2  z = {z:.1f} {zu}")
            return [im]

        anim = FuncAnimation(fig, update, frames=nfr, blit=False)
        anim.save(f"{prefix}_prop.gif", writer=PillowWriter(fps=args.fps))
        plt.close(fig)
        print(f"wrote {prefix}_prop.gif ({nfr} frames)")
    else:
        print("note: /field/frames がないため GIF は生成しません"
              " (入力に `frames = <interval>` を指定すると記録されます)")

    # --- 位相の伝搬アニメーション (複素電界がある場合) ---
    if frames_c is not None:
        from matplotlib.animation import FuncAnimation, PillowWriter
        nfr = frames_c.shape[0]
        # 強度の低い領域の位相は無意味なので、振幅で不透明度を落として表示する
        amp = np.abs(frames_c)
        amax = amp.max() if amp.max() > 0 else 1.0
        fig, ax = plt.subplots(figsize=(6, 5))
        im = ax.imshow(np.angle(frames_c[0]), origin="lower", extent=extent_xy(),
                       cmap="twilight", vmin=-np.pi, vmax=np.pi)
        im.set_alpha(np.clip(amp[0] / amax * 4.0, 0, 1))
        ax.set_xlabel(f"x [{unit}]")
        ax.set_ylabel(f"y [{unit}]")
        cb = fig.colorbar(im, ax=ax)
        cb.set_label("arg(E) [rad]")

        def update_ph(i):
            im.set_data(np.angle(frames_c[i]))
            im.set_alpha(np.clip(amp[i] / amax * 4.0, 0, 1))
            if frames_z is not None:
                z, zu = frames_z[i] * 1e6, unit
            elif dz and frame_interval:
                z, zu = i * frame_interval * dz * 1e6, unit
            else:
                z, zu = i, "frame"
            ax.set_title(f"arg(E(x, y))  z = {z:.1f} {zu}")
            return [im]

        anim = FuncAnimation(fig, update_ph, frames=nfr, blit=False)
        anim.save(f"{prefix}_phase.gif", writer=PillowWriter(fps=args.fps))
        plt.close(fig)
        print(f"wrote {prefix}_phase.gif ({nfr} frames)")


if __name__ == "__main__":
    main()
