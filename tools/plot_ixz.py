#!/usr/bin/env python3
"""OpenBPM の出力 (time_series_data.h5) を可視化する後処理スクリプト。

生成物 (入力 HDF5 と同じディレクトリに出力):
  - <prefix>_ixz.png    : 伝搬マップ |E(x, y=Ny/2, z)|^2 (/field/Ixz)
  - <prefix>_final.png  : 最終電界 |E(x,y)|^2 と屈折率分布 (/field/Efinal_*, n_out_r)
  - <prefix>_prop.gif   : |E(x,y)|^2 の伝搬アニメーション (/field/frames がある場合のみ。
                          ソルバ入力で `frames = <interval>` を指定して記録する)

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
        dx = read_scalar(f, "/metadata/grid_dx")
        dz = read_scalar(f, "/metadata/grid_dz")
        dy = read_scalar(f, "/metadata/grid_dy")
        frame_interval = read_scalar(f, "/metadata/frame_interval", 0)

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

        def update(i):
            im.set_data(scale(frames[i]))
            z = (i * frame_interval * dz * 1e6) if (dz and frame_interval) else i
            ax.set_title(f"|E(x, y)|^2  z = {z:.1f} {unit if dz else 'frame'}")
            return [im]

        anim = FuncAnimation(fig, update, frames=nfr, blit=False)
        anim.save(f"{prefix}_prop.gif", writer=PillowWriter(fps=args.fps))
        plt.close(fig)
        print(f"wrote {prefix}_prop.gif ({nfr} frames)")
    else:
        print("note: /field/frames がないため GIF は生成しません"
              " (入力に `frames = <interval>` を指定すると記録されます)")


if __name__ == "__main__":
    main()
