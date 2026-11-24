#!/usr/bin/env python3
"""Compose the README example image from the PNG previews written by `dgm3d.py view --mode all`.

Usage:
  tools/readme_image.py PREFIX [OUT] [--title TEXT]

  PREFIX   the --out prefix of the view run, e.g. fjh  (reads fjh_view.png, fjh_wedge.png, …)
  OUT      default docs/examples.png

Example:
  ./dgm3d.py view --lat 47.0748 --lon 12.7506 --az 268 --fov 60 --dist 12 --mode all --out fjh
  tools/readme_image.py fjh --title "Großglockner from the Kaiser-Franz-Josefs-Höhe"

Needs matplotlib (Fedora: sudo dnf install python3-matplotlib).
"""
import argparse, os, sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

PANELS = [
    ("view",    "The view from the standpoint (rendered)"),
    ("wedge",   "wedge — terrain in the field of view"),
    ("visible", "visible — orange: seen from the standpoint"),
    ("relief",  "relief — panorama plaque, near ridges stand out"),
    ("litho",   "litho — lithophane, backlit"),
]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("prefix")
    p.add_argument("out", nargs="?", default="docs/examples.png")
    p.add_argument("--title", default="")
    p.add_argument("--credit", default="Elevation data: © BEV – Bundesamt für Eich- und Vermessungswesen, "
                                       "CC BY 4.0")
    p.add_argument("--dpi", type=int, default=90)
    a = p.parse_args(argv)

    imgs = {}
    for key, _ in PANELS:
        f = f"{a.prefix}_{key}.png"
        if not os.path.exists(f):
            sys.exit(f"{f} not found — run dgm3d.py view ... --mode all --out {a.prefix} first")
        imgs[key] = mpimg.imread(f)

    fig = plt.figure(figsize=(16, 13))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.1, 1.3, 0.9], hspace=0.18, wspace=0.05)
    slots = {"view": gs[0, :], "wedge": gs[1, 0], "visible": gs[1, 1],
             "relief": gs[2, 0], "litho": gs[2, 1]}
    for key, label in PANELS:
        ax = fig.add_subplot(slots[key])
        img = imgs[key]
        ax.imshow(img, cmap="gray" if img.ndim == 2 else None, interpolation="lanczos")
        ax.set_title(label, fontsize=13, loc="left")
        ax.axis("off")
    if a.title:
        fig.suptitle(a.title, fontsize=17, y=0.985)
    fig.text(0.5, 0.008, a.credit, ha="center", fontsize=10, color="#666")
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=a.dpi, bbox_inches="tight", facecolor="white")
    print(f"→ {a.out}")


if __name__ == "__main__":
    main()
