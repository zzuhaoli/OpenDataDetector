#!/usr/bin/env python3

import uproot
import matplotlib.pyplot as plt
import numpy as np
import os
import hist
import mplhep
import argparse
from pathlib import Path
import json

plt.rcParams["ytick.right"] = plt.rcParams["xtick.top"] = True
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["font.size"] = 12.0
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["legend.frameon"] = False
plt.rcParams["legend.columnspacing"] = 0.2
plt.rcParams["legend.handletextpad"] = 0.2
plt.rcParams["legend.labelspacing"] = 0.2
plt.rcParams["legend.borderpad"] = 0
plt.rcParams["legend.handlelength"] = 1.0


parser = argparse.ArgumentParser(description="Make material composition plots")
parser.add_argument("input", help="Input root file with histograms")
parser.add_argument("config", type=Path, help="Input config file with colours")
parser.add_argument("output", type=Path, help="Output directory for plots")
args = parser.parse_args()

args.output.mkdir(parents=True, exist_ok=True)

rf = uproot.open(args.input)

tracker = {
    "beampipe": "Beam pipe",
    "pixel": "Pixel",
    "sstrips": "Short Strips",
    "lstrips": "Long Strips",
    "solenoid": "Solenoid",
}

calo = {
    "ecalbarrel": "ECal barrel",
    "ecalendcap": "ECal endcap",
    "hcalbarrel": "HCal barrel",
    "hcalendcap": "HCal endcap",
}

with args.config.open() as f:
    config = json.load(f)
configcolours = {det: config[det][0]["colour"] for det in config.keys() }

for group, names in [
    ("tracker_", {**tracker}),
    ("calo_", {**calo}),
    ("", {**tracker, **calo}),
]:
    for q in (
        "x0", 
        "l0",
    ):
        for qq in ("phi", "eta"):
            h = None
            hists = {}
            labels = []
            detector = None
            for k in rf:
                name, _ = k.split(";", 1)
                if not name.endswith("all"): continue
                if not q in name: continue
                if not qq in name: continue
                if name.startswith("detector"):
                    detector = values.copy()
                    continue

                o = rf[k].to_hist()
                o.axes[0].label = qq
                l, _ = k.split("_", 1)
                if not l in names:
                    continue
                hists[l] = o
                l = names[l]

            labels = [v for k, v in names.items() if k in hists]
            colours  = [configcolours[k] for k in names.keys() if k in hists]
            hists = [hists[k] for k in names.keys() if k in hists]
                
            fig, ax = plt.subplots()
            mplhep.histplot(hists, ax=ax, stack=True, histtype="fill", label=labels, color=colours)
            ymin, ymax = ax.get_ylim()
            ax.set_ylim(top=1.2*ymax)
            ax.legend(ncol=3)

            ylab = {"l0": "\lambda_i", "x0": "X_0"}[q]
            ax.set_ylabel(rf"${ylab}$")
            ax.set_xlim(hists[0].axes[0].edges[0], hists[0].axes[0].edges[-1])
            fig.savefig(args.output / f"{group}{q}_vs_{qq}.pdf")

