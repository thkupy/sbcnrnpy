#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Nov 2020
Function to analyze the additional data (per reviewer request)
for Figure8a.
These results show the stimulus dependence of the simulations in Fig. 8a.
Analysis will only be presented in text of the revised manuscript.
@author: Thomas Kuenzel (kuenzel{at}bio2.rwth-aachen.de)
"""

import os
import matplotlib
import sys

if os.environ.get("DISPLAY", "") == "":
    print("no display found. Using non-interactive Agg backend")
    matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.stats import linregress

# some global plot settings
plt.rc("font", family="serif", serif="times")
plt.rcParams["pdf.fonttype"] = 42
plt.rc("text", usetex=True)
plt.rc("xtick", labelsize="x-small")
plt.rc("ytick", labelsize="x-small")
plt.rc("axes", labelsize="small")


def makeplot(axh, xd, yd, infot):
    print("-----------------")
    print(infot)
    print(str(np.mean(yd)) + "+/-" + str(np.std(yd)) + ", N=" + str(len(yd)))
    alpha = 0.05
    x = np.linspace(np.min(xd), np.max(xd), 100)
    axh.plot(xd, yd, "ko-")
    n = np.size(xd)
    a, b, r, p, err = linregress(xd, yd)
    axh.plot(x, np.polyval((a, b), x), "r--")
    y = a * x + b
    sd = 1.0 / (n - 2.0) * np.sum((yd - a * xd - b) ** 2)
    sd = np.sqrt(sd)
    sxd = np.sum((xd - xd.mean()) ** 2)
    sx = (x - xd.mean()) ** 2
    q = scipy.stats.t.ppf(1.0 - alpha / 2, n - 2)
    dy = q * sd * np.sqrt(1.0 / n + sx / sxd)
    yl = y - dy
    yu = y + dy
    axh.fill_between(x, yl, yu, alpha=0.3, facecolor="blue", edgecolor="none")
    yl = plt.ylim()
    xl = plt.xlim()
    axh.annotate(
        "r="
        + str(np.round(r, 2))
        + ", r²="
        + str(np.round(r ** 2, 2))
        + ", p="
        + str(np.round(p, 2)),
        (xl[0], yl[1]),
    )
    axh.annotate(
        str(np.round(np.mean(yd), 2)) + " +/- " + str(np.round(np.std(yd), 2)),
        (xl[0], yl[0]),
    )
    print("y = " + str(a) + " * x + " + str(b))
    print("r = " + str(np.round(r, 2)))
    print("r²= " + str(np.round(r ** 2, 2)))
    if p < 0.05:
        print("*** p = " + str(np.round(p, 5)) + " ***")
    else:
        print("p = " + str(np.round(p, 5)))


if __name__ == "__main__":
    df = pd.read_csv("./results/Figure8a_data.csv", delimiter=" ", index_col="N")
    df1 = df[df.rlf == 1].copy()
    df1.sort_values(by=["db"], inplace=True)
    df2 = df[df.atcf == 1].copy()
    df2.sort_values(by=["cf"], inplace=True)
    df3 = df[df.at1500 == 1].copy()
    df3.sort_values(by=["freq"], inplace=True)

    fh = plt.figure()
    ax1 = fh.add_subplot(6, 3, 1)
    ax1.set_ylabel("LmaxAP (µm)")
    ax1.set_title("Levels @1500Hz CF")
    x = df1.db
    y = df1.Lmaxap
    makeplot(ax1, x, y, "db,Lmaxap")

    ax2 = fh.add_subplot(6, 3, 2)
    ax2.set_title("CFs, at CF")
    x = df2.cf
    y = df2.Lmaxap
    makeplot(ax2, x, y, "cf,Lmaxap")

    ax3 = fh.add_subplot(6, 3, 3)
    ax3.set_title("Frequencies @1500Hz CF")
    x = df3.freq
    y = df3.Lmaxap
    makeplot(ax3, x, y, "freq,Lmaxap")

    ax4 = fh.add_subplot(6, 3, 4)
    makeplot(ax4, df1.db, df1.Dmaxap, "db,Dmaxap")
    ax4.set_ylabel("DmaxAP (µm)")

    ax5 = fh.add_subplot(6, 3, 5)
    makeplot(ax5, df2.cf, df2.Dmaxap, "cf,Dmaxap")

    ax6 = fh.add_subplot(6, 3, 6)
    makeplot(ax6, df3.freq, df3.Dmaxap, "freq,Dmaxap")

    ax7 = fh.add_subplot(6, 3, 7)
    makeplot(ax7, df1.db, df1.Lmaxvs, "db,Lmaxvs")
    ax7.set_ylabel("Lmaxvs (µm)")

    ax8 = fh.add_subplot(6, 3, 8)
    makeplot(ax8, df2.cf, df2.Lmaxvs, "cf,Lmaxvs")

    ax9 = fh.add_subplot(6, 3, 9)
    makeplot(ax9, df3.freq, df3.Lmaxvs, "freq,Lmaxvs")

    ax10 = fh.add_subplot(6, 3, 10)
    makeplot(ax10, df1.db, df1.Dmaxvs, "db,Dmaxvs")
    ax10.set_ylabel("Dmaxvs (µm)")

    ax11 = fh.add_subplot(6, 3, 11)
    makeplot(ax11, df2.cf, df2.Dmaxvs, "cf,Dmaxvs")

    ax12 = fh.add_subplot(6, 3, 12)
    makeplot(ax12, df3.freq, df3.Dmaxvs, "freq,Dmaxvs")

    ax13 = fh.add_subplot(6, 3, 13)
    makeplot(ax13, df1.db, df1.dapmaxap, "db,apchange")
    ax13.set_ylabel("AP change")

    ax14 = fh.add_subplot(6, 3, 14)
    makeplot(ax14, df2.cf, df2.dapmaxap, "cf,apchange")

    ax15 = fh.add_subplot(6, 3, 15)
    makeplot(ax15, df3.freq, df3.dapmaxap, "freq,apchange")

    ax16 = fh.add_subplot(6, 3, 16)
    makeplot(ax16, df1.db, df1.dvsmaxvs, "db,vschange")
    ax16.set_ylabel("Dmaxvs (µm)")
    ax16.set_xlabel("Level (dB SPL)")

    ax17 = fh.add_subplot(6, 3, 17)
    makeplot(ax17, df2.cf, df2.dvsmaxvs, "cf,vschange")
    ax17.set_xlabel("CF (Hz)")

    ax18 = fh.add_subplot(6, 3, 18)
    makeplot(ax18, df3.freq, df3.dvsmaxvs, "freq,vschange")
    ax18.set_xlabel("Freq (Hz)")

    print("----------------------------")
    print("Grand Average...")
    print("Dmaxap:")
    print(
        str(np.mean(df.Dmaxap))
        + "+/-"
        + str(np.std(df.Dmaxap))
        + ", N="
        + str(len(df.index))
    )
    print("Lmaxap:")
    print(
        str(np.mean(df.Lmaxap))
        + "+/-"
        + str(np.std(df.Lmaxap))
        + ", N="
        + str(len(df.index))
    )
    print("Dmaxvs:")
    print(
        str(np.mean(df.Dmaxvs))
        + "+/-"
        + str(np.std(df.Dmaxvs))
        + ", N="
        + str(len(df.index))
    )
    print("Lmaxvs:")
    print(
        str(np.mean(df.Lmaxvs))
        + "+/-"
        + str(np.std(df.Lmaxvs))
        + ", N="
        + str(len(df.index))
    )
    print("apchange:")
    print(
        str(np.mean(df.dapmaxap))
        + "+/-"
        + str(np.std(df.dapmaxap))
        + ", N="
        + str(len(df.index))
    )
    print("vschange:")
    print(
        str(np.mean(df.dvsmaxvs))
        + "+/-"
        + str(np.std(df.dvsmaxvs))
        + ", N="
        + str(len(df.index))
    )

    plt.tight_layout()
    plt.show()
