#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This function creates summary figures for the impedance of the
realistic 3D-dendrite. They will be presented as figure 10F of the
revised manuscript.
This is related to a reviewer request from the submission at JNP.

Created October 2020
@author: Thomas Kuenzel (kuenzel{at}bio2.rwth-aachen.de)
"""

import os

if os.environ.get("DISPLAY", "") == "":
    print("no display found. Using non-interactive Agg backend")
    matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np
import semodels as sem
import sehelper as seh
import sys
import multiprocessing
from functools import partial
import scipy.signal as sig
import time
from scipy.optimize.minpack import curve_fit
from scipy.stats import chisquare


if __name__ == "__main__":
    R1 = np.load("./results/impedance_res_3d_summarydata_cell1.npy").item()
    R2 = np.load("./results/impedance_res_3d_summarydata_cell2.npy").item()
    R3 = np.load("./results/impedance_res_3d_summarydata_cell3.npy").item()

    allPL = np.concatenate(
        (
            R1["PL"],
            R1["PLp"],
            R2["PL"],
            R2["PLp"],
            R3["PL"],
            R3["PLp"],
        )
    )
    allZM = np.concatenate(
        (
            R1["ZM"],
            R1["ZMp"],
            R2["ZM"],
            R2["ZMp"],
            R3["ZM"],
            R3["ZMp"],
        )
    )
    allFR = np.concatenate(
        (
            R1["FR"],
            R1["FRp"],
            R2["FR"],
            R2["FRp"],
            R3["FR"],
            R3["FRp"],
        )
    )
    allZR = np.concatenate(
        (
            R1["ZR"],
            R1["ZRp"],
            R2["ZR"],
            R2["ZRp"],
            R3["ZR"],
            R3["ZRp"],
        )
    )
    decfun = lambda x, y0, tau, off: y0 * np.exp(-x / tau) + off
    p0 = [15.0, 50.0, 2.0]
    pzm, covzm = curve_fit(decfun, allPL, allZM, p0=p0)
    pzr, covzr = curve_fit(decfun, allPL, allZR, p0=p0)
    fitx = np.linspace(np.min(allPL), np.max(allPL), 500)
    lfp = np.polyfit(allPL, allFR, 1)
    lfy = lfp[0] * fitx + lfp[1]

    ms = 3
    fs = 7
    fh = plt.figure(figsize=(4, 7))
    ax1 = fh.add_subplot(311)
    ax1.plot(R1["PL"], R1["ZM"], "o", color="k", markersize=ms)
    ax1.plot(R1["PLp"], R1["ZMp"], "o", color=(0.7, 0.7, 0.7), markersize=ms)
    ax1.plot(R2["PL"], R2["ZM"], "x", color="k", markersize=ms)
    ax1.plot(R2["PLp"], R2["ZMp"], "x", color=(0.7, 0.7, 0.7), markersize=ms)
    ax1.plot(R3["PL"], R3["ZM"], "s", color="k", markersize=ms)
    ax1.plot(R3["PLp"], R3["ZMp"], "s", color=(0.7, 0.7, 0.7), markersize=ms)
    ax1.set_xlabel("path length ($\mu m$)")
    ax1.set_ylabel("Zmax ($M\Omega$)")
    y0, tau, off = pzm
    resfun = lambda x: y0 * np.exp(-x / tau) + off
    ax1.plot(fitx, resfun(fitx), "r-")
    res = allZM - decfun(allPL, *pzm)
    ss_r = np.sum(res ** 2)
    ss_t = np.sum((allZM - np.mean(allZM)) ** 2)
    r_s = 1 - (ss_r / ss_t)
    annotxt1 = (
        r"$y = "
        + str(np.round(y0, 2))
        + r"* e^{\frac{-x}{"
        + str(np.round(tau, 2))
        + r"}}+$"
        + str(np.round(off, 2))
    )
    annotxt2 = r"$r^{2}=$" + str(np.round(r_s, 2))
    ax1.annotate(annotxt1, (175, 13), fontsize=fs)
    ax1.annotate(annotxt2, (175, 12), fontsize=fs)
    chi2, p = chisquare(allZM, f_exp=decfun(allPL, *pzm), ddof=3)
    annotxt3 = r"$\chi^2$=" + str(round(chi2, 4)) + ", p=" + str(round(p, 4))
    ax1.annotate(annotxt3, (175, 11), fontsize=fs)

    ax2 = fh.add_subplot(312)
    ax2.plot(R1["PL"], R1["FR"], "o", color="k", markersize=ms)
    ax2.plot(R1["PLp"], R1["FRp"], "o", color=(0.7, 0.7, 0.7), markersize=ms)
    ax2.plot(R2["PL"], R2["FR"], "x", color="k", markersize=ms)
    ax2.plot(R2["PLp"], R2["FRp"], "x", color=(0.7, 0.7, 0.7), markersize=ms)
    ax2.plot(R3["PL"], R3["FR"], "s", color="k", markersize=ms)
    ax2.plot(R3["PLp"], R3["FRp"], "s", color=(0.7, 0.7, 0.7), markersize=ms)
    ax2.set_xlabel("path length ($\mu m$)")
    ax2.set_ylabel("$f_{res} (Hz)$")
    ax2.plot(fitx, lfy, "r-")
    res = allFR - (lfp[0] * allPL + lfp[1])
    ss_r = np.sum(res ** 2)
    ss_t = np.sum((allFR - np.mean(allFR)) ** 2)
    r_s = 1 - (ss_r / ss_t)
    annotxt1 = (
        r"$y = " + str(np.round(lfp[0], 2)) + r"*x+" + str(np.round(lfp[1], 2)) + r"$"
    )
    annotxt2 = r"$r^{2}=$" + str(np.round(r_s, 2))
    ax2.annotate(annotxt1, (175, 90), fontsize=fs)
    ax2.annotate(annotxt2, (175, 80), fontsize=fs)
    chi2, p = chisquare(allFR, f_exp=(lfp[0] * allPL + lfp[1]), ddof=2)
    annotxt3 = r"$\chi^2$=" + str(round(chi2, 4)) + ", p=" + str(round(p, 4))
    ax2.annotate(annotxt3, (175, 70), fontsize=fs)

    #    anno1 = "y=" + str(round(fy[0], 3)) + "*x+" + str(round(fy[1], 3))
    #    anno2 = "y=" + str(round(fya[0], 3)) + "*x+" + str(round(fya[1], 3))
    #    plt.annotate(anno1, xy=(250, 16), color="r")

    ax3 = fh.add_subplot(313)
    ax3.plot(R1["PL"], R1["ZR"], "o", color="k", markersize=ms)
    ax3.plot(R1["PLp"], R1["ZRp"], "o", color=(0.7, 0.7, 0.7), markersize=ms)
    ax3.plot(R2["PL"], R2["ZR"], "x", color="k", markersize=ms)
    ax3.plot(R2["PLp"], R2["ZRp"], "x", color=(0.7, 0.7, 0.7), markersize=ms)
    ax3.plot(R3["PL"], R3["ZR"], "s", color="k", markersize=ms)
    ax3.plot(R3["PLp"], R3["ZRp"], "s", color=(0.7, 0.7, 0.7), markersize=ms)
    ax3.set_xlabel("path length ($\mu m$)")
    ax3.set_ylabel("$Z_{max} - Z_{0} (M\Omega)$")
    y0, tau, off = pzr
    resfun = lambda x: y0 * np.exp(-x / tau) + off
    ax3.plot(fitx, resfun(fitx), "r-")
    res = allZR - decfun(allPL, *pzr)
    ss_r = np.sum(res ** 2)
    ss_t = np.sum((allZR - np.mean(allZR)) ** 2)
    r_s = 1 - (ss_r / ss_t)
    annotxt1 = (
        r"$y = "
        + str(np.round(y0, 2))
        + r"* e^{\frac{-x}{"
        + str(np.round(tau, 2))
        + r"}}+$"
        + str(np.round(off, 2))
    )
    annotxt2 = r"$r^{2}=$" + str(np.round(r_s, 2))
    ax3.annotate(annotxt1, (175, 10), fontsize=fs)
    ax3.annotate(annotxt2, (175, 9), fontsize=fs)
    chi2, p = chisquare(allZR, f_exp=decfun(allPL, *pzr), ddof=3)
    annotxt3 = r"$\chi^2$=" + str(round(chi2, 4)) + ", p=" + str(round(p, 4))
    ax3.annotate(annotxt3, (175, 8), fontsize=fs)
    plt.tight_layout()
    with PdfPages("./figs/imp3d_fres_zmax_vs_path_summary.pdf") as pp:
        pp.savefig(fh)
        plt.close(fh)
