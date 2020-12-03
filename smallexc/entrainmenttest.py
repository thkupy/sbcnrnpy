#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Additional simulations to show the influence of dendritic synapses
on entrainment for the stylized dendrite and one realistic 3D-dendrite.
Entrainment index is calculated as in Joris et al., 1994.
This will be presented in Fig.7 and 11 of the revised manuscript.

This was a reviewer request from the JNP submission.
"""
import sys
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from neuron import h
import semodels as sem
import sehelper as seh


def runonce(x, cf=450.0, freq=500.0, dbspl=65.0, dur=1000.0):
    thisdur = (dur - 5.0) / 1000.0
    N = 32
    anf_num = (N + 1, 0, 0)
    SO = seh.create_sound(
        fs=100e3,
        freq=freq,
        duration=thisdur,
        dbspl=dbspl,
    )
    SP = seh.create_spk(SO, fs=100e3, N=N + 1, cf=cf, seed=65564 + x, anf_num=anf_num)
    S = np.array(SP["spikes"] * 1000.0)
    R0 = sem.SE_BS(
        S,
        Simdur=dur,
        dt=0.01,
        G_EoH=0.055,
        StochasticEoH=True,
        EoHseed=917634 + x,
        N=32,
        G_Dend=0.0,
    )
    print(".")
    RD = sem.SE_BS(
        S,
        Simdur=dur,
        dt=0.01,
        G_EoH=0.055,
        StochasticEoH=True,
        EoHseed=917634 + x,
        N=32,
        G_Dend=0.016,
    )
    print(".")
    R03 = sem.SE_3D(
        S,
        Simdur=dur,
        dt=0.01,
        G_EoH=0.055,
        StochasticEoH=True,
        EoHseed=917634 + x,
        N=32,
        G_Dend=0.0,
        cell=1,
    )
    print(".")
    RD3 = sem.SE_3D(
        S,
        Simdur=dur,
        dt=0.01,
        G_EoH=0.055,
        StochasticEoH=True,
        EoHseed=917634 + x,
        N=32,
        G_Dend=0.064,
        cell=1,
    )
    print(".")
    Ev0 = seh.SimpleDetectAP(R0["Vm"], thr=-100, dt=0.01, LM=-20, RM=10)
    EvD = seh.SimpleDetectAP(RD["Vm"], thr=-100, dt=0.01, LM=-20, RM=10)
    Ev03 = seh.SimpleDetectAP(R03["Vm"], thr=-100, dt=0.01, LM=-20, RM=10)
    EvD3 = seh.SimpleDetectAP(RD3["Vm"], thr=-100, dt=0.01, LM=-20, RM=10)
    VS0, phi0, R0, phases0 = seh.vectorstrength(Ev0["PeakT"], freq, [0.0, dur])
    VSD, phiD, RD, phasesD = seh.vectorstrength(EvD["PeakT"], freq, [0.0, dur])
    VS03, phi03, R03, phases03 = seh.vectorstrength(Ev03["PeakT"], freq, [0.0, dur])
    VSD3, phiD3, RD3, phasesD3 = seh.vectorstrength(EvD3["PeakT"], freq, [0.0, dur])
    nAP0 = len(Ev0["PeakT"])
    nAPD = len(EvD["PeakT"])
    nAP03 = len(Ev03["PeakT"])
    nAPD3 = len(EvD3["PeakT"])
    ISI0 = np.diff(Ev0["PeakT"])
    ISID = np.diff(EvD["PeakT"])
    ISI03 = np.diff(Ev03["PeakT"])
    ISID3 = np.diff(EvD3["PeakT"])
    cycdurms = 1000.0 / freq
    ncycles = (thisdur * 1000.0) / cycdurms
    bins = np.arange(0.0, 25.0, cycdurms)
    hist0, bin_edges0 = np.histogram(ISI0, bins)
    histD, bin_edgesD = np.histogram(ISID, bins)
    hist03, bin_edges03 = np.histogram(ISI03, bins)
    histD3, bin_edgesD3 = np.histogram(ISID3, bins)
    EI0 = hist0[0] / ncycles
    EID = histD[0] / ncycles
    EI03 = hist03[0] / ncycles
    EID3 = histD3[0] / ncycles
    return (nAP0, nAPD, VS0, VSD, EI0, EID, nAP03, nAPD3, VS03, VSD3, EI03, EID3)


def makeplots(R, F, CF):
    from matplotlib.ticker import ScalarFormatter
    from scipy.ndimage import gaussian_filter1d

    filtsig = 1.0  # .5 gaussian_filter(AP0, sigma=filtsig, order=0)(as Figure 7/11)
    for idat in range(12):
        R[:, idat] = gaussian_filter1d(
            R[:, idat],
            sigma=filtsig,
            order=0,
        )
    plt.rc("font", family="serif", serif="Times")
    plt.rc("text", usetex=True)
    plt.rc("xtick", labelsize=7)
    plt.rc("ytick", labelsize=7)
    plt.rc("axes", labelsize=8)
    X = np.log2(F / mycf)  # x-axis in octaves re cf
    fh = plt.figure(figsize=(12.0 / 2.54, 5.0 / 2.54))
    ax = fh.subplots(ncols=3, sharex="all")  # ,subplot_kw=dict(xscale="log")
    fh.suptitle("Influence of dendritic inputs on rate, VS and EI @" + str(int(CF)))
    ax[0].plot(X, R[:, 0], "r-", linewidth=0.5)
    ax[0].plot(X, R[:, 1], "k-", linewidth=0.5)
    ax[0].plot(X, R[:, 6], "r--", linewidth=0.5)
    ax[0].plot(X, R[:, 7], "k--", linewidth=0.5)
    ax[0].set_xlim(left=-4, right=2)
    ax[0].get_xaxis().set_major_formatter(ScalarFormatter())
    ax[0].set_xlabel("Frequency (Octaves re CF)")
    ax[0].set_ylabel("Response rate (AP/s)")
    ax[1].plot(X, R[:, 2], "r-", linewidth=0.5)
    ax[1].plot(X, R[:, 3], "k-", linewidth=0.5)
    ax[1].plot(X, R[:, 8], "r--", linewidth=0.5)
    ax[1].plot(X, R[:, 9], "k--", linewidth=0.5)
    ax[1].set_xlim(left=-4, right=2)
    ax[1].get_xaxis().set_major_formatter(ScalarFormatter())
    ax[1].set_xlabel("Frequency (Octaves re CF)")
    ax[1].set_ylabel("Vectorstrength")
    ax[2].plot(X, R[:, 4], "r-", linewidth=0.5)
    ax[2].plot(X, R[:, 5], "k-", linewidth=0.5)
    ax[2].plot(X, R[:, 10], "r--", linewidth=0.5)
    ax[2].plot(X, R[:, 11], "k--", linewidth=0.5)
    ax[2].set_xlim(left=-4, right=2)
    ax[2].set_xticks((-4, -2, 0, 2))
    ax[2].set_ylim((0.0, 0.35))
    ax[2].yaxis.tick_right()
    ax[2].yaxis.set_label_position("right")
    ax[2].get_xaxis().set_major_formatter(ScalarFormatter())
    ax[2].set_xlabel("Frequency (Octaves re CF)")
    ax[2].set_ylabel("Entrainment Index")
    mxI = np.argmax(np.sqrt((R[:, 11] - R[:, 10]) ** 2))
    fmaxdiff = X[mxI]
    if fmaxdiff < 4.0:
        mxI = np.argmin(np.abs(np.abs(X) - 4.0))
        fmaxdiff = X[mxI]
    mxCF = np.argmin(np.abs(X))
    v_atmaxdiff = R[mxI, 11]
    v_atcf = R[mxCF, 11]
    maxdiff = R[mxI, 11] - R[mxI, 10]
    maxdiff_cf = R[mxCF, 11] - R[mxCF, 10]
    ax[2].annotate(
        "+" + str(np.round(maxdiff, 3)),
        xy=(fmaxdiff, v_atmaxdiff),
        xytext=(0.0, 0.25),
        arrowprops=dict(facecolor="black", width=0.5, headwidth=2, headlength=3),
        fontsize=5,
    )
    ax[2].annotate(
        "+" + str(np.round(maxdiff_cf, 3)),
        xy=(0.0, v_atcf),
        xytext=(0.5, 0.15),
        arrowprops=dict(facecolor="black", width=0.5, headwidth=2, headlength=3),
        fontsize=5,
    )
    plt.tight_layout()
    #
    fh2 = plt.figure(figsize=(12.0 / 2.54, 5.0 / 2.54))
    ax = fh2.subplots(ncols=3, sharex="all")  # subplot_kw=dict(xscale="log")
    fh2.suptitle("Change caused by dendritic inputs @" + str(int(CF)))
    ax[0].plot(X, R[:, 1] - R[:, 0], "g-", linewidth=0.5)
    ax[0].plot(X, R[:, 7] - R[:, 6], "g--", linewidth=0.5)
    ax[0].set_xlim(left=-4, right=2)
    ax[0].get_xaxis().set_major_formatter(ScalarFormatter())
    ax[0].set_xlabel("Frequency (Octaves re CF)")
    ax[0].set_ylabel("Rate change (dAP/s)")
    ax[1].plot(X, R[:, 3] - R[:, 2], "g-", linewidth=0.5)
    ax[1].plot(X, R[:, 9] - R[:, 8], "g--", linewidth=0.5)
    ax[1].set_xlim(left=-4, right=2)
    ax[1].get_xaxis().set_major_formatter(ScalarFormatter())
    ax[1].set_xlabel("Frequency (Octaves re CF)")
    ax[1].set_ylabel("VS change")
    ax[2].plot(X, R[:, 5] - R[:, 4], "g-", linewidth=0.5)
    ax[2].plot(X, R[:, 11] - R[:, 10], "g--", linewidth=0.5)
    ax[2].set_xlim(left=-4, right=2)
    ax[2].get_xaxis().set_major_formatter(ScalarFormatter())
    ax[2].set_xlabel("Frequency (Octaves re CF)")
    ax[2].set_ylabel("EI change")
    plt.tight_layout()
    return (fh, fh2)


if __name__ == "__main__":
    # Checking for arguments and/or setting defaults
    addedargs = sys.argv[1:]  #
    myargs = [1, 25, 1500.1, 77]  # weload, nfreqs, cf, db
    for iarg, thisarg in enumerate(addedargs):
        myargs[iarg] = float(thisarg)
    weload = bool(myargs[0])
    nfreqs = int(myargs[1])
    mycf = float(myargs[2])
    dbspl = float(myargs[3])
    filesuffix = str(int(mycf))
    #
    if weload:
        saves = np.load("./results/entrainment_" + filesuffix + ".npy").tolist()
        result = saves["R"]
        mycf = saves["cf"]
        nfreqs = saves["nfreqs"]
        freqs = np.round(np.geomspace(125.0, 2.5 * mycf, nfreqs))
    else:
        result = np.zeros((nfreqs, 12))
        freqs = np.round(np.geomspace(125.0, 2.5 * mycf, nfreqs))
        for icond, thisfreq in enumerate(freqs):
            print(thisfreq)
            result[icond, :] = runonce(icond, cf=mycf, freq=thisfreq)
        saves = {}
        saves["R"] = result
        saves["nfreqs"] = nfreqs
        saves["cf"] = mycf
        np.save(
            "./results/entrainment_" + filesuffix + ".npy", saves, allow_pickle=True
        )
    fighandles = makeplots(result, freqs, mycf)
    pp = PdfPages("./figs/Entrainmenttest_" + filesuffix + ".pdf")
    pp.savefig(fighandles[0])
    pp.savefig(fighandles[1])
    pp.close()
