#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This function explores the frequency-dependence of the dendrite-to-soma transfer 
resistance and phase. This will be presented as Fig.1E/F, Fig.2F/G and Fig.5H/J.
This was a reviewer request from the submission at JNP.

You can pass arguments to this function as in:
    python3 impedance_testing_stylized.py weload nconds ncores

    weload: 1 = load data instead of calculating new (use 0 or anything to recalc)
    nconds: Nr. of conditions varied for L, D, gKLT_d, IH_d and Ra
    ncores: how many processes you'd like to spawn

Created September 2020
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

ncores = 1

plt.rc("font", family="serif", serif="times")
plt.rcParams["pdf.fonttype"] = 42
plt.rc("text", usetex=True)
plt.rc("xtick", labelsize="x-small")
plt.rc("ytick", labelsize="x-small")
plt.rc("axes", labelsize="small")


def average_over_cyle(D, freq, dt, dur):
    IPI = 1000.0 / freq
    allc = []
    cycsamps = int(np.round(IPI / dt))
    c = list(range(0, len(D) - cycsamps, cycsamps))
    for ii, cc in enumerate(c):
        allc.append(D[cc : cc + cycsamps])
    av = np.mean(np.array(allc), axis=0)
    st = np.std(np.array(allc), axis=0)
    meanD = np.mean(av)
    stdD = np.std(av)
    modamp = np.max(av) - np.min(av)
    return (av, st, meanD, stdD, modamp)


def runonce(
    x,
    impfreq,
    L,
    D,
    gKLT_d,
    gIH_d,
    Ra,
):
    simdur = 1000.0  # 1000
    pretime = 100.0
    dt = 0.01
    N = 1
    S = np.array(
        [
            [
                0.0,
            ],
            [
                0.0,
            ],
        ]
    )
    T = np.linspace(0, simdur - dt, int(round(simdur / dt)))
    impamp = 0.1
    impphase = 0.0
    IS = impamp * np.sin(impfreq[x] / 1000.0 * 2 * np.pi * T + impphase)
    R = sem.SE_BS(
        S,
        Simdur=simdur,
        dt=dt,
        L=L[x],
        D=D[x],
        N=N,
        gKLT_d=gKLT_d[x],
        gIH_d=gIH_d[x],
        EoHseed=8347236,
        impedancetesting=1,
        impstim=(impamp, impfreq[x], impphase),
        Ra=Ra[x],
    )
    A = average_over_cyle(R["Vm"], impfreq[x], dt, simdur)
    cytax = np.linspace(0, 1, A[0].size)
    vmax = np.max(A[0])
    imax = np.argmax(A[0])
    vmin = np.min(A[0])
    imin = np.argmin(A[0])
    r_megohm = (vmax - A[2]) / impamp
    phi = cytax[imax]
    return (r_megohm, phi)


def runcurve(
    L=50.0,
    D=3.0,
    gKLT_d=0.0085,
    gIH_d=0.001,
    Ra=150.0,
    nconds=50,
    fmin=10.0,
    fmax=250.0,
):
    R = []
    phi = []
    L = np.repeat(L, nconds)
    D = np.repeat(D, nconds)
    gKLT_d = np.repeat(gKLT_d, nconds)
    gIH_d = np.repeat(gIH_d, nconds)
    print(Ra)
    Ra = np.repeat(Ra, nconds)
    freqs = np.round(np.linspace(fmin, fmax, nconds), 1)
    poolmess = partial(
        runonce,
        impfreq=freqs,
        L=L,
        D=D,
        gKLT_d=gKLT_d,
        gIH_d=gIH_d,
        Ra=Ra,
    )
    with multiprocessing.Pool(ncores) as activepool:
        r = activepool.map(poolmess, range(nconds))
    return (freqs, r)


def makefig(R, X, xlabeltext, suptitletext, ylimset):
    from scipy.ndimage import gaussian_filter1d

    fh = plt.figure(figsize=(16.0 / 2.54, 10.0 / 2.54))
    fh.suptitle(suptitletext, fontsize=12)
    sp1 = fh.add_subplot(2, 3, 1)
    sp1.set_xscale("log")
    sp1.set_xlabel("Input Freq. (Hz)")
    sp1.set_ylabel(r"$Z_{T} (M\Omega)$")
    sp1.set_xticks((10, 20, 50, 100, 200))
    sp1.set_ylim(ylimset[0])
    sp2 = fh.add_subplot(2, 3, 2)
    sp2.set_xscale("log")
    sp2.set_xlabel("Input Freq. (Hz)")
    sp2.set_ylabel("Phaseshift (rad)")
    sp2.set_xticks((10, 20, 50, 100, 200))
    sp2.set_ylim(ylimset[1])
    sp3 = fh.add_subplot(2, 3, 4)
    sp3.set_xlabel(xlabeltext)
    sp3.set_ylabel(r"Res. Amp. ($M\Omega$)")
    sp3.xaxis.set_major_locator(plt.MaxNLocator(4))
    sp3.set_ylim(ylimset[2])
    sp4 = fh.add_subplot(2, 3, 5)
    sp4.set_ylabel("Res. Freq. (Hz)")
    sp4.set_xlabel(xlabeltext)
    sp4.xaxis.set_major_locator(plt.MaxNLocator(4))
    sp4.set_ylim(ylimset[3])
    sp5 = fh.add_subplot(2, 3, 6)
    sp5.set_ylabel("Phaseshift @100Hz (rad)")
    sp5.set_xlabel(xlabeltext)
    sp5.xaxis.set_major_locator(plt.MaxNLocator(4))
    sp5.set_ylim(ylimset[4])
    relres = []
    maxres = []
    resfreq = []
    phaseatfreq = []
    col = np.linspace(0.0, 0.8, len(X))
    for ix, x in enumerate(X):
        dat = np.array(R[ix][1])
        mycol = (col[ix], col[ix], col[ix])
        normr = np.array(dat[:, 0])
        sp1.plot(R[ix][0], normr, color=mycol)
        sp1.text(R[ix][0][-1], normr[-1], str(round(x, 4)), fontsize=7)
        phasedat = np.unwrap((0.25 - dat[:, 1]))
        sp2.plot(R[ix][0], phasedat, color=mycol)
        maxres.append(np.max(normr))
        relres.append(np.max(normr - dat[0, 0]))
        imaxres = np.argmax(normr)
        phaseatfreq.append(phasedat[18])
        resfreq.append(R[ix][0][imaxres])
        sp1.plot(resfreq[-1], maxres[-1], "kd", markersize=5)
    # sp3.plot(X, relres, "k.")
    sp3.plot(
        X, gaussian_filter1d(relres, sigma=1.0), "k"
    )  # mild filter (4 freq spacing)
    # sp4.plot(X, resfreq, "k.")
    sp4.plot(
        X, gaussian_filter1d(resfreq, sigma=1.0), "k"
    )  # mild filter (4 freq spacing)
    # sp5.plot(X, phaseatfreq, "k.")
    sp5.plot(
        X, gaussian_filter1d(phaseatfreq, sigma=1.0), "k"
    )  # mild filter (4 freq spacing)
    plt.tight_layout()
    return fh


if __name__ == "__main__":
    addedargs = sys.argv[1:]  # args are (fname, weload, nconds, ncores)
    myargs = [1, 5, 4]
    for iarg, thisarg in enumerate(addedargs):
        myargs[iarg] = int(thisarg)
    nconds = myargs[1]
    ncores = myargs[2]
    if myargs[0] == 1:
        weload = True
    else:
        weload = False
    #
    if weload:
        resL = np.load("./results/impedance_res_stylized_L.npy")
        resL = resL.tolist()
        datashape_L = np.shape(resL)
        allL = np.linspace(5.0, 500.0, datashape_L[0])
        #
        resD = np.load("./results/impedance_res_stylized_D.npy")
        resD = resD.tolist()
        datashape_D = np.shape(resD)
        allD = np.linspace(0.1, 4.0, datashape_D[0])
        #
        resK = np.load("./results/impedance_res_stylized_K.npy")
        resK = resK.tolist()
        datashape_K = np.shape(resK)
        allK = np.linspace(0.0, 0.04, datashape_K[0])
        #
        resH = np.load("./results/impedance_res_stylized_H.npy")
        resH = resH.tolist()
        datashape_H = np.shape(resH)
        allH = np.linspace(0.0, 0.005, datashape_H[0])
        #
        resR = np.load("./results/impedance_res_stylized_R.npy")
        resR = resR.tolist()
        datashape_R = np.shape(resR)
        allRa = np.linspace(25.0, 300.0, datashape_R[0])
    else:
        allL = np.linspace(5.0, 500.0, nconds)
        allD = np.linspace(0.1, 4.0, nconds)
        allK = np.linspace(0.0, 0.04, nconds)
        allH = np.linspace(0.0, 0.005, nconds)
        allRa = np.linspace(10.0, 250.0, nconds)
        resL = []
        resD = []
        resK = []
        resH = []
        resR = []
        for icond in range(nconds):
            resL.append(runcurve(L=allL[icond]))
        np.save("./results/impedance_res_stylized_L.npy", resL)
        for icond in range(nconds):
            resD.append(runcurve(D=allD[icond]))
        np.save("./results/impedance_res_stylized_D.npy", resD)
        for icond in range(nconds):
            resK.append(runcurve(gKLT_d=allK[icond]))
        np.save("./results/impedance_res_stylized_K.npy", resK)
        for icond in range(nconds):
            resH.append(runcurve(gIH_d=allH[icond]))
        np.save("./results/impedance_res_stylized_H.npy", resH)
        for icond in range(nconds):
            resR.append(runcurve(Ra=allRa[icond]))
        np.save("./results/impedance_res_stylized_R.npy", resR)
    #
    fhandles = []
    ylimset1 = ((0.0, 28.0), (-0.4, 0.05), (0.0, 4.0), (85.0, 160.0), (-0.11, -0.03))
    ylimset2 = (
        (13.5, 22.1),
        (-0.25, 0.05),
        (2.0, 4.5),
        (110.0, 165.0),
        (-0.05, -0.019),
    )
    ylimset3 = ((15.9, 23.0), (-0.25, 0.05), (0.0, 4.5), (0.0, 170.0), (-0.1, 0.0))
    fhandles.append(makefig(resL, allL, "L(dend1) (µm)", "Length of Dend1", ylimset1))
    fhandles.append(makefig(resD, allD, "D(dend1) (µm)", "Diameter of Dend1", ylimset1))
    fhandles.append(makefig(resK, allK, "gKLT(dend1) (µS)", "gKLT of Dend1", ylimset2))
    fhandles.append(makefig(resH, allH, "gH(dend1) (µS)", "gH of Dend1", ylimset2))
    fhandles.append(
        makefig(resR, allRa, r"$Ra (M\Omega)$", "Axial Resistance", ylimset3)
    )
    with PdfPages("./figs/Impedance_stylized.pdf") as pp:
        for fhandle in fhandles:
            pp.savefig(fhandle)
            plt.close(fhandle)
