#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This function explores the frequency-dependence of the dendrite-to-soma transfer 
resistance and phase. This will be presented as Fig.10D/E/G/H.
This was a reviewer request from the submission at JNP.

You can pass arguments to this function as in:
    python3 impedance_testing_3d.py weload nconds ncores cell

    weload: 1 = load data instead of calculating new (use 0 or anything to recalc)
    nconds: Nr. of conditions varied for Ra (processing time)
    ncores: how many processes you'd like to spaw
    cell: which cell (1,2 or 3 at the moment (see semodels.BS_3D))

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
import scipy.signal as sig
import time

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
    Ra,
    imptesttype,
    impdendtype,
    impdendid,
    cell,
):
    time.sleep(np.random.uniform())  # testing whether this helps parallel processing
    # errors caused by reading the sloc files..-
    print(str(x))
    simdur = 1000.0  # 1000
    pretime = 100.0
    dt = 0.01
    N = 32
    S = np.zeros((N, 1))
    # T = np.linspace(0, simdur - dt, int(round(simdur / dt)))
    impamp = 0.1
    impphase = 0.0
    # IS = impamp * np.sin(impfreq[x]/1000.0 * 2 * np.pi * T + impphase)
    if imptesttype[x] == 1:
        thisfreq = impfreq[x]
    else:
        thisfreq = impfreq
    R = sem.SE_3D(
        S,
        Simdur=simdur,
        dt=dt,
        cell=cell[x],
        Ra=Ra[x],
        EoHseed=8347236,
        impedancetesting=imptesttype[x],
        impstim=(impamp, thisfreq, impphase),
        impdendtype=impdendtype[x],
        impdendid=impdendid[x],
    )
    if imptesttype[x] == 1:
        A = average_over_cyle(R["Vm"], impfreq[x], dt, simdur)
        cytax = np.linspace(0, 1, A[0].size)
        vmax = np.max(A[0])
        imax = np.argmax(A[0])
        vmin = np.min(A[0])
        imin = np.argmin(A[0])
        r_megohm = (vmax - A[2]) / impamp
        phi = cytax[imax]
        return (r_megohm, phi)
    elif imptesttype[x] == -1:
        fres, zmax, Pxx, f = analyze_chirpdata(R["Vm"], dt, simdur)
        return (fres, zmax, Pxx, f)


def runcurve(
    Ra=150.0,
    nconds=41,
    fmin=10.0,
    fmax=250.0,
    impdendtype=0,
    impdendid=-1,
    cell=1,
):
    R = []
    phi = []
    print("Ra: " + str(Ra))  # DEBUG
    Ra = np.repeat(Ra, nconds)
    freqs = np.round(np.linspace(fmin, fmax, nconds), 1)
    imptesttype = np.repeat(1, nconds)
    impdendtype = np.repeat(impdendtype, nconds)
    impdendid = np.repeat(impdendid, nconds)
    cell = np.repeat(cell, nconds)
    poolmess = partial(
        runonce,
        impfreq=freqs,
        Ra=Ra,
        imptesttype=imptesttype,
        impdendtype=impdendtype,
        impdendid=impdendid,
        cell=cell,
    )
    with multiprocessing.Pool(ncores) as activepool:
        r = activepool.map(poolmess, range(nconds))
    return (freqs, r)


def runmapping(
    impdendtype=0,
    ndends=5,
    cell=1,
):
    impfreqs = [10.0, 250.0]
    Ra = np.repeat(150.0, ndends)
    imptesttype = np.repeat(-1, ndends)
    impdendtype = np.repeat(impdendtype, ndends)
    impdendid = np.arange(ndends)
    cell = np.repeat(cell, ndends)
    poolmess = partial(
        runonce,
        impfreq=impfreqs,
        Ra=Ra,
        imptesttype=imptesttype,
        impdendtype=impdendtype,
        impdendid=impdendid,
        cell=cell,
    )
    with multiprocessing.Pool(ncores) as activepool:
        r = activepool.map(poolmess, range(ndends))
    return r


def mapfig(N, R, dendtype, cellnr):
    M = sem.SE_3D(
        S=np.zeros((32, 1)),
        Simdur=10.0,
        dt=0.1,
        cell=cellnr,
        impedancetesting=1,
        impdendtype=dendtype,
        impdendid=N,
        debug=1,
    )
    thishand = plt.figure()
    plt.semilogx(R[3], R[2], "k-")
    plt.plot(R[0], R[1], "ro")
    plt.plot(R[3][0], R[2][0], "mo")
    plt.text(R[0], R[1], str(R[0]) + "Hz/" + str(R[1]) + "$(M\Omega)$")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("$Z (M\Omega)$")
    if dendtype == 0:
        plt.title("Cell " + str(cellnr) + ", distal dend " + str(N))
    else:
        plt.title("Cell " + str(cellnr) + ", proximal dend " + str(N))
    plt.ylim((0, 15))
    plt.tight_layout()
    return thishand, M["impathlen"], R[0], R[1], R[2][0]


def makefig(R, X, xlabeltext, suptitletext):
    fh = plt.figure(figsize=(16.0 / 2.54, 10.0 / 2.54))
    fh.suptitle(suptitletext, fontsize=12)
    sp1 = fh.add_subplot(2, 3, 1)
    sp1.set_xscale("log")
    sp1.set_xlabel("Input Freq. (Hz)")
    sp1.set_ylabel(r"$Z_{T} (M\Omega)$")
    sp1.set_xticks((10, 20, 50, 100, 200))
    sp2 = fh.add_subplot(2, 3, 2)
    sp2.set_xscale("log")
    sp2.set_xlabel("Input Freq. (Hz)")
    sp2.set_ylabel("Phaseshift (rad)")
    sp2.set_xticks((10, 20, 50, 100, 200))
    sp3 = fh.add_subplot(2, 3, 4)
    sp3.set_xlabel(xlabeltext)
    sp3.set_ylabel(r"Res. Amp. ($M\Omega$)")
    sp3.xaxis.set_major_locator(plt.MaxNLocator(4))
    sp4 = fh.add_subplot(2, 3, 5)
    sp4.set_ylabel("Res. Freq. (Hz)")
    sp4.set_xlabel(xlabeltext)
    sp4.xaxis.set_major_locator(plt.MaxNLocator(4))
    sp5 = fh.add_subplot(2, 3, 6)
    sp5.set_ylabel("Phaseshift @100Hz (rad)")
    sp5.set_xlabel(xlabeltext)
    sp5.xaxis.set_major_locator(plt.MaxNLocator(4))
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
    sp3.plot(X, relres, "k")
    sp4.plot(X, resfreq, "k")
    sp5.plot(X, phaseatfreq, "k")
    plt.tight_layout()
    return fh


def FilteredSignal(signal, fs, cutoff):
    B, A = sig.butter(1, cutoff / (fs / 2), btype="low")
    filtered_signal = sig.filtfilt(B, A, signal, axis=0)
    return filtered_signal


def analyze_chirpdata(R, dt, dur):
    t = np.linspace(dt, dur, int(dur / dt)) / 1000.0  # in s
    w = 0.1 * sig.chirp(t, f0=10.0, f1=250.0, t1=dur / 1000.0)  # in Hz and s
    fax = np.linspace(10.0, 250.0, w.size)
    Zraw = (R - np.mean(R)) / 0.1
    Zu = sig.hilbert(Zraw - np.mean(Zraw))
    fZu = FilteredSignal(np.abs(Zu), 1000.0 / dt, 5.0)
    imax = np.argmax(fZu)
    fres = fax[imax]
    zmax = fZu[imax]
    print("fres=" + str(fres))
    return (fres, zmax, fZu, fax)


if __name__ == "__main__":
    addedargs = sys.argv[1:]  # args are (fname, weload, nconds, ncores, cell)
    myargs = [1, 5, 4, 1]
    for iarg, thisarg in enumerate(addedargs):
        myargs[iarg] = int(thisarg)
    nconds = myargs[1]
    ncores = myargs[2]
    cell = myargs[3]
    if myargs[0] == 1:
        weload = True
    else:
        weload = False
    #
    nprox, ndist = sem.SE_3D(np.array((0,)), cell=cell, listinfo=True)
    allprox = np.arange(nprox)
    alldist = np.arange(ndist)
    print("cell" + str(cell) + ": NProx=" + str(nprox) + ", NDist=" + str(ndist))
    if weload and os.path.isfile("./results/impedance_res_3d_R_" + str(cell) + ".npy"):
        resR = np.load("./results/impedance_res_3d_R_" + str(cell) + ".npy")
        resR = resR.tolist()
        datashape_R = np.shape(resR)
        allRa = np.linspace(25.0, 300.0, datashape_R[0])
        #
        resD = np.load("./results/impedance_res_3d_D_" + str(cell) + ".npy")
        resD = resD.tolist()
        resP = np.load("./results/impedance_res_3d_P_" + str(cell) + ".npy")
        resP = resP.tolist()
    else:
        allRa = np.linspace(10.0, 250.0, nconds)
        resR = []
        for icond in range(nconds):
            resR.append(runcurve(Ra=allRa[icond], cell=cell))
        np.save("./results/impedance_res_3d_R_" + str(cell) + ".npy", resR)
        resD = runmapping(impdendtype=0, ndends=ndist, cell=cell)
        np.save("./results/impedance_res_3d_D_" + str(cell) + ".npy", resD)
        resP = runmapping(impdendtype=1, ndends=nprox, cell=cell)
        np.save("./results/impedance_res_3d_P_" + str(cell) + ".npy", resP)
    #
    fhandle1 = makefig(resR, allRa, "Ra (Ohm cm)", "Axial Resistance")
    fhandles2 = []
    PL = []
    FR = []
    ZM = []
    ZR = []
    PLp = []
    FRp = []
    ZMp = []
    ZRp = []
    fhandles3 = []
    #
    with PdfPages("./figs/Impedance_3d_Ra_" + str(cell) + ".pdf") as pp:
        pp.savefig(fhandle1)
    for ndend, dendres in enumerate(resD):
        fh, pathlen, fres, zmax, z0 = mapfig(ndend, dendres, 0, cell)
        fhandles2.append(fh)
        PL.append(pathlen)
        FR.append(fres)
        ZM.append(zmax)
        ZR.append(zmax - z0)
    with PdfPages("./figs/Impedance_3d_AllDist_" + str(cell) + ".pdf") as pp:
        for thishandle in fhandles2:
            pp.savefig(thishandle)
            plt.close(thishandle)
    for ndend, dendres in enumerate(resP):
        fh, pathlen, fres, zmax, z0 = mapfig(ndend, dendres, 1, cell)
        PLp.append(pathlen)
        FRp.append(fres)
        ZMp.append(zmax)
        ZRp.append(zmax - z0)
        fhandles3.append(fh)
    with PdfPages("./figs/Impedance_3d_AllProx_" + str(cell) + ".pdf") as pp:
        for thishandle in fhandles3:
            pp.savefig(thishandle)
            plt.close(thishandle)
    #
    allPL = np.concatenate((PL, PLp))
    allZM = np.concatenate((ZM, ZMp))
    allFR = np.concatenate((FR, FRp))
    allZR = np.concatenate((ZR, ZRp))
    fy = np.polyfit(PL, ZM, 1)
    fya = np.polyfit(allPL, allZM, 1)
    fy2 = np.polyfit(PL, FR, 1)
    fy2a = np.polyfit(allPL, allFR, 1)
    fy3 = np.polyfit(PL, ZR, 1)
    fy3a = np.polyfit(allPL, allZR, 1)
    fx = np.linspace(np.min(allPL), np.max(allPL), 100)
    #
    fh = plt.figure()
    plt.subplot(311)
    plt.plot(fx, fy[0] * fx + fy[1], "r-")
    plt.plot(fx, fya[0] * fx + fya[1], "k-")
    anno1 = "y=" + str(round(fy[0], 3)) + "*x+" + str(round(fy[1], 3))
    anno2 = "y=" + str(round(fya[0], 3)) + "*x+" + str(round(fya[1], 3))
    plt.annotate(anno1, xy=(250, 16), color="r")
    plt.annotate(anno2, xy=(250, 14), color="k")
    plt.plot(PL, ZM, "r.")
    plt.plot(PLp, ZMp, "m.")
    plt.xlim((0, 300))
    plt.ylim((0, 17))
    plt.xlabel("path length ($\mu m$)")
    plt.ylabel("Zmax ($M\Omega$)")
    plt.subplot(312)
    plt.plot(fx, fy2[0] * fx + fy2[1], "r-")
    plt.plot(fx, fy2a[0] * fx + fy2a[1], "k-")
    anno3 = "y=" + str(round(fy2[0], 3)) + "*x+" + str(round(fy2[1], 3))
    anno4 = "y=" + str(round(fy2a[0], 3)) + "*x+" + str(round(fy2a[1], 3))
    plt.plot(PL, FR, "r.")
    plt.plot(PLp, FRp, "m.")
    plt.annotate(anno3, xy=(250, 140), color="r")
    plt.annotate(anno4, xy=(250, 120), color="k")
    plt.xlabel("path length ($\mu m$)")
    plt.ylabel("$f_{res} (Hz)$")
    plt.xlim((0, 300))
    plt.ylim((0, 150))
    plt.subplot(313)
    plt.plot(fx, fy3[0] * fx + fy3[1], "r-")
    plt.plot(fx, fy3a[0] * fx + fy3a[1], "k-")
    anno5 = "y=" + str(round(fy3[0], 3)) + "*x+" + str(round(fy3[1], 3))
    anno6 = "y=" + str(round(fy3a[0], 3)) + "*x+" + str(round(fy3a[1], 3))
    plt.plot(PL, ZR, "r.")
    plt.plot(PLp, ZRp, "m.")
    plt.annotate(anno5, xy=(250, 11), color="r")
    plt.annotate(anno6, xy=(250, 9), color="k")
    plt.xlabel("path length ($\mu m$)")
    plt.ylabel("$Z_{max} - Z_{0} (M\Omega)$")
    plt.xlim((0, 300))
    plt.ylim((-5, 12))

    plt.tight_layout()
    with PdfPages("./figs/imp3d_fres_zmax_vs_path_" + str(cell) + ".pdf") as pp:
        pp.savefig(fh)
        plt.close(fh)
    # Generate and save per cell result
    SD = {}
    SD["fx"] = fx
    SD["fy"] = fy
    SD["fya"] = fya
    SD["fy2"] = fy2
    SD["fy2a"] = fy2a
    SD["fy3"] = fy3
    SD["fy3a"] = fy3a
    SD["PL"] = PL
    SD["PLp"] = PLp
    SD["FR"] = FR
    SD["FRp"] = FRp
    SD["ZM"] = ZM
    SD["ZMp"] = ZMp
    SD["ZR"] = ZR
    SD["ZRp"] = ZRp
    np.save("./results/impedance_res_3d_summarydata_cell" + str(cell) + ".npy", SD)
