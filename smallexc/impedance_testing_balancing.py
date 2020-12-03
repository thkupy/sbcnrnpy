#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This function explores whether impedance can be "balanced" by specific values
of gKLT and gH. To speed up things this uses frequency chirps instead of individual
frequencies to analyze impedance.
Part of these results are presented as Fig.8D.
This was a reviewer request from the submission at JNP.

Control by input parameters: load, conds, cores, dur, dt, mp
load: 1 loads existing data, 0 reruns no matter what
conds: how many nxn conditions to test
cores: how many jobs to spawn
dur: duration of chirp
dt: sampling interval
mp: 1 uses parallel processing

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
import scipy.signal as sig
import semodels as sem
import sehelper as seh
import sys
import multiprocessing
from functools import partial

plt.rc("font", family="serif", serif="times")
plt.rcParams["pdf.fonttype"] = 42
plt.rc("text", usetex=True)
plt.rc("xtick", labelsize="x-small")
plt.rc("ytick", labelsize="x-small")
plt.rc("axes", labelsize="small")

hspan = (0.0, 0.015)


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
    fZu = FilteredSignal(np.abs(Zu), 1000.0 / dt, 10.0)
    imax = np.argmax(fZu)
    fres = fax[imax]
    zmax = fZu[imax]
    return (fres, zmax, fZu, fax)


def runonce(x, P):
    print(str(x) + ": " + str(P["gKLT_d"][x]) + "/" + str(P["gIH_d"][x]))
    simdur = P["dur"][x]
    dt = P["dt"][x]
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
    impamp = 0.1
    impphase = 0.0
    impfreq = (10.0, 250.0)
    R = sem.SE_BS(
        S,
        Simdur=P["dur"][x],
        dt=dt,
        N=N,
        gKLT_d=P["gKLT_d"][x],
        gIH_d=P["gIH_d"][x],
        EoHseed=23947543,
        impedancetesting=-1,  # chirp mode in semodels.SE_BS
        impstim=(impamp, impfreq, impphase),
    )
    fres, zmax, Pxx, f = analyze_chirpdata(R["Vm"], dt, simdur)
    return (fres, zmax)


def mphandler(P):
    poolmess = partial(runonce, P=P)
    if P["mp"] == 1:
        with multiprocessing.Pool(P["cores"]) as activepool:
            r = activepool.map(poolmess, P["Number"])
    else:
        r = list(map(poolmess, P["Number"]))  # for debug
    return r


def makefig(R, P):
    from matplotlib.ticker import ScalarFormatter, AutoMinorLocator, MaxNLocator
    from scipy import interpolate

    #
    nticks = 5
    frestemp = R[0][:, 0]  # return [SAPr,Smeanv,Smodv]
    zmaxtemp = R[0][:, 1]
    fres = np.reshape(frestemp, (P["N"], P["N"]))
    zmax = np.reshape(zmaxtemp, (P["N"], P["N"]))
    fh = plt.figure(figsize=(16.0 / 2.54, 10.0 / 2.54))
    #
    ax1 = plt.subplot(221)
    CS1 = ax1.contourf(
        np.unique(P["gIH_d"]),
        np.unique(P["gKLT_d"]),
        fres,
        21,
    )  # repeated = y, tiled = x!!
    cs = ax1.contour(
        np.unique(P["gIH_d"]),
        np.unique(P["gKLT_d"]),
        fres,
        levels=[
            150.0,
        ],
    )
    isofrespath = cs.collections[0].get_paths()
    isofresverts = isofrespath[0].vertices.copy()
    isofres_h = isofresverts[:, 0]
    isofres_klt = isofresverts[:, 1]
    h_hr = np.linspace(hspan[0], hspan[1], 1000)
    interpfun = interpolate.interp1d(isofres_h, isofres_klt, fill_value="extrapolate")
    klt_hr = interpfun(h_hr)
    plt.title("fres")
    plt.ylabel("$gKLT_d$")
    plt.xlabel("$gH_d$")
    ax1.yaxis.set_major_locator(MaxNLocator(nticks))
    ax1.xaxis.set_major_locator(MaxNLocator(nticks))
    cbar1 = plt.colorbar(CS1, use_gridspec=True)
    cbar1.ax.set_ylabel("f(Hz)")
    tl = MaxNLocator(nbins=5)
    cbar1.locator = tl
    cbar1.update_ticks()
    #
    ax2 = plt.subplot(222)
    CS2 = ax2.contourf(
        np.unique(P["gIH_d"]),
        np.unique(P["gKLT_d"]),
        zmax,
        21,
    )  # repeated = y, tiled = x!!
    ax2.plot(h_hr, klt_hr)
    plt.title("Zmax")
    plt.ylabel("$gKLT_d$")
    plt.xlabel("$gH_d$")
    ax2.yaxis.set_major_locator(MaxNLocator(nticks))
    ax2.xaxis.set_major_locator(MaxNLocator(nticks))
    cbar2 = plt.colorbar(CS2, use_gridspec=True)
    cbar2.ax.set_ylabel("Z(Megohm)")
    tl = MaxNLocator(nbins=5)
    cbar2.locator = tl
    cbar2.update_ticks()
    #
    ax3 = plt.subplot(223)
    CS3 = ax3.contourf(
        np.unique(P["gIH_d"]),
        np.unique(P["gKLT_d"]),
        fres,
        21,
    )  # repeated = y, tiled = x!!
    plt.title("fres")
    plt.ylabel("$gKLT_d$")
    plt.xlabel("$gH_d$")
    ax3.yaxis.set_major_locator(MaxNLocator(nticks))
    ax3.xaxis.set_major_locator(MaxNLocator(nticks))
    cbar3 = plt.colorbar(CS3, use_gridspec=True)
    cbar3.ax.set_ylabel("f(Hz)")
    tl = MaxNLocator(nbins=5)
    cbar3.locator = tl
    cbar3.update_ticks()
    #
    ax4 = plt.subplot(224)
    h_hr = np.linspace(hspan[0], hspan[1], 1000)
    interpfun = interpolate.interp1d(isofres_h, isofres_klt, fill_value="extrapolate")
    klt_hr = interpfun(h_hr)
    klt_hr[klt_hr < 0.0] = 0.0
    CS4 = ax4.contourf(
        np.unique(P["gIH_d"]),
        np.unique(P["gKLT_d"]),
        zmax,
        21,
    )  # repeated = y, tiled = x!!
    cs = ax4.contour(
        np.unique(P["gIH_d"]),
        np.unique(P["gKLT_d"]),
        zmax,
        levels=[
            18.33,
        ],
    )
    isozmaxpath = cs.collections[0].get_paths()
    isozmaxverts = isozmaxpath[0].vertices.copy()
    isofres_h2 = isozmaxverts[:, 0]
    isofres_klt2 = isozmaxverts[:, 1]
    h_hr2 = np.linspace(hspan[0], hspan[1], 1000)
    interpfun = interpolate.interp1d(isofres_h2, isofres_klt2, fill_value="extrapolate")
    klt_hr2 = interpfun(h_hr2)
    klt_hr2[klt_hr2 < 0.0] = 0.0
    ax3.plot(h_hr2, klt_hr2, "r-")
    plt.title("Zmax")
    plt.ylabel("$gKLT_d$")
    plt.xlabel("$gH_d$")
    ax4.yaxis.set_major_locator(MaxNLocator(nticks))
    ax4.xaxis.set_major_locator(MaxNLocator(nticks))
    cbar4 = plt.colorbar(CS4, use_gridspec=True)
    cbar4.ax.set_ylabel("Z(Megohm)")
    tl = MaxNLocator(nbins=5)
    cbar4.locator = tl
    cbar4.update_ticks()
    plt.tight_layout()
    out = {"h": h_hr, "klt": klt_hr, "h2": h_hr2, "klt2": klt_hr2}
    np.save("./results/fres_comp_klt_h.npy", out, allow_pickle=True)
    return fh


if __name__ == "__main__":
    addedargs = sys.argv[1:]
    myargs = [1, 3, 4, 3000.0, 0.01, 1]
    for iarg, thisarg in enumerate(addedargs):
        myargs[iarg] = thisarg
    nconds = int(myargs[1])
    ncores = int(myargs[2])
    dur = float(myargs[3])
    dt = float(myargs[4])
    mp = int(myargs[5])
    if int(myargs[0]) == 1:
        weload = True
    else:
        weload = False
    #
    if os.path.isfile("./results/impedance_balancing.npy") and weload:
        print("Data for impedance_balancing found... loading!")
        output = np.load("./results/impedance_balancing.npy", allow_pickle=True)
        P = np.load("./results/impedance_balancing_P.npy", allow_pickle=True)
        P = P.tolist()
    else:
        # Model Parameters (all in a linearly aranged fashion, so that a minimal
        # amount of programming is required to change the experiment).
        P = {}
        P["N"] = int(nconds)  # just to be sure
        P["cores"] = int(ncores)  # just to be sure
        P["mp"] = mp
        P["TotalN"] = int(P["N"] ** 2)
        P["Number"] = list(range(P["TotalN"]))
        P["dur"] = np.repeat(dur, P["TotalN"])
        P["dt"] = np.repeat(dt, P["TotalN"])
        P["gLEAK_d"] = np.repeat(0.001, P["TotalN"])
        P["gKLT_s"] = np.repeat(0.017, P["TotalN"])
        P["gIH_s"] = np.repeat(0.002, P["TotalN"])
        # Now define the two variable parameters. The repeated = y, the tiled = x!!
        allKLT = np.linspace(0.0, 0.04, P["N"])
        allH = np.linspace(hspan[0], hspan[1], P["N"])
        P["gKLT_d"] = np.repeat(allKLT, P["N"])
        P["gIH_d"] = np.tile(allH, P["N"])
        # make go!
        output = []
        output.append(mphandler(P))
        output = np.array(output)
        np.save("./results/impedance_balancing.npy", output, allow_pickle=True)
        np.save("./results/impedance_balancing_P.npy", P, allow_pickle=True)
    #
    ###Plot results
    fhandle = makefig(
        R=output,
        P=P,
    )
    pp = PdfPages("./figs/impedance_balancing.pdf")
    pp.savefig()
    pp.close()
