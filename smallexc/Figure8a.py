#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 28 Aug 2019
This file contains all code necessary to generate Figure 8, part a:
    L vs. D on spiking
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
import semodels as sem
import sehelper as seh
import multiprocessing
from functools import partial

# some global plot settings
plt.rc("font", family="serif", serif="times")
plt.rcParams["pdf.fonttype"] = 42
plt.rc("text", usetex=True)
plt.rc("xtick", labelsize="x-small")
plt.rc("ytick", labelsize="x-small")
plt.rc("axes", labelsize="small")


def runonce(x, P):
    ####GENERATE SPIKE INPUTS (the spiketimes are always the same to improve comparability)
    thisdur = (P["dur"][x] - 5.0) / 1000.0
    N = P["Nsyn"][x]
    anf_num = (N + 1, 0, 0)
    if P["AdvZseed"]:
        thisZseed = P["Zseed"] + x
    else:
        thisZseed = P["Zseed"]
    if P["AdvEoHseed"]:
        thisEoHseed = P["EoHseed"] + x
    else:
        thisEoHseed = P["EoHseed"]
    SO = seh.create_sound(
        fs=100e3, freq=P["freq"][x], duration=thisdur, dbspl=P["dB"][x]
    )
    SP = seh.create_spk(
        SO, fs=100e3, N=P["Nsyn"][x] + 1, cf=P["cf"][x], seed=thisZseed, anf_num=anf_num
    )  # generate enough spiketrains
    S = np.array(SP["spikes"] * 1000.0)
    R0 = sem.SE_BS(
        S,
        Simdur=P["dur"][x],
        dt=P["dt"][x],
        G_EoH=P["G_EoH"][x],
        StochasticEoH=P["Stochastic_EoH"][x],
        N=P["Nsyn"][x],
        G_Dend=0.0,
        L=P["L"][x],
        D=P["D"][x],
        constR=P["constR"][x],
        gKLT_d=P["gKLT_d"][x],
        gIH_d=P["gIH_d"][x],
        gLEAK_d=P["gLEAK_d"][x],
        EoHseed=thisEoHseed,
    )
    RD = sem.SE_BS(
        S,
        Simdur=P["dur"][x],
        dt=P["dt"][x],
        G_EoH=P["G_EoH"][x],
        StochasticEoH=P["Stochastic_EoH"][x],
        N=P["Nsyn"][x],
        G_Dend=P["gsyn"][x],
        L=P["L"][x],
        D=P["D"][x],
        constR=P["constR"][x],
        gKLT_d=P["gKLT_d"][x],
        gIH_d=P["gIH_d"][x],
        gLEAK_d=P["gLEAK_d"][x],
        EoHseed=thisEoHseed,
    )
    Ev0 = seh.SimpleDetectAP(R0["Vm"], thr=-100, dt=P["dt"][x], LM=-20, RM=10)
    EvD = seh.SimpleDetectAP(RD["Vm"], thr=-100, dt=P["dt"][x], LM=-20, RM=10)
    APr0 = len(Ev0["PeakT"]) / (P["dur"][x] / 1000.0)
    APrD = len(EvD["PeakT"]) / (P["dur"][x] / 1000.0)
    VS0, phi0, Ray0, phases0 = seh.vectorstrength(
        Ev0["PeakT"], P["freq"][x], [0.0, P["dur"][x]]
    )
    VSD, phiD, RayD, phasesD = seh.vectorstrength(
        EvD["PeakT"], P["freq"][x], [0.0, P["dur"][x]]
    )
    print((str(x)))
    return [VS0, phi0, Ray0, APr0, VSD, phiD, RayD, APrD]


def myMPhandler(P):
    p = multiprocessing.Pool(P["cores"])
    poolmess = partial(runonce, P=P)
    if P["mp"]:
        r = p.map(poolmess, P["Number"])
    else:  # for debug
        r = list(map(poolmess, P["Number"]))  # for debug
    return r


def plotres(output, P, x, y, xlabs, ylabs):
    from matplotlib.ticker import ScalarFormatter, AutoMinorLocator, MaxNLocator
    import scipy.ndimage as ndimage

    titsize = 9
    mycmap = "bone"
    nticks = 5
    r_crit = 0.001  # 0.001
    outVS0 = output[0][:, 0]
    outPH0 = output[0][:, 1]
    outRC0 = output[0][:, 2]
    outAP0 = output[0][:, 3]
    outVSD = output[0][:, 4]
    outPHD = output[0][:, 5]
    outRCD = output[0][:, 6]
    outAPD = output[0][:, 7]

    outVS0[
        np.logical_or(outRC0 > r_crit, outRCD > r_crit)
    ] = 0.0  # set VS of conditions to zero that fail the rayleigh test
    outPH0[
        np.logical_or(outRC0 > r_crit, outRCD > r_crit)
    ] = 0.0  # set phi of conditions to zero that fail the rayleigh test
    outVSD[np.logical_or(outRC0 > r_crit, outRCD > r_crit)] = 0.0
    outPHD[np.logical_or(outRC0 > r_crit, outRCD > r_crit)] = 0.0

    pshift = seh.get_anglediff(outPH0, outPHD)
    pshift = np.reshape(pshift, (P["N"], P["N"]))
    APdiff = outAPD - outAP0
    APdiff = np.reshape(APdiff, (P["N"], P["N"]))
    VSdiff = outVSD - outVS0
    VSdiff = np.reshape(VSdiff, (P["N"], P["N"]))

    filtsig = 0.5
    APdiff = ndimage.gaussian_filter(APdiff, sigma=filtsig, order=0)
    VSdiff = ndimage.gaussian_filter(VSdiff, sigma=filtsig, order=0)
    #
    # debug fig
    fdebug = plt.figure()
    plt.subplot(221)
    c1 = plt.contourf(
        x,
        y,
        np.reshape(outAP0, (P["N"], P["N"])),
        21,
    )
    plt.contour(x, y, np.reshape(outAP0, (P["N"], P["N"])), levels=(50,), colors="w")
    plt.colorbar(c1, use_gridspec=True)
    plt.title("AP0")
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    plt.subplot(222)
    c2 = plt.contourf(
        x,
        y,
        np.reshape(outVS0, (P["N"], P["N"])),
        21,
    )
    plt.contour(x, y, np.reshape(outVS0, (P["N"], P["N"])), levels=(0.01,), colors="w")
    plt.colorbar(c2, use_gridspec=True)
    plt.title("VS0")
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    plt.subplot(223)
    c3 = plt.contourf(
        x,
        y,
        np.reshape(outAPD, (P["N"], P["N"])),
        21,
    )
    plt.contour(x, y, np.reshape(outAPD, (P["N"], P["N"])), levels=(50,), colors="w")
    plt.colorbar(c3, use_gridspec=True)
    plt.title("APD")
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    plt.subplot(224)
    c4 = plt.contourf(
        x,
        y,
        np.reshape(outVSD, (P["N"], P["N"])),
        21,
    )
    plt.contour(x, y, np.reshape(outVSD, (P["N"], P["N"])), levels=(0.01,), colors="w")
    plt.colorbar(c4, use_gridspec=True)
    plt.title("VSD")
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    plt.tight_layout()
    #
    # Publication Figure -- Differences Only
    fwidth = 5.5  # cm
    fhandle = plt.figure(figsize=(fwidth / 2.54, (fwidth * 2.0) / 2.54), dpi=600)
    #
    ax1 = plt.subplot(311)
    CS1 = plt.contourf(
        x, y, APdiff, np.linspace(-30, 30, 21), extend="both", cmap="seismic"
    )  # repeated = y, tiled = x!!
    plt.contour(
        x,
        y,
        np.reshape(outAP0, (P["N"], P["N"])),
        levels=(50.0,),
        colors="k",
        linewidths=(0.5,),
    )
    plt.title("AP Rate change")
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    ax1.yaxis.set_major_locator(MaxNLocator(nticks))
    ax1.xaxis.set_major_locator(MaxNLocator(nticks))
    cbar1 = plt.colorbar(CS1, use_gridspec=True)
    cbar1.ax.set_ylabel(r"$\Delta$ Rate (Hz)")
    tl = MaxNLocator(nbins=5)
    cbar1.locator = tl
    cbar1.update_ticks()
    #
    ax2 = plt.subplot(312)
    CS2 = plt.contourf(
        x, y, VSdiff, np.linspace(-0.1, 0.1, 21), extend="both", cmap="seismic"
    )
    plt.contour(
        x,
        y,
        np.reshape(outAP0, (P["N"], P["N"])),
        levels=(50.0,),
        colors="k",
        linestyles=("dashed"),
        linewidths=(0.5,),
    )
    plt.contour(
        x,
        y,
        np.reshape(outVS0, (P["N"], P["N"])),
        levels=(0.01,),
        colors="k",
        linewidths=(0.5,),
    )
    plt.title(r"Vectorstrength change")
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    ax2.yaxis.set_major_locator(MaxNLocator(nticks))
    ax2.xaxis.set_major_locator(MaxNLocator(nticks))
    cbar2 = plt.colorbar(CS2, use_gridspec=True)
    cbar2.ax.set_ylabel(r"$\Delta$ Vectorstrength")
    tl = MaxNLocator(nbins=5)
    cbar2.locator = tl
    cbar2.update_ticks()
    ax2.plot(5.35, 291.0, "wx", markersize=3)  # position of example traces
    ax2.plot(3.0, 50.0, "wx", markersize=3)  # position of example traces
    #
    ax3 = plt.subplot(313)
    CS3 = plt.contourf(
        x, y, pshift, np.linspace(-0.1, 0.1, 21), extend="both", cmap="seismic"
    )  # repeated = y, tiled = x!!
    plt.title("Phase shift")
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    ax3.yaxis.set_major_locator(MaxNLocator(nticks))
    ax3.xaxis.set_major_locator(MaxNLocator(nticks))
    cbar3 = plt.colorbar(CS3, use_gridspec=True)
    cbar3.ax.set_ylabel(r"$\Delta \varphi$ (cycles)")
    tl = MaxNLocator(nbins=5)
    cbar3.locator = tl
    cbar3.update_ticks()
    plt.tight_layout()
    #
    # output
    cf = P["cf"][0]
    freq = P["freq"][0]
    db = P["dB"][0]
    conds = P["N"]
    with open("./results/Figure8.txt", "a") as resfile:
        resfile.write("#################################################\n")
        resfile.write("CF: " + str(cf) + "Hz\n")
        resfile.write("Frequency: " + str(freq) + "Hz\n")
        resfile.write("Level: " + str(db) + "dB SPL\n")
        resfile.write("Resolution: " + str(conds) + "x" + str(conds) + " conditions\n")
        val1 = np.max(APdiff, axis=0)
        arg1 = np.argmax(APdiff, axis=0)
        val2 = np.max(val1)
        arg1 = arg1[np.argmax(val1)]
        arg2 = np.argmax(val1)
        resfile.write(
            "MAX AP: L="
            + str(y[arg1])
            + ", D="
            + str(x[arg2])
            + " => dAP="
            + str(val2)
            + " spikes\n"
        )
        val1 = np.min(APdiff, axis=0)
        arg1 = np.argmin(APdiff, axis=0)
        val2 = np.min(val1)
        arg1 = arg1[np.argmin(val1)]
        arg2 = np.argmin(val1)
        resfile.write(
            "MIN AP: L="
            + str(y[arg1])
            + ", D="
            + str(x[arg2])
            + " => dAP="
            + str(val2)
            + " spikes\n"
        )
        val1 = np.max(VSdiff, axis=0)
        arg1 = np.argmax(VSdiff, axis=0)
        val2 = np.max(val1)
        arg1 = arg1[np.argmax(val1)]
        arg2 = np.argmax(val1)
        resfile.write(
            "MAX VS: L="
            + str(y[arg1])
            + ", D="
            + str(x[arg2])
            + " => dVS="
            + str(val2)
            + "\n"
        )
        #
        val1 = np.min(VSdiff, axis=0)
        arg1 = np.argmin(VSdiff, axis=0)
        val2 = np.min(val1)
        arg1 = arg1[np.argmin(val1)]
        arg2 = np.argmin(val1)
        resfile.write(
            "MIN VS: L="
            + str(y[arg1])
            + ", D="
            + str(x[arg2])
            + " => dVS="
            + str(val2)
            + "\n"
        )
        resfile.write(".\n")
    return (fhandle, fdebug)


def mkexampletraces(minL, minD, maxL, maxD):
    """
    As per reviewer request we generate example traces for extreme conditions
    """
    dur = 5000.0
    dur_s = dur / 1000.0
    dt = 0.01
    freq = 1389.0
    cf = 1501.0
    dB = 65.0
    zseed = 45453
    eohseed = 34823
    stocheoh = True
    geoh = 0.055
    N = 32
    gsyn = 0.016
    anf_num = (N + 1, 0, 0)
    SO = seh.create_sound(fs=100e3, freq=freq, duration=dur_s, dbspl=dB)
    SP = seh.create_spk(SO, fs=100e3, N=N + 1, cf=cf, seed=zseed, anf_num=anf_num)
    S = np.array(SP["spikes"] * 1000.0)
    print(".")
    RX1 = sem.SE_BS(
        S,
        Simdur=dur,
        dt=dt,
        G_EoH=geoh,
        StochasticEoH=stocheoh,
        N=N,
        G_Dend=0.0,
        L=minL,
        D=minD,
        constR=False,
        EoHseed=eohseed,
    )
    print(".")
    RX2 = sem.SE_BS(
        S,
        Simdur=dur,
        dt=dt,
        G_EoH=geoh,
        StochasticEoH=stocheoh,
        N=N,
        G_Dend=gsyn,
        L=minL,
        D=minD,
        constR=False,
        EoHseed=eohseed,
    )
    print(".")
    RX3 = sem.SE_BS(
        S,
        Simdur=dur,
        dt=dt,
        G_EoH=geoh,
        StochasticEoH=stocheoh,
        N=N,
        G_Dend=0.0,
        L=maxL,
        D=maxD,
        constR=False,
        EoHseed=eohseed,
    )
    print(".")
    RX4 = sem.SE_BS(
        S,
        Simdur=dur,
        dt=dt,
        G_EoH=geoh,
        StochasticEoH=stocheoh,
        N=N,
        G_Dend=gsyn,
        L=maxL,
        D=maxD,
        constR=False,
        EoHseed=eohseed,
    )
    Ev1 = seh.SimpleDetectAP(RX1["Vm"], thr=-100, dt=dt, LM=-20, RM=10)
    Ev2 = seh.SimpleDetectAP(RX2["Vm"], thr=-100, dt=dt, LM=-20, RM=10)
    Ev3 = seh.SimpleDetectAP(RX3["Vm"], thr=-100, dt=dt, LM=-20, RM=10)
    Ev4 = seh.SimpleDetectAP(RX4["Vm"], thr=-100, dt=dt, LM=-20, RM=10)
    APr1 = np.round(len(Ev1["PeakT"]) / dur_s)
    APr2 = np.round(len(Ev2["PeakT"]) / dur_s)
    APr3 = np.round(len(Ev3["PeakT"]) / dur_s)
    APr4 = np.round(len(Ev4["PeakT"]) / dur_s)
    VS1, phi1, Ray1, phases1 = seh.vectorstrength(Ev1["PeakT"], freq, [0.0, dur])
    VS2, phi2, Ray2, phases2 = seh.vectorstrength(Ev2["PeakT"], freq, [0.0, dur])
    VS3, phi3, Ray3, phases3 = seh.vectorstrength(Ev3["PeakT"], freq, [0.0, dur])
    VS4, phi4, Ray4, phases4 = seh.vectorstrength(Ev4["PeakT"], freq, [0.0, dur])
    tx = np.linspace(0, dur - dt, int(dur / dt))
    fhand = plt.figure()
    ax1 = plt.subplot(211)
    ax1.set_title(str(APr1) + ", " + str(VS1) + " / " + str(APr2) + ", " + str(VS2))
    ax1.plot(tx, RX1["Vm"], "k")
    ax1.plot(Ev1["PeakT"], Ev1["PeakV"], color="g", marker=".", linestyle="None")
    ax1.plot(tx, RX2["Vm"], "r")
    ax1.plot(Ev2["PeakT"], Ev2["PeakV"], color="r", marker=".", linestyle="None")
    ax1.set_xlim((300.0, 500.0))
    ax2 = plt.subplot(212)
    ax2.set_title(str(APr3) + ", " + str(VS3) + " / " + str(APr4) + ", " + str(VS4))
    ax2.plot(tx, RX3["Vm"], color="k", linewidth=0.5)
    ax2.plot(Ev3["PeakT"], Ev3["PeakV"], color="g", marker=".", linestyle="None")
    ax2.plot(tx, RX4["Vm"], color="r", linewidth=0.5)
    ax2.plot(Ev4["PeakT"], Ev4["PeakV"], color="r", marker=".", linestyle="None")
    ax2.set_xlim((300.0, 500.0))
    plt.tight_layout()
    return fhand


if __name__ == "__main__":
    # Checking for arguments and/or setting defaults
    addedargs = sys.argv[
        1:
    ]  # fname, weload, nconds, ncores, freq, cf, level, mkexamples
    myargs = [1, 5, 4, 1389.0, 1501.0, 65.0, 1]  # as above starting with weload
    for iarg, thisarg in enumerate(addedargs):
        myargs[iarg] = float(thisarg)
    weload = bool(myargs[0])
    nconds = int(myargs[1])
    ncores = int(myargs[2])
    freq = myargs[3]
    cf = myargs[4]
    db = myargs[5]
    makeexamples = bool(myargs[6])
    fsuffix = str(int(freq)) + "_" + str(int(cf)) + "_" + str(int(db))
    ####LOAD DATA FOR Figure8a (IF IT EXISTS)
    if weload and os.path.isfile("./results/Figure8a_" + fsuffix + ".npy"):
        print("Data for Figure8a found... loading!")
        output = np.load("./results/Figure8a_" + fsuffix + ".npy", allow_pickle=True)
        P = np.load("./results/Figure8a_P_" + fsuffix + ".npy", allow_pickle=True)
        P = P.tolist()
    else:
        conditions = nconds  # 32
        cores = ncores  # 20
        dur = 5000.0
        dt = 0.01
        output = []
        # Model Parameters (all in a linearly aranged fashion, so that a minimal
        # amount of programming is required to change the experiment).
        P = {}
        P["N"] = conditions
        P["cores"] = cores
        P["TotalN"] = int(P["N"] ** 2)
        P["Number"] = list(range(P["TotalN"]))
        P["mp"] = True
        P["Zseed"] = 45453
        P["AdvZseed"] = False
        P["EoHseed"] = 34823
        P["AdvEoHseed"] = False
        #########################################
        P["dur"] = np.repeat(dur, P["TotalN"])
        P["dt"] = np.repeat(dt, P["TotalN"])
        P["G_EoH"] = np.repeat(0.055, P["TotalN"])
        P["Stochastic_EoH"] = np.repeat(True, P["TotalN"])
        P["freq"] = np.repeat(freq, P["TotalN"])
        P["cf"] = np.repeat(cf, P["TotalN"])
        P["dB"] = np.repeat(db, P["TotalN"])
        P["Nsyn"] = np.repeat(32, P["TotalN"])
        P["gsyn"] = np.repeat(0.016, P["TotalN"])
        P["constR"] = np.repeat(False, P["TotalN"])
        P["gLEAK_d"] = np.repeat(0.001, P["TotalN"])
        P["gKLT_d"] = np.repeat(0.0085, P["TotalN"])
        P["gIH_d"] = np.repeat(0.001, P["TotalN"])
        P["gKLT_s"] = np.repeat(0.017, P["TotalN"])
        P["gIH_s"] = np.repeat(0.002, P["TotalN"])

        # Now define the two variable parameters. The repeated = y, the tiled = x!!
        allL = np.linspace(1.0, 501.0, P["N"])
        allD = np.linspace(1.0, 10.0, P["N"])
        P["L"] = np.repeat(allL, P["N"])
        P["D"] = np.tile(allD, P["N"])

        # make go!
        output.append(myMPhandler(P))
        output = np.array(output)
        np.save("./results/Figure8a_" + fsuffix + ".npy", output, allow_pickle=True)
        np.save("./results/Figure8a_P_" + fsuffix + ".npy", P, allow_pickle=True)

    ###Plot results
    fhandles = plotres(
        output=output,
        P=P,
        x=np.unique(P["D"]),
        y=np.unique(P["L"]),
        xlabs=r"Diameter $(\mu m)$",
        ylabs="Length $(\mu m)$",
    )
    with PdfPages("./figs/Figure8a_" + fsuffix + ".pdf") as pp:
        pp.savefig(fhandles[0])
        pp.savefig(fhandles[1])
    if makeexamples:
        fhandle2 = mkexampletraces(291.0, 5.35, 50.0, 3.0)
        pp = PdfPages("./figs/Figure8a_" + fsuffix + "_Examples.pdf")
        pp.savefig()
        pp.close()
