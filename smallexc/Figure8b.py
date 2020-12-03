#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 28 Aug 2019
This file contains all code necessary to generate Figure 8, part b:
    Nsyn vs. gsyn on spiking
@author: Thomas Kuenzel (kuenzel{at}bio2.rwth-aachen.de)
"""

import os
import matplotlib

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

    fwidth = 5.5  # cm
    fhandle = plt.figure(figsize=(fwidth / 2.54, (fwidth * 2.0) / 2.54), dpi=600)

    ax1 = plt.subplot(311)
    CS1 = plt.contourf(
        x, y, APdiff, np.linspace(-30, 30, 21), extend="both", cmap="seismic"
    )  # repeated = y, tiled = x!!
    plt.title("AP Rate change")
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    ax1.yaxis.set_major_locator(MaxNLocator(nticks))
    ax1.xaxis.set_major_locator(MaxNLocator(nticks - 1))
    cbar1 = plt.colorbar(CS1, use_gridspec=True)
    cbar1.ax.set_ylabel(r"$\Delta$ Rate (Hz)")
    tl = MaxNLocator(nbins=5)
    cbar1.locator = tl
    cbar1.update_ticks()
    ax1.plot(0.016, 32.0, "wx", markersize=3)  # position of example traces
    ax1.plot(0.048, 32.0, "wx", markersize=3)  # position of example traces
    #
    ax2 = plt.subplot(312)
    CS2 = plt.contourf(
        x, y, VSdiff, np.linspace(-0.1, 0.1, 21), extend="both", cmap="seismic"
    )
    plt.title(r"Vectorstrength change")
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    ax2.yaxis.set_major_locator(MaxNLocator(nticks))
    ax2.xaxis.set_major_locator(MaxNLocator(nticks - 1))
    cbar2 = plt.colorbar(CS2, use_gridspec=True)
    cbar2.ax.set_ylabel(r"$\Delta$ Vectorstrength")
    tl = MaxNLocator(nbins=5)
    cbar2.locator = tl
    cbar2.update_ticks()

    ax3 = plt.subplot(313)
    CS3 = plt.contourf(
        x, y, pshift, np.linspace(-0.1, 0.1, 21), extend="both", cmap="seismic"
    )  # repeated = y, tiled = x!!
    plt.title("Phase shift")
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    ax3.yaxis.set_major_locator(MaxNLocator(nticks))
    ax3.xaxis.set_major_locator(MaxNLocator(nticks - 1))
    cbar3 = plt.colorbar(CS3, use_gridspec=True)
    cbar3.ax.set_ylabel(r"$\Delta \varphi$ (cycles)")
    tl = MaxNLocator(nbins=5)
    cbar3.locator = tl
    cbar3.update_ticks()
    plt.tight_layout()
    return fhandle


def mkexampletraces():
    """
    As per reviewer request we generate example traces for two extreme conditions,
    one at N=32 and g=0.06 (negative effect on rate) and one at N=32 and g=0.03
    (positive effect on rate).
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
    geoh = 0.55
    L = 50.0
    D = 3.0
    N = 32
    anf_num = (N + 1, 0, 0)
    kltd = 0.0085
    hd = 0.001
    ld = 0.001
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
        G_Dend=0.016,
        L=L,
        D=D,
        constR=False,
        gKLT_d=kltd,
        gIH_d=hd,
        gLEAK_d=ld,
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
        G_Dend=0.048,
        L=L,
        D=D,
        constR=False,
        gKLT_d=kltd,
        gIH_d=hd,
        gLEAK_d=ld,
        EoHseed=eohseed,
    )
    Ev1 = seh.SimpleDetectAP(RX1["Vm"], thr=-100, dt=dt, LM=-20, RM=10)
    Ev2 = seh.SimpleDetectAP(RX2["Vm"], thr=-100, dt=dt, LM=-20, RM=10)
    APr1 = len(Ev1["PeakT"]) / dur_s
    APr2 = len(Ev2["PeakT"]) / dur_s
    VS1, phi1, Ray1, phases1 = seh.vectorstrength(Ev1["PeakT"], freq, [0.0, dur])
    VS2, phi2, Ray2, phases2 = seh.vectorstrength(Ev2["PeakT"], freq, [0.0, dur])
    tx = np.linspace(0, dur - dt, int(dur / dt))
    fhand = plt.figure()
    ax1 = plt.subplot(211)
    ax1.set_title("Rate=" + str(APr1) + ", VS=" + str(VS1))
    ax2 = plt.subplot(212)
    ax2.set_title("Rate=" + str(APr2) + ", VS=" + str(VS2))
    ax1.plot(tx, RX1["Vm"], "k")
    ax1.set_xlim((800.0, 1000.0))
    ax2.plot(tx, RX2["Vm"], color="r", linewidth=0.5)
    ax2.set_xlim((800.0, 1000.0))
    ax1.plot(Ev1["PeakT"], Ev1["PeakV"], color="g", marker=".", linestyle="None")
    ax2.plot(Ev2["PeakT"], Ev2["PeakV"], color="r", marker=".", linestyle="None")
    plt.tight_layout()
    return fhand


if __name__ == "__main__":
    ####LOAD DATA FOR Figure8b (IF IT EXISTS)
    if os.path.isfile("./results/Figure8b.npy"):
        print("Data for Figure8b found... loading!")
        output = np.load("./results/Figure8b.npy", allow_pickle=True)
        P = np.load("./results/Figure8b_P.npy", allow_pickle=True)
        P = P.tolist()
    else:
        conditions = 32
        cores = 20
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
        P["freq"] = np.repeat(1389.0, P["TotalN"])
        P["cf"] = np.repeat(1501.0, P["TotalN"])
        P["dB"] = np.repeat(65.0, P["TotalN"])
        P["tau_dendsyn"] = np.repeat(2.0, P["TotalN"])
        P["L"] = np.repeat(50.0, P["TotalN"])
        P["D"] = np.repeat(3.0, P["TotalN"])
        P["constR"] = np.repeat(False, P["TotalN"])
        P["IPI"] = np.repeat(1000.0 / P["freq"], P["TotalN"])
        P["gLEAK_d"] = np.repeat(0.001, P["TotalN"])
        P["gKLT_s"] = np.repeat(0.017, P["TotalN"])
        P["gIH_s"] = np.repeat(0.002, P["TotalN"])
        P["gKLT_d"] = np.repeat(0.0085, P["TotalN"])
        P["gIH_d"] = np.repeat(0.001, P["TotalN"])

        # Now define the two variable parameters. The repeated = y, the tiled = x!!
        allg = np.linspace(0.001, 0.064, P["N"])
        alln = np.linspace(1, 128, P["N"], dtype=int)
        P["Nsyn"] = np.repeat(alln, P["N"])
        P["gsyn"] = np.tile(allg, P["N"])

        # make go!
        output.append(myMPhandler(P))
        output = np.array(output)
        np.save("./results/Figure8b.npy", output, allow_pickle=True)
        np.save("./results/Figure8b_P.npy", P, allow_pickle=True)

    ###Plot results
    fhandle = plotres(
        output=output,
        P=P,
        x=np.unique(P["gsyn"]),
        y=np.unique(P["Nsyn"]),
        xlabs=r"$g_{syn} (\mu S)$",
        ylabs="$N_{syn}$",
    )
    pp = PdfPages("./figs/Figure8b.pdf")
    pp.savefig()
    pp.close()
    fhandle2 = mkexampletraces()
    pp = PdfPages("./figs/Figure8b_Examples.pdf")
    pp.savefig()
    pp.close()
