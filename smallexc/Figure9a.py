#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This file contains all code necessary to generate Figure 9, part a:
    freq vs. cf on spiking

Created on 28 Aug 2019
Modified sept 2020
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
    from scipy import stats

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
        x,
        y,
        APdiff,
        np.linspace(-30, 30, 21),
        cmap="seismic",
        extend="both",
    )  # repeated = y, tiled = x!!
    ax1.set_title("AP Rate change")
    ax1.set_xlabel(xlabs)
    ax1.set_ylabel(ylabs)
    ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_xticks((0.2, 0.5, 1.0, 2.0))
    ax1.set_yticks((0.2, 0.5, 1.0, 2.0))
    ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    #
    # Linear regression: Reviewer request
    yl = ax1.get_ylim()
    xl = ax1.get_xlim()
    ax1.plot([0.125, 2.5], [0.125, 2.5], "k--", linewidth=0.5)
    ax1.plot(x, y[np.argmax(APdiff, 0)], "k.", markersize=2)
    s, intc, r, p, std = stats.linregress(x, y[np.argmax(APdiff, 0)])
    ax1.plot(x, s * x + intc, "y-", linewidth=0.5)
    ax1.set_ylim(yl)
    ax1.set_xlim(xl)
    outtext = (
        "bestAPdiff regression, y="
        + str(round(s, 3))
        + "*x+"
        + str(round(intc, 3))
        + ", r="
        + str(round(r, 3))
        + ", p="
        + str(round(p, 3))
    )
    print(outtext)
    #
    cbar1 = plt.colorbar(CS1, use_gridspec=True)
    cbar1.ax.set_ylabel(r"$\Delta$ Rate (Hz)")
    tl = MaxNLocator(nbins=5)
    cbar1.locator = tl
    cbar1.update_ticks()
    #
    #
    ax2 = plt.subplot(312)
    CS2 = plt.contourf(
        x,
        y,
        VSdiff,
        np.linspace(-0.1, 0.1, 21),
        cmap="seismic",
        extend="both",
    )
    ax2.set_title(r"Vectorstrength change")
    ax2.set_xlabel(xlabs)
    ax2.set_ylabel(ylabs)
    ax2.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_xticks((0.2, 0.5, 1.0, 2.0))
    ax2.set_yticks((0.2, 0.5, 1.0, 2.0))
    ax2.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax2.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    #
    # Linear regression: Reviewer request
    yl = ax2.get_ylim()
    xl = ax2.get_xlim()
    ax2.plot([0.125, 2.5], [0.125, 2.5], "k--", linewidth=0.5)
    ax2.plot(x, y[np.argmax(VSdiff, 0)], "k.", markersize=2)
    s, intc, r, p, std = stats.linregress(x, y[np.argmax(VSdiff, 0)])
    ax2.plot(x, s * x + intc, "y-", linewidth=0.5)
    ax2.set_ylim(yl)
    ax2.set_xlim(xl)
    outtext = (
        "best VSdiff regression, y="
        + str(round(s, 3))
        + "*x+"
        + str(round(intc, 3))
        + ", r="
        + str(round(r, 3))
        + ", p="
        + str(round(p, 3))
    )
    print(outtext)
    cbar2 = plt.colorbar(CS2, use_gridspec=True)
    cbar2.ax.set_ylabel(r"$\Delta$ Vectorstrength")
    tl = MaxNLocator(nbins=5)
    cbar2.locator = tl
    cbar2.update_ticks()
    #
    #
    ax3 = plt.subplot(313)
    CS3 = plt.contourf(
        x,
        y,
        pshift,
        np.linspace(-0.1, 0.1, 21),
        cmap="seismic",
        extend="both",
    )  # repeated = y, tiled = x!!
    ax3.set_title("Phase shift")
    ax3.set_xlabel(xlabs)
    ax3.set_ylabel(ylabs)
    ax3.set_yscale("log")
    ax3.set_xscale("log")
    ax3.set_xticks((0.2, 0.5, 1.0, 2.0))
    ax3.set_yticks((0.2, 0.5, 1.0, 2.0))
    ax3.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax3.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    #
    # Linear regression: Reviewer request
    yl = ax3.get_ylim()
    xl = ax3.get_xlim()
    ax3.plot([0.0, 2.5], [0.0, 2.5], "k--", linewidth=0.5)
    ax3.plot(x, y[np.argmax(np.abs(VSdiff), 0)], "k.", markersize=2)
    s, intc, r, p, std = stats.linregress(x, y[np.argmax(np.abs(VSdiff), 0)])
    ax3.plot(x, s * x + intc, "y-", linewidth=0.5)
    ax3.set_ylim(yl)
    ax3.set_xlim(xl)
    outtext = (
        "best phidiff regression, y="
        + str(round(s, 3))
        + "*x+"
        + str(round(intc, 3))
        + ", r="
        + str(round(r, 3))
        + ", p="
        + str(round(p, 3))
    )
    print(outtext)
    #
    cbar3 = plt.colorbar(CS3, use_gridspec=True)
    cbar3.ax.set_ylabel(r"$\Delta \varphi$ (cycles)")
    tl = MaxNLocator(nbins=5)
    cbar3.locator = tl
    cbar3.update_ticks()

    plt.tight_layout()
    return fhandle


if __name__ == "__main__":
    ####LOAD DATA FOR Figure9a (IF IT EXISTS)
    if os.path.isfile("./results/Figure9a.npy"):
        print("Data for Figure9a found... loading!")
        output = np.load("./results/Figure9a.npy", allow_pickle=True)
        P = np.load("./results/Figure9a_P.npy", allow_pickle=True)
        P = P.tolist()
    else:
        conditions = 25
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
        P["AdvZseed"] = True
        P["EoHseed"] = 34823
        P["AdvEoHseed"] = True
        #########################################
        P["dur"] = np.repeat(dur, P["TotalN"])
        P["dt"] = np.repeat(dt, P["TotalN"])
        P["G_EoH"] = np.repeat(0.055, P["TotalN"])
        P["Stochastic_EoH"] = np.repeat(True, P["TotalN"])
        P["dB"] = np.repeat(65.0, P["TotalN"])
        P["Nsyn"] = np.repeat(int(32), P["TotalN"])
        P["gsyn"] = np.repeat(0.016, P["TotalN"])
        P["tau_dendsyn"] = np.repeat(2.0, P["TotalN"])
        P["L"] = np.repeat(50.0, P["TotalN"])
        P["D"] = np.repeat(3.0, P["TotalN"])
        P["constR"] = np.repeat(False, P["TotalN"])
        P["gLEAK_d"] = np.repeat(0.001, P["TotalN"])
        P["gKLT_d"] = np.repeat(0.0085, P["TotalN"])
        P["gKLT_s"] = np.repeat(0.017, P["TotalN"])
        P["gIH_s"] = np.repeat(0.002, P["TotalN"])
        P["gIH_d"] = np.repeat(0.001, P["TotalN"])

        # Now define the two variable parameters. The repeated = y, the tiled = x!!
        allcf = np.round(np.geomspace(125.1, 2501.0, P["N"]))
        allfreq = np.round(np.geomspace(125.1, 2501.0, P["N"]))
        P["freq"] = np.repeat(allfreq, P["N"])
        P["cf"] = np.tile(allcf, P["N"])

        # make go!
        output.append(myMPhandler(P))
        output = np.array(output)
        np.save("./results/Figure9a.npy", output, allow_pickle=True)
        np.save("./results/Figure9a_P.npy", P, allow_pickle=True)

    ###Plot results
    fhandle = plotres(
        output=output,
        P=P,
        x=np.unique(P["cf"]) / 1000.0,
        y=np.unique(P["freq"]) / 1000.0,
        xlabs=r"$CF (kHz)$",
        ylabs=r"$Frequency (kHz)$",
    )
    pp = PdfPages("./figs/Figure9a.pdf")
    pp.savefig()
    pp.close()
