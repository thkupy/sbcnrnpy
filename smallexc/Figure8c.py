#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 28 Aug 2019

revised sept 2020: upon reviewer request you can also "balance" gklt with gh with#
respect to the resonance frequency of the dendrite.

This file contains all code necessary to generate Figure 8, part c:
    KLT vs. freq/cf on spiking
@author: Thomas Kuenzel (kuenzel{at}bio2.rwth-aachen.de)
"""

import os
import sys
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


def lookup_balancing_gh(gklt, luwhat):
    try:
        LUT = np.load("./results/fres_comp_klt_h.npy")
    except:
        msg = "Lookup table not found. Run smallexc/impedance_testing_balancing.py"
        print(msg)
        return -1
    LUT = LUT.tolist()
    if luwhat == 1:  # keep fres constant
        myidx = np.argmin(np.abs(LUT["klt"] - gklt))
        return LUT["h"][myidx]
    else:  # keep zmax constant
        myidx = np.argmin(np.abs(LUT["klt2"] - gklt))
        return LUT["h2"][myidx]


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
    ax1.xaxis.set_major_locator(MaxNLocator(nticks))
    cbar1 = plt.colorbar(CS1, use_gridspec=True)
    cbar1.ax.set_ylabel(r"$\Delta$ Rate (Hz)")
    tl = MaxNLocator(nbins=5)
    cbar1.locator = tl
    cbar1.update_ticks()

    ax2 = plt.subplot(312)
    CS2 = plt.contourf(
        x, y, VSdiff, np.linspace(-0.1, 0.1, 21), extend="both", cmap="seismic"
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
    return fhandle


if __name__ == "__main__":
    # Checking for arguments and/or setting defaults
    addedargs = sys.argv[1:]
    if len(addedargs) == 0:
        compensate_klt = 0
    else:
        compensate_klt = int(addedargs[0])
    if compensate_klt == 1:
        filesuffix = "_kltcomp_f.npy"
        figsuffix = "_kltcomp_f.pdf"
    elif compensate_klt == 2:
        filesuffix = "_kltcomp_z.npy"
        figsuffix = "_kltcomp_z.pdf"
    else:
        filesuffix = ".npy"
        figsuffix = ".pdf"

    ####LOAD DATA FOR Figure8c (IF IT EXISTS)
    if os.path.isfile("./results/Figure8c" + filesuffix):
        print("Data for Figure8c found... loading!")
        output = np.load("./results/Figure8c" + filesuffix, allow_pickle=True)
        P = np.load("./results/Figure8c_P" + filesuffix, allow_pickle=True)
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
        P["dB"] = np.repeat(65.0, P["TotalN"])
        P["Nsyn"] = np.repeat(int(32), P["TotalN"])
        P["gsyn"] = np.repeat(0.016, P["TotalN"])
        P["tau_dendsyn"] = np.repeat(2.0, P["TotalN"])
        P["L"] = np.repeat(50.0, P["TotalN"])
        P["D"] = np.repeat(3.0, P["TotalN"])
        P["constR"] = np.repeat(False, P["TotalN"])
        P["gLEAK_d"] = np.repeat(0.001, P["TotalN"])
        P["gKLT_s"] = np.repeat(0.017, P["TotalN"])
        P["gIH_s"] = np.repeat(0.002, P["TotalN"])

        # Now define the two variable parameters. The repeated = y, the tiled = x!!
        allklt = np.linspace(0.00001, 0.02, P["N"])
        allcf = np.round(np.geomspace(125.1, 2501.0, P["N"]))
        allfreq = np.round(allcf * 0.875)
        P["gKLT_d"] = np.repeat(allklt, P["N"])
        P["freq"] = np.tile(allfreq, P["N"])
        P["cf"] = np.tile(allcf, P["N"])
        if compensate_klt > 0.5:
            allh = np.array(())
            for thisk in allklt:
                allh = np.append(allh, lookup_balancing_gh(thisk, compensate_klt))
            P["gIH_d"] = np.repeat(allh, P["N"])
        else:
            P["gIH_d"] = np.repeat(0.001, P["TotalN"])
        # make go!
        output.append(myMPhandler(P))
        output = np.array(output)
        np.save("./results/Figure8c" + filesuffix, output, allow_pickle=True)
        np.save("./results/Figure8c_P" + filesuffix, P, allow_pickle=True)

    ###Plot results
    fhandle = plotres(
        output=output,
        P=P,
        x=np.unique(P["cf"]) / 1000.0,
        y=np.unique(P["gKLT_d"]) * 1000.0,
        xlabs=r"$CF (kHz)$",
        ylabs=r"$g_{KLTdend} (nS)$",
    )
    pp = PdfPages("./figs/Figure8c" + figsuffix)
    pp.savefig()
    pp.close()
