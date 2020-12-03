#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This file contains all code necessary to generate Figure 6 panel G5/G6. This is a
reviewer request from the JNP submission: "show the experiment in Figure 6G1-4 with
different number of dendritic total synaptic conductance".
Here, profiles at 0.75 cycles will be calculated for a "conditions" range of different gsyntotals and endbulb g.

To test this rapidly, change dt (line 183) to coarser values and the and number of
conditions (line 180) to lower (5x5 works well).
Adjust number of available cores in line 181.

Created sept 2020
@author: kuenzel{at}bio2.rwth-aachen.de
"""

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import semodels as sem
import sehelper as seh
import multiprocessing
import os
from functools import partial
from math import ceil


def runonce(x, P):
    ####GENERATE SPIKE INPUTS (the spiketimes are always the same)
    # generate dendrite spiketimes (always the same)
    thisdur = (P["dur"][x] - 5.0) / 1000.0
    SO = seh.create_sound(
        fs=100e3, freq=P["freq"][x], duration=thisdur, dbspl=P["dB"][x]
    )
    SP = seh.create_spk(
        SO, fs=100e3, N=P["Nsyn"][x], cf=P["cf"][x], seed=3433833
    )  # generate enough spiketrains
    cycledur = 1.0 / P["freq"][x]
    shift = cycledur * 100.0  # we wait 100 cycle durations to start with the EoH
    currlen = len(SP)
    # SP.loc[currlen] = [
    #    SP.loc[currlen - 1][0],
    #    SP.loc[currlen - 1][1],
    #    np.array((shift + P["EoH_Time"][x],)),
    #    SP.loc[currlen - 1][3],
    # ]
    SP.loc[currlen] = [
        np.array((shift + P["EoH_Time"][x],)),
        SP.loc[currlen - 1][1],
        SP.loc[currlen - 1][2],
        SP.loc[currlen - 1][3],
    ]
    S = np.array(SP["spikes"] * 1000.0)
    R0 = sem.SE_BS(
        S,
        Simdur=P["dur"][x],
        dt=P["dt"][x],
        G_EoH=P["G_EoH"][x],
        StochasticEoH=P["Stochastic_EoH"][x],
        N=P["Nsyn"][x],
        G_Dend=0.0,
        tau_Dend=P["tau_dendsyn"][x],
        L=P["L"][x],
        D=P["D"][x],
        constR=P["constR"][x],
        gKLT_d=P["gKLT_d"][x],
        gIH_d=P["gIH_d"][x],
        gLEAK_d=P["gLEAK_d"][x],
        EoHseed=65232,
    )
    R = sem.SE_BS(
        S,
        Simdur=P["dur"][x],
        dt=P["dt"][x],
        G_EoH=P["G_EoH"][x],
        StochasticEoH=P["Stochastic_EoH"][x],
        N=P["Nsyn"][x],
        G_Dend=P["gsyn"][x],
        tau_Dend=P["tau_dendsyn"][x],
        L=P["L"][x],
        D=P["D"][x],
        constR=P["constR"][x],
        gKLT_d=P["gKLT_d"][x],
        gIH_d=P["gIH_d"][x],
        gLEAK_d=P["gLEAK_d"][x],
        EoHseed=65232,
    )
    print((str(P["G_EoH"][x]) + " " + str(P["gsyn"][x]) + " " + str(x)))
    return [R0["Vm"], R["Vm"]]


def myMPhandler(P):
    p = multiprocessing.Pool(P["cores"])
    poolmess = partial(runonce, P=P)
    if P["mp"]:
        r = p.map(poolmess, P["Number"])
    else:  # for debug
        r = list(map(poolmess, P["Number"]))  # for debug
    return r


def partialcmap(thiscmap="bone", minval=0.0, maxval=0.75, n=256):
    import matplotlib.colors as colors

    cmapIn = plt.get_cmap(thiscmap)
    new_cmap = colors.LinearSegmentedColormap.from_list(
        "trunc({n},{a:.2f},{b:.2f})".format(n=cmapIn.name, a=minval, b=maxval),
        cmapIn(np.linspace(minval, maxval, n)),
    )
    return new_cmap


def plotres(output, P, x, y):
    from matplotlib.ticker import ScalarFormatter, AutoMinorLocator, MaxNLocator
    import scipy.ndimage as ndimage

    plt.rc("font", family="serif", serif="Times")
    plt.rc("text", usetex=True)
    plt.rc("xtick", labelsize=7)
    plt.rc("ytick", labelsize=7)
    plt.rc("axes", labelsize=8)
    RT0 = output[0][:, 0]
    RT = output[0][:, 1]
    cycledur = (1.0 / P["freq"][0]) * 1000.0
    shift = cycledur * 100.0  # we wait 100 cycle durations to start with the EoH
    Ishift = int(round(shift))
    spax = int(ceil(np.sqrt(P["TotalN"])))
    fail = []
    delay = []
    amp = []
    slope = []
    for iii in range(P["TotalN"]):
        EoHTime = P["EoH_Time"][iii] * 1000.0
        IEoH = int(np.round((shift + EoHTime) / P["dt"][0]))
        Ileft = int(round(shift / P["dt"][0]))
        Iright = Ileft + (int(round((1.5 * cycledur) / P["dt"][0])))
        Vmax = np.max(RT[iii][Ileft:Iright])
        Imax = np.argmax(RT[iii][Ileft:Iright]) + Ileft
        Tmax = Imax * P["dt"][0]
        amplitude = Vmax - RT[iii][IEoH]
        amp.append(amplitude)
        if Vmax > -30.0:
            fail.append(False)
        else:
            fail.append(True)
        tdel = Tmax - (shift + P["EoH_Time"][0])
        delay.append(tdel)
        slope.append(amplitude / tdel)
    amp = np.array(amp)
    delay = np.array(delay)
    fail = np.array(fail)
    delay[fail] = np.nan
    F = np.reshape(fail, (P["N"], P["N"]))
    A = np.reshape(amp, (P["N"], P["N"]))
    D = np.reshape(delay, (P["N"], P["N"]))
    cbkwargs = {"format": "%.0f"}
    #
    fh1 = plt.figure(figsize=(8.5 / 2.54, (8.5 / 2.54) / 1.25), dpi=600)
    ax1 = plt.subplot(2, 2, 2)
    ax2 = plt.subplot(2, 2, 4, sharex=ax1)
    mycmap = partialcmap()
    mycols = mycmap(np.linspace(0.0, 1.0, P["N"]))
    for ig, thisg in enumerate(x):
        ax1.plot(y * 1000.0, A[:, ig], color=mycols[ig])
        ax2.plot(y * 1000.0, D[:, ig], color=mycols[ig])
    ax1.get_xaxis().set_major_locator(MaxNLocator(nbins=3, min_n_ticks=3))
    ax1.get_yaxis().set_major_locator(MaxNLocator(nbins=3, min_n_ticks=4))
    ax2.get_xaxis().set_major_locator(MaxNLocator(nbins=3, min_n_ticks=3))
    ax2.get_yaxis().set_major_locator(MaxNLocator(nbins=3, min_n_ticks=3))
    ax1.set_xlabel(r"$g_{EoH}$ (nS)")
    ax2.set_xlabel(r"$g_{EoH}$ (nS)")
    ax1.set_ylabel(r"Amplitude (mV)")
    ax2.set_ylabel(r"Delay (ms)")
    sm = plt.cm.ScalarMappable(cmap=mycmap, norm=plt.Normalize(vmin=0.0, vmax=16.0))
    plt.colorbar(sm, ax=ax1, **cbkwargs)
    cbar = plt.colorbar(sm, use_gridspec=True, ax=ax2, **cbkwargs)
    cbar.ax.set_ylabel(r"$g_{syn}$ (nS)")
    tl = MaxNLocator(nbins=4)
    cbar.locator = tl
    cbar.update_ticks()
    return fh1


if __name__ == "__main__":
    ####LOAD DATA FOR Figure6C (IF IT EXISTS)
    if os.path.isfile("./results/Figure6C.npy"):
        print("Data for Figure6C found... loading!")
        output = np.load("./results/Figure6C.npy", allow_pickle=True)
        P = np.load("./results/Figure6C_P.npy", allow_pickle=True)
        P = P.tolist()
    else:
        conditions = 20
        cores = 20
        dur = 400.0
        dt = 0.01
        output = []

        # Model Parameters (all in a linearly aranged fashion, so that a minimal
        # amount of programming is required to change the experiment).
        P = {}
        P["N"] = conditions
        P["cores"] = cores
        P["TotalN"] = int(P["N"] ** 2)
        P["Number"] = list(range(P["TotalN"]))
        P["distritype"] = "normal"
        P["DistParam"] = 1.0
        P["mp"] = True
        #########################################
        P["dur"] = np.repeat(dur, P["TotalN"])
        P["dt"] = np.repeat(dt, P["TotalN"])
        P["Stochastic_EoH"] = np.repeat(False, P["TotalN"])
        P["L"] = np.repeat(50.0, P["TotalN"])
        P["D"] = np.repeat(3.0, P["TotalN"])
        P["Nsyn"] = np.repeat(32, P["TotalN"])
        P["tau_dendsyn"] = np.repeat(2.0, P["TotalN"])
        P["constR"] = np.repeat(False, P["TotalN"])
        P["freq"] = np.repeat(500.0, P["TotalN"])
        P["cf"] = np.repeat(751.0, P["TotalN"])
        P["dB"] = np.repeat(65.0, P["TotalN"])
        P["IPI"] = np.repeat(1000.0 / P["freq"], P["TotalN"])
        P["gLEAK_d"] = np.repeat(0.001, P["TotalN"])
        P["gKLT_s"] = np.repeat(0.017, P["TotalN"])
        P["gIH_s"] = np.repeat(0.002, P["TotalN"])
        P["gKLT_d"] = np.repeat(0.0085, P["TotalN"])
        P["gIH_d"] = np.repeat(0.001, P["TotalN"])
        # set EoH onset time:
        thistime = 0.75 / P["freq"][0]  # in seconds
        P["EoH_Time"] = np.repeat(thistime, P["TotalN"])
        # Now define the two variable parameters. The repeated = y, the tiled = x!!:
        allgEoH = np.linspace(0.025, 0.04, P["N"])  # in muS
        allgDend = np.linspace(0.0, 0.016, P["N"])  # in muS
        P["gsyn"] = np.tile(allgDend, P["N"])
        P["G_EoH"] = np.repeat(allgEoH, P["N"])
        # make go!
        output.append(myMPhandler(P))
        output = np.array(output)
        np.save("./results/Figure6C.npy", output, allow_pickle=True)
        np.save("./results/Figure6C_P.npy", P, allow_pickle=True)

    ###Plot results
    with PdfPages("./figs/Figure6C.pdf") as pp:  # Exp-specific file/variable names here
        fhandle = plotres(
            output=output,
            P=P,
            x=np.unique(P["gsyn"]),
            y=np.unique(P["G_EoH"]),
        )
        plt.tight_layout()
        pp.savefig(fhandle)
        plt.close()
