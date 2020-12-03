#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This file contains all code necessary to generate Figure 6 panel F & G.

The example traces in F were taken from a run with only 5x5 conditions (see line 283),
the 2D-plots were produced from 25x25 conditions.

To test this rapidly, change dt (line 286) to coarser values and the and number of
conditions (line 283) to lower (5x5 works well).
Adjust number of available cores in line 284.

Created on 21 Aug 2019
Revised sept 2020
@author: kuenzel
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
    print((str(P["G_EoH"][x]) + " " + str(P["EoH_Time"][x]) + " " + str(x)))
    # thistrace=R['Vm'][int(round(shift-cycledur)):int(round(shift+(2.0*cycledur)))]
    # return thistrace
    return [R0["Vm"], R["Vm"]]


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

    plt.rc("font", family="serif", serif="Times")
    plt.rc("text", usetex=True)
    plt.rc("xtick", labelsize=7)
    plt.rc("ytick", labelsize=7)
    plt.rc("axes", labelsize=8)
    titsize = 8
    mycmap = "bone"

    RT0 = output[0][:, 0]
    RT = output[0][:, 1]
    #    T=np.reshape(RT,(P['N'],P['N']))
    #    T0=np.reshape(RT0,(P['N'],P['N']))
    # fwidth=7.0#cm
    # fhandle=plt.figure(figsize=(fwidth/2.54, (fwidth*2)/2.54), dpi=300)
    cycledur = (1.0 / P["freq"][0]) * 1000.0
    shift = cycledur * 100.0  # we wait 100 cycle durations to start with the EoH
    Ishift = int(round(shift))

    fh1 = plt.figure()
    eP = P
    eP["G_EoH"] = np.repeat(0.0, P["TotalN"])
    noEoH = runonce(0, eP)
    noEoH = noEoH[1]
    Tax = np.linspace(
        0.0, P["dur"][0] - P["dt"][0], int(round(P["dur"][0] / P["dt"][0]))
    )
    spax = int(ceil(np.sqrt(P["TotalN"])))
    fail = []
    fail0 = []
    delay0 = []
    delay = []
    amp0 = []
    amp = []
    slope = []
    slope0 = []

    for iii in range(P["TotalN"]):
        EoHTime = P["EoH_Time"][iii] * 1000.0
        fh1.add_subplot(spax, spax, iii + 1)
        plt.plot(Tax, noEoH, "r")
        plt.plot(Tax, RT0[iii], "k-")
        plt.plot(Tax, RT[iii], "g-")
        IEoH = int(np.round((shift + EoHTime) / P["dt"][0]))
        plt.plot(Tax[IEoH], RT[iii][IEoH], "m^")

        plt.xlim((shift, shift + (1.5 * cycledur)))
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
        Vmax0 = np.max(RT0[iii][Ileft:Iright])
        Imax0 = np.argmax(RT0[iii][Ileft:Iright]) + Ileft
        Tmax0 = Imax0 * P["dt"][0]
        amplitude0 = Vmax0 - RT0[iii][IEoH]
        amp0.append(amplitude0)
        if Vmax0 > -30.0:
            fail0.append(False)
        else:
            fail0.append(True)
        tdel = Tmax - (shift + EoHTime)
        tdel0 = Tmax0 - (shift + EoHTime)
        delay.append(tdel)
        delay0.append(tdel0)
        slope.append(amplitude / tdel)
        slope0.append(amplitude0 / tdel0)
        plt.plot(Tmax, Vmax, "gx")
        plt.plot(Tmax0, Vmax0, "kx")

    fh2 = plt.figure(figsize=(8.5 / 2.54, (8.5 / 2.54) / 1.25), dpi=600)
    F0 = np.reshape(fail0, (P["N"], P["N"]))
    F = np.reshape(fail, (P["N"], P["N"]))
    A0 = np.reshape(amp0, (P["N"], P["N"]))
    A = np.reshape(amp, (P["N"], P["N"]))
    D0 = np.reshape(delay0, (P["N"], P["N"]))
    D = np.reshape(delay, (P["N"], P["N"]))

    cbkwargs = {"format": "%.1f"}
    ax1 = plt.subplot(2, 2, 1)
    # absomin = np.min(np.array([np.min(A0),np.min(A)]))
    # absomax = np.max(np.array([np.max(A0),np.max(A)]))
    conts = np.linspace(30.0, 75.0, 17)
    CS1 = plt.contourf((x * 1000.0) / 2.0, y * 1000.0, A0, conts, cmap=mycmap)
    plt.contour(
        (x * 1000.0) / 2.0, y * 1000.0, F0, (0.001,), colors="yellow"
    )  # repeated = y, tiled = x!!
    plt.title(r"AP amplitude, $g_{dendsyn} = 0nS$", fontsize=titsize)
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    ax1.get_xaxis().set_major_formatter(ScalarFormatter())
    ax1.get_yaxis().set_major_formatter(ScalarFormatter())
    cbar1 = plt.colorbar(CS1, use_gridspec=True, **cbkwargs)
    cbar1.ax.set_ylabel(r"Amplitude (mV)")
    tl = MaxNLocator(nbins=5)
    cbar1.locator = tl
    cbar1.update_ticks()

    ax2 = plt.subplot(2, 2, 2)
    CS2 = plt.contourf((x * 1000.0) / 2.0, y * 1000.0, A, conts, cmap=mycmap)
    plt.contour(
        (x * 1000.0) / 2.0, y * 1000.0, F, (0.001,), colors="yellow"
    )  # repeated = y, tiled = x!!
    plt.title(r"AP amplitude, $g_{dendsyn} = 16nS$", fontsize=titsize)
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    ax2.get_xaxis().set_major_formatter(ScalarFormatter())
    ax2.get_yaxis().set_major_formatter(ScalarFormatter())
    cbar2 = plt.colorbar(CS2, use_gridspec=True, **cbkwargs)
    cbar2.ax.set_ylabel(r"Amplitude (mV)")
    tl = MaxNLocator(nbins=5)
    cbar2.locator = tl
    cbar2.update_ticks()

    # absomin = np.min(np.array([np.min(D0),np.min(D)]))
    # absomax = np.max(np.array([np.max(D0),np.max(D)]))
    conts = np.linspace(0.3, 1.2, 17)
    ax3 = plt.subplot(2, 2, 3)
    CS3 = plt.contourf((x * 1000.0) / 2.0, y * 1000.0, D0, conts, cmap=mycmap)
    plt.contour(
        (x * 1000.0) / 2.0, y * 1000.0, F0, (0.001,), colors="yellow"
    )  # repeated = y, tiled = x!!
    plt.title(r"AP delay, $g_{dendsyn} = 0nS$", fontsize=titsize)
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    ax3.get_xaxis().set_major_formatter(ScalarFormatter())
    ax3.get_yaxis().set_major_formatter(ScalarFormatter())
    cbar3 = plt.colorbar(CS3, use_gridspec=True, **cbkwargs)
    cbar3.ax.set_ylabel(r"Delay (ms)")
    tl = MaxNLocator(nbins=5)
    cbar3.locator = tl
    cbar3.update_ticks()

    ax4 = plt.subplot(2, 2, 4)
    CS4 = plt.contourf((x * 1000.0) / 2.0, y * 1000.0, D, conts, cmap=mycmap)
    plt.contour(
        (x * 1000.0) / 2.0, y * 1000.0, F, (0.001,), colors="yellow"
    )  # repeated = y, tiled = x!!
    plt.title(r"AP delay, $g_{dendsyn} = 16nS$", fontsize=titsize)
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    ax4.get_xaxis().set_major_formatter(ScalarFormatter())
    ax4.get_yaxis().set_major_formatter(ScalarFormatter())
    cbar4 = plt.colorbar(CS4, use_gridspec=True, **cbkwargs)
    cbar4.ax.set_ylabel(r"Delay (ms)")
    tl = MaxNLocator(nbins=5)
    cbar4.locator = tl
    cbar4.update_ticks()
    return (fh1, fh2)


if __name__ == "__main__":
    ####LOAD DATA FOR Figure6B (IF IT EXISTS)
    if os.path.isfile("./results/Figure6B.npy"):
        print("Data for Figure6B found... loading!")
        output = np.load("./results/Figure6B.npy", allow_pickle=True)
        P = np.load("./results/Figure6B_P.npy", allow_pickle=True)
        P = P.tolist()
    else:
        conditions = 25
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
        P["gsyn"] = np.repeat(0.016, P["TotalN"])
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

        # Now define the two variable parameters. The repeated = y, the tiled = x!!
        alltimes = np.linspace(0.0, ((1.0 / P["freq"][0])), P["N"])  # in seconds
        allgs = np.linspace(0.025, 0.04, P["N"])  # in muS
        P["EoH_Time"] = np.tile(alltimes, P["N"])
        P["G_EoH"] = np.repeat(allgs, P["N"])

        # make go!
        output.append(myMPhandler(P))
        output = np.array(output)
        np.save("./results/Figure6B.npy", output, allow_pickle=True)
        np.save("./results/Figure6B_P.npy", P, allow_pickle=True)

    ###Plot results
    with PdfPages("./figs/Figure6B.pdf") as pp:  # Exp-specific file/variable names here
        fhandle = plotres(
            output=output,
            P=P,
            x=np.unique(P["EoH_Time"]),
            y=np.unique(P["G_EoH"]),
            xlabs=r"EoH time (ms, re cycle)",
            ylabs=r"gEoH (nS)",
        )
        plt.tight_layout()
        pp.savefig(fhandle[0])
        pp.savefig(fhandle[1])
        plt.close()
