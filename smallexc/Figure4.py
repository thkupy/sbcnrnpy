#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jul 12 2019
This file contains all code necessary to generate Figure 4.
@author: Thomas Kuenzel (kuenzel{at}bio2.rwth-aachen.de)
"""
import os
import matplotlib

# if os.environ.get('DISPLAY','') == '':
#    print('no display found. Using non-interactive Agg backend')#
matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import semodels as sem
import sehelper as seh
import multiprocessing
from functools import partial


def average_over_cyle(D, freq, dt, dur):
    IPI = 1000.0 / freq
    allc = []
    cycsamps = int(np.round(IPI / dt))
    c = list(
        range(0, len(D) - cycsamps, cycsamps)
    )  # len(D)-cycsamps to make sure that we do not try to access samples > len(D)
    for ii, cc in enumerate(c):
        allc.append(D[cc : cc + cycsamps])
    av = np.mean(np.array(allc), axis=0)
    st = np.std(np.array(allc), axis=0)
    meanD = np.mean(av)
    stdD = np.std(av)
    modamp = np.max(av) - np.min(av)
    return (av, st, meanD, stdD, modamp)


def runonce(x, P):
    print(x)
    ####MODEL INPUT = SOUND
    thisdur = (P["dur"][x] - 5.0) / 1000.0
    SO = seh.create_sound(
        fs=100e3, freq=P["freq"][x], duration=thisdur, dbspl=P["dB"][x]
    )
    SP = seh.create_spk(
        SO, fs=100e3, N=P["Nsyn"][x], cf=P["cf"][x], seed=2121977
    )  # generate enough spiketrains
    S = np.array(SP["spikes"] * 1000.0)
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
    )
    Ev = seh.SimpleDetectAP(R["Vm"], thr=-100, dt=P["dt"][x], LM=-20, RM=10)
    SAPr = len(Ev["PeakT"]) / (P["dur"][x] / 1000.0)
    sav, sst, Smeanv, Stdv, Smodv = average_over_cyle(
        D=R["Vm"], freq=P["freq"][x], dt=P["dt"][x], dur=P["dur"][x]
    )
    #
    # save example data (if required)
    if P["save_examples"] and np.any(P["which_example_to_save"] == x):
        XT = {}
        XT["Sound"] = R["Vm"]
        XT["Av_Sound"] = sav
        XT["St_Sound"] = sst
        thisfn = P["Example_PN"] + P["Example_FN"] + str(x) + ".npy"
        np.save(thisfn, XT, allow_pickle=True)
    #
    return [SAPr, Smeanv, Smodv]


def myMPhandler(P):
    p = multiprocessing.Pool(P["cores"])
    poolmess = partial(runonce, P=P)
    if P["mp"]:
        r = p.map(poolmess, P["Number"])
    else:  # for debug
        r = list(map(poolmess, P["Number"]))  # for debug
    return r


def basic_parameters_1d(N=21, cores=4, dur=500.0, dt=0.025):
    P = {}
    P["mp"] = True
    P["N"] = int(conditions)  # just to be sure
    P["cores"] = int(cores)  # just to be sure
    P["Number"] = list(range(P["N"]))
    P["distritype"] = ["normal", "alpha"]
    P["DistParam"] = [0.505, 1.72]
    P["save_examples"] = True
    P["N-examples"] = 5
    P["Example_PN"] = "./results/tmp/"
    P["Example_FN"] = "ex_"
    P["which_example_to_save"] = np.linspace(0, P["N"] - 1, P["N-examples"], dtype=int)
    P["var_key"] = " "
    P["titletext"] = " "
    P["labeltext"] = " "
    #########################################
    P["dur"] = np.repeat(dur, P["N"])
    P["dt"] = np.repeat(dt, P["N"])
    P["G_EoH"] = np.repeat(0.0, P["N"])
    P["Stochastic_EoH"] = np.repeat(False, P["N"])
    P["Nsyn"] = np.repeat(32, P["N"])
    P["gsyn"] = np.repeat(0.016, P["N"])
    P["tau_dendsyn"] = np.repeat(2.0, P["N"])
    P["L"] = np.repeat(50.0, P["N"])
    P["D"] = np.repeat(3.0, P["N"])
    P["constR"] = np.repeat(False, P["N"])
    P["freq"] = np.repeat(200.0, P["N"])
    P["cf"] = np.repeat(200.0, P["N"])
    P["dB"] = np.repeat(60.0, P["N"])
    P["IPI"] = np.repeat(5.0, P["N"])
    P["gLEAK_d"] = np.repeat(0.001, P["N"])
    P["gKLT_s"] = np.repeat(0.017, P["N"])
    P["gIH_s"] = np.repeat(0.002, P["N"])
    P["gKLT_d"] = np.repeat(0.0085, P["N"])
    P["gIH_d"] = np.repeat(0.001, P["N"])
    return P


def process_examplefiles(P):
    exfilelist = []
    ST = []
    AvS = []
    StS = []
    for file in os.listdir(P["Example_PN"]):
        if file.endswith(".npy"):
            exfilelist.append(os.path.join(P["Example_PN"], file))
    for thisfile in sorted(exfilelist):
        tx = np.load(thisfile, allow_pickle=True)
        tx = tx.tolist()
        ST.append(tx["Sound"])
        AvS.append(tx["Av_Sound"])
        StS.append(tx["St_Sound"])
    XT = {}
    XT["Sound"] = ST
    XT["Av_Sound"] = AvS
    XT["St_Sound"] = StS
    # remove the files now
    for thisfile in exfilelist:
        os.remove(thisfile)
    return XT


def plotres_2d(output, P, x, y, xlabs, ylabs):
    import matplotlib.font_manager
    from matplotlib.ticker import ScalarFormatter, AutoMinorLocator, MaxNLocator
    import scipy.ndimage as ndimage

    plt.rc("font", family="serif", serif="Times")
    plt.rc("text", usetex=True)
    plt.rc("xtick", labelsize=7)
    plt.rc("ytick", labelsize=7)
    plt.rc("axes", labelsize=8)
    nticks = 5
    mycmap = "bone"
    filtsig = 0.25

    outV = output[0][:, 1]  # return [SAPr,Smeanv,Smodv]
    outM = output[0][:, 2]
    outA = output[0][:, 0]
    V = ndimage.gaussian_filter(
        np.reshape(outV, (P["N"], P["N"])), sigma=filtsig, order=0
    )
    M = ndimage.gaussian_filter(
        np.reshape(outM, (P["N"], P["N"])), sigma=filtsig, order=0
    )
    A = ndimage.gaussian_filter(
        np.reshape(outA, (P["N"], P["N"])), sigma=filtsig, order=0
    )

    fwidth = 7.0  # cm
    fhandle = plt.figure(figsize=(fwidth / 2.54, (fwidth * 2) / 2.54), dpi=300)

    ax1 = plt.subplot(311)
    CS1 = plt.contourf(x, y, V, 21, cmap=mycmap)  # repeated = y, tiled = x!!
    # CS1=plt.pcolormesh(x,y,V,shading='gouraud',cmap=mycmap)#repeated = y, tiled = x!!
    plt.title("Mean Vm")
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    ax1.yaxis.set_major_locator(MaxNLocator(nticks))
    ax1.xaxis.set_major_locator(MaxNLocator(nticks))
    cbar1 = plt.colorbar(CS1, use_gridspec=True)
    cbar1.ax.set_ylabel("Vm (mV)")
    tl = MaxNLocator(nbins=5)
    cbar1.locator = tl
    cbar1.update_ticks()

    ax2 = plt.subplot(312)
    CS2 = plt.contourf(x, y, M, 21, cmap=mycmap)  # repeated = y, tiled = x!!
    # CS2=plt.pcolormesh(x,y,M,shading='gouraud',cmap=mycmap)#repeated = y, tiled = x!!
    plt.title("Vm Modulation")
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    ax2.yaxis.set_major_locator(MaxNLocator(nticks))
    ax2.xaxis.set_major_locator(MaxNLocator(nticks))
    cbar2 = plt.colorbar(CS2, use_gridspec=True)
    cbar2.ax.set_ylabel(r"Vm mod. ($\pm$ mV)")
    tl = MaxNLocator(nbins=5)
    cbar2.locator = tl
    cbar2.update_ticks()
    tl = MaxNLocator(nbins=5)
    cbar2.locator = tl
    cbar2.update_ticks()

    ax3 = plt.subplot(313)
    CS3 = plt.contourf(x, y, A, 21, cmap=mycmap)  # repeated = y, tiled = x!!
    # CS3=plt.pcolormesh(x,y,A,shading='gouraud',cmap=mycmap)#repeated = y, tiled = x!!
    plt.title("APs")
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    ax3.yaxis.set_major_locator(MaxNLocator(nticks))
    ax3.xaxis.set_major_locator(MaxNLocator(nticks))
    cbar3 = plt.colorbar(CS3, use_gridspec=True)
    cbar3.ax.set_ylabel("AP rate (Hz)")

    plt.tight_layout()
    return fhandle


def plotresult_1d(output, example, P):
    # Imports and Settings
    from matplotlib.ticker import ScalarFormatter, AutoMinorLocator, MaxNLocator
    import scipy.ndimage as ndimage

    plt.rc("font", family="serif", serif="Times")
    plt.rc("text", usetex=True)
    plt.rc("xtick", labelsize=7)
    plt.rc("ytick", labelsize=7)
    plt.rc("axes", labelsize=8)
    titsize = 9
    mycmap = "bone"
    nticks = 5
    filtsig = 0.25
    r_crit = 0.001
    #
    Scol = "k"
    Scmap = "bone"
    #
    NUM_COLORS = P["N-examples"] + 3
    Scm = plt.get_cmap(Scmap)
    unique_Scolors = [Scm(1.0 * (i + 1) / NUM_COLORS) for i in range(NUM_COLORS)]
    #
    # Prep data [-->SAPr,Smeanv,Smodv]
    SAPr = output[0][:, 0]
    Smeanv = output[0][:, 1]
    Smodv = output[0][:, 2]
    XD = P[P["var_key"]]
    #
    # Create Figures
    fwidth = 8.5  # cm
    fwidth_full = 17.0  # cm
    fhandle1 = plt.figure(figsize=(5 / 2.54, 11 / 2.54), dpi=300)

    ax1 = plt.subplot(311)
    # To specify the number of ticks on both or any single axes
    plt.plot(XD, Smeanv, color=Scol)
    plt.ylabel("Mean Vm (mV)")
    plt.xlabel(P["labeltext"], labelpad=15)
    plt.title(P["titletext"], fontsize=titsize)
    ax1.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax1.yaxis.set_major_locator(plt.MaxNLocator(5))
    #
    ax2 = plt.subplot(312)
    plt.plot(XD, Smodv, color=Scol)
    plt.ylabel("Vm Modulation (mV)", labelpad=15)
    plt.xlabel(P["labeltext"])
    ax2.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax2.yaxis.set_major_locator(plt.MaxNLocator(5))
    plt.tight_layout()
    #
    if P["save_examples"]:
        fhandle2 = plt.figure(figsize=(5.5 / 2.54, (fwidth) / 2.54), dpi=300)
        s = []
        s.append(plt.subplot(211))
        s.append(plt.subplot(212))
        taxis = np.linspace(0, P["dur"][0] - P["dt"][0], P["dur"][0] / P["dt"][0])
        cycaxis = np.linspace(
            0, P["IPI"][0] - P["dt"][0], int(np.round(P["IPI"][0] / P["dt"][0]))
        )
        for ix in range(P["N-examples"]):
            if P["var_key"] == "IPI":
                cycaxis = np.linspace(0.0, 1.0, len(example["Av_Sound"][ix]))
                xlabtext2 = "Cycle (norm.)"
            else:
                xlabtext2 = "IPI (ms)"
            s[0].plot(taxis, example["Sound"][ix], color=unique_Scolors[ix])
            s[1].plot(cycaxis, example["Av_Sound"][ix], color=unique_Scolors[ix])
        s[0].set_xlim(100.0, 150.0)
        s[0].set_title("Sound / Example trace", fontsize=titsize)
        s[0].set_xlabel("Time (ms)")
        legstring = [
            "%.3f" % number for number in P[P["var_key"]][P["which_example_to_save"]]
        ]
        s[0].legend(legstring, fontsize=5)
        s[1].set_title("Sound / averaged trace", fontsize=titsize)
        s[1].set_xlabel(xlabtext2)
        s[1].legend(legstring, fontsize=6)
        for thish in s:
            thish.xaxis.set_major_locator(plt.MaxNLocator(5))
            thish.yaxis.set_major_locator(plt.MaxNLocator(5))
    else:
        fhandle2 = []
    plt.tight_layout()
    return ((fhandle1, fhandle2), (ax1, ax2, s))


if __name__ == "__main__":
    Fig4_loadifexists = [True, True, True]
    conditions = 64  # production = 64
    conditions_2d = 25  # production = 25
    cores = 20  # production = 20
    dur = 5000.0  # production = 5000.0
    dt = 0.01  # production = 0.01

    myylim1 = (-64.0, -52.5)
    myylim2 = (0, 4.5)
    #####################################
    #####Fig4a:         #################
    #####---------------#################
    ##### Variable -> tau dendsyn #######
    ##### R_dend = free #################
    ##### No Endbulb    #################
    #####################################
    ####LOAD DATA FOR Fig4a (IF IT EXISTS)
    if (
        os.path.isfile("./results/Fig4a.npy") and Fig4_loadifexists[0]
    ):  # Exp-specific file/variable names here
        print("Data for Fig4a found... loading!")
        output4a = np.load(
            "./results/Fig4a.npy", allow_pickle=True
        )  # Exp-specific file/variable names here
        Ex4a = np.load(
            "./results/Fig4a_Ex.npy", allow_pickle=True
        )  # Exp-specific file/variable names here
        Ex4a = Ex4a.tolist()
        P4a = np.load(
            "./results/Fig4a_P.npy", allow_pickle=True
        )  # Exp-specific file/variable names here
        P4a = P4a.tolist()  # Exp-specific file/variable names here
    else:
        output4a = []  # Exp-specific file/variable names here
        P4a = basic_parameters_1d(
            conditions, cores, dur, dt
        )  # Exp-specific file/variable names here
        ####VARIABLE Parameter for EXP1a####
        P4a["mp"] = True
        P4a["var_key"] = "tau_dendsyn"
        P4a["titletext"] = r"Synaptic decay tau (dendritic)"
        P4a["labeltext"] = r"Decay tau (ms)"
        P4a["tau_dendsyn"] = np.linspace(0.2, 20.0, P4a["N"])  # Exp-specific code here
        # make go!
        output4a.append(myMPhandler(P4a))  # Exp-specific file/variable names here
        output4a = np.array(output4a)  # Exp-specific file/variable names here
        if P4a["save_examples"]:
            Ex4a = process_examplefiles(P4a)
        np.save(
            "./results/Fig4a.npy", output4a, allow_pickle=True
        )  # Exp-specific file/variable names here
        np.save(
            "./results/Fig4a_P.npy", P4a, allow_pickle=True
        )  # Exp-specific file/variable names here
        np.save(
            "./results/Fig4a_Ex.npy", Ex4a, allow_pickle=True
        )  # Exp-specific file/variable names here
    ###Plot results
    with PdfPages("./figs/Fig4a.pdf") as pp:  # Exp-specific file/variable names here
        fhandles, axhandles = plotresult_1d(output4a, Ex4a, P4a)
        # panel specific axis limits and other artwork
        axhandles[0].set_ylim(myylim1)
        axhandles[1].set_ylim(myylim2)
        # axhandles[2][0].set_ylim((-65.0,-51.0))
        # axhandles[2][1].set_ylim((-65.0,-51.0))
        #
        pp.savefig(fhandles[0], bbox_inches="tight")
        pp.savefig(fhandles[1], bbox_inches="tight")
        plt.close()

    with PdfPages(
        "./figs/Fig4a_Insets.pdf"
    ) as pp:  # Exp-specific file/variable names here
        fhandles, axhandles = plotresult_1d(output4a, Ex4a, P4a)
        # panel specific axis limits and other artwork
        axhandles[0].set_ylim((-64.5, -61.0))
        axhandles[1].set_ylim(myylim2)
        axhandles[0].set_xlim((0.2, 1.2))
        axhandles[1].set_xlim((0.2, 1.2))
        #
        pp.savefig(fhandles[0], bbox_inches="tight")
        plt.close()

    #####################################
    #####Fig4b:         #################
    #####---------------#################
    ##### Variable -> freq ##############
    ##### R_dend = free #################
    ##### No Endbulb    #################
    #####################################
    ####LOAD DATA FOR Fig4b (IF IT EXISTS)
    if (
        os.path.isfile("./results/Fig4b.npy") and Fig4_loadifexists[1]
    ):  # Exp-specific file/variable names here
        print("Data for Fig4b found... loading!")
        output4b = np.load(
            "./results/Fig4b.npy", allow_pickle=True
        )  # Exp-specific file/variable names here
        Ex4b = np.load(
            "./results/Fig4b_Ex.npy", allow_pickle=True
        )  # Exp-specific file/variable names here
        Ex4b = Ex4b.tolist()
        P4b = np.load(
            "./results/Fig4b_P.npy", allow_pickle=True
        )  # Exp-specific file/variable names here
        P4b = P4b.tolist()  # Exp-specific file/variable names here
    else:
        output4b = []  # Exp-specific file/variable names here
        P4b = basic_parameters_1d(
            conditions, cores, dur, dt
        )  # Exp-specific file/variable names here
        ####VARIABLE Parameter for EXP1b###
        P4b["mp"] = True
        P4b["var_key"] = "IPI"
        P4b["titletext"] = r"Input frequency"
        P4b["labeltext"] = r"Cycle duration (ms)"
        P4b["IPI"] = np.round(
            np.linspace(2.0, 7.5, P4b["N"]), 2
        )  # Exp-specific code here
        P4b["freq"] = 1000.0 / P4b["IPI"]
        P4b["cf"] = P4b["freq"]
        # make go!
        output4b.append(myMPhandler(P4b))  # Exp-specific file/variable names here
        output4b = np.array(output4b)  # Exp-specific file/variable names here
        if P4b["save_examples"]:
            Ex4b = process_examplefiles(P4b)
        np.save(
            "./results/Fig4b.npy", output4b, allow_pickle=True
        )  # Exp-specific file/variable names here
        np.save(
            "./results/Fig4b_P.npy", P4b, allow_pickle=True
        )  # Exp-specific file/variable names here
        np.save(
            "./results/Fig4b_Ex.npy", Ex4b, allow_pickle=True
        )  # Exp-specific file/variable names here
    ###Plot results
    with PdfPages("./figs/Fig4b.pdf") as pp:  # Exp-specific file/variable names here
        fhandle, axhandles = plotresult_1d(output4b, Ex4b, P4b)
        # panel specific axis limits and other artwork
        axhandles[0].set_ylim(myylim1)
        axhandles[1].set_ylim(myylim2)
        # axhandles[2][0].set_ylim((-65.0,-51.0))
        # axhandles[2][1].set_ylim((-65.0,-51.0))
        pp.savefig(fhandle[0], bbox_inches="tight")
        pp.savefig(fhandle[1], bbox_inches="tight")
        plt.close()

    #####################################
    #####Fig4e:         #################
    #####---------------#################
    ##### Variable -> tau dend & freq ###
    ##### R_dend = fix  #################
    ##### No Endbulb    #################
    #####################################
    ####LOAD DATA FOR Fig4e (IF IT EXISTS)
    if os.path.isfile("./results/Fig4e.npy") and Fig4_loadifexists[2]:
        print("Data for Fig4e found... loading!")
        output4e = np.load("./results/Fig4e.npy", allow_pickle=True)
        P4e = np.load("./results/Fig4e_P.npy", allow_pickle=True)
        P4e = P4e.tolist()
    else:
        # Model Parameters (all in a linearly aranged fashion, so that a minimal
        # amount of programming is required to change the experiment).
        P4e = {}
        P4e["mp"] = True
        P4e["N"] = int(conditions_2d)  # just to be sure
        P4e["cores"] = int(cores)  # just to be sure
        P4e["TotalN"] = int(P4e["N"] ** 2)
        P4e["Number"] = list(range(P4e["TotalN"]))
        P4e["save_examples"] = False
        #########################################
        P4e["dur"] = np.repeat(dur, P4e["TotalN"])
        P4e["dt"] = np.repeat(dt, P4e["TotalN"])
        P4e["L"] = np.repeat(50.0, P4e["TotalN"])
        P4e["D"] = np.tile(3.0, P4e["TotalN"])
        P4e["G_EoH"] = np.repeat(0.0, P4e["TotalN"])
        P4e["Stochastic_EoH"] = np.repeat(False, P4e["TotalN"])
        P4e["Nsyn"] = np.repeat(32, P4e["TotalN"])
        P4e["gsyn"] = np.repeat(0.016, P4e["TotalN"])
        P4e["constR"] = np.repeat(False, P4e["TotalN"])
        P4e["dB"] = np.repeat(60.0, P4e["TotalN"])
        P4e["gKLT_d"] = np.repeat(0.0085, P4e["TotalN"])
        P4e["gIH_d"] = np.repeat(0.001, P4e["TotalN"])
        P4e["gLEAK_d"] = np.repeat(0.001, P4e["TotalN"])
        P4e["gKLT_s"] = np.repeat(0.017, P4e["TotalN"])
        P4e["gIH_s"] = np.repeat(0.002, P4e["TotalN"])

        # Now define the two variable parameters. The repeated = y, the tiled = x!!
        tau = np.linspace(0.2, 20.0, P4e["N"])
        ipi = np.round(np.linspace(2.0, 7.5, P4e["N"]), 2)  # Exp-specific code here
        freq = 1000.0 / ipi
        P4e["IPI"] = np.tile(ipi, P4e["N"])
        P4e["freq"] = np.tile(freq, P4e["N"])
        P4e["cf"] = P4e["freq"]
        P4e["tau_dendsyn"] = np.repeat(tau, P4e["N"])

        # make go!
        output4e = []
        output4e.append(myMPhandler(P4e))
        output4e = np.array(output4e)
        np.save("./results/Fig4e.npy", output4e, allow_pickle=True)
        np.save("./results/Fig4e_P.npy", P4e, allow_pickle=True)

    ###Plot results
    fhandle = plotres_2d(
        output=output4e,
        P=P4e,
        x=np.unique(P4e["IPI"]),
        y=np.unique(P4e["tau_dendsyn"]),
        xlabs=r"IPI (ms)",
        ylabs=r"Dendritic decay tau (ms)",
    )  # x=np.unique(P4e['freq']); xlabs=r'Input frequency (Hz)'

    pp = PdfPages("./figs/Fig4e.pdf")
    pp.savefig(bbox_inches="tight")
    pp.close()
