#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Jul 12 2019
This file contains all code necessary to generate Figure 5.
@author: Thomas Kuenzel (kuenzel{at}bio2.rwth-aachen.de)
"""
import os
import matplotlib

# if os.environ.get('DISPLAY','') == '':
#    print('no display found. Using non-interactive Agg backend')
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
    Fig5_loadifexists = [True, True, True, True]
    conditions = 64  # production = 64
    conditions_2d = 25  # production = 25
    cores = 20  # production = 20
    dur = 2000.0  # production = 5000.0
    dt = 0.01  # production = 0.01

    # myxlim=(0.0,)
    myylim1 = (-63.0, -61.0)
    myylim2 = (0.0, 4.5)
    #####################################
    #####Fig5a:         #################
    #####---------------#################
    ##### Variable -> gKLT ##############
    ##### R_dend = free #################
    ##### No Endbulb    #################
    #####################################
    ####LOAD DATA FOR Fig5a (IF IT EXISTS)
    if (
        os.path.isfile("./results/Fig5a.npy") and Fig5_loadifexists[0]
    ):  # Exp-specific file/variable names here
        print("Data for Fig5a found... loading!")
        output5a = np.load(
            "./results/Fig5a.npy", allow_pickle=True
        )  # Exp-specific file/variable names here
        Ex5a = np.load(
            "./results/Fig5a_Ex.npy", allow_pickle=True
        )  # Exp-specific file/variable names here
        Ex5a = Ex5a.tolist()
        P5a = np.load(
            "./results/Fig5a_P.npy", allow_pickle=True
        )  # Exp-specific file/variable names here
        P5a = P5a.tolist()  # Exp-specific file/variable names here
    else:
        output5a = []  # Exp-specific file/variable names here
        P5a = basic_parameters_1d(
            conditions, cores, dur, dt
        )  # Exp-specific file/variable names here
        ####VARIABLE Parameter for EXP1a####
        P5a["mp"] = True
        P5a["var_key"] = "gKLT_d"
        P5a["titletext"] = r"Dendritic $g_{KLT}$"
        P5a["labeltext"] = r"$g_{KLT}$ ($\mu$S)"
        P5a["gKLT_d"] = np.linspace(0.0001, 0.035, P5a["N"])  # Exp-specific code here
        # make go!
        output5a.append(myMPhandler(P5a))  # Exp-specific file/variable names here
        output5a = np.array(output5a)  # Exp-specific file/variable names here
        if P5a["save_examples"]:
            Ex5a = process_examplefiles(P5a)
        np.save(
            "./results/Fig5a.npy", output5a, allow_pickle=True
        )  # Exp-specific file/variable names here
        np.save(
            "./results/Fig5a_P.npy", P5a, allow_pickle=True
        )  # Exp-specific file/variable names here
        np.save(
            "./results/Fig5a_Ex.npy", Ex5a, allow_pickle=True
        )  # Exp-specific file/variable names here
    ###Plot results
    with PdfPages("./figs/Fig5a.pdf") as pp:  # Exp-specific file/variable names here
        fhandles, axhandles = plotresult_1d(output5a, Ex5a, P5a)
        # panel specific axis limits and other artwork
        # axhandles[0].set_xlim(myxlim)
        # axhandles[0].set_xlim(myxlim)
        axhandles[0].set_ylim(myylim1)
        axhandles[1].set_ylim(myylim2)
        # axhandles[2][0].set_ylim((-65.0,-51.0))
        # axhandles[2][1].set_ylim((-65.0,-51.0))
        #
        pp.savefig(fhandles[0], bbox_inches="tight")
        pp.savefig(fhandles[1], bbox_inches="tight")
        plt.close()

    #####################################
    #####Fig5b:         #################
    #####---------------#################
    ##### Variable -> gH ################
    ##### R_dend = free #################
    ##### No Endbulb    #################
    #####################################
    ####LOAD DATA FOR Fig5b (IF IT EXISTS)
    if (
        os.path.isfile("./results/Fig5b.npy") and Fig5_loadifexists[1]
    ):  # Exp-specific file/variable names here
        print("Data for Fig5b found... loading!")
        output5b = np.load(
            "./results/Fig5b.npy", allow_pickle=True
        )  # Exp-specific file/variable names here
        Ex5b = np.load(
            "./results/Fig5b_Ex.npy", allow_pickle=True
        )  # Exp-specific file/variable names here
        Ex5b = Ex5b.tolist()
        P5b = np.load(
            "./results/Fig5b_P.npy", allow_pickle=True
        )  # Exp-specific file/variable names here
        P5b = P5b.tolist()  # Exp-specific file/variable names here
    else:
        output5b = []  # Exp-specific file/variable names here
        P5b = basic_parameters_1d(
            conditions, cores, dur, dt
        )  # Exp-specific file/variable names here
        ####VARIABLE Parameter for EXP1b###
        P5b["mp"] = True
        P5b["var_key"] = "gIH_d"
        P5b["titletext"] = r"Dendritic $g_{H}$"
        P5b["labeltext"] = r"$g_{H}$ ($\mu$S)"
        P5b["gIH_d"] = np.linspace(0.0001, 0.005, P5b["N"])  # Exp-specific code here
        # make go!
        output5b.append(myMPhandler(P5b))  # Exp-specific file/variable names here
        output5b = np.array(output5b)  # Exp-specific file/variable names here
        if P5b["save_examples"]:
            Ex5b = process_examplefiles(P5b)
        np.save(
            "./results/Fig5b.npy", output5b, allow_pickle=True
        )  # Exp-specific file/variable names here
        np.save(
            "./results/Fig5b_P.npy", P5b, allow_pickle=True
        )  # Exp-specific file/variable names here
        np.save(
            "./results/Fig5b_Ex.npy", Ex5b, allow_pickle=True
        )  # Exp-specific file/variable names here
    ###Plot results
    with PdfPages("./figs/Fig5b.pdf") as pp:  # Exp-specific file/variable names here
        fhandle, axhandles = plotresult_1d(output5b, Ex5b, P5b)
        # panel specific axis limits and other artwork
        axhandles[0].set_ylim(myylim1)
        axhandles[1].set_ylim(myylim2)
        # axhandles[2][0].set_ylim((-65.0,-51.0))
        # axhandles[2][1].set_ylim((-65.0,-51.0))
        pp.savefig(fhandle[0], bbox_inches="tight")
        pp.savefig(fhandle[1], bbox_inches="tight")
        plt.close()

    #####################################
    #####Fig5c:         #################
    #####---------------#################
    ##### Variable -> gLEAK #############
    ##### R_dend = free #################
    ##### No Endbulb    #################
    #####################################
    ####LOAD DATA FOR Fig5c (IF IT EXISTS)
    if (
        os.path.isfile("./results/Fig5c.npy") and Fig5_loadifexists[2]
    ):  # Exp-specific file/variable names here
        print("Data for Fig5c found... loading!")
        output5c = np.load(
            "./results/Fig5c.npy", allow_pickle=True
        )  # Exp-specific file/variable names here
        Ex5c = np.load(
            "./results/Fig5c_Ex.npy", allow_pickle=True
        )  # Exp-specific file/variable names here
        Ex5c = Ex5c.tolist()
        P5c = np.load(
            "./results/Fig5c_P.npy", allow_pickle=True
        )  # Exp-specific file/variable names here
        P5c = P5c.tolist()  # Exp-specific file/variable names here
    else:
        output5c = []  # Exp-specific file/variable names here
        P5c = basic_parameters_1d(
            conditions, cores, dur, dt
        )  # Exp-specific file/variable names here
        ####VARIABLE Parameter for EXP1b###
        P5c["mp"] = True
        P5c["var_key"] = "gLEAK_d"
        P5c["titletext"] = r"Dendritic $g_{Leak}$"
        P5c["labeltext"] = r"$g_{Leak}$ ($\mu$S)"
        P5c["gLEAK_d"] = np.linspace(0.0001, 0.01, P5c["N"])  # Exp-specific code here
        # make go!
        output5c.append(myMPhandler(P5c))  # Exp-specific file/variable names here
        output5c = np.array(output5c)  # Exp-specific file/variable names here
        if P5c["save_examples"]:
            Ex5c = process_examplefiles(P5c)
        np.save(
            "./results/Fig5c.npy", output5c, allow_pickle=True
        )  # Exp-specific file/variable names here
        np.save(
            "./results/Fig5c_P.npy", P5c, allow_pickle=True
        )  # Exp-specific file/variable names here
        np.save(
            "./results/Fig5c_Ex.npy", Ex5c, allow_pickle=True
        )  # Exp-specific file/variable names here
    ###Plot results
    with PdfPages("./figs/Fig5c.pdf") as pp:  # Exp-specific file/variable names here
        fhandle, axhandles = plotresult_1d(output5c, Ex5c, P5c)
        # panel specific axis limits and other artwork
        axhandles[0].set_ylim(myylim1)
        axhandles[1].set_ylim(myylim2)
        # axhandles[2][0].set_ylim((-65.0,-51.0))
        # axhandles[2][1].set_ylim((-65.0,-51.0))
        pp.savefig(fhandle[0], bbox_inches="tight")
        pp.savefig(fhandle[1], bbox_inches="tight")
        plt.close()

    #####################################
    #####Fig5e:         #################
    #####---------------#################
    ##### Variable -> tau dend & freq ###
    ##### R_dend = fix  #################
    ##### No Endbulb    #################
    #####################################
    ####LOAD DATA FOR Fig5e (IF IT EXISTS)
    if os.path.isfile("./results/Fig5e.npy") and Fig5_loadifexists[3]:
        print("Data for Fig5e found... loading!")
        output5e = np.load("./results/Fig5e.npy", allow_pickle=True)
        P5e = np.load("./results/Fig5e_P.npy", allow_pickle=True)
        P5e = P5e.tolist()
    else:
        # Model Parameters (all in a linearly aranged fashion, so that a minimal
        # amount of programming is required to change the experiment).
        P5e = {}
        P5e["mp"] = True
        P5e["N"] = int(conditions_2d)  # just to be sure
        P5e["cores"] = int(cores)  # just to be sure
        P5e["TotalN"] = int(P5e["N"] ** 2)
        P5e["Number"] = list(range(P5e["TotalN"]))
        P5e["save_examples"] = False
        #########################################
        P5e["dur"] = np.repeat(dur, P5e["TotalN"])
        P5e["dt"] = np.repeat(dt, P5e["TotalN"])
        P5e["L"] = np.repeat(50.0, P5e["TotalN"])
        P5e["D"] = np.tile(3.0, P5e["TotalN"])
        P5e["G_EoH"] = np.repeat(0.0, P5e["TotalN"])
        P5e["Stochastic_EoH"] = np.repeat(False, P5e["TotalN"])
        P5e["Nsyn"] = np.repeat(32, P5e["TotalN"])
        P5e["gsyn"] = np.repeat(0.016, P5e["TotalN"])
        P5e["constR"] = np.repeat(False, P5e["TotalN"])
        P5e["dB"] = np.repeat(60.0, P5e["TotalN"])
        P5e["IPI"] = np.repeat(5.0, P5e["TotalN"])
        P5e["freq"] = np.repeat(200.0, P5e["TotalN"])
        P5e["cf"] = P5e["freq"]
        P5e["tau_dendsyn"] = np.repeat(2.0, P5e["TotalN"])
        P5e["gLEAK_d"] = np.repeat(0.001, P5e["TotalN"])
        P5e["gKLT_s"] = np.repeat(0.017, P5e["TotalN"])
        P5e["gIH_s"] = np.repeat(0.002, P5e["TotalN"])

        # Now define the two variable parameters. The repeated = y, the tiled = x!!
        allKLT = np.linspace(1.0e-9, 0.02, P5e["N"])
        allIH = np.linspace(1.0e-9, 0.01, P5e["N"])
        P5e["gKLT_d"] = np.repeat(allKLT, P5e["N"])
        P5e["gIH_d"] = np.tile(allIH, P5e["N"])

        # make go!
        output5e = []
        output5e.append(myMPhandler(P5e))
        output5e = np.array(output5e)
        np.save("./results/Fig5e.npy", output5e, allow_pickle=True)
        np.save("./results/Fig5e_P.npy", P5e, allow_pickle=True)

    ###Plot results
    fhandle = plotres_2d(
        output=output5e,
        P=P5e,
        x=np.unique(P5e["gIH_d"]),
        y=np.unique(P5e["gKLT_d"]),
        xlabs=r"gH ($\mu$S)",
        ylabs=r"gKLT ($\mu$S)",
    )

    pp = PdfPages("./figs/Fig5e.pdf")
    pp.savefig(bbox_inches="tight")
    pp.close()
