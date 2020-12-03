#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This file contains all code necessary to generate the 3D surface plot included in
Figure 4. This was a reviewer request from the submission at JNP.

Created sept 2020
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


def runonce(x, P):
    # print("Round " + str(P4F["freq"][x]) + ", cond: " + str(x))
    thisdur = (P["dur"][x] - 5.0) / 1000.0
    SO = seh.create_sound(
        fs=100e3, freq=P["freq"][x], duration=thisdur, dbspl=P["dB"][x]
    )
    SP = seh.create_spk(SO, fs=100e3, N=P["Nsyn"][x], cf=P["cf"][x], seed=2121977)
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
    return [SAPr, Smeanv, Smodv]


def myMPhandler(P):
    # p = multiprocessing.Pool(P["cores"])
    poolmess = partial(runonce, P=P)
    if P["mp"]:
        with multiprocessing.Pool(P["cores"]) as activepool:
            r = activepool.map(poolmess, range(P["Number"]))
        # r = p.map(poolmess, P["Number"])
    else:  # for debug
        r = list(map(poolmess, P["Number"]))  # for debug
    return r


def plotres(R, P):
    import matplotlib.font_manager
    from matplotlib.ticker import ScalarFormatter, AutoMinorLocator, MaxNLocator
    import scipy.ndimage as ndimage
    from mpl_toolkits.mplot3d import Axes3D

    plt.rc("font", family="serif", serif="Times")
    plt.rc("text", usetex=True)
    plt.rc("xtick", labelsize=7)
    plt.rc("ytick", labelsize=7)
    plt.rc("axes", labelsize=8)
    nticks = 5
    titsize = 9
    mycmap = "bone"
    filtsig = 0.25
    fhandle, ax = plt.subplots(P4F["N3"], 3)  # just for debug
    x = np.unique(P4F["D"])
    y = np.unique(P4F["L"])
    m = []
    for nfreq, freq in enumerate(P4F["allFreq"]):
        V = ndimage.gaussian_filter(R["allmeanv"][:, :, nfreq], sigma=filtsig, order=0)
        M = ndimage.gaussian_filter(R["allmodv"][:, :, nfreq], sigma=filtsig, order=0)
        A = ndimage.gaussian_filter(R["allnap"][:, :, nfreq], sigma=filtsig, order=0)
        ax[nfreq, 0].contourf(x, y, V, 21, cmap=mycmap)  # repeated = y, tiled = x!!
        ax[nfreq, 0].set_title("@:" + str(freq) + "Hz")
        ax[nfreq, 0].yaxis.set_major_locator(MaxNLocator(nticks))
        ax[nfreq, 0].xaxis.set_major_locator(MaxNLocator(nticks))
        ax[nfreq, 1].contourf(x, y, M, 21, cmap=mycmap)  # repeated = y, tiled = x!!
        ax[nfreq, 1].yaxis.set_major_locator(MaxNLocator(nticks))
        ax[nfreq, 1].xaxis.set_major_locator(MaxNLocator(nticks))
        ax[nfreq, 2].contourf(x, y, A, 21, cmap=mycmap)  # repeated = y, tiled = x!!
        ax[nfreq, 2].yaxis.set_major_locator(MaxNLocator(nticks))
        ax[nfreq, 2].xaxis.set_major_locator(MaxNLocator(nticks))
        voltcont = 0.25  # this determines which modulation contour in mV is drawn
        m.append(
            np.take(y, np.argmin((R["allmodv"][:, :, nfreq] - voltcont) ** 2, axis=0))
        )
    plt.tight_layout()
    fhandle2 = plt.figure()
    ax = fhandle2.add_subplot(111, projection="3d")
    cycdur = 1000.0 / P4F["allFreq"]
    m = np.array(m)
    X = np.tile(x, P4F["N"])
    X = np.reshape(X, (P4F["N"], P4F["N"]))
    Y = np.repeat(cycdur, P4F["N"])
    Y = np.reshape(Y, (P4F["N"], P4F["N"]))
    ax.plot_surface(
        X,
        Y,
        m,
        edgecolor="k",
        linewidth=0.25,
        rstride=1,
        cstride=1,
        cmap="bone",
        antialiased=True,
        alpha=0.75,
    )
    ax.view_init(elev=45, azim=200)  # 45/300?
    ax.set_title(str(voltcont) + "mV modulation surface")
    ax.set_xlabel("Diameter (µm)")
    ax.set_ylabel("Cycle dur. (ms)")
    ax.set_zlabel("Length (Hz)")
    # ax.figure.canvas.draw()
    return fhandle, fhandle2


if __name__ == "__main__":
    conditions_3d = 15
    conditions_2d = 15
    cores = 20
    dur = 5000.0  # production = 5000.0
    dt = 0.01  # production = 0.01

    myylim1 = (-64.5, -60.5)
    myylim2 = (0.0, 5.5)
    ################################
    #####Fig2 Part B:         ######
    #####-----------------------####
    ##### Variables -> LxDxFreq ####
    ##### No Endbulb    ############
    ################################
    ####LOAD DATA (IF IT EXISTS)
    if os.path.isfile("./results/Figure4F.npy"):
        print("Data for Figure4F found... loading!")
        R = np.load("./results/Figure4F.npy", allow_pickle=True)
        R = R.tolist()
        P4F = np.load("./results/Figure4F_P.npy", allow_pickle=True)
        P4F = P4F.tolist()
    else:
        # Model Parameters (all in a linearly aranged fashion, so that a minimal
        # amount of programming is required to change the experiment).
        P4F = {}
        P4F["mp"] = True
        P4F["N"] = int(conditions_2d)
        P4F["N3"] = int(conditions_3d)
        P4F["cores"] = int(cores)  # just to be sure
        P4F["TotalN"] = int(P4F["N"] ** 2)
        P4F["Number"] = list(range(P4F["TotalN"]))
        P4F["save_examples"] = False
        #########################################
        P4F["dur"] = np.repeat(dur, P4F["TotalN"])
        P4F["dt"] = np.repeat(dt, P4F["TotalN"])
        P4F["G_EoH"] = np.repeat(0.0, P4F["TotalN"])
        P4F["Stochastic_EoH"] = np.repeat(False, P4F["TotalN"])
        P4F["Nsyn"] = np.repeat(32, P4F["TotalN"])
        P4F["gsyn"] = np.repeat(0.016, P4F["TotalN"])
        P4F["tau_dendsyn"] = np.repeat(2.0, P4F["TotalN"])
        P4F["constR"] = np.repeat(False, P4F["TotalN"])
        P4F["dB"] = np.repeat(60.0, P4F["TotalN"])
        P4F["gKLT_d"] = np.repeat(0.0085, P4F["TotalN"])
        P4F["gIH_d"] = np.repeat(0.001, P4F["TotalN"])
        P4F["gLEAK_d"] = np.repeat(0.001, P4F["TotalN"])
        P4F["gKLT_s"] = np.repeat(0.017, P4F["TotalN"])
        P4F["gIH_s"] = np.repeat(0.002, P4F["TotalN"])
        #
        # Now define the two variable parameters. The repeated = y, the tiled = x!!
        allL = np.linspace(10.0, 250.0, P4F["N"])
        allD = np.linspace(0.1, 10.0, P4F["N"])
        P4F["L"] = np.repeat(allL, P4F["N"])
        P4F["D"] = np.tile(allD, P4F["N"])
        # Now define the third variable: CF
        allFreq = np.round(np.linspace(125.0, 666.6, conditions_3d))
        P4F["allFreq"] = allFreq
        P4F["which3dcond"] = 0
        #
        allmodv = np.empty((conditions_2d, conditions_2d, conditions_3d), dtype=float)
        allmeanv = np.empty((conditions_2d, conditions_2d, conditions_3d), dtype=float)
        allnap = np.empty((conditions_2d, conditions_2d, conditions_3d), dtype=float)
        # Now run the model conditions_3d x conditions_2d x conditions_2d times
        for nfreq, thisfreq in enumerate(allFreq):
            output = []
            P4F["freq"] = np.repeat(thisfreq, P4F["TotalN"])
            P4F["IPI"] = np.repeat(1000.0 / thisfreq, P4F["TotalN"])
            P4F["cf"] = np.repeat(thisfreq, P4F["TotalN"])
            P4F["which3dcond"] = nfreq
            output.append(myMPhandler(P4F))
            output = np.array(output)
            SAPr = output[0][:, 0]
            Smeanv = output[0][:, 1]
            Smodv = output[0][:, 2]
            allnap[:, :, nfreq] = np.reshape(SAPr, (P4F["N"], P4F["N"]))
            allmeanv[:, :, nfreq] = np.reshape(Smeanv, (P4F["N"], P4F["N"]))
            allmodv[:, :, nfreq] = np.reshape(Smodv, (P4F["N"], P4F["N"]))
        R = {"allnap": allnap, "allmeanv": allmeanv, "allmodv": allmodv}
        np.save("./results/Figure4F.npy", R, allow_pickle=True)
        np.save("./results/Figure4F_P.npy", P4F, allow_pickle=True)

    ###Plot results
    fhandle, fhandle2 = plotres(
        R=R,
        P=P4F,
    )
    pp = PdfPages("./figs/Figure4F.pdf")
    pp.savefig(fhandle)
    pp.savefig(fhandle2)
    pp.close()
