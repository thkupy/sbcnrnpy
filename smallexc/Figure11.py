#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This file contains all code necessary to generate Figure 11.
-->Realistic dendrite frequency response areas, VS and preferred phase...
Code modified after original code by Elisabeth Koert.
@author: Thomas Kuenzel (kuenzel{at}bio2.rwth-aachen.de)
"""
##TK add my dirs to path
import sys
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import semodels as sem
import sehelper as seh
import multiprocessing
import os
from functools import partial


def runonce(x, P):
    ####hard coded settings, rarely change these:
    APdetect_thr = -75
    ####GENERATE SPIKE INPUTS (the spiketimes are always the same to improve comparability)
    thisdur = (P["dur"][x] - 5.0) / 1000.0
    N = P["Nsyn"][x] + 1
    anf_num = (N, 0, 0)
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
        SO, fs=100e3, N=N, cf=P["cf"][x], seed=thisZseed, anf_num=anf_num
    )  # generate enough spiketrains
    S = np.array(SP["spikes"] * 1000.0)
    R0 = sem.SE_3D(
        S,
        Simdur=P["dur"][x],
        dt=P["dt"][x],
        G_EoH=P["G_EoH"][x],
        StochasticEoH=P["Stochastic_EoH"][x],
        N=P["Nsyn"][x],
        G_Dend=P["gsyn"][x] * 0.0,  # set dendritic synapses to zero
        tau_Dend=P["tau_dendsyn"][x],
        gKLT_d=P["gKLT_d"][x],
        gIH_d=P["gIH_d"][x],
        gLEAK_d=P["gLEAK_d"][x],
        cell=P["cell"][x],
        EoHseed=thisEoHseed,
        somaticsynapsesonly=P["synapselocationtype"][x],
    )
    #
    RD = sem.SE_3D(
        S,
        Simdur=P["dur"][x],
        dt=P["dt"][x],
        G_EoH=P["G_EoH"][x],
        StochasticEoH=P["Stochastic_EoH"][x],
        N=P["Nsyn"][x],
        G_Dend=P["gsyn"][x],
        tau_Dend=P["tau_dendsyn"][x],
        gKLT_d=P["gKLT_d"][x],
        gIH_d=P["gIH_d"][x],
        gLEAK_d=P["gLEAK_d"][x],
        cell=P["cell"][x],
        EoHseed=thisEoHseed,
        somaticsynapsesonly=P["synapselocationtype"][x],
    )
    Ev0 = seh.SimpleDetectAP(R0["Vm"], thr=APdetect_thr, dt=P["dt"][x], LM=-20, RM=10)
    EvD = seh.SimpleDetectAP(RD["Vm"], thr=APdetect_thr, dt=P["dt"][x], LM=-20, RM=10)
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
    poolmess = partial(runonce, P=P)
    if P["mp"]:
        with multiprocessing.Pool(P["cores"]) as activepool:
            r = activepool.map(poolmess, P["Number"])
    else:  # for debug
        r = list(map(poolmess, P["Number"]))  # for debug
    return r


def plotres(output, P, x, y, xlabs, ylabs):
    from matplotlib.ticker import ScalarFormatter, AutoMinorLocator, MaxNLocator
    from matplotlib import cm
    from matplotlib.colors import ListedColormap
    import scipy.ndimage as ndimage

    plt.rc("font", family="serif", serif="Times")
    plt.rc("text", usetex=True)
    plt.rc("xtick", labelsize=7)
    plt.rc("ytick", labelsize=7)
    plt.rc("axes", labelsize=8)
    titsize = 9
    mycmap = "bone"

    r_crit = 0.001  # 0.001
    outVS0 = output[0][:, 0]
    outPH0 = output[0][:, 1]
    outRC0 = output[0][:, 2]
    outAP0 = output[0][:, 3]
    outVSD = output[0][:, 4]
    outPHD = output[0][:, 5]
    outRCD = output[0][:, 6]
    outAPD = output[0][:, 7]
    #
    outVS0nodiff = outVS0.copy()  # for the somaticsynapse condition only
    outVS0nodiff[outRC0 > r_crit] = 0.0  # for the somaticsynapse condition only
    outPH0nodiff = outPH0.copy()  # for the somaticsynapse condition only
    outPH0nodiff[outRC0 > r_crit] = 0.0  # for the somaticsynapse condition only
    VS0nodiff = np.reshape(
        outVS0nodiff, (P["N"], P["N"])
    )  # somaticsynapse condition only
    PH0nodiff = np.reshape(
        outPH0nodiff, (P["N"], P["N"])
    )  # somaticsynapse condition only
    #
    outVS0[
        np.logical_or(outRC0 > r_crit, outRCD > r_crit)
    ] = 0.0  # set VS of conditions to zero that fail the rayleigh test
    outPH0[
        np.logical_or(outRC0 > r_crit, outRCD > r_crit)
    ] = 0.0  # set phi of conditions to zero that fail the rayleigh test
    outVSD[
        np.logical_or(outRC0 > r_crit, outRCD > r_crit)
    ] = 0.0  # set VS of conditions to zero that fail the rayleigh test
    outPHD[
        np.logical_or(outRC0 > r_crit, outRCD > r_crit)
    ] = 0.0  # set VS of conditions to zero that fail the rayleigh test
    pshift = seh.get_anglediff(outPH0, outPHD)  # new test
    pshift = np.reshape(pshift, (P["N"], P["N"]))  # new test
    APdiff = outAPD - outAP0
    APdiff = np.reshape(APdiff, (P["N"], P["N"]))  # new test
    VSdiff = outVSD - outVS0
    VSdiff = np.reshape(VSdiff, (P["N"], P["N"]))  # new test
    VS0 = np.reshape(outVS0, (P["N"], P["N"]))
    PH0 = np.reshape(outPH0, (P["N"], P["N"]))
    RC0 = np.reshape(outRC0, (P["N"], P["N"]))
    AP0 = np.reshape(outAP0, (P["N"], P["N"]))
    VSD = np.reshape(outVSD, (P["N"], P["N"]))
    PHD = np.reshape(outPHD, (P["N"], P["N"]))
    RCD = np.reshape(outRCD, (P["N"], P["N"]))
    APD = np.reshape(outAPD, (P["N"], P["N"]))

    filtsig = 0.5
    AP0 = ndimage.gaussian_filter(AP0, sigma=filtsig, order=0)
    APD = ndimage.gaussian_filter(APD, sigma=filtsig, order=0)
    APdiff = ndimage.gaussian_filter(APdiff, sigma=filtsig, order=0)
    VS0 = ndimage.gaussian_filter(VS0, sigma=filtsig, order=0)
    VSD = ndimage.gaussian_filter(VSD, sigma=filtsig, order=0)
    VSdiff = ndimage.gaussian_filter(VSdiff, sigma=filtsig, order=0)

    fwidth = 17.0  # cm
    fhandle1 = plt.figure(figsize=(fwidth / 2.54, (fwidth / 1.618) / 2.54), dpi=600)

    # row 2: VS
    cbkwargs = {"format": "%.1f"}
    absomin = np.min(np.array([np.min(VS0), np.min(VSD)]))
    absomax = np.max(np.array([np.max(VS0), np.max(VSD)]))
    conts = np.linspace(absomin, absomax, 21)
    ax1 = plt.subplot(334)
    if P["synapselocationtype"][0] == 1:
        VS0nodiff = ndimage.gaussian_filter(VS0nodiff, sigma=filtsig, order=0)
        CS1b = plt.contour(
            x, y, VS0nodiff, (0.6,), colors="yellow"
        )  # repeated=y,tiled=x
        CS1 = plt.contourf(
            x,
            y,
            VS0nodiff,
            conts,
            cmap=mycmap,
        )  # repeated = y, tiled = x!!
    else:
        CS1b = plt.contour(x, y, VS0, (0.6,), colors="yellow")  # repeated=y, tiled=x!!
        CS1 = plt.contourf(
            x,
            y,
            VS0,
            conts,
            cmap=mycmap,
        )  # repeated = y, tiled = x!!
    plt.title(r"VS, $g_{dendsyn} = 0nS$", fontsize=titsize)
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    ax1.set_xscale("log")
    # now set good xticks (remember, we are from 125 to P["cf"] * 2.501 Hz)
    xl = ax1.get_xlim()
    if xl[1] > 5000:
        goodticks = [500, 1000, 2000, 5000]
    elif (xl[1] < 5000) and (xl[1] > 2000):
        goodticks = [200, 500, 1000, 2000]
    else:
        goodticks = [200, 500, 1000]
    ax1.set_xticks(goodticks)
    ax1.set_yticks([0, 20, 40, 60, 80])
    ax1.get_xaxis().set_major_formatter(ScalarFormatter())
    ax1.get_yaxis().set_major_formatter(ScalarFormatter())
    cbar1 = plt.colorbar(CS1, use_gridspec=True, **cbkwargs)
    cbar1.ax.set_ylabel("Vectorstrength")
    tl = MaxNLocator(nbins=5)
    cbar1.locator = tl
    cbar1.update_ticks()
    ax1.annotate(
        "",
        xy=(P["cf"][0], 10.0),
        xytext=(P["cf"][0], -2.0),
        arrowprops=dict(color="white", arrowstyle="wedge"),
    )
    #
    ax2 = plt.subplot(335)
    CS2 = plt.contourf(x, y, VSD, conts, cmap=mycmap)  # repeated = y, tiled = x!!
    CS2b = plt.contour(x, y, VSD, (0.6,), colors="yellow")  # repeated = y, tiled = x!!
    # CS2=plt.pcolormesh(x,y,VSD,shading='gouraud',cmap=mycmap)#repeated = y, tiled = x!!
    plt.title(r"VS, $g_{dendsyn} = 16nS$", fontsize=titsize)
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    ax2.set_xscale("log")
    ax2.set_xticks(goodticks)
    ax2.set_yticks([0, 20, 40, 60, 80])
    ax2.get_xaxis().set_major_formatter(ScalarFormatter())
    ax2.get_yaxis().set_major_formatter(ScalarFormatter())
    cbar2 = plt.colorbar(CS2, use_gridspec=True, **cbkwargs)
    cbar2.ax.set_ylabel("Vectorstrength")
    tl = MaxNLocator(nbins=5)
    cbar2.locator = tl
    cbar2.update_ticks()
    ax2.annotate(
        "",
        xy=(P["cf"][0], 10.0),
        xytext=(P["cf"][0], -2.0),
        arrowprops=dict(color="white", arrowstyle="wedge"),
    )
    #
    ax3 = plt.subplot(336)
    CS3 = plt.contourf(
        x,
        y,
        VSdiff,
        np.linspace(-0.1, 0.1, 21),
        cmap="seismic",
        extend="both",
    )  # repeated = y, tiled = x!!
    plt.title("VS change", fontsize=titsize)
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    plt.xscale("log")
    ax3.set_xscale("log")
    ax3.set_xticks(goodticks)
    ax3.set_yticks([0, 20, 40, 60, 80])
    ax3.get_xaxis().set_major_formatter(ScalarFormatter())
    ax3.get_yaxis().set_major_formatter(ScalarFormatter())
    cbar3 = plt.colorbar(CS3, use_gridspec=True)
    cbar3.ax.set_ylabel(r"$\Delta$ Vectorstrength")
    tl = MaxNLocator(nbins=5)
    cbar3.locator = tl
    cbar3.update_ticks()
    ax3.annotate(
        "",
        xy=(P["cf"][0], 10.0),
        xytext=(P["cf"][0], -2.0),
        arrowprops=dict(color="black", arrowstyle="wedge"),
    )
    #
    # create cyclic colormap
    cmapN = 256
    linmap = cm.get_cmap("bone")
    cmapvals = np.ones((cmapN, 4))
    cmapvals[0 : int(cmapN / 2), :] = linmap(np.linspace(0.0, 1.0, cmapN / 2))
    cmapvals[int(cmapN / 2) :, :] = linmap(np.linspace(1.0, 0.0, cmapN / 2))
    cycliccmp = ListedColormap(cmapvals)
    #
    # row 3: preferred phase
    absomin = np.min(np.array([np.min(PH0), np.min(PHD)]))
    absomax = np.max(np.array([np.max(PH0), np.max(PHD)]))
    conts = np.linspace(absomin, absomax, 21)
    ax4 = plt.subplot(337)
    if P["synapselocationtype"][0] == 1:
        CS4 = plt.contourf(
            x, y, PH0nodiff, conts, cmap=cycliccmp
        )  # repeated=y, tiled=x!
    else:
        CS4 = plt.contourf(x, y, PH0, conts, cmap=cycliccmp)  # repeated =y, tiled =x!!
    plt.title(r"$\varphi$ , $g_{dendsyn} = 0nS$", fontsize=titsize)
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    plt.xscale("log")
    ax4.set_xscale("log")
    ax4.set_xticks(goodticks)
    ax4.set_yticks([0, 20, 40, 60, 80])
    ax4.get_xaxis().set_major_formatter(ScalarFormatter())
    ax4.get_yaxis().set_major_formatter(ScalarFormatter())
    cbar4 = plt.colorbar(CS4, use_gridspec=True, **cbkwargs)
    cbar4.ax.set_ylabel(r"$\varphi$ (cycles)")
    tl = MaxNLocator(nbins=5)
    cbar4.locator = tl
    cbar4.update_ticks()
    ax4.annotate(
        "",
        xy=(P["cf"][0], 10.0),
        xytext=(P["cf"][0], -2.0),
        arrowprops=dict(color="white", arrowstyle="wedge"),
    )
    #
    ax5 = plt.subplot(338)
    CS5 = plt.contourf(x, y, PHD, conts, cmap=cycliccmp)  # repeated = y, tiled = x!!
    # CS5=plt.pcolormesh(x,y,PHD,shading='gouraud',cmap=mycmap)#repeated = y, tiled = x!!
    plt.title(r"$\varphi$ , $g_{dendsyn} = 16nS$", fontsize=titsize)
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    plt.xscale("log")
    ax5.set_xscale("log")
    ax5.set_xticks(goodticks)
    ax5.set_yticks([0, 20, 40, 60, 80])
    ax5.get_xaxis().set_major_formatter(ScalarFormatter())
    ax5.get_yaxis().set_major_formatter(ScalarFormatter())
    cbar5 = plt.colorbar(CS5, use_gridspec=True, **cbkwargs)
    cbar5.ax.set_ylabel(r"$\phi$ (cycles)")
    tl = MaxNLocator(nbins=5)
    cbar5.locator = tl
    cbar5.update_ticks()
    ax5.annotate(
        "",
        xy=(P["cf"][0], 10.0),
        xytext=(P["cf"][0], -2.0),
        arrowprops=dict(color="white", arrowstyle="wedge"),
    )
    #
    ax6 = plt.subplot(339)
    # pshift = seh.get_anglediff(PH0,PHD)
    # print pshift
    CS6 = plt.contourf(
        x,
        y,
        pshift,
        np.linspace(-0.15, 0.15, 21),
        cmap="seismic",
        extend="both",
    )  # repeated = y, tiled = x!!
    plt.title("Phase shift", fontsize=titsize)
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    plt.xscale("log")
    ax6.set_xscale("log")
    ax6.set_xticks(goodticks)
    ax6.set_yticks([0, 20, 40, 60, 80])
    ax6.get_xaxis().set_major_formatter(ScalarFormatter())
    ax6.get_yaxis().set_major_formatter(ScalarFormatter())
    cbar6 = plt.colorbar(CS6, use_gridspec=True)
    cbar6.ax.set_ylabel(r"$\Delta \varphi$ (cycles)")
    tl = MaxNLocator(nbins=5)
    cbar6.locator = tl
    cbar6.update_ticks()
    ax6.annotate(
        "",
        xy=(P["cf"][0], 10.0),
        xytext=(P["cf"][0], -2.0),
        arrowprops=dict(color="black", arrowstyle="wedge"),
    )

    # row 1: AP rate
    cbkwargs = {"format": "%.0f"}
    absomin = np.min(np.array([np.min(AP0), np.min(APD)]))
    absomax = np.max(np.array([np.max(AP0), np.max(APD)]))
    conts = np.linspace(absomin, absomax, 21)
    ax7 = plt.subplot(331)
    CS7 = plt.contourf(x, y, AP0, conts, cmap=mycmap)  # repeated = y, tiled = x!!
    CS7b = plt.contour(
        x,
        y,
        AP0,
        (
            80.0,
            170.0,
        ),
        colors="yellow",
    )  # repeated = y, tiled = x!!
    # CS7=plt.pcolormesh(x,y,AP0,shading='gouraud',cmap=mycmap)#repeated = y, tiled = x!!
    plt.title("AP rate, $g_{dendsyn} = 0nS$", fontsize=titsize)
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    plt.xscale("log")
    ax7.set_xscale("log")
    ax7.set_xticks(goodticks)
    ax7.set_yticks([0, 20, 40, 60, 80])
    ax7.get_xaxis().set_major_formatter(ScalarFormatter())
    ax7.get_yaxis().set_major_formatter(ScalarFormatter())
    cbar7 = plt.colorbar(CS7, use_gridspec=True, **cbkwargs)
    cbar7.ax.set_ylabel("Rate (Hz)")
    tl = MaxNLocator(nbins=5)
    cbar7.locator = tl
    cbar7.update_ticks()
    ax7.annotate(
        "",
        xy=(P["cf"][0], 10.0),
        xytext=(P["cf"][0], -2.0),
        arrowprops=dict(color="white", arrowstyle="wedge"),
    )
    #
    ax8 = plt.subplot(332)
    CS8 = plt.contourf(x, y, APD, conts, cmap=mycmap)  # repeated = y, tiled = x!!
    CS8b = plt.contour(
        x,
        y,
        APD,
        (
            80.0,
            170.0,
        ),
        colors="yellow",
    )  # repeated = y, tiled = x!!
    # CS8=plt.pcolormesh(x,y,APD,shading='gouraud',cmap=mycmap)#repeated = y, tiled = x!!
    plt.title("AP rate, $g_{dendsyn} = 16nS$", fontsize=titsize)
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    plt.xscale("log")
    ax8.set_xscale("log")
    ax8.set_xticks(goodticks)
    ax8.set_yticks([0, 20, 40, 60, 80])
    ax8.get_xaxis().set_major_formatter(ScalarFormatter())
    ax8.get_yaxis().set_major_formatter(ScalarFormatter())
    cbar8 = plt.colorbar(CS8, use_gridspec=True, **cbkwargs)
    cbar8.ax.set_ylabel("Rate (Hz)")
    tl = MaxNLocator(nbins=5)
    cbar8.locator = tl
    cbar8.update_ticks()
    ax8.annotate(
        "",
        xy=(P["cf"][0], 10.0),
        xytext=(P["cf"][0], -2.0),
        arrowprops=dict(color="white", arrowstyle="wedge"),
    )
    #
    ax9 = plt.subplot(333)
    CS9 = plt.contourf(
        x,
        y,
        APdiff,
        np.linspace(-30, 30, 21),
        cmap="seismic",
        extend="both",
    )  # repeated = y, tiled = x!!
    plt.title("Rate change", fontsize=titsize)
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    plt.xscale("log")
    ax9.set_xscale("log")
    ax9.set_xticks(goodticks)
    ax9.set_yticks([0, 20, 40, 60, 80])
    ax9.get_xaxis().set_major_formatter(ScalarFormatter())
    ax9.get_yaxis().set_major_formatter(ScalarFormatter())
    cbar9 = plt.colorbar(CS9, use_gridspec=True)
    cbar9.ax.set_ylabel(r"$\Delta$ Rate (Hz)")
    tl = MaxNLocator(nbins=5)
    cbar9.locator = tl
    cbar9.update_ticks()
    ax9.annotate(
        "",
        xy=(P["cf"][0], 10.0),
        xytext=(P["cf"][0], -2.0),
        arrowprops=dict(color="black", arrowstyle="wedge"),
    )
    plt.tight_layout()
    #
    # Second Figure: Lineplot at one Level (Reviewer Request)
    # Axes are fixed so that different versions (500,1000,1500,2000 Hz) can be combined
    # into one plot later.
    x = np.log2(x / P["cf"][0])
    fhandle2 = plt.figure(figsize=(fwidth / 2.54, (fwidth / 1.618) / 2.54), dpi=600)
    ax20 = plt.subplot(3, 4, 4)
    ax20.plot(x, AP0[-2, :], "k")
    ax20.plot(x, APD[-2, :], "g")
    ax20.set_title("@" + str(round(y[-2])) + "dB SPL")
    ax20.set_xlim(left=-4, right=2)
    ax20.get_xaxis().set_major_formatter(ScalarFormatter())
    ax20.set_ylabel("AP Rate (Hz)")
    ax20.set_ylim((0, 250))
    ax21 = plt.subplot(3, 4, 8)
    if P["synapselocationtype"][0] == 1:
        ax21.plot(x, VS0nodiff[-2, :], "k")
    else:
        ax21.plot(x, VS0[-2, :], "k")
    ax21.plot(x, VSD[-2, :], "g")
    ax21.set_xlim(left=-4, right=2)
    ax21.get_xaxis().set_major_formatter(ScalarFormatter())
    ax21.set_ylabel("VS")
    ax21.set_ylim((0, 1))
    ax22 = plt.subplot(3, 4, 12)
    if P["synapselocationtype"][0] == 1:
        ax22.plot(x, PH0[-2, :], "k")
    else:
        ax22.plot(x, PH0[-2, :], "k")
    ax22.plot(x, PHD[-2, :], "g")
    ax22.set_xlim(left=-4, right=2)
    ax22.get_xaxis().set_major_formatter(ScalarFormatter())
    ax22.set_ylabel("Phase (cycles)")
    ax22.set_ylim((0, 1))
    #
    plt.tight_layout()
    return (fhandle1, fhandle2)


if __name__ == "__main__":
    # Checking for arguments and/or setting defaults
    addedargs = sys.argv[1:]  # fname, weload, nconds, ncores, cell, cf, synmodel, gsyn
    myargs = [1, 5, 4, 1, 1500.1, 0, 0.064]  # as above starting with weload
    for iarg, thisarg in enumerate(addedargs):
        myargs[iarg] = float(thisarg)
    weload = bool(myargs[0])
    nconds = int(myargs[1])
    ncores = int(myargs[2])
    thiscell = int(myargs[3])
    thiscf = myargs[4]
    usesomaticsynapsemodel = int(myargs[5])
    gsyn = myargs[6]
    if usesomaticsynapsemodel:
        filesuffix = "cell" + str(thiscell) + "_" + str(round(thiscf)) + "Hz_somsyn"
    else:
        filesuffix = "cell" + str(thiscell) + "_" + str(round(thiscf)) + "Hz"
    if not gsyn == 0.064:
        filesuffix = filesuffix + "_" + str(1000.0 * gsyn) + "nS"
    # debug
    print("Frequency: " + str(round(thiscf)))
    print("Synapse location type:" + str(usesomaticsynapsemodel))
    # debug
    if os.path.isfile("./results/Figure11_" + filesuffix + ".npy"):
        print("Data for Figure11 (@" + filesuffix + ") found... loading!")
        output = np.load("./results/Figure11_" + filesuffix + ".npy", allow_pickle=True)
        P = np.load("./results/Figure11_P_" + filesuffix + ".npy", allow_pickle=True)
        P = P.tolist()
    else:
        dur = 5000.0
        dt = 0.01
        output = []
        # Model Parameters (all in a linearly aranged fashion, so that a minimal
        # amount of programming is required to change the experiment).
        P = {}
        P["N"] = nconds
        P["cores"] = ncores
        P["TotalN"] = int(P["N"] ** 2)
        P["Number"] = list(range(P["TotalN"]))
        P["mp"] = True
        P["Zseed"] = 45453
        P["AdvZseed"] = True
        P["EoHseed"] = 34823
        P["AdvEoHseed"] = True
        #########################################
        P["cell"] = np.repeat(thiscell, P["TotalN"])  # which model neuron is used
        #########################################
        P["dur"] = np.repeat(dur, P["TotalN"])
        P["dt"] = np.repeat(dt, P["TotalN"])
        P["G_EoH"] = np.repeat(0.055, P["TotalN"])
        P["Stochastic_EoH"] = np.repeat(True, P["TotalN"])
        P["Nsyn"] = np.repeat(32, P["TotalN"])  # 32
        P["gsyn"] = np.repeat(0.064, P["TotalN"])  # 0.048
        P["tau_dendsyn"] = np.repeat(2.0, P["TotalN"])
        P["gLEAK_d"] = np.repeat(0.001, P["TotalN"])
        P["gKLT_d"] = np.repeat(0.0085, P["TotalN"])
        P["gIH_d"] = np.repeat(0.001, P["TotalN"])
        P["gKLT_s"] = np.repeat(0.017, P["TotalN"])
        P["gIH_s"] = np.repeat(0.002, P["TotalN"])
        # Now define the two variable parameters. The repeated = y, the tiled = x!!
        allfreq = np.round(np.geomspace(125.0, thiscf * 2.5, P["N"]))
        alllev = np.round(np.linspace(0.0, 80.0, P["N"]), 1)
        P["freq"] = np.tile(allfreq, P["N"])
        P["cf"] = np.repeat(thiscf, P["TotalN"])
        P["dB"] = np.repeat(alllev, P["N"])
        #
        # State whether the default model or the somatic synapse model is used
        P["synapselocationtype"] = np.repeat(usesomaticsynapsemodel, P["TotalN"])
        #
        # make go!
        output.append(myMPhandler(P))
        output = np.array(output)
        np.save("./results/Figure11_" + filesuffix + ".npy", output, allow_pickle=True)
        np.save("./results/Figure11_P_" + filesuffix + ".npy", P, allow_pickle=True)

    ###Plot results
    fhandle = plotres(
        output=output,
        P=P,
        x=np.unique(P["freq"]),
        y=np.unique(P["dB"]),
        xlabs="Frequency (Hz)",
        ylabs="Level (dB SPL)",
    )
    pp = PdfPages("./figs/Figure11_" + filesuffix + ".pdf")
    pp.savefig(fhandle[0])
    pp.savefig(fhandle[1])
    pp.close()
