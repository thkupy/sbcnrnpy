#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 21 Aug 2019
Modified September 2020 upon reviewer comments at JNP.
This file contains code necessary to generate Figure 6B-E.
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
from scipy import stats

plt.rc("font", family="serif", serif="times")
plt.rcParams["pdf.fonttype"] = 42
plt.rc("text", usetex=True)
plt.rc("xtick", labelsize="x-small")
plt.rc("ytick", labelsize="x-small")
plt.rc("axes", labelsize="small")


def runmain():
    rep = 10
    simdur = 1000.0
    dt = 0.01
    N = 33
    freq = 1251.0  # 1251.0
    cf = 1001.0  # 1501.0
    db = 70
    SO = seh.create_sound(
        fs=100e3, freq=freq, duration=(simdur * 0.6) / 1000.0, dbspl=db
    )
    Silence = seh.create_sound(
        fs=100e3, freq=freq, duration=(simdur * 0.4) / 1000.0, dbspl=-10
    )
    R32 = []
    R0 = []
    Ev0 = []
    Fail0 = []
    Ev32 = []
    Fail32 = []
    APr0 = []
    APr32 = []
    allVS0 = []
    allVS32 = []
    allPhi0 = []
    allPhi32 = []
    #
    SIN = []
    DIN = []
    SR0 = []
    DR0 = []
    SF0 = []
    DF0 = []
    SR32 = []
    DR32 = []
    SF32 = []
    DF32 = []
    for thisrep in range(rep):
        SP = seh.create_spk(
            np.concatenate((SO, Silence)),
            fs=100e3,
            N=N,
            cf=cf,
            seed=65562 + thisrep,
            anf_num=(N, 0, 0),
        )
        S = np.array(SP["spikes"] * 1000.0)
        R32.append(
            sem.SE_BS(
                S,
                Simdur=simdur,
                dt=dt,
                G_EoH=0.055,
                StochasticEoH=True,
                N=N,
                G_Dend=0.016,
                EoHseed=8347236 + thisrep,
            )
        )
        R0.append(
            sem.SE_BS(
                S,
                Simdur=simdur,
                dt=dt,
                G_EoH=0.055,
                StochasticEoH=True,
                N=N,
                G_Dend=0.0,
                EoHseed=8347236 + thisrep,
            )
        )
        Ev0.append(seh.SimpleDetectAP(R0[-1]["Vm"], thr=-100, dt=dt, LM=-20, RM=10))
        Fail0.append(seh.SimpleGetFails(S[-1], Ev0[-1], R0[-1], simdur, dt))
        Ev32.append(seh.SimpleDetectAP(R32[-1]["Vm"], thr=-100, dt=dt, LM=-20, RM=10))
        Fail32.append(seh.SimpleGetFails(S[-1], Ev32[-1], R32[-1], simdur, dt))
        APr0.append(len(Ev0[-1]["PeakT"]) / (simdur / 1000.0))
        APr32.append(len(Ev32[-1]["PeakT"]) / (simdur / 1000.0))
        #
        tVS0, tphi0, tRay0, tphases0 = seh.vectorstrength(
            Ev0[-1]["PeakT"], freq, [0.0, simdur * 0.6]
        )
        tVS32, tphi32, tRay32, tphases32 = seh.vectorstrength(
            Ev32[-1]["PeakT"], freq, [0.0, simdur * 0.6]
        )
        allVS0.append(tVS0)
        allVS32.append(tVS32)
        allPhi0.append(tphi0)
        allPhi32.append(tphi32)
        #
        # generate more results for writing...
        stimlen = simdur * 0.6
        stimlen_s = (simdur * 0.6) / 1000.0
        spontlen_s = (simdur - stimlen) / 1000.0
        myEv0 = np.array(Ev0[-1]["PeakT"])
        myFail0 = np.array(Fail0[-1]["PeakT"])
        myEv32 = np.array(Ev32[-1]["PeakT"])
        myFail32 = np.array(Fail32[-1]["PeakT"])
        SIN.append(np.sum(S[-1] > stimlen) / spontlen_s)
        DIN.append(np.sum(S[-1] < stimlen) / stimlen_s)
        SR0.append(np.sum(myEv0 > stimlen) / spontlen_s)
        DR0.append(np.sum(myEv0 < stimlen) / stimlen_s)
        SF0.append(np.sum(myFail0 > stimlen) / spontlen_s)
        DF0.append(np.sum(myFail0 < stimlen) / stimlen_s)
        SR32.append(np.sum(myEv32 > stimlen) / spontlen_s)
        DR32.append(np.sum(myEv32 < stimlen) / stimlen_s)
        SF32.append(np.sum(myFail32 > stimlen) / spontlen_s)
        DF32.append(np.sum(myFail32 < stimlen) / stimlen_s)
    SIN = np.array(SIN)
    DIN = np.array(DIN)
    SR0 = np.array(SR0)
    DR0 = np.array(DR0)
    SF0 = np.array(SF0)
    DF0 = np.array(DF0)
    SR32 = np.array(SR32)
    DR32 = np.array(DR32)
    SF32 = np.array(SF32)
    DF32 = np.array(DF32)
    # Pool spiketimes and construct smooth cyclehist
    Spt0 = np.array(())
    for s0 in Ev0:
        Spt0 = np.append(Spt0, s0["PeakT"])
    Spt32 = np.array(())
    for s32 in Ev32:
        Spt32 = np.append(Spt32, s32["PeakT"])
    VS0, phi0, Ray0, phases0 = seh.vectorstrength(Spt0, freq, [0.0, simdur * 0.6])
    VS32, phi32, Ray32, phases32 = seh.vectorstrength(Spt32, freq, [0.0, simdur * 0.6])
    T = np.linspace(0, simdur - dt, int(round(simdur / dt)))
    sT = np.linspace(0, (simdur * 0.8) - dt, int(round((simdur * 0.8) / dt)))
    print("######Statistical comparison#########")
    print("Total number of AP: " + str(len(Spt0)) + " vs. " + str(len(Spt32)))
    outtext1 = str(np.mean(APr0)) + "+/-" + str(round(np.std(APr0), 2))
    outtext2 = str(np.mean(APr32)) + "+/-" + str(round(np.std(APr32), 2))
    print("Mean AP rate: " + outtext1 + " vs. " + outtext2)
    Tval, p = stats.ttest_rel(APr0, APr32)
    print("T=" + str(round(Tval, 3)) + ", p=" + str(p) + ", N=" + str(rep))
    outtext4 = str(np.mean(allVS0)) + "+/-" + str(round(np.std(allVS0), 2))
    outtext5 = str(np.mean(allVS32)) + "+/-" + str(round(np.std(allVS32), 2))
    print("Mean VS: " + outtext4 + " vs. " + outtext5)
    Tval, p = stats.ttest_rel(allVS0, allVS32)
    print("T=" + str(round(Tval, 3)) + ", p=" + str(p) + ", N=" + str(rep))
    outtext6 = str(np.mean(allPhi0)) + "+/-" + str(round(np.std(allPhi0), 2))
    outtext7 = str(np.mean(allPhi32)) + "+/-" + str(round(np.std(allPhi32), 2))
    print("Mean Phi: " + outtext6 + " vs. " + outtext7)
    Tval, p = stats.ttest_rel(allPhi0, allPhi32)
    print("T=" + str(round(Tval, 3)) + ", p=" + str(p) + ", N=" + str(rep))
    #
    print(
        "Driven input rate: "
        + str(round(np.mean(DIN), 2))
        + "+/-"
        + str(round(np.std(DIN), 2))
    )
    outtext1 = str(round(np.mean(DR0), 2)) + "+/-" + str(round(np.std(DR0), 2))
    outtext2 = str(round(np.mean(DR32), 2)) + "+/-" + str(round(np.std(DR32), 2))
    print("Mean driven AP rate: " + outtext1 + " vs. " + outtext2)
    Tval, p = stats.ttest_rel(DR0, DR32)
    print("T=" + str(round(Tval, 3)) + ", p=" + str(p) + ", N=" + str(rep))
    #
    print(
        "Spont input rate: "
        + str(round(np.mean(SIN), 2))
        + "+/-"
        + str(round(np.std(SIN), 2))
    )
    outtext1 = str(round(np.mean(SR0), 2)) + "+/-" + str(round(np.std(SR0), 2))
    outtext2 = str(round(np.mean(SR32), 2)) + "+/-" + str(round(np.std(SR32), 2))
    print("Mean spont AP rate: " + outtext1 + " vs. " + outtext2)
    Tval, p = stats.ttest_rel(SR0, SR32)
    print("T=" + str(round(Tval, 3)) + ", p=" + str(p) + ", N=" + str(rep))
    #
    outtext1 = str(round(np.mean(DF0), 2)) + "+/-" + str(round(np.std(DF0), 2))
    outtext2 = str(round(np.mean(DF32), 2)) + "+/-" + str(round(np.std(DF32), 2))
    print("Mean driven Fail rate: " + outtext1 + " vs. " + outtext2)
    Tval, p = stats.ttest_rel(DF0, DF32)
    print("T=" + str(round(Tval, 3)) + ", p=" + str(p) + ", N=" + str(rep))
    #
    outtext1 = str(round(np.mean(SF0), 2)) + "+/-" + str(round(np.std(SF0), 2))
    outtext2 = str(round(np.mean(SF32), 2)) + "+/-" + str(round(np.std(SF32), 2))
    print("Mean spont Fail rate: " + outtext1 + " vs. " + outtext2)
    Tval, p = stats.ttest_rel(SF0, SF32)
    print("T=" + str(round(Tval, 3)) + ", p=" + str(p) + ", N=" + str(rep))
    print("Exact values for traces shown:")
    print("Driven AP rate: " + str(DR0[0]) + " vs. " + str(DR32[0]))
    print("Driven Fail rate: " + str(DF0[0]) + " vs. " + str(DF32[0]))
    print("Spont AP rate: " + str(SR0[0]) + " vs. " + str(SR32[0]))
    print("Spont Fail rate: " + str(SF0[0]) + " vs. " + str(SF32[0]))
    # draw panels
    fwidth = 8.5
    fh = plt.figure(figsize=(fwidth / 2.54, fwidth / 2.54))
    fh.subplots_adjust(
        left=0.15, bottom=0.125, right=0.95, top=0.95, wspace=0.55, hspace=0.55
    )
    sp1 = fh.add_subplot(2, 3, (1, 2))
    sp1.plot(T, R0[0]["Vm"], "k-", linewidth=0.5)
    sp1.plot(Ev0[0]["PeakT"], Ev0[0]["PeakV"], "g.", markersize=2)
    sp1.plot(Fail0[0]["PeakT"], Fail0[0]["PeakV"], "r.", markersize=2)
    sp1.plot(S[-1], np.zeros((len(S[-1]),)) + 20.0, "b.", markersize=2)
    sp1.set_ylabel("Vm (mV)")
    sp1.set_xlabel("Time (ms)")
    sp1.spines["top"].set_visible(False)
    sp1.spines["right"].set_visible(False)
    sp1.plot(T, (np.concatenate((SO, Silence)) * 75) - 75.0, color="g", linewidth=0.5)
    #
    sp1b = fh.add_subplot(2, 3, 3)
    H0 = np.histogram(phases0, np.linspace(0, 1, 33))
    sp1b.plot(H0[1][1:], H0[0], "k")
    sp1b.set_ylabel("AP/bin")
    sp1b.set_xlabel("Phase (cycles)")
    sp1b.set_xticks((0, 0.5, 1))
    sp1b.set_yticks((0, 60, 120))
    sp1b.plot((phi0, phi0), (0, sp1b.get_ylim()[1]), "r-")
    annotext = r"VS=" + str(round(VS0, 2)) + "\n" + r"$\phi$=" + str(round(phi0, 2))
    sp1b.annotate(annotext, xy=(phi0, sp1b.get_ylim()[1]), fontsize=5)
    sp1b.spines["top"].set_visible(False)
    sp1b.spines["right"].set_visible(False)
    #
    sp2 = fh.add_subplot(2, 3, (4, 5))
    sp2.plot(T, R32[0]["Vm"], "k-", linewidth=0.5)
    sp2.plot(Ev32[0]["PeakT"], Ev32[0]["PeakV"], "g.", markersize=2)
    sp2.plot(Fail32[0]["PeakT"], Fail32[0]["PeakV"], "r.", markersize=2)
    sp2.plot(S[-1], np.zeros((len(S[-1]),)) + 20.0, "b.", markersize=2)
    sp2.plot(T, (np.concatenate((SO, Silence)) * 75) - 75.0, color="g", linewidth=0.5)
    sp2.set_ylabel("Vm (mV)")
    sp2.set_xlabel("Time (ms)")
    sp2.spines["top"].set_visible(False)
    sp2.spines["right"].set_visible(False)
    #
    sp2b = fh.add_subplot(2, 3, 6)
    H32 = np.histogram(phases32, np.linspace(0, 1, 33))
    sp2b.plot(H32[1][1:], H32[0], "k")
    sp2b.set_ylabel("AP/bin")
    sp2b.set_xlabel("Phase (cycles)")
    sp2b.set_xticks((0, 0.5, 1))
    sp2b.set_yticks((0, 60, 120))
    sp2b.plot((phi32, phi32), (0, sp2b.get_ylim()[1]), "r-")
    annotext = r"VS=" + str(round(VS32, 2)) + "\n" + r"$\phi$=" + str(round(phi32, 2))
    sp2b.annotate(annotext, xy=(phi32, sp2b.get_ylim()[1]), fontsize=5)
    sp2b.spines["top"].set_visible(False)
    sp2b.spines["right"].set_visible(False)
    #
    pp = PdfPages("./figs/Figure6_raw.pdf")
    pp.savefig(dpi=600)
    pp.close()


if __name__ == "__main__":
    runmain()
