#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created 2019-07-11
Modified 2020-09-09 upon reviewer comments
Creates raw figures that were used to prepare 
    figure 1 of the Koert & Kuenzel Paper.
@author: Thomas Kuenzel (kuenzel{at}bio2.rwth-aachen.de)
"""

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import ScalarFormatter, AutoMinorLocator, MaxNLocator
import numpy as np
import smallexc.semodels as sem
import smallexc.sehelper as seh
import time

plt.rc("font", family="serif", serif="times")
plt.rcParams["pdf.fonttype"] = 42
plt.rc("text", usetex=True)
plt.rc("xtick", labelsize="x-small")
plt.rc("ytick", labelsize="x-small")
plt.rc("axes", labelsize="small")


def mycol(irep, nrep):
    C = ((1.0 / nrep) * irep, 0, (1 - ((1.0 / nrep) * irep)))
    return C


def runmain():
    simdur = 150.0
    dt = 0.01
    N = 32
    freq = 125
    SO = seh.create_sound(
        fs=100e3, freq=freq, duration=(simdur * 0.8) / 1000.0, dbspl=65.0
    )
    SP = seh.create_spk(SO, fs=100e3, N=N, cf=2.0 * freq, seed=65562)
    S = np.array(SP["spikes"] * 1000.0)
    R32 = sem.SE_BS(
        S, Simdur=simdur, dt=dt, G_EoH=0.0, StochasticEoH=False, N=N, G_Dend=0.016
    )
    R1 = sem.SE_BS(
        S, Simdur=simdur, dt=dt, G_EoH=0.0, StochasticEoH=False, N=1, G_Dend=0.0005
    )
    T = np.linspace(0, simdur - dt, int(round(simdur / dt)))
    sT = np.linspace(0, (simdur * 0.8) - dt, int(round((simdur * 0.8) / dt)))

    fwidth = 8.5
    fh = plt.figure(figsize=(fwidth / 2.54, (fwidth / 1.618) / 2.54))
    sp1 = fh.add_subplot(232)
    sp1.plot(T, R1["Vm"], linewidth=0.5)
    sp1.plot(sT, (SO * 6) - 65.0, color="r", linewidth=0.5)
    sp1.set_ylim((-66.0, -58.0))
    sp1.set_ylabel("Vm (mV)")
    sp1.set_xticks((0, 50, 100, 150))
    sp1.set_yticks((-65, -63, -61, -59))
    sp1.spines["top"].set_visible(False)
    sp1.spines["right"].set_visible(False)

    sp2 = fh.add_subplot(235)
    sp2.plot(T, R32["Vm"], linewidth=0.5)
    sp2.plot(sT, (SO * 6) - 65.0, color="r", linewidth=0.5)
    sp2.set_ylim((-66.0, -58.0))
    sp2.set_ylabel("Vm (mV)")
    sp2.set_xlabel("Time (ms)")
    sp2.set_xticks((0, 50, 100, 150))
    sp2.set_yticks((-65, -63, -61, -59))
    sp2.spines["top"].set_visible(False)
    sp2.spines["right"].set_visible(False)

    sp3 = fh.add_subplot(233)
    sp3.plot(S[0], np.ones(len(S[0])) * 0, "k.", markersize=1)
    sp3.plot(sT, ((SO * 3.0) - 0.5), color="r", linewidth=0.5)
    sp3.set_xlim((0, simdur))
    sp3.set_ylim((-1.1, N / 10.0))
    sp3.axes.get_yaxis().set_visible(False)
    sp3.set_xticks((0, 50, 100, 150))
    sp3.spines["top"].set_visible(False)
    sp3.spines["right"].set_visible(False)

    sp4 = fh.add_subplot(236)
    for s in range(N):
        sp4.plot(S[s], np.ones(len(S[s])) * (s / 10.0), "k.", markersize=1)
    sp4.set_xlim((0, simdur))
    sp4.plot(sT, ((SO * 3.0) - 0.5), color="r", linewidth=0.5)
    sp4.set_ylim((-1.1, N / 10.0))
    sp4.set_xlabel("Time (ms)")
    sp4.axes.get_yaxis().set_visible(False)
    sp4.set_xticks((0, 50, 100, 150))
    sp4.spines["top"].set_visible(False)
    sp4.spines["right"].set_visible(False)

    plt.tight_layout()
    # plt.show()
    pp = PdfPages("./figs/Figure1_raw.pdf")
    pp.savefig(dpi=600)
    pp.close()


if __name__ == "__main__":
    runmain()
