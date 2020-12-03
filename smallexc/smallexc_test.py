#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Run this first after installing the smallexc model, if this works, all
other simulations should run fine.
Created on Thu Jun  4 10:33:53 2020
@author: Thomas Kuenzel (kuenzel{at}bio2.rwth-aachen.de)
"""

import matplotlib.pyplot as plt
import numpy as np
import semodels as sem
import sehelper as seh


def runmain():
    simdur = 400.0  # 1000
    dt = 0.01
    N = 32
    freq = 1450.0
    cf = 1501.0
    db = 70  # 70?
    SO = seh.create_sound(
        fs=100e3, freq=freq, duration=(simdur * 0.6) / 1000.0, dbspl=db
    )
    Silence = seh.create_sound(
        fs=100e3, freq=freq, duration=(simdur * 0.4) / 1000.0, dbspl=-10
    )
    SP = seh.create_spk(
        np.concatenate((SO, Silence)),
        fs=100e3,
        N=N + 1,
        cf=cf,
        seed=65562,
        anf_num=(N, 0, 0),
    )
    S = np.array(SP["spikes"] * 1000.0)
    R32Som = sem.SE_BS(
        S,
        Simdur=simdur,
        dt=dt,
        G_EoH=0.055,
        StochasticEoH=True,
        N=N,
        G_Dend=0.016,
        EoHseed=8347236,
        somaticsynapsesonly=1,
    )
    R32 = sem.SE_BS(
        S,
        Simdur=simdur,
        dt=dt,
        G_EoH=0.055,
        StochasticEoH=True,
        N=N,
        G_Dend=0.016,
        EoHseed=8347236,
    )
    R0 = sem.SE_BS(
        S,
        Simdur=simdur,
        dt=dt,
        G_EoH=0.055,
        StochasticEoH=True,
        N=N,
        G_Dend=0.0,
        EoHseed=8347236,
    )

    R3d_32 = sem.SE_3D(
        S,
        Simdur=simdur,
        dt=dt,
        G_EoH=0.055,
        StochasticEoH=True,
        N=N,
        G_Dend=0.064,
        EoHseed=8347236,
    )
    R3d_32som = sem.SE_3D(
        S,
        Simdur=simdur,
        dt=dt,
        G_EoH=0.055,
        StochasticEoH=True,
        N=N,
        G_Dend=0.064,
        EoHseed=8347236,
        somaticsynapsesonly=1,
    )
    R3d_0 = sem.SE_3D(
        S,
        Simdur=simdur,
        dt=dt,
        G_EoH=0.055,
        StochasticEoH=True,
        N=N,
        G_Dend=0.0,
        EoHseed=8347236,
    )
    T = np.linspace(0, simdur - dt, int(round(simdur / dt)))

    fh = plt.figure()
    sp1 = fh.add_subplot(4, 1, 1)
    sp1.set_title("Stylized Dendrite")
    sp1.plot(T, R0["Vm"], "b-", linewidth=0.5, label="0nS")
    sp1.plot(T, R32["Vm"], "g-", linewidth=0.5, label="16nS")
    sp1.set_ylabel("Vm (mV)")
    sp1.set_xlabel("Time (ms)")
    sp1.plot(T, (np.concatenate((SO, Silence)) * 75) - 75.0, color="g", linewidth=0.5)
    sp1.legend()

    sp2 = fh.add_subplot(4, 1, 2)
    sp2.plot(T, R32["Vm"], "g-", linewidth=0.5, label="Dendritic")
    sp2.plot(T, R32Som["Vm"], "r-", linewidth=0.5, label="Somatic")
    sp2.plot(T, (np.concatenate((SO, Silence)) * 75) - 75.0, color="g", linewidth=0.5)
    sp2.set_ylabel("Vm (mV)")
    sp2.set_xlabel("Time (ms)")
    sp2.spines["top"].set_visible(False)
    sp2.spines["right"].set_visible(False)
    sp2.legend()

    sp3 = fh.add_subplot(4, 1, 3)
    sp3.set_title("3D Dendrite")
    sp3.plot(T, R3d_0["Vm"], "b-", linewidth=0.5, label="0nS")
    sp3.plot(T, R3d_32["Vm"], "g-", linewidth=0.5, label="32nS")
    sp3.plot(T, (np.concatenate((SO, Silence)) * 75) - 75.0, color="g", linewidth=0.5)
    sp3.set_ylabel("Vm (mV)")
    sp3.set_xlabel("Time (ms)")
    sp3.spines["top"].set_visible(False)
    sp3.spines["right"].set_visible(False)
    sp3.legend()

    sp4 = fh.add_subplot(4, 1, 4)
    sp4.plot(T, R3d_32["Vm"], "g-", linewidth=0.5, label="Dendritic")
    sp4.plot(T, R3d_32som["Vm"], "r-", linewidth=0.5, label="Somatic")
    sp4.plot(T, (np.concatenate((SO, Silence)) * 75) - 75.0, color="g", linewidth=0.5)
    sp4.set_ylabel("Vm (mV)")
    sp4.set_xlabel("Time (ms)")
    sp4.spines["top"].set_visible(False)
    sp4.spines["right"].set_visible(False)
    sp4.legend()

    plt.tight_layout()
    plt.show()


def check_req_directories():
    import os.path
    print("results directory exists: " + str(os.path.exists("./results")))
    print("figs directory exists: " + str(os.path.exists("./figs")))


if __name__ == "__main__":
    check_req_directories()
    runmain()
