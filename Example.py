#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 10:33:53 2020
@author: Thomas Kuenzel (kuenzel{at}bio2.rwth-aachen.de)
"""

import matplotlib.pyplot as plt
import numpy as np
import smallexc.semodels as sem
import smallexc.sehelper as seh

def runmain():
    simdur = 400.0#1000
    dt=0.01
    N=33
    freq = 1250.0
    cf = 1501.0
    db = 70
    SO=seh.create_sound(fs=100e3,freq=freq,duration=(simdur*0.6)/1000.0,dbspl=db)
    Silence=seh.create_sound(fs=100e3,freq=freq,duration=(simdur*0.4)/1000.0,dbspl=-10)
    SP=seh.create_spk(np.concatenate((SO,Silence)),fs=100e3,N=N,cf=cf,seed=65562,anf_num=(N,0,0))
    S=np.array(SP['spikes']*1000.0)
    R32=sem.SE_BS(S,Simdur=simdur,dt=dt,G_EoH=0.055,
                  StochasticEoH=True,N=N,G_Dend=0.016,EoHseed=8347236)
    R0=sem.SE_BS(S,Simdur=simdur,dt=dt,G_EoH=0.055,
                 StochasticEoH=True,N=N,G_Dend=0.0,EoHseed=8347236)
    T=np.linspace(0,simdur-dt,int(round(simdur/dt)))
    
    fh=plt.figure()
    sp1=fh.add_subplot(2,1,1)
    sp1.plot(T,R0['Vm'],'k-',linewidth=0.5)
    sp1.set_ylabel('Vm (mV)')
    sp1.set_xlabel('Time (ms)')
    sp1.plot(T,(np.concatenate((SO,Silence))*75)-75.0,color='g',linewidth=0.5)
    
    sp2=fh.add_subplot(2,1,2)
    sp2.plot(T,R32['Vm'],'k-',linewidth=0.5)
    sp2.plot(T,(np.concatenate((SO,Silence))*75)-75.0,color='g',linewidth=0.5)
    sp2.set_ylabel('Vm (mV)')
    sp2.set_xlabel('Time (ms)')
    sp2.spines['top'].set_visible(False)
    sp2.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.show()
    
if __name__ == '__main__':
    runmain()