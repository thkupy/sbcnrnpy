#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 21 Aug 2019
This file contains all code necessary to generate Figure 6.
@author: Thomas Kuenzel (kuenzel{at}bio2.rwth-aachen.de)
"""
import os
import matplotlib
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import smallexc.semodels as sem
import smallexc.sehelper as seh
import multiprocessing
from functools import partial

plt.rc('font', family='serif',serif='times')
plt.rcParams['pdf.fonttype'] = 42
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')
plt.rc('axes', labelsize='small')

def runmain():
    simdur = 1000.0#1000
    dt=0.01
    N=33
    freq = 1251.0#1389.0#501.0
    cf = 1501.0#1501.0#751.0
    db = 70#65.0#60
    SO=seh.create_sound(fs=100e3,freq=freq,duration=(simdur*0.6)/1000.0,dbspl=db)
    Silence=seh.create_sound(fs=100e3,freq=freq,duration=(simdur*0.4)/1000.0,dbspl=-10)
    SP=seh.create_spk(np.concatenate((SO,Silence)),fs=100e3,N=N,cf=cf,seed=65562,anf_num=(N,0,0))
    S=np.array(SP['spikes']*1000.0)
    R32=sem.SE_BS(S,Simdur=simdur,dt=dt,G_EoH=0.055,
                  StochasticEoH=True,N=N,G_Dend=0.016,EoHseed=8347236)
    R0=sem.SE_BS(S,Simdur=simdur,dt=dt,G_EoH=0.055,
                 StochasticEoH=True,N=N,G_Dend=0.0,EoHseed=8347236)
    T=np.linspace(0,simdur-dt,int(round(simdur/dt)))
    sT=np.linspace(0,(simdur*0.8)-dt,int(round((simdur*0.8)/dt)))
    Ev0=seh.SimpleDetectAP(R0['Vm'],thr=-100,dt=dt,LM=-20,RM=10)
    Fail0=seh.SimpleGetFails(S[-1],Ev0,R0,simdur,dt)
    Ev32=seh.SimpleDetectAP(R32['Vm'],thr=-100,dt=dt,LM=-20,RM=10)
    Fail32=seh.SimpleGetFails(S[-1],Ev32,R32,simdur,dt)
    APr0=len(Ev0['PeakT'])/(simdur/1000.0)
    APr32=len(Ev32['PeakT'])/(simdur/1000.0)
    print str(APr0) + 'vs' + str(APr32)
    VS0,phi0,Ray0,phases0=seh.vectorstrength(Ev0['PeakT'],freq,[0.0,simdur*0.6])
    VS32,phi32,Ray32,phases32=seh.vectorstrength(Ev32['PeakT'],freq,[0.0,simdur*0.6])
    
    fwidth=8.5
    fh=plt.figure(figsize=(fwidth/2.54, fwidth/2.54))
    fh.subplots_adjust(left=0.15, bottom=0.125, right=0.95, top=0.95, wspace=0.55, hspace=0.55)
    sp1=fh.add_subplot(2,3,(1,2))
    sp1.plot(T,R0['Vm'],'k-',linewidth=0.5)
    sp1.plot(Ev0['PeakT'],Ev0['PeakV'],'g.',Markersize=2)
    sp1.plot(Fail0['PeakT'],Fail0['PeakV'],'r.',Markersize=2)
    sp1.plot(S[-1],np.zeros((len(S[-1]),))+20.0,'b.',Markersize=2)
    #sp1.plot(sT,(SO*6)-65.0,color='r',linewidth=0.5)
#    sp1.set_ylim((-66.0,-58.0))
    sp1.set_ylabel('Vm (mV)')
    sp1.set_xlabel('Time (ms)')
#    sp1.set_xticks((0,50,100,150))
#    sp1.set_yticks((-65,-63,-61,-59))
    #sp1.set_xlim((50.0,450.0))
    sp1.spines['top'].set_visible(False)
    sp1.spines['right'].set_visible(False)
    sp1.plot(T,(np.concatenate((SO,Silence))*75)-75.0,color='g',linewidth=0.5)
    
    sp1b=fh.add_subplot(2,3,3)
    H0=np.histogram(phases0,np.linspace(0,1,33))
    sp1b.plot(H0[1][1:],H0[0],'k')
    sp1b.set_ylabel('AP/bin')
    sp1b.set_xlabel('Phase (cycles)')
    sp1b.set_xticks((0,0.5,1))
    sp1b.set_yticks((0,5,10))
    sp1b.plot((phi0,phi0),(0,10),'r-')
    annotext=r'VS = '+str(round(VS0,2))
    sp1b.annotate(annotext,xy=(phi0,10),fontsize=5)
    sp1b.spines['top'].set_visible(False)
    sp1b.spines['right'].set_visible(False)
    
    sp2=fh.add_subplot(2,3,(4,5))
    sp2.plot(T,R32['Vm'],'k-',linewidth=0.5)
    sp2.plot(Ev32['PeakT'],Ev32['PeakV'],'g.',Markersize=2)
    sp2.plot(Fail32['PeakT'],Fail32['PeakV'],'r.',Markersize=2)
    sp2.plot(S[-1],np.zeros((len(S[-1]),))+20.0,'b.',Markersize=2)
    #sp2.set_xlim((0.0,50.0))
    sp2.plot(T,(np.concatenate((SO,Silence))*75)-75.0,color='g',linewidth=0.5)
    #sp2.set_ylim((-66.0,-58.0))
    sp2.set_ylabel('Vm (mV)')
    sp2.set_xlabel('Time (ms)')
    #sp2.set_xticks((0,50,100,150))
    #sp2.set_yticks((-65,-63,-61,-59))
    sp2.spines['top'].set_visible(False)
    sp2.spines['right'].set_visible(False)
#
    sp2b=fh.add_subplot(2,3,6)
    H32=np.histogram(phases32,np.linspace(0,1,33))
    sp2b.plot(H32[1][1:],H32[0],'k')
    #sp2b.plot(phases32,np.ones((len(phases32))),'b.',markersize=1)
    #for iph,phi in enumerate(phases32):
    #    sp2b.plot((phi,phi),(0,1),'-',color=(.8,.8,.8),linewidth=0.25)
    #sp2b.plot((phi32,phi32),(0,VS32),'b-',linewidth=0.5)
    sp2b.set_ylabel('AP/bin')
    sp2b.set_xlabel('Phase (cycles)')
    sp2b.set_xticks((0,0.5,1))
    sp2b.set_yticks((0,5,10))
    sp2b.plot((phi32,phi32),(0,10),'r-')
    annotext=r'VS = '+str(round(VS32,2))
    sp2b.annotate(annotext,xy=(phi32,10),fontsize=5)
    sp2b.spines['top'].set_visible(False)
    sp2b.spines['right'].set_visible(False)
    
#    sp3=fh.add_subplot(233)
#    sp3.plot(S[0],np.ones(len(S[0]))*0,'k.',markersize=1)
#    sp3.plot(sT,((SO*3.0)-0.5),color='r',linewidth=0.5)
#    sp3.set_xlim((0,simdur))
#    sp3.set_ylim((-1.1,N/10.0))
#    sp3.axes.get_yaxis().set_visible(False)
#    sp3.set_xticks((0,50,100,150))
#    sp3.spines['top'].set_visible(False)
#    sp3.spines['right'].set_visible(False)
#              
#    sp4=fh.add_subplot(236)
#    for s in range(N):
#        sp4.plot(S[s],np.ones(len(S[s]))*(s/10.0),'k.',markersize=1)
#    sp4.set_xlim((0,simdur))
#    sp4.plot(sT,((SO*3.0)-0.5),color='r',linewidth=0.5)
#    sp4.set_ylim((-1.1,N/10.0))
#    sp4.set_xlabel('Time (ms)')
#    sp4.axes.get_yaxis().set_visible(False) 
#    sp4.set_xticks((0,50,100,150))
#    sp4.spines['top'].set_visible(False)
#    sp4.spines['right'].set_visible(False)
              
    #plt.tight_layout()
    #plt.show()
    pp = PdfPages('./Figure6/Figure6_raw.pdf')
    pp.savefig(dpi=600)
    pp.close()
    
if __name__ == '__main__':
    runmain()
