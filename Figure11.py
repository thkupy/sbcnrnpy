#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This file contains all code necessary to generate Figure 11.
-->Realistic dendrite frequency response areas, VS and preferred phase...
Code modified after original code by Elisabeth Koert.
@author: Thomas Kuenzel (kuenzel{at}bio2.rwth-aachen.de)
"""
##TK add my dirs to path
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import smallexc.semodels as sem
import smallexc.sehelper as seh
import multiprocessing
import os
from functools import partial


Exp4_loadifexists = True#do not run sim if data already exists

def runonce(x,P):
    ####hard coded settings, rarely change these:
    APdetect_thr = -75
    ####GENERATE SPIKE INPUTS (the spiketimes are always the same to improve comparability)
    thisdur=(P['dur'][x]-5.0)/1000.0
    N=P['Nsyn'][x]+1
    anf_num = (N,0,0)
    if P['AdvZseed']:
        thisZseed = P['Zseed']+x
    else:
        thisZseed = P['Zseed']
    if P['AdvEoHseed']:
        thisEoHseed = P['EoHseed']+x
    else:
        thisEoHseed = P['EoHseed']
    SO=seh.create_sound(fs=100e3,freq=P['freq'][x],duration=thisdur,dbspl=P['dB'][x])
    SP=seh.create_spk(SO,fs=100e3,N=N,cf=P['cf'][x],seed=thisZseed,anf_num=anf_num)#generate enough spiketrains
    S=np.array(SP['spikes']*1000.0)
    R0=sem.SE_3D(S,Simdur=P['dur'][x],dt=P['dt'][x],
                 G_EoH=P['G_EoH'][x],StochasticEoH=P['Stochastic_EoH'][x],
                 N=P['Nsyn'][x],G_Dend=0.0,tau_Dend=P['tau_dendsyn'][x],
                 gKLT_d=P['gKLT_d'][x],gIH_d=P['gIH_d'][x],
                 gLEAK_d=P['gLEAK_d'][x],
                 cell=P['cell'][x],EoHseed=thisEoHseed,)

    RD=sem.SE_3D(S,Simdur=P['dur'][x],dt=P['dt'][x],
                 G_EoH=P['G_EoH'][x],StochasticEoH=P['Stochastic_EoH'][x],
                 N=P['Nsyn'][x],G_Dend=P['gsyn'][x],tau_Dend=P['tau_dendsyn'][x],
                 gKLT_d=P['gKLT_d'][x],gIH_d=P['gIH_d'][x],
                 gLEAK_d=P['gLEAK_d'][x],
                 cell=P['cell'][x],EoHseed=thisEoHseed)

    Ev0=seh.SimpleDetectAP(R0['Vm'],thr=APdetect_thr,dt=P['dt'][x],LM=-20,RM=10)
    EvD=seh.SimpleDetectAP(RD['Vm'],thr=APdetect_thr,dt=P['dt'][x],LM=-20,RM=10)

    APr0=len(Ev0['PeakT'])/(P['dur'][x]/1000.0)
    APrD=len(EvD['PeakT'])/(P['dur'][x]/1000.0)

    VS0,phi0,Ray0,phases0=seh.vectorstrength(Ev0['PeakT'],P['freq'][x],[0.0,P['dur'][x]])
    VSD,phiD,RayD,phasesD=seh.vectorstrength(EvD['PeakT'],P['freq'][x],[0.0,P['dur'][x]])
    print(str(x))
    return [VS0,phi0,Ray0,APr0,VSD,phiD,RayD,APrD]

def myMPhandler(P):
    p = multiprocessing.Pool(P['cores'])
    poolmess=partial(runonce, P=P)
    if P['mp']:
        r = p.map(poolmess, P['Number'])
    else:#for debug
        r = list(map(poolmess, P['Number']))#for debug
    return r

def plotres(output,P,x,y,xlabs,ylabs):
    from matplotlib.ticker import ScalarFormatter,AutoMinorLocator,MaxNLocator
    import scipy.ndimage as ndimage

    plt.rc('font', family='serif', serif='Times')
    plt.rc('text', usetex=True)
    plt.rc('xtick', labelsize=7)
    plt.rc('ytick', labelsize=7)
    plt.rc('axes', labelsize=8)
    titsize=9
    mycmap = 'bone'
    
    r_crit = 0.001#0.001
    outVS0=output[0][:,0]
    outPH0=output[0][:,1]
    outRC0=output[0][:,2]
    outAP0=output[0][:,3]
    outVSD=output[0][:,4]
    outPHD=output[0][:,5]
    outRCD=output[0][:,6]
    outAPD=output[0][:,7]
    
    outVS0[np.logical_or(outRC0>r_crit,outRCD>r_crit)]=0.0#set VS of conditions to zero that fail the rayleigh test
    outPH0[np.logical_or(outRC0>r_crit,outRCD>r_crit)]=0.0#set phi of conditions to zero that fail the rayleigh test
    outVSD[np.logical_or(outRC0>r_crit,outRCD>r_crit)]=0.0#set VS of conditions to zero that fail the rayleigh test
    outPHD[np.logical_or(outRC0>r_crit,outRCD>r_crit)]=0.0#set VS of conditions to zero that fail the rayleigh test
    pshift=seh.get_anglediff(outPH0,outPHD)#new test
    pshift=np.reshape(pshift,(P['N'],P['N']))#new test
    APdiff = outAPD-outAP0
    APdiff = np.reshape(APdiff,(P['N'],P['N']))#new test
    VSdiff = outVSD-outVS0
    VSdiff = np.reshape(VSdiff,(P['N'],P['N']))#new test
    VS0=np.reshape(outVS0,(P['N'],P['N']))
    PH0=np.reshape(outPH0,(P['N'],P['N']))
    RC0=np.reshape(outRC0,(P['N'],P['N']))
    AP0=np.reshape(outAP0,(P['N'],P['N']))
    VSD=np.reshape(outVSD,(P['N'],P['N']))
    PHD=np.reshape(outPHD,(P['N'],P['N']))
    RCD=np.reshape(outRCD,(P['N'],P['N']))
    APD=np.reshape(outAPD,(P['N'],P['N']))
    
    filtsig = 0.5
    AP0=ndimage.gaussian_filter(AP0, sigma=filtsig, order=0)
    APD=ndimage.gaussian_filter(APD, sigma=filtsig, order=0)
    APdiff=ndimage.gaussian_filter(APdiff, sigma=filtsig, order=0)
    VS0=ndimage.gaussian_filter(VS0, sigma=filtsig, order=0)
    VSD=ndimage.gaussian_filter(VSD, sigma=filtsig, order=0)
    VSdiff=ndimage.gaussian_filter(VSdiff, sigma=filtsig, order=0)

    fwidth=17.0#cm
    fhandle=plt.figure(figsize=(fwidth/2.54, (fwidth/1.618)/2.54), dpi=600)

    #row 2: VS
    cbkwargs = {'format': '%.1f'}
    absomin = np.min(np.array([np.min(VS0),np.min(VSD)]))
    absomax = np.max(np.array([np.max(VS0),np.max(VSD)]))
    conts = np.linspace(absomin,absomax,21)
    ax1=plt.subplot(334)
    CS1=plt.contourf(x,y,VS0,conts,cmap=mycmap)#repeated = y, tiled = x!!
    CS1b=plt.contour(x,y,VS0,(0.5,),colors='yellow')#repeated = y, tiled = x!!
    #CS1=plt.pcolormesh(x,y,VS0,shading='gouraud',cmap=mycmap)#repeated = y, tiled = x!!
    plt.title(r'VS, $g_{dendsyn} = 0nS$',fontsize=titsize)
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    ax1.set_xscale('log')
    ax1.set_xticks([200,500,1000,2000])
    ax1.set_yticks([0,20,40,60,80])
    ax1.get_xaxis().set_major_formatter(ScalarFormatter())
    ax1.get_yaxis().set_major_formatter(ScalarFormatter())  
    cbar1 = plt.colorbar(CS1,use_gridspec=True,extend='both', **cbkwargs)
    cbar1.ax.set_ylabel('Vectorstrength')
    tl=MaxNLocator(nbins=5)
    cbar1.locator = tl
    cbar1.update_ticks()
    ax1.annotate('', xy=(1500.0,10.0),xytext=(1500.0,-2.0),
            arrowprops=dict(color='white',arrowstyle="wedge"))
    #  
    ax2=plt.subplot(335)
    CS2=plt.contourf(x,y,VSD,conts,cmap=mycmap)#repeated = y, tiled = x!!
    CS2b=plt.contour(x,y,VSD,(0.5,),colors='yellow')#repeated = y, tiled = x!!
    #CS2=plt.pcolormesh(x,y,VSD,shading='gouraud',cmap=mycmap)#repeated = y, tiled = x!!
    plt.title(r'VS, $g_{dendsyn} = 16nS$',fontsize=titsize)
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    ax2.set_xscale('log')
    ax2.set_xticks([200,500,1000,2000])
    ax2.set_yticks([0,20,40,60,80])
    ax2.get_xaxis().set_major_formatter(ScalarFormatter())
    ax2.get_yaxis().set_major_formatter(ScalarFormatter())  
    cbar2 = plt.colorbar(CS2,use_gridspec=True,extend='both', **cbkwargs)
    cbar2.ax.set_ylabel('Vectorstrength')
    tl=MaxNLocator(nbins=5)
    cbar2.locator = tl
    cbar2.update_ticks()
    ax2.annotate('', xy=(1500.0,10.0),xytext=(1500.0,-2.0),
            arrowprops=dict(color='white',arrowstyle="wedge"))
    #
    ax3=plt.subplot(336)
    CS3=plt.contourf(x,y,VSdiff,np.linspace(-0.05,0.05,21),extend='both',cmap='seismic')#repeated = y, tiled = x!!
    plt.title('VS change',fontsize=titsize)
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    plt.xscale('log')
    ax3.set_xscale('log')
    ax3.set_xticks([200,500,1000,2000])
    ax3.set_yticks([0,20,40,60,80])
    ax3.get_xaxis().set_major_formatter(ScalarFormatter())
    ax3.get_yaxis().set_major_formatter(ScalarFormatter())   
    cbar3 = plt.colorbar(CS3,use_gridspec=True,extend='both',cmap='seismic')
    cbar3.ax.set_ylabel(r'$\Delta$ Vectorstrength')
    tl=MaxNLocator(nbins=5)
    cbar3.locator = tl
    cbar3.update_ticks()
    ax3.annotate('', xy=(1500.0,10.0),xytext=(1500.0,-2.0),
            arrowprops=dict(color='black',arrowstyle="wedge"))

    #row 3: preferred phase
    absomin = np.min(np.array([np.min(PH0),np.min(PHD)]))
    absomax = np.max(np.array([np.max(PH0),np.max(PHD)]))
    conts = np.linspace(absomin,absomax,21)
    ax4=plt.subplot(337)
    CS4=plt.contourf(x,y,PH0,conts,cmap=mycmap)#repeated = y, tiled = x!!
    #CS4=plt.pcolormesh(x,y,PH0,shading='gouraud',cmap=mycmap)#repeated = y, tiled = x!!
    plt.title(r'$\varphi$ , $g_{dendsyn} = 0nS$',fontsize=titsize)
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    plt.xscale('log')
    ax4.set_xscale('log')
    ax4.set_xticks([200,500,1000,2000])
    ax4.set_yticks([0,20,40,60,80])
    ax4.get_xaxis().set_major_formatter(ScalarFormatter())
    ax4.get_yaxis().set_major_formatter(ScalarFormatter())  
    cbar4 = plt.colorbar(CS4,use_gridspec=True,extend='both', **cbkwargs)
    cbar4.ax.set_ylabel(r'$\varphi$ (cycles)')
    tl=MaxNLocator(nbins=5)
    cbar4.locator = tl
    cbar4.update_ticks()
    ax4.annotate('', xy=(1500.0,10.0),xytext=(1500.0,-2.0),
            arrowprops=dict(color='white',arrowstyle="wedge"))
    #
    ax5=plt.subplot(338)
    CS5=plt.contourf(x,y,PHD,conts,cmap=mycmap)#repeated = y, tiled = x!!
    #CS5=plt.pcolormesh(x,y,PHD,shading='gouraud',cmap=mycmap)#repeated = y, tiled = x!!
    plt.title(r'$\varphi$ , $g_{dendsyn} = 16nS$',fontsize=titsize)
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    plt.xscale('log')
    ax5.set_xscale('log')
    ax5.set_xticks([200,500,1000,2000])
    ax5.set_yticks([0,20,40,60,80])
    ax5.get_xaxis().set_major_formatter(ScalarFormatter())
    ax5.get_yaxis().set_major_formatter(ScalarFormatter())  
    cbar5 = plt.colorbar(CS5,use_gridspec=True,extend='both', **cbkwargs)
    cbar5.ax.set_ylabel(r'$\phi$ (cycles)')
    tl=MaxNLocator(nbins=5)
    cbar5.locator = tl
    cbar5.update_ticks()
    ax5.annotate('', xy=(1500.0,10.0),xytext=(1500.0,-2.0),
            arrowprops=dict(color='white',arrowstyle="wedge"))
    #
    ax6=plt.subplot(339)
    #pshift = seh.get_anglediff(PH0,PHD)
    #print pshift
    CS6=plt.contourf(x,y,pshift,np.linspace(-0.08,0.08,21),extend='both',cmap='seismic')#repeated = y, tiled = x!!
    plt.title('Phase shift',fontsize=titsize)
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    plt.xscale('log')
    ax6.set_xscale('log')
    ax6.set_xticks([200,500,1000,2000])
    ax6.set_yticks([0,20,40,60,80])
    ax6.get_xaxis().set_major_formatter(ScalarFormatter())
    ax6.get_yaxis().set_major_formatter(ScalarFormatter())  
    cbar6 = plt.colorbar(CS6,use_gridspec=True,extend=mycmap)
    cbar6.ax.set_ylabel(r'$\Delta \varphi$ (cycles)')
    tl=MaxNLocator(nbins=5)
    cbar6.locator = tl
    cbar6.update_ticks()
    ax6.annotate('', xy=(1500.0,10.0),xytext=(1500.0,-2.0),
            arrowprops=dict(color='black',arrowstyle="wedge"))
    
    #row 1: AP rate
    cbkwargs = {'format': '%.0f'}
    absomin = np.min(np.array([np.min(AP0),np.min(APD)]))
    absomax = np.max(np.array([np.max(AP0),np.max(APD)]))
    conts = np.linspace(absomin,absomax,21)
    ax7=plt.subplot(331)
    CS7=plt.contourf(x,y,AP0,conts,cmap=mycmap)#repeated = y, tiled = x!!
    CS7b=plt.contour(x,y,AP0,(90.0,200.0,),colors='yellow')#repeated = y, tiled = x!!
    #CS7=plt.pcolormesh(x,y,AP0,shading='gouraud',cmap=mycmap)#repeated = y, tiled = x!!
    plt.title('AP rate, $g_{dendsyn} = 0nS$',fontsize=titsize)
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    plt.xscale('log')
    ax7.set_xscale('log')
    ax7.set_xticks([200,500,1000,2000])
    ax7.set_yticks([0,20,40,60,80])
    ax7.get_xaxis().set_major_formatter(ScalarFormatter())
    ax7.get_yaxis().set_major_formatter(ScalarFormatter())  
    cbar7 = plt.colorbar(CS7,use_gridspec=True,extend='both', **cbkwargs)
    cbar7.ax.set_ylabel('Rate (Hz)')
    tl=MaxNLocator(nbins=5)
    cbar7.locator = tl
    cbar7.update_ticks()
    ax7.annotate('', xy=(1500.0,10.0),xytext=(1500.0,-2.0),
            arrowprops=dict(color='white',arrowstyle="wedge"))
    #
    ax8=plt.subplot(332)
    CS8=plt.contourf(x,y,APD,conts,cmap=mycmap)#repeated = y, tiled = x!!
    CS8b=plt.contour(x,y,APD,(90.0,200.0),colors='yellow')#repeated = y, tiled = x!!
    #CS8=plt.pcolormesh(x,y,APD,shading='gouraud',cmap=mycmap)#repeated = y, tiled = x!!
    plt.title('AP rate, $g_{dendsyn} = 16nS$',fontsize=titsize)
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    plt.xscale('log')
    ax8.set_xscale('log')
    ax8.set_xticks([200,500,1000,2000])
    ax8.set_yticks([0,20,40,60,80])
    ax8.get_xaxis().set_major_formatter(ScalarFormatter())
    ax8.get_yaxis().set_major_formatter(ScalarFormatter())  
    cbar8 = plt.colorbar(CS8,use_gridspec=True,extend='both', **cbkwargs)
    cbar8.ax.set_ylabel('Rate (Hz)')
    tl=MaxNLocator(nbins=5)
    cbar8.locator = tl
    cbar8.update_ticks()
    ax8.annotate('', xy=(1500.0,10.0),xytext=(1500.0,-2.0),
            arrowprops=dict(color='white',arrowstyle="wedge"))
    #
    ax9=plt.subplot(333)
    CS9=plt.contourf(x,y,APdiff,np.linspace(-10,10,21),extend='both',cmap='seismic')#repeated = y, tiled = x!!
    plt.title('Rate change',fontsize=titsize)
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    plt.xscale('log')
    ax9.set_xscale('log')
    ax9.set_xticks([200,500,1000,2000])
    ax9.set_yticks([0,20,40,60,80])
    ax9.get_xaxis().set_major_formatter(ScalarFormatter())
    ax9.get_yaxis().set_major_formatter(ScalarFormatter())  
    cbar9 = plt.colorbar(CS9,use_gridspec=True,extend='both')
    cbar9.ax.set_ylabel(r'$\Delta$ Rate (Hz)')
    tl=MaxNLocator(nbins=5)
    cbar9.locator = tl
    cbar9.update_ticks()
    ax9.annotate('', xy=(1500.0,10.0),xytext=(1500.0,-2.0),
            arrowprops=dict(color='black',arrowstyle="wedge"))
    plt.tight_layout()
    return fhandle

if __name__ == '__main__':
    if os.path.isfile('./results/Figure11.npy') and Exp4_loadifexists:
        print("Data for Figure11 found... loading!")
        output = np.load('./results/Figure11.npy', allow_pickle=True)
        P = np.load('./results/Figure11_P.npy', allow_pickle=True)
        P = P.tolist()
    else:
        conditions=32 #number of conditions per varied parameter (here only 2d experiments)
        cores=20
        dur=5000.0
        dt=0.01
        output=[]
        #Model Parameters (all in a linearly aranged fashion, so that a minimal
        #amount of programming is required to change the experiment).
        P={}
        P['N']=conditions
        P['cores']=cores
        P['TotalN'] = int(P['N']**2)
        P['Number'] = range(P['TotalN'])
        P['mp']=True
        P['Zseed']=45453
        P['AdvZseed']=True
        P['EoHseed']=34823
        P['AdvEoHseed']=True
        #########################################
        P['cell'] = np.repeat(1,P['TotalN'])#which model neuron is used
        #########################################
        P['dur'] = np.repeat(dur,P['TotalN'])
        P['dt'] = np.repeat(dt,P['TotalN'])
        P['G_EoH'] = np.repeat(0.055,P['TotalN'])
        P['Stochastic_EoH'] = np.repeat(True,P['TotalN'])
        P['Nsyn'] = np.repeat(32,P['TotalN'])#32
        P['gsyn'] = np.repeat(0.064,P['TotalN'])#0.048
        P['tau_dendsyn'] = np.repeat(2.0,P['TotalN'])
        P['IPI'] = np.repeat(0.0,P['TotalN'])#irrelevant
        P['gLEAK_d'] = np.repeat(0.001,P['TotalN'])
        P['gKLT_d'] = np.repeat(0.0085,P['TotalN'])
        P['gIH_d'] = np.repeat(0.001,P['TotalN'])
        P['gKLT_s'] = np.repeat(0.017,P['TotalN'])
        P['gIH_s'] = np.repeat(0.002,P['TotalN'])
        #Now define the two variable parameters. The repeated = y, the tiled = x!!
        thiscf = 1500.0
        allfreq= np.round(np.geomspace(125.0,thiscf*2.5,P['N']))
        alllev = np.round(np.linspace(0.0,80.0,P['N']),1)
        P['freq'] = np.tile(allfreq,P['N'])
        P['cf'] = np.repeat(thiscf,P['TotalN'])
        P['dB'] = np.repeat(alllev,P['N'])

        #make go!
        output.append(myMPhandler(P))
        output = np.array(output)
        np.save('./results/Figure11.npy',output, allow_pickle=True)
        np.save('./results/Figure11_P.npy',P, allow_pickle=True)

    ###Plot results
    fhandle=plotres(output=output,P=P,x=np.unique(P['freq']),y=np.unique(P['dB']),
            xlabs=u'Frequency (Hz)',ylabs=u'Level (dB SPL)')
    pp = PdfPages('./figs/Figure11.pdf')
    pp.savefig()
    pp.close()
