#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Jul 11 2019
This file contains all code necessary to generate Figure 2.
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

def average_over_cyle(D,freq,dt,dur):
    IPI = 1000.0/freq
    allc=[]
    cycsamps=int(np.round(IPI/dt))
    c=range(0,len(D)-cycsamps,cycsamps)
    for ii,cc in enumerate(c):
        allc.append(D[cc:cc+cycsamps])
    av = np.mean(np.array(allc),axis=0)
    st = np.std(np.array(allc),axis=0)
    meanD = np.mean(av)
    stdD = np.std(av)
    modamp = np.max(av)-np.min(av)
    return(av,st,meanD,stdD,modamp)

def runonce(x,P):
    print(x)
    thisdur=(P['dur'][x]-5.0)/1000.0
    SO=seh.create_sound(fs=100e3,freq=P['freq'][x],duration=thisdur,dbspl=P['dB'][x])
    SP=seh.create_spk(SO,fs=100e3,N=P['Nsyn'][x],cf=P['cf'][x],seed=2121977)#generate enough spiketrains
    S=np.array(SP['spikes']*1000.0)
    R=sem.SE_BS(S,Simdur=P['dur'][x],dt=P['dt'][x],G_EoH=P['G_EoH'][x],
              StochasticEoH=P['Stochastic_EoH'][x],
              N=P['Nsyn'][x],G_Dend=P['gsyn'][x],tau_Dend=P['tau_dendsyn'][x],
              L=P['L'][x],D=P['D'][x],constR=P['constR'][x],
              gKLT_d=P['gKLT_d'][x],gIH_d=P['gIH_d'][x],gLEAK_d=P['gLEAK_d'][x])
    Ev=seh.SimpleDetectAP(R['Vm'],thr=-100,dt=P['dt'][x],LM=-20,RM=10)
    SAPr=len(Ev['PeakT'])/(P['dur'][x]/1000.0)  
    sav,sst,Smeanv,Stdv,Smodv=average_over_cyle(D=R['Vm'],freq=P['freq'][x],dt=P['dt'][x],dur=P['dur'][x])
    #save example data (if required)
    if P['save_examples'] and np.any(P['which_example_to_save']==x):
        XT={}
        XT['Sound']=R['Vm']
        XT['Av_Sound']=sav
        XT['St_Sound']=sst
        thisfn = P['Example_PN'] + P['Example_FN'] + str(x) + '.npy'
        np.save(thisfn,XT,allow_pickle=True)       
    #
    return [SAPr,Smeanv,Smodv]

def myMPhandler(P):
    p = multiprocessing.Pool(P['cores'])
    poolmess=partial(runonce, P=P)
    if P['mp']:
        r = p.map(poolmess, P['Number'])
    else:#for debug
        r = list(map(poolmess, P['Number']))#for debug
    return r
 
def basic_parameters_1d(N=21,cores=4,dur=500.0,dt=0.025):
    P={}
    P['mp'] = True
    P['N']=int(conditions)#just to be sure
    P['cores']=int(cores)#just to be sure
    P['Number'] = range(P['N'])
    P['distritype']=['normal','alpha']
    P['DistParam']=[0.505,1.72]
    P['save_examples']=True
    P['N-examples']=5
    P['Example_PN']='./results/tmp/'
    P['Example_FN']='ex_'
    P['which_example_to_save']=np.linspace(0,P['N']-1,P['N-examples'],dtype=int)
    P['var_key']=' '
    P['titletext']=' '
    P['labeltext']=' '
    #########################################
    P['dur'] = np.repeat(dur,P['N'])
    P['dt'] = np.repeat(dt,P['N'])
    P['G_EoH'] = np.repeat(0.0,P['N'])
    P['Stochastic_EoH'] = np.repeat(False,P['N'])
    P['Nsyn'] = np.repeat(32,P['N'])
    P['gsyn'] = np.repeat(0.016,P['N'])
    P['tau_dendsyn'] = np.repeat(2.0,P['N'])
    P['L'] = np.repeat(50.0,P['N'])
    P['D'] = np.repeat(3.0,P['N'])      
    P['constR'] = np.repeat(False,P['N'])
    P['freq'] = np.repeat(200.0,P['N'])
    P['cf'] = np.repeat(200.0,P['N'])
    P['dB'] = np.repeat(60.0,P['N'])
    P['IPI'] = np.repeat(5.0,P['N'])        
    P['gLEAK_d'] = np.repeat(0.001,P['N'])
    P['gKLT_s'] = np.repeat(0.017,P['N'])
    P['gIH_s'] = np.repeat(0.002,P['N'])
    P['gKLT_d'] = np.repeat(0.0085,P['N'])
    P['gIH_d'] = np.repeat(0.001,P['N'])
    return(P)

def process_examplefiles(P):
    exfilelist = []
    NT=[]
    AvN=[]
    StN=[]
    AT=[]
    AvA=[]
    StA=[]
    ST=[]
    AvS=[]
    StS=[]
    for file in os.listdir(P['Example_PN']):
        if file.endswith(".npy"):
            exfilelist.append(os.path.join(P['Example_PN'], file))
    for thisfile in sorted(exfilelist):
        tx=np.load(thisfile,allow_pickle=True)
        tx=tx.tolist()
        ST.append(tx['Sound'])
        AvS.append(tx['Av_Sound'])
        StS.append(tx['St_Sound'])
    XT={}
    XT['Sound']=ST
    XT['Av_Sound']=AvS
    XT['St_Sound']=StS
    #remove the files now
    for thisfile in exfilelist:
        os.remove(thisfile)
    return(XT)

def plotres_2d(output,P,x,y,xlabs,ylabs):
    import matplotlib.font_manager
    from matplotlib.ticker import ScalarFormatter,AutoMinorLocator,MaxNLocator
    import scipy.ndimage as ndimage

    plt.rc('font', family='serif', serif='Times')
    plt.rc('text', usetex=True)
    plt.rc('xtick', labelsize=7)
    plt.rc('ytick', labelsize=7)
    plt.rc('axes', labelsize=8)
    nticks = 5
    titsize=9
    mycmap = 'bone'
    filtsig = 0.25
    
    outV=output[0][:,1]#return [SAPr,Smeanv,Smodv]
    outM=output[0][:,2]
    outA=output[0][:,0]
    V=ndimage.gaussian_filter(np.reshape(outV,(P['N'],P['N'])), sigma=filtsig, order=0)
    M=ndimage.gaussian_filter(np.reshape(outM,(P['N'],P['N'])), sigma=filtsig, order=0)
    A=ndimage.gaussian_filter(np.reshape(outA,(P['N'],P['N'])), sigma=filtsig, order=0)
    
    fwidth=7.0#cm
    fhandle=plt.figure(figsize=(fwidth/2.54, (fwidth*2)/2.54), dpi=300)
    
    ax1=plt.subplot(311)
    CS1=plt.contourf(x,y,V,21,cmap=mycmap)#repeated = y, tiled = x!!
    #CS1=plt.pcolormesh(x,y,V,shading='gouraud',cmap=mycmap)#repeated = y, tiled = x!!
    plt.title('Mean Vm')
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    ax1.yaxis.set_major_locator(MaxNLocator(nticks))
    ax1.xaxis.set_major_locator(MaxNLocator(nticks))
    cbar1 = plt.colorbar(CS1,use_gridspec=True)
    cbar1.ax.set_ylabel('Vm (mV)')
    tl=MaxNLocator(nbins=5)
    cbar1.locator = tl
    cbar1.update_ticks()

    ax2=plt.subplot(312)
    CS2=plt.contourf(x,y,M,21,cmap=mycmap)#repeated = y, tiled = x!!
    #CS2=plt.pcolormesh(x,y,M,shading='gouraud',cmap=mycmap)#repeated = y, tiled = x!!
    plt.title('Vm Modulation')
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    ax2.yaxis.set_major_locator(MaxNLocator(nticks))
    ax2.xaxis.set_major_locator(MaxNLocator(nticks))
    cbar2 = plt.colorbar(CS2,use_gridspec=True)
    cbar2.ax.set_ylabel(r'Vm mod. ($\pm$ mV)')
    tl=MaxNLocator(nbins=5)
    cbar2.locator = tl
    cbar2.update_ticks()
    tl=MaxNLocator(nbins=5)
    cbar2.locator = tl
    cbar2.update_ticks()

    ax3=plt.subplot(313)
    CS3=plt.contourf(x,y,A,21,cmap=mycmap)#repeated = y, tiled = x!!
    #CS3=plt.pcolormesh(x,y,A,shading='gouraud',cmap=mycmap)#repeated = y, tiled = x!!
    plt.title('APs')
    plt.xlabel(xlabs)
    plt.ylabel(ylabs)
    ax3.yaxis.set_major_locator(MaxNLocator(nticks))
    ax3.xaxis.set_major_locator(MaxNLocator(nticks))
    cbar3 = plt.colorbar(CS3,use_gridspec=True)
    cbar3.ax.set_ylabel('AP rate (Hz)')

    plt.tight_layout()
    return fhandle
   
def plotresult_1d(output,example,P):
    #Imports and Settings
    from matplotlib.ticker import ScalarFormatter,AutoMinorLocator,MaxNLocator
    import scipy.ndimage as ndimage
    plt.rc('font', family='serif', serif='Times')
    plt.rc('text', usetex=True)
    plt.rc('xtick', labelsize=7)
    plt.rc('ytick', labelsize=7)
    plt.rc('axes', labelsize=8)
    titsize=9
    mycmap = 'bone'
    nticks = 5
    filtsig = 0.25
    r_crit = 0.001
    #
    Ncol = 'b'
    Ncmap = 'Blues'
    Acol = 'g'
    Acmap = 'Greens'
    Scol = 'k'
    Scmap = 'bone'
    #
    NUM_COLORS = P['N-examples']+3
    Scm = plt.get_cmap(Scmap)
    unique_Scolors=[ Scm(1.*(i+1)/NUM_COLORS) for i in range(NUM_COLORS)]
    #
    #Prep data [-->SAPr,Smeanv,Smodv]
    SAPr=output[0][:,0]
    Smeanv=output[0][:,1]
    Smodv=output[0][:,2]
    XD=P[P['var_key']]
    #
    #Create Figures
    fwidth = 8.5#cm
    fwidth_full = 17.0#cm
    
    fhandle1=plt.figure(figsize=(5/2.54, 11/2.54), dpi=300) 
    ax1=plt.subplot(311)
    # To specify the number of ticks on both or any single axes
    plt.plot(XD,Smeanv,color=Scol)
    plt.ylabel('Mean Vm (mV)')
    plt.xlabel(P['labeltext'])
    plt.title(P['titletext'],fontsize=titsize)
    ax1.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax1.yaxis.set_major_locator(plt.MaxNLocator(5))
    #
    ax2=plt.subplot(312)
    plt.plot(XD,Smodv,color=Scol)
    plt.ylabel('Vm Modulation (mV)')
    plt.xlabel(P['labeltext'])
    ax2.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax2.yaxis.set_major_locator(plt.MaxNLocator(5))
    plt.tight_layout()
    #
    if P['save_examples']:
        fhandle2=plt.figure(figsize=(5.5/2.54, (fwidth)/2.54), dpi=300) 
        s=[]
        s.append(plt.subplot(211))
        s.append(plt.subplot(212))
        taxis = np.linspace(0,P['dur'][0]-P['dt'][0],P['dur'][0]/P['dt'][0])
        cycaxis = np.linspace(0,P['IPI'][0]-P['dt'][0],int(np.round(P['IPI'][0]/P['dt'][0])))
        for ix in range(P['N-examples']):
            if P['var_key']=='IPI':
                cycaxis = np.linspace(0.0,1.0,len(example['Av_Sound'][ix]))
                xlabtext2='Cycle (norm.)'
            else:
                xlabtext2='IPI (ms)'
            s[0].plot(taxis,example['Sound'][ix],color=unique_Scolors[ix])
            s[1].plot(cycaxis,example['Av_Sound'][ix],color=unique_Scolors[ix])
        s[0].set_xlim(100.0, 150.0)
        s[0].set_title('Sound / Example trace',fontsize=titsize)
        s[0].set_xlabel('Time (ms)')
        legstring = ["%.3f" % number for number in P[P['var_key']][P['which_example_to_save']]]
        s[0].legend(legstring,fontsize=5)
        s[1].set_title('Sound / averaged trace',fontsize=titsize)
        s[1].set_xlabel(xlabtext2)
        s[1].legend(legstring,fontsize=6)
        for thish in s:
            thish.xaxis.set_major_locator(plt.MaxNLocator(5))
            thish.yaxis.set_major_locator(plt.MaxNLocator(5))
    else:
        fhandle2=[]
    plt.tight_layout()
    return((fhandle1,fhandle2),(ax1,ax2,s))
    
if __name__ == '__main__':
    Fig2_loadifexists=[True,True,True,True,True]
    conditions=64
    conditions_2d = 25
    cores=20
    dur = 2000.0#production = 5000.0
    dt = 0.01#production = 0.01
    
    myylim1=(-64.5,-60.5)
    myylim2=(0.0,5.5)
    #####################################
    #####Fig2a:         #################
    #####---------------#################
    ##### Variable -> L #################
    ##### R_dend = free #################
    ##### No Endbulb    #################
    #####################################
    ####LOAD DATA FOR Fig2a (IF IT EXISTS)
    if os.path.isfile('./results/Fig2a.npy') and Fig2_loadifexists[0]:#Exp-specific file/variable names here
        print("Data for Fig2a found... loading!")
        output2a = np.load('./results/Fig2a.npy', allow_pickle=True)#Exp-specific file/variable names here
        Ex2a = np.load('./results/Fig2a_Ex.npy', allow_pickle=True)#Exp-specific file/variable names here
        Ex2a = Ex2a.tolist()
        P2a = np.load('./results/Fig2a_P.npy', allow_pickle=True)#Exp-specific file/variable names here
        P2a = P2a.tolist()#Exp-specific file/variable names here
    else:
        output2a=[]#Exp-specific file/variable names here
        P2a = basic_parameters_1d(conditions,cores,dur,dt)#Exp-specific file/variable names here
        ####VARIABLE Parameter for EXP1a####
        P2a['var_key']='L'
        P2a['titletext']=r'Dendrite 1 Length'
        P2a['labeltext']=r'L(Dend1) ($\mu$m)'
        P2a['L']=np.linspace(1.0,501.0,P2a['N'])#Exp-specific code here      
        #make go!
        output2a.append(myMPhandler(P2a))#Exp-specific file/variable names here
        output2a = np.array(output2a)#Exp-specific file/variable names here
        if P2a['save_examples']:
            Ex2a=process_examplefiles(P2a)
        np.save('./results/Fig2a.npy',output2a, allow_pickle=True)#Exp-specific file/variable names here
        np.save('./results/Fig2a_P.npy',P2a, allow_pickle=True)#Exp-specific file/variable names here
        np.save('./results/Fig2a_Ex.npy',Ex2a, allow_pickle=True)#Exp-specific file/variable names here
    ###Plot results 
    with PdfPages('./figs/Fig2a.pdf') as pp:#Exp-specific file/variable names here
        fhandles,axhandles=plotresult_1d(output2a,Ex2a,P2a)
        #panel specific axis limits and other artwork
	axhandles[0].set_ylim(myylim1)
	axhandles[1].set_ylim(myylim2)
	axhandles[2][0].set_ylim((-65.0,-51.0))
        axhandles[2][1].set_ylim((-65.0,-51.0))
        #
        pp.savefig(fhandles[0])
        pp.savefig(fhandles[1])
        plt.close()

    #####################################
    #####Fig2b:         #################
    #####---------------#################
    ##### Variable -> D #################
    ##### R_dend = free #################
    ##### No Endbulb    #################
    #####################################
    ####LOAD DATA FOR Fig2b (IF IT EXISTS)
    if os.path.isfile('./results/Fig2b.npy') and Fig2_loadifexists[1]:#Exp-specific file/variable names here
        print("Data for Fig2b found... loading!")
        output2b = np.load('./results/Fig2b.npy', allow_pickle=True)#Exp-specific file/variable names here
        Ex2b = np.load('./results/Fig2b_Ex.npy', allow_pickle=True)#Exp-specific file/variable names here
        Ex2b = Ex2b.tolist()
        P2b = np.load('./results/Fig2b_P.npy', allow_pickle=True)#Exp-specific file/variable names here
        P2b = P2b.tolist()#Exp-specific file/variable names here
    else:
        output2b=[]#Exp-specific file/variable names here
        P2b = basic_parameters_1d(conditions,cores,dur,dt)#Exp-specific file/variable names here
        ####VARIABLE Parameter for EXP1b###
        P2b['var_key']='D'
        P2b['titletext']=r'Dendrite 1 Diameter'
        P2b['labeltext']=r'D(Dend1) ($\mu$m)'
        P2b['D']=np.linspace(1.0,10.0,P2b['N'])#Exp-specific code here      
        #make go!
        output2b.append(myMPhandler(P2b))#Exp-specific file/variable names here
        output2b = np.array(output2b)#Exp-specific file/variable names here
        if P2b['save_examples']:
            Ex2b=process_examplefiles(P2b)
        np.save('./results/Fig2b.npy',output2b, allow_pickle=True)#Exp-specific file/variable names here
        np.save('./results/Fig2b_P.npy',P2b, allow_pickle=True)#Exp-specific file/variable names here
        np.save('./results/Fig2b_Ex.npy',Ex2b, allow_pickle=True)#Exp-specific file/variable names here
    ###Plot results 
    with PdfPages('./figs/Fig2b.pdf') as pp:#Exp-specific file/variable names here
        fhandle,axhandles=plotresult_1d(output2b,Ex2b,P2b)
        #panel specific axis limits and other artwork
	axhandles[0].set_ylim(myylim1)
        axhandles[1].set_ylim(myylim2)
        axhandles[2][0].set_ylim((-65.0,-51.0))
        axhandles[2][1].set_ylim((-65.0,-51.0))
        pp.savefig(fhandle[0])
        pp.savefig(fhandle[1])
        plt.close()
        
    #####################################
    #####Fig2c:         #################
    #####---------------#################
    ##### Variable -> L #################
    ##### R_dend = fix  #################
    ##### No Endbulb    #################
    #####################################
    ####LOAD DATA FOR Fig2c (IF IT EXISTS)
    if os.path.isfile('./results/Fig2c.npy') and Fig2_loadifexists[2]:#Exp-specific file/variable names here
        print("Data for Fig2c found... loading!")
        output2c = np.load('./results/Fig2c.npy', allow_pickle=True)#Exp-specific file/variable names here
        Ex2c = np.load('./results/Fig2c_Ex.npy', allow_pickle=True)#Exp-specific file/variable names here
        Ex2c = Ex2c.tolist()
        P2c = np.load('./results/Fig2c_P.npy', allow_pickle=True)#Exp-specific file/variable names here
        P2c = P2c.tolist()#Exp-specific file/variable names here
    else:
        output2c=[]#Exp-specific file/variable names here
        P2c = basic_parameters_1d(conditions,cores,dur,dt)#Exp-specific file/variable names here
        ####VARIABLE Parameter for EXP1c####
        P2c['var_key']='L'
        P2c['titletext']=r'Dendrite 1 Length ($R_{m}$=fixed)'
        P2c['labeltext']=r'L(Dend1) ($\mu$m)'
        P2c['L']=np.linspace(1.0,501.0,P2c['N'])#Exp-specific code here     
        P2c['constR'] = np.repeat(True,P2c['N'])
        #make go!
        output2c.append(myMPhandler(P2c))#Exp-specific file/variable names here
        output2c = np.array(output2c)#Exp-specific file/variable names here
        if P2c['save_examples']:
            Ex2c=process_examplefiles(P2c)
        np.save('./results/Fig2c.npy',output2c, allow_pickle=True)#Exp-specific file/variable names here
        np.save('./results/Fig2c_P.npy',P2c, allow_pickle=True)#Exp-specific file/variable names here
        np.save('./results/Fig2c_Ex.npy',Ex2c, allow_pickle=True)#Exp-specific file/variable names here
    ###Plot results 
    with PdfPages('./figs/Fig2c.pdf') as pp:#Exp-specific file/variable names here
        fhandle,axhandles=plotresult_1d(output2c,Ex2c,P2c)
        axhandles[0].set_ylim(myylim1)
        axhandles[1].set_ylim(myylim2)
	axhandles[2][0].set_ylim((-65.0,-51.0))
        axhandles[2][1].set_ylim((-65.0,-51.0))
        pp.savefig(fhandle[0])
        pp.savefig(fhandle[1])
        plt.close()

    #####################################
    #####Fig2d:         #################
    #####---------------#################
    ##### Variable -> D #################
    ##### R_dend = fix  #################
    ##### No Endbulb    #################
    #####################################
    ####LOAD DATA FOR Fig2d (IF IT EXISTS)
    if os.path.isfile('./results/Fig2d.npy') and Fig2_loadifexists[3]:#Exp-specific file/variable names here
        print("Data for Fig2d found... loading!")
        output2d = np.load('./results/Fig2d.npy', allow_pickle=True)#Exp-specific file/variable names here
        Ex2d = np.load('./results/Fig2d_Ex.npy', allow_pickle=True)#Exp-specific file/variable names here
        Ex2d = Ex2d.tolist()
        P2d = np.load('./results/Fig2d_P.npy', allow_pickle=True)#Exp-specific file/variable names here
        P2d = P2d.tolist()#Exp-specific file/variable names here
    else:
        output2d=[]#Exp-specific file/variable names here
        P2d = basic_parameters_1d(conditions,cores,dur,dt)#Exp-specific file/variable names here
        ####VARIABLE Parameter for EXP1d####
        P2d['var_key']='D'
        P2d['titletext']=r'Dendrite 1 Diameter ($R_{m}$=fixed)'
        P2d['labeltext']=r'D(Dend1) ($\mu$m)'
        P2d['D']=np.linspace(1.0,10.0,P2d['N'])#Exp-specific code here      
        P2d['constR'] = np.repeat(True,P2d['N'])
        #make go!
        output2d.append(myMPhandler(P2d))#Exp-specific file/variable names here
        output2d = np.array(output2d)#Exp-specific file/variable names here
        if P2d['save_examples']:
            Ex2d=process_examplefiles(P2d)
        np.save('./results/Fig2d.npy',output2d, allow_pickle=True)#Exp-specific file/variable names here
        np.save('./results/Fig2d_P.npy',P2d, allow_pickle=True)#Exp-specific file/variable names here
        np.save('./results/Fig2d_Ex.npy',Ex2d, allow_pickle=True)#Exp-specific file/variable names here
    ###Plot results 
    with PdfPages('./figs/Fig2d.pdf') as pp:#Exp-specific file/variable names here
        fhandle,axhandles=plotresult_1d(output2d,Ex2d,P2d)
        axhandles[0].set_ylim(myylim1)
        axhandles[1].set_ylim(myylim2)
	axhandles[2][0].set_ylim((-65.0,-51.0))
        axhandles[2][1].set_ylim((-65.0,-51.0))
        pp.savefig(fhandle[0])
        pp.savefig(fhandle[1])
        plt.close()



    #####################################
    #####Fig2e:         #################
    #####---------------#################
    ##### Variable -> L & D #############
    ##### R_dend = fix  #################
    ##### No Endbulb    #################
    #####################################
    ####LOAD DATA FOR Fig2e (IF IT EXISTS)
    if os.path.isfile('./results/Fig2e.npy') and Fig2_loadifexists[4]:
        print("Data for Fig2e found... loading!")
        output2e = np.load('./results/Fig2e.npy', allow_pickle=True)
        P2e = np.load('./results/Fig2e_P.npy', allow_pickle=True)
        P2e = P2e.tolist()
    else:
        #Model Parameters (all in a linearly aranged fashion, so that a minimal
        #amount of programming is required to change the experiment).
        P2e={}
        P2e['mp'] = True
        P2e['N']=int(conditions_2d)#just to be sure
        P2e['cores']=int(cores)#just to be sure
        P2e['TotalN'] = int(P2e['N']**2)
        P2e['Number'] = range(P2e['TotalN'])
        P2e['save_examples']=False
       #########################################
        P2e['dur'] = np.repeat(dur,P2e['TotalN'])
        P2e['dt'] = np.repeat(dt,P2e['TotalN'])
        P2e['G_EoH'] = np.repeat(0.0,P2e['TotalN'])
        P2e['Stochastic_EoH'] = np.repeat(False,P2e['TotalN'])
        P2e['Nsyn'] = np.repeat(32,P2e['TotalN'])
        P2e['gsyn'] = np.repeat(0.016,P2e['TotalN'])
        P2e['tau_dendsyn'] = np.repeat(2.0,P2e['TotalN'])
        P2e['constR'] = np.repeat(False,P2e['TotalN'])
        P2e['freq'] = np.repeat(200.0,P2e['TotalN'])
        P2e['IPI'] = np.repeat(5.0,P2e['TotalN'])
        P2e['cf'] = np.repeat(200.0,P2e['TotalN'])
        P2e['dB'] = np.repeat(60.0,P2e['TotalN'])
        P2e['gKLT_d'] = np.repeat(0.0085,P2e['TotalN'])
        P2e['gIH_d'] = np.repeat(0.001,P2e['TotalN'])
        P2e['gLEAK_d'] = np.repeat(0.001,P2e['TotalN'])
        P2e['gKLT_s'] = np.repeat(0.017,P2e['TotalN'])
        P2e['gIH_s'] = np.repeat(0.002,P2e['TotalN'])

        #Now define the two variable parameters. The repeated = y, the tiled = x!!
        allL = np.linspace(10.0,250.0,P2e['N'])
        allD = np.linspace(0.1,10.0,P2e['N'])
        P2e['L'] = np.repeat(allL,P2e['N'])
        P2e['D'] = np.tile(allD,P2e['N'])

        #make go!
        output2e=[]
        output2e.append(myMPhandler(P2e))
        output2e = np.array(output2e)
        np.save('./results/Fig2e.npy',output2e, allow_pickle=True)
        np.save('./results/Fig2e_P.npy',P2e, allow_pickle=True)

    ###Plot results
    fhandle=plotres_2d(output=output2e,P=P2e,x=np.unique(P2e['D']),y=np.unique(P2e['L']),
            xlabs=r'Dendrite1 diameter ($\mu$m)',ylabs=r'Dendrite1 length ($\mu$m)')
    pp = PdfPages('./figs/Fig2e.pdf')
    pp.savefig()
    pp.close()
