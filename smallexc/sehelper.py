#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 10:09:28 2019
This module contains various helper and analysis functions
@author: Thomas Kuenzel (kuenzel{at}bio2.rwth-aachen.de)
"""

import numpy as np

def SimpleDetectAP(V,thr=-100,dt=0.01,LM=-20,RM=10):
    """
    Detect spikes in simulated Vm without knowing the Spt or with many inputs.
    Using a dV/dt threshold of -100mV/ms usually is robust.
    """
    T = np.linspace(0,(len(V)*dt)-dt,len(V))
    dV = np.diff(V)/dt
    Ilow=np.where(dV<thr)[0]
    Ilow = np.concatenate(([0],Ilow))
    dIlow=np.diff(Ilow)
    firstIlow=np.where(dIlow>1.1)[0]
    DetectI=Ilow[firstIlow+1]
    DetectT = T[DetectI]
    PeakI = []
    PeakT = []
    PeakV = []
    for nEv,IEv in enumerate(DetectI):
        if IEv+LM < 0:
            localI=V[0:IEv+RM].argmax()-IEv
            PeakV.append(V[0:IEv+RM].max())
        elif IEv+RM > len(V):
            localI=V[IEv+LM:len(V)].argmax()+LM
            PeakV.append(V[IEv+LM:len(V)].max())
        else:
            localI=V[IEv+LM:IEv+RM].argmax()+LM
            PeakV.append(V[IEv+LM:IEv+RM].max())
        PeakI.append(IEv+localI)
        PeakT.append(T[PeakI[-1]])
            
    Res = {}
    Res['PeakI']=PeakI
    Res['PeakT']=PeakT
    Res['PeakV']=PeakV
    Res['DetectI']=DetectI
    Res['DetectT']=DetectT
    return(Res)

def SimpleGetFails(IN,OUT,R,dur=500.0,dt=0.01,windowsize=1.0):
    isAP = []
    isFail = []
    FV=[]
    FI=[]
    FT=[]
    for isp,sp in enumerate(IN):
        #print isp #DEBUG
        if isp == len(IN)-1:
            ISI_in = dur-sp
        else:
            ISI_in = abs(IN[isp+1]-sp)
        sp_I = int(round(sp/dt))
        winend_I = int(round((sp+ISI_in)/dt))
        #print ISI_in #DEBUG
        #print sp_I#DEBUG
        #print winend_I#DEBUG
        alldelays=OUT['PeakT']-sp
        alldelays[alldelays<0]=999.9
        thisdelay=np.min(alldelays)
        #print thisdelay#DEBUG
        if thisdelay < ISI_in and thisdelay < windowsize:
            #print 'isAP!'#DEBUG
            isAP.append(True)
            isFail.append(False)
        else:
            #print 'isFail!'#DEBUG
            isAP.append(False)
            isFail.append(True)
            FV.append(np.max(R['Vm'][sp_I:winend_I]))
            FI.append(np.argmax(R['Vm'][sp_I:winend_I]))
            FT.append(sp+(FI[-1]*dt))
    R={}
    R['isAP']=isAP
    R['isFail']=isFail
    R['PeakI']=FI
    R['PeakV']=FV
    R['PeakT']=FT
    return R

def get_anglediff(PH0,PHD):
    """
    Calculate correct difference and sign of difference between two angles
    """
    result = np.zeros(len(PHD))
    for iph in range(len(PHD)):
        a = (PH0[iph] - PHD[iph]) % 1.0
        b = (PHD[iph] - PH0[iph]) % 1.0
        if a<b:
            result[iph]=-a
        else:
            result[iph]=b
    return result


def g_template(style='e',dur=2,dt=0.01, 
    tau_rise=0.15, tau_decay=0.2, gmax=0.055,cutshort=False):
    """
    This function returns either excitatory (style = 'e') or inhibitory (style='i') 
    conduction templates of duration dur at sampling interval dt.
    The conduction has an (approximate!) rise 
    and decay tau of tau_rise and tau_decay and a maximal g of gmax. Please note: 
    the tau_rise of the excitatory template is a somewhat complicate issue, 
    since it is not a standard exponential.
    Nevertheless, shorter tau_rise means shorter rise time.
    """
    if dur < tau_decay * 5:
            dur = tau_decay * 5
    if style == 'e':
        MagicNumb = 25
        x = np.linspace(0,dur,dur/dt)
        y_rise = 1/(1+np.exp(-(MagicNumb*x-(MagicNumb*tau_rise))))
        y_rise = y_rise - min(y_rise)
        y_rise = y_rise / max(y_rise)
        y_decay = np.exp((x/ -tau_decay))
        y = y_rise * y_decay
        y = y/max(y)
        y = y*gmax
        if cutshort:
            if len(y) > dur/dt:
                y = y[0:dur/dt]
        return(y)
    elif style == 'i':
        x = np.linspace(0,dur,dur/dt);
        y_rise = (1-np.exp(x/-tau_rise))
        y_decay = np.exp( (x/ -tau_decay))
        y = y_rise * y_decay
        y = y/max(y)
        y = y*gmax
        if cutshort:
            if len(y) > dur/dt:
                y = y[0:dur/dt]
        return(y)
    else:
        print('choose either e or i for template style (give as string!)')

def mysmooth(y,b,wind=21):
    """
    smoothing by convolution with a Kaiser window
    """
    #wind=21
    # extending the data at beginning and at the end
    s = np.r_[y[wind-1:0:-1],y,y[-1:-wind:-1]]
    w = np.kaiser(wind,b)
    newy = np.convolve(w/w.sum(),s,mode='valid')
    return newy[int(np.floor(wind/2)):len(newy)-int(np.floor(wind/2))]

def make_g_trace(dur,dt,spt,template,rnd=False,gceil=0.1,rseed=-1):
    """
    make_g_trace(dur,dt,spt,template,rnd=False,gceil=0.1):
    ...
    This function will return a trace of dur duration (in ms) at sampling 
    interval dt in fractions of ms (=dur/dt samples) where at spt times 
    (spike arrival times in ms) a template shape has been convolved with 
    the trace, i.e. a shape has been inserted.
    If template is a list of numpy.array instead of a single numpy.array, 
    it must contain as many templates as are events in spt -- then members 
    of template will be convoluted with the trace individually when it is 
    their turn.
    """
    if rseed != -1:
        np.random.seed(rseed)
    rawtrace = np.zeros(int(np.round(dur/dt))) 
    if type(template) is np.ndarray:
        for nspt in spt:
            if rnd:
                rawtrace[int(round(nspt/dt))] = 0.2*np.random.randn()+1.0
                #this has a "fixed" SD of ~0.2
            else:
                rawtrace[int(round(nspt/dt))-1] = 1
        trace = np.convolve(rawtrace,template)
        if trace.size > dur/dt:
            trace = trace[0:int(np.round(dur/dt))]
    else:
        print('template must be a single numpy-array! Returning zeros!')
        trace = np.zeros(dur/dt)
    trace[trace > gceil] = gceil
    trace = mysmooth(trace,32)
    return(trace)

def randspiketrain(tstop=400,deadt=1.1,scale=15):
    """
    Returns random spiketimes from drawn from shifted exponential distribution of ISI.
    tstop: largest possible spiketime
    dt: temporal resolution
    deadt: smalles ISI possible
    scale: scale of the exponential function, roughly larger numbers give more long ISI
    """
    spt = []
    lastspt = 0
    while lastspt < tstop:
        ISI = np.random.exponential(scale)+deadt
        thisspt = lastspt+ISI
        lastspt = thisspt
        if thisspt < tstop:
            spt.append(thisspt)
    return(np.array(spt))

def stochastic_pulsetrain(N=16,distritype='normal',DistParam=1.0,
                          IPI=10.0,dur=400.0,dt=0.01,debug=False):
    from scipy.stats import alpha
    #
    #step 1 calc the perfect unjittered train
    sp_raw_t = np.array(range(int(round(IPI/dt)),int(round(dur/dt)),int(round(IPI/dt))))*dt
    allspt=[]
    for ntrain in range(N):
        spt=np.zeros((len(sp_raw_t,)))
        #step 2 jitter every spike according to distritype
        for nspt,sptr in enumerate(sp_raw_t):
            if distritype=='normal':
                jspt=np.random.normal(sptr,DistParam)
                if jspt >= sptr+IPI:
                    spt[nspt] = sptr+((jspt-sptr)%IPI)
                elif jspt <= sptr-IPI:
                    spt[nspt] = sptr-(np.abs(jspt-sptr)%IPI)
                else:
                    spt[nspt]=jspt
            elif distritype=='alpha':
                jspt=alpha.rvs(DistParam,loc=sptr)
                if jspt >= sptr+IPI:
                    spt[nspt]=sptr+(np.random.rand()*IPI)
                else:
                    spt[nspt]=jspt    
            else:
                spt[nspt]=sptr
        allspt.append(np.sort(spt))
    if debug:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.hist(allspt[0]-sp_raw_t)
    return(allspt,sp_raw_t)


def create_sound(fs=100e3,freq=500,duration=0.4,dbspl=60):
    '''
    Create stimulus
    '''
    import thorns.waves as wv
    S = wv.ramped_tone(fs=fs,freq=freq,duration=duration,dbspl=dbspl)
    return(S)
    
def create_spk(S,fs=100e3,N=100,cf=500,seed=0,species='cat',anf_num=[]):
    '''
    Create ANF spiketrains using the cochlea package
    '''
    if len(anf_num)==0:
        anf_num=(0,N,0)
    import cochlea
    T = cochlea.run_zilany2014(S,fs=fs,anf_num=anf_num,cf=cf,species=species,seed=seed)
    return(T)

def spt2cyc(spt,freq,TW):
    """
    helper for vectorstrength
    """
    if len(TW)==2:
        inTW = (spt>TW[0]) * (spt<TW[1])
        spt = spt[inTW]
    wspt = np.exp(2*np.pi*1j*1e-3*spt*freq)
    if len(spt)<2:
        R = 0
        alpha = 1
    Nspt = len(spt)
    cycs = np.angle(wspt)/2/np.pi
    return(wspt,cycs)
    
def vectorstrength(spt,freq,TW):
    """
    This function calculates vectorstrength and Rayleigh-Statistics for an array
    of spiketimes (in ms) given a stimulus frequency in Hz! A timewindow, given
    as TW = (begin,end) in ms can be applied so spiketimes outside this timewindow
    are disregarded.
    """
    spt = np.array(spt)
    if spt.size == 0:
        return(0.0,0.0,1.0,[])
    wspt,cycs = spt2cyc(spt,freq,TW)
    N = len(cycs)
    R = np.mean(wspt)
    VS = np.abs(R)
    phi = np.angle(R)/2/np.pi
    phi = 1+phi-np.round(phi)
    phi = phi%1
    R = simple_rayleighsign(VS,phi,N)
    phases = cycs = 1+cycs - np.min(np.round(cycs))
    phases = phases%1
    return(VS,phi,R,phases)

def simple_rayleighsign(VS,phi,N):
    """
    This follows N.I. Fisher, 'Statistical analysis of circular data', p.70
    """
    Rbar = VS
    phibar = phi*2*np.pi
    Z = N*(Rbar**2)
    P = np.exp(-Z)*(1+(2*Z - Z**2)/(4*N) - (24*Z - 132*Z**2 + 76*Z**3 - 9*Z**4)/(288*N**2))
    return(P)
    
def cyclehist(spt,freq,TW,nbins=33,doplot=False,normed=False):
    """
    Calculate and possibly plot a cylcehistogram.
    
    cyclehist(spt,freq,TW,nbins=33,doplot=False,normed=False):
    """
    import matplotlib.pyplot as plt
    
    spt = np.array(spt)
    wspt,cycs = spt2cyc(spt,freq,TW) 
    cycs = 1+cycs - np.min(np.round(cycs))
    cycs = cycs%1
    N = len(cycs)
    thebins = np.linspace(1.0/nbins,1,nbins)
    vals,edges=np.histogram(cycs,bins=thebins)
    edges = edges[1:]-np.diff(edges)
    vals = np.array(vals)
    if normed:
        vals = vals/np.float(np.max(vals))
    if doplot:
        R = np.mean(wspt)
        VS = np.abs(R)
        phi = np.angle(R)/2/np.pi
        phi = 1+phi-np.round(phi)
        phi = phi%1
        P = simple_rayleighsign(VS,phi,N)
        plt.figure()
        plt.bar(edges,vals,1.0/nbins,align='center',edgecolor='b')
        plt.vlines(phi,0,1,color='r',linewidth=2)
        plt.title( 'Cyclehistogram @' + str(np.round(freq)) + 'Hz, VS=' +str("%.2f" %VS) + 'P=' + str("%.2f" %P) )
    return(edges,vals)
