#! Imports
import math
import numpy
import pylab
import os
import sys


# Add directories to python path
sys.path.append(os.getcwd())                            
parent_dir='/'.join(os.getcwd().split('/')[0:-1])       
                   
model_dir=   '/'.join(os.getcwd().split('/')[0:-1])        
code_dir=  '/'.join(os.getcwd().split('/')[0:-2]) 
picture_dir=  '/'.join(os.getcwd().split('/')[0:-3]) + '/pictures'     
                
sys.path.append(model_dir) 
sys.path.append(code_dir+'/nest_toolbox') 
spath  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]
siter=0

from model_params import models                               # Then import models     
from src import misc, my_nest, my_signals, plot_settings
from src.my_population import MyGroup 
from src.my_axes import MyAxes 

neuronModels=['SNR_aeif']
synapseModels=[ 'MSN_SNR_gaba_p0', 'MSN_SNR_gaba_p1','MSN_SNR_gaba_p2']

# Simulate or use stored data
LOAD=True

def plot_steady_state_freq(ax, frequencies, relativeFacilitation ):
    colors = misc.make_N_colors('Blues', 5)
    ax.plot(frequencies,relativeFacilitation[0,:],**{'label':'Data 1',
                                                     'color':colors[1]})  
    ax.plot(frequencies,relativeFacilitation[1,:],**{'label':'Data 1+2',
                                                     'color':colors[2]})
    ax.plot(frequencies,relativeFacilitation[2,:],**{'label':'Data 2',
                                                     'color':colors[3]})
    ax.set_ylabel('Pss/P1') 
    ax.set_xlabel('Firing rate (spikes/s)')
    ax.my_set_no_ticks( yticks=4, xticks = 4 ) 
    ax.text( 0.75, 0.7,'Set 1' , transform=ax.transAxes, **{ 'color' : colors[1] }) 
    ax.text( 0.4, 0.4, 'Set 1+2' , transform=ax.transAxes, **{ 'color' : colors[2] }) 
    ax.text( 0.2, 0.2, 'Set 2' , transform=ax.transAxes, **{ 'color' : colors[3] }) 

def plot_steady_state_freq_zoom(ax, frequencies, relativeFacilitation ):
    colors = misc.make_N_colors('Blues', 5)
    ax.plot(frequencies, relativeFacilitation[0,:],**{'label':'Data 1',
                                                     'color':colors[1]})  
    ax.plot(frequencies, relativeFacilitation[1,:],**{'label':'Data 1+2',
                                                     'color':colors[2]})
    ax.plot(frequencies, relativeFacilitation[2,:],**{'label':'Data 2',
                                                     'color':colors[3]})
    ax.set_ylabel('Pss/P1') 
    ax.set_xlabel('Firing rate (spikes/s)')
    ax.my_set_no_ticks( yticks=4, xticks = 4 ) 
    ax.text( 0.7, 0.08,'Set 1' , transform=ax.transAxes, **{ 'color' : colors[1] }) 
    ax.text( 0.6, 0.4, 'Set 1+2' , transform=ax.transAxes, **{ 'color' : colors[2] }) 
    ax.text( 0.5, 0.7, 'Set 2' , transform=ax.transAxes, **{ 'color' : colors[3] })  


def plot_max_freq(ax, frequencies, relativeFacilitation ):
    colors = misc.make_N_colors('Blues', 5)
    ax.plot(frequencies,relativeFacilitation[0,:],**{'label':'Data 1',
                                                     'color':colors[1]})  
    ax.plot(frequencies,relativeFacilitation[1,:],**{'label':'Data 1+2',
                                                     'color':colors[2]})
    ax.plot(frequencies,relativeFacilitation[2,:],**{'label':'Data 2',
                                                     'color':colors[3]})
    ax.set_ylabel('Pmax/P1') 
    ax.set_xlabel('Firing rate (spikes/s)')
    ax.my_set_no_ticks( yticks=4, xticks = 4 ) 
    ax.text( 0.4, 0.4,'Set 1' , transform=ax.transAxes, **{ 'color' : colors[1] }) 
    ax.text( 0.4, 0.3, 'Set 1+2' , transform=ax.transAxes, **{ 'color' : colors[2] }) 
    ax.text( 0.4, 0.2, 'Set 2' , transform=ax.transAxes, **{ 'color' : colors[3] }) 

def plot_recovery(ax, revoceryTimes, relativeRecovery ):
    colors = misc.make_N_colors('Blues', 5)
    ax.plot(revoceryTimes,relativeRecovery[0,:],**{'label':'Data 1',
                                                 'color':colors[1]})  
    ax.plot(revoceryTimes,relativeRecovery[1,:],**{'label':'Data 1+2',
                                                 'color':colors[2]})
    ax.plot(revoceryTimes,relativeRecovery[2,:],**{'label':'Data 2',
                                                 'color':colors[3]})
    ax.set_ylabel('Prec/P1') 
    ax.set_xlabel('Time (ms)')
    ax.my_set_no_ticks( yticks=4, xticks = 4 ) 
    ax.text( 0.4, 0.5,'Set 1' , transform=ax.transAxes, **{ 'color' : colors[1] }) 
    ax.text( 0.4, 0.4, 'Set 1+2' , transform=ax.transAxes, **{ 'color' : colors[2] }) 
    ax.text( 0.4, 0.3, 'Set 2' , transform=ax.transAxes, **{ 'color' : colors[3] }) 
    
def plot_text(ax, infoString=''):
    
    my_nest.ResetKernel()
    model_list=models()
    my_nest.MyLoadModels( model_list, neuronModels )
    my_nest.MyLoadModels( model_list, synapseModels )
    
    SNR = MyGroup( neuronModels[0], 2, mm_dt = 0.1)
    statusSNR = my_nest.GetStatus( SNR[:] )[0]
    
    statusSynapes=[]
    for s in synapseModels:
        statusSynapes.append( my_nest.GetDefaults(s) )
    tb = ''     
    tb = tb + infoString
    
    for ss in statusSynapes:
        tb = tb + '\n'  
        tb = tb + ' %s %10s\n' % ( 'Synapse', ss['synapsemodel'])   
        tb = tb + ' %s %5s %3s\n' % ( 'Weight', 
                                      str( round( ss['weight'], 1) ), 'nS')
        
        if 'U' in ss.keys():                             
            tb = tb + ' %s %5s %3s \n' % ( 'U:', str ( ss['U'] ), '--' ) 
            tb = tb + ' %s %5s %3s \n' % ( 'tau_fac:', 
                                           str ( ss['tau_fac'] ), 'ms' )                                 
            tb = tb + ' %s %5s %3s \n' % ( 'tau_rec:', 
                                           str ( ss['tau_rec'] ), 'ms' )
            tb = tb + ' %s %5s %3s \n' % ( 'tau_psc:', 
                                           str ( ss['tau_psc'] ), 'ms' ) 

    ax.text( 0.85, 0.5, tb , fontsize= font_size_text,
             horizontalalignment='right',
             verticalalignment='center',
             transform=ax.transAxes,     # to define coordinates in right scale
             **{ 'fontname' : 'monospace' })                           
    
    ax.my_remove_axis( xaxis=True, yaxis=True )
    ax.my_remove_spine(left=True,  bottom=True, right=True, top=True)

def simulate_steady_state_freq(frequencies, flag='ss'):
    global siter
    
    relativeFacilitation=[]
    model_list=models()
    data={}
    n=len(frequencies)
    
    for syn in synapseModels:
        my_nest.ResetKernel()       
        my_nest.MyLoadModels( model_list, neuronModels )
        my_nest.MyLoadModels( model_list, [syn])
        
        ss=my_nest.GetDefaults(syn)       
        synapticEficacy = ss['weight']*ss['U'] 
            
        SNR = MyGroup( neuronModels[0], n, mm_dt = .1, params={'I_e':-150.},
                       record_from=['g_GABAA_1'], spath=spath, 
                       siter=str(siter)  )
        siter+=1
        
        tSim=3*1000/frequencies[0]  
        spikeTimes=[]        
        for f in frequencies :
            isi  = 1000./f
            spikeTimes.append(numpy.arange(1,tSim,isi))
        
        
        if not LOAD:
            for target, st in zip(SNR, spikeTimes ) :
                source = my_nest.Create('spike_generator', 
                                    params={'spike_times':st} )
                my_nest.SetDefaults(syn, params={'delay':1.})
                my_nest.Connect(source, [target], model=syn)
        
            my_nest.MySimulate(tSim)
            SNR.get_signal( 'g','g_GABAA_1', stop=tSim ) # retrieve signal
            SNR.save_signal( 'g','g_GABAA_1', stop=tSim )
        
        elif LOAD: 
            SNR.load_signal( 'g','g_GABAA_1')

        signal=SNR.signals['g_GABAA_1']
        
        tmpSteadyState=[]
        for i, st in enumerate(spikeTimes, start=1):
            
            if SNR.mm_dt==0.1:  indecies=numpy.int64(numpy.ceil(st*10))+9
            elif SNR.mm_dt==1.: indecies=numpy.int64(numpy.ceil(st))
            
            values=signal[i].signal[indecies]-signal[i].signal[indecies-1]
            
            if flag=='ss':  tmpSteadyState.append(values[-1]/synapticEficacy)
            if flag=='max': tmpSteadyState.append(max(values)/synapticEficacy)
            
        relativeFacilitation.append(tmpSteadyState)
        
    relativeFacilitation=numpy.array(relativeFacilitation)
        
    return frequencies, relativeFacilitation

def simulate_recovery(revoceryTimes):
    global siter
    
    relativeRecovery=[]
    model_list=models()
    data={}
    n=len(revoceryTimes)
    
    for syn in synapseModels:
        my_nest.ResetKernel()       
        my_nest.MyLoadModels( model_list, neuronModels )
        my_nest.MyLoadModels( model_list, [syn])
        
        ss=my_nest.GetDefaults(syn)       
        synapticEficacy = ss['weight']*ss['U'] 

        SNR = MyGroup( neuronModels[0], n, mm_dt = .1, params={'I_e':-150.},
                       record_from=['g_GABAA_1'], spath=spath, 
                       siter=str(siter) )
        siter+=1
        
        tSim=5000
        spikeTimes=[]
        for rt in revoceryTimes:
            spikeTimes.append(numpy.array([1.,11.,21.,31.,41.,41+rt]))
            
        if not LOAD:
            for target, st in zip(SNR, spikeTimes ) :
       
                source = my_nest.Create('spike_generator', 
                                    params={'spike_times':st} )
                my_nest.SetDefaults(syn, params={'delay':1.})
                my_nest.Connect(source, [target], model=syn)
        
            my_nest.MySimulate(tSim)
            SNR.get_signal( 'g','g_GABAA_1', stop=tSim ) # retrieve signal
            SNR.save_signal( 'g','g_GABAA_1', stop=tSim )
        
        elif LOAD: 
            SNR.load_signal( 'g','g_GABAA_1')
        
        signal=SNR.signals['g_GABAA_1']
        
        tmpSteadyState=[]
        for i, st in enumerate(spikeTimes, start=1):
            
            if SNR.mm_dt==0.1:  indecies=numpy.int64(numpy.ceil(st*10))+9
            elif SNR.mm_dt==1.: indecies=numpy.int64(numpy.ceil(st))
            
            values=signal[i].signal[indecies]-signal[i].signal[indecies-1]
            
            tmpSteadyState.append(values[-1]/synapticEficacy)
            #tmpSteadyState.append(max(values)/synapticEficacy)
            
        relativeRecovery.append(tmpSteadyState)
        
    relativeRecovery=numpy.array(relativeRecovery)
        
    return revoceryTimes, relativeRecovery

print 'Simulation'

# SIMULATION

# Steady state conductance all frequencies
frequenciesAll, steadyStateFacilitationAll = simulate_steady_state_freq(numpy.arange(1.,100.,1.))

# Max conductance all frequencies
if LOAD:
    siter=0
    frequenciesAll, maxFacilitationAll = simulate_steady_state_freq(numpy.arange(1.,100.,1.), 
                                                                    flag='max')

# Recovery spike
revoceryTimes, relativeRecovery = simulate_recovery(numpy.arange(100,5000,200))

# Steady state conductance 0-1 frequencies
frequenciesZoom, relativeFacilitationZoom = simulate_steady_state_freq(numpy.arange(0.1,1.,.1))

# DISPLAY
plot_settings.set_mode(mode='by_fontsize', w = 750.0, h = 400.0, fontsize=12)
font_size_text = 9
fig = pylab.figure( facecolor = 'w' )

ax_list = []
ax_list.append( MyAxes(fig, [ .1, .4, .18,  .26 ] ) )    # text box
ax_list.append( MyAxes(fig, [ .35, .6, .24, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .73, .6, .24, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .35, .1, .24, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .73, .1, .24, .34 ] ) )    # 

# Text
ax=ax_list[0]
plot_text(ax)

# Example steady state all
ax=ax_list[1]
plot_steady_state_freq(ax, frequenciesAll, steadyStateFacilitationAll)

# Example steady state zoom
ax=ax_list[2]
plot_steady_state_freq_zoom(ax, frequenciesZoom, relativeFacilitationZoom)

# Max
ax=ax_list[3]
plot_max_freq(ax, frequenciesAll, maxFacilitationAll)

# Example recovery spike
ax=ax_list[4]
plot_recovery(ax, revoceryTimes, relativeRecovery)
#ax.legend(numpoints=1, loc='best')

pylab.show()

name = sys.argv[0].split('/')[-1].split('.')[0]
fig.savefig( picture_dir + '/' + name  + '.svg', dpi = 500, format = 'svg')
