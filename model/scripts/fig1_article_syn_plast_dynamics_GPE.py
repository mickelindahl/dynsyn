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

sys.path.append(model_dir) 
sys.path.append(code_dir+'/nest_toolbox') 

from model_params import models                               # Then import models     
from src import misc, my_nest, my_signals, plot_settings
from src.my_population import MyGroup 
from src.my_axes import MyAxes 

### GLOBALS ###

# Paths and naming for saving data and picutes
FILE_NAME = sys.argv[0].split('/')[-1].split('.')[0]
model_name=os.getcwd().split('/')[-2]
PICTURE_PATH='/'.join(os.getcwd().split('/')[0:-3]) + '/pictures/'+model_name    
OUTPUT_PATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]

# Models
NEURON_MODELS=['SNR_aeif']
SYNAPSE_MODELS=['GPE_SNR_gaba_p']

HEADER_SIMULATION_SETUP=( '**** BEGINNING GENERAL SCRIPT SETUP ****\n'+
                          'FILE_NAME:'+str(FILE_NAME)+'\n'+                         
                          'PICTURE_PATH:'+str(PICTURE_PATH)+'\n'+  
                          'OUTPUT_PATH:'+str(OUTPUT_PATH)+'\n\n'+    
                          'NEURON_MODEL:'+str(NEURON_MODELS)+'\n'+
                          'SYNAPSE_MODELS:'+str(SYNAPSE_MODELS)+'\n'+
                          '**** END GENERAL SCRIPT SETUP ****\n')

### END GLOBALS ###

def plot_steady_state_freq(ax, frequencies, relativeFacilitation ):
    colors = ['r','c']
    linestyles=['--','-']
    labels=[r'$ref_{30 Hz}^{GPe}$', r'$dep^{GPe}$']
    coords=[[0.35, 0.2], [0.1, 0.7]]
    syn_static=['GPE_SNR_gaba_s_ref']

    #print relativeFacilitation[0, 23:26]
    #print frequencies[23:26]

    ax.plot(frequencies,relativeFacilitation[0,:],**{'color':colors[1]})
    my_nest.ResetKernel()
    model_list, model_dict=models() 
    syns=[SYNAPSE_MODELS[0]]
    syns.extend(syn_static)
    
    my_nest.MyLoadModels( model_list, syns )
    xdata=ax.lines[0].get_xdata()
    ss=my_nest.GetDefaults(SYNAPSE_MODELS[0])  
    synapticEficacy = ss['weight']*ss['U'] 
    sw=my_nest.GetDefaults('GPE_SNR_gaba_s_ref')      
    ax.plot([min(xdata), max(xdata)],[sw['weight']/synapticEficacy,
                                      sw['weight']/synapticEficacy],
                                      **{'color':colors[0],'linestyle':linestyles[0]})
    t_a = ax.transAxes
    
    for coord, label, color in zip(coords, labels, colors):
        ax.text( coord[0], coord[1], label , transform=ax.transAxes, 
                 fontsize=pylab.rcParams['font.size']+2,
                 **{ 'color' : color})  
    
    
    ytext=ax.set_ylabel(r'$p_{ss}$/$p_1$') 
    pylab.setp(ytext, fontsize=16.)
    ax.set_xlabel('Firing rate (spikes/s)')
    ax.my_set_no_ticks( yticks=6, xticks = 6 ) 
      
    ax.set_xlim(misc.adjust_limit([0,100]))
    ax.set_ylim(misc.adjust_limit([0,1]))
    pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend
    
def plot_recovery(ax, revoceryTimes, relativeRecovery ):
    colors = ['c']
    ax.plot(revoceryTimes,relativeRecovery[0,:],**{'color':colors[0]})
    ax.text( 0.2, 0.35, r'$dep^{GPe}$',
            transform=ax.transAxes, 
                 fontsize=pylab.rcParams['font.size']+2,
                 **{ 'color' : colors[0]})


    ytext=ax.set_ylabel(r'Pair pulse ratio') 
    #pylab.setp(ytext, fontsize=16.)
    
    ax.set_xlabel('Inter stimulus interval (ms)')

    ax.set_xlim(misc.adjust_limit([0,5000]))
    ax.set_ylim(misc.adjust_limit([0,4]))
    ax.my_set_no_ticks( yticks=6, xticks = 5 )
    pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend
    
def plot_text(ax, infoString=''):
    
    my_nest.ResetKernel()
    model_list, model_dict=models()
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    my_nest.MyLoadModels( model_list, SYNAPSE_MODELS )
    my_nest.MyLoadModels( model_list, ['GPE_SNR_gaba_s_ref','GPE_SNR_gaba_s_max'] )
    
    SNR = MyGroup( NEURON_MODELS[0], 2, mm_dt = 0.1)
    statusSNR = my_nest.GetStatus( SNR[:] )[0]
    
    statusSynapes=[]
    for s in SYNAPSE_MODELS:
        statusSynapes.append( my_nest.GetDefaults(s) )
    tb = ''     
    tb = tb + infoString
    
    for ss in statusSynapes:
        tb = tb + '\n'  
        tb = tb + ' %s %10s\n' % ( 'Synapse', ss['synapsemodel'])   
        tb = tb + ' %s %5s %3s\n' % ( 'Weight', 
                                      str( round( ss['weight'], 1) ), 'nS')
                                    
        tb = tb + ' %s %5s %3s \n' % ( 'U:', str ( ss['U'] ), '--' ) 
        tb = tb + ' %s %5s %3s \n' % ( 'tau_fac:', 
                                           str ( ss['tau_fac'] ), 'ms' )                                 
        tb = tb + ' %s %5s %3s \n' % ( 'tau_rec:', 
                                           str ( ss['tau_rec'] ), 'ms' )
        tb = tb + ' %s %5s %3s \n' % ( 'tau_psc:', 
                                           str ( ss['tau_psc'] ), 'ms' ) 
        synapticEficacy = ss['weight']*ss['U'] 
        tb = tb + ' %s %5s %3s \n' % ( 'P1:',synapticEficacy, 'nS' )
    
    sw=my_nest.GetDefaults('GPE_SNR_gaba_s_ref')   
    tb = tb + ' %s %5s %3s \n' % ( r'$ref_{32 Hz}^{GPe}$ fraction', 
                                           str ( sw['weight']/synapticEficacy), '--' ) 
    sw=my_nest.GetDefaults('GPE_SNR_gaba_s_max')
    tb = tb + ' %s %5s %3s \n' % ( r'$dep^{GPe}$ max', 
                                           str ( sw['weight']), '--' )
    ax.text( 0.85, 0.5, tb , fontsize= font_size_text,
             horizontalalignment='right',
             verticalalignment='center',
             transform=ax.transAxes,     # to define coordinates in right scale
             **{ 'fontname' : 'monospace' })                           
    
    ax.my_remove_axis( xaxis=True, yaxis=True )
    ax.my_remove_spine(left=True,  bottom=True, right=True, top=True)

def simulate_steady_state_freq(frequencies, flag='ss', load=True):
    
     # Path were raw data is saved. For example the spike trains.
    save_result_at=OUTPUT_PATH+'/simulate_steady_state_freq.pkl'
    save_header_at=OUTPUT_PATH+'/simulate_steady_state_freq_header'   
    
    relativeFacilitation=[]
    n=len(frequencies)
    if not load:    
        for syn in SYNAPSE_MODELS:
            my_nest.ResetKernel()   
            model_list, model_dict=models()    
            my_nest.MyLoadModels( model_list, NEURON_MODELS )
            my_nest.MyLoadModels( model_list, [syn])
            
                    
            SNR = MyGroup( NEURON_MODELS[0], n, mm=True, mm_dt = .1, 
                               params={'I_e':-150.},
                           record_from=['g_GABAA_2'] )

            tSim=5*1000/frequencies[0]  
            spikeTimes=[]    
            tmpSteadyState=[]    
            for f in frequencies :

                isi  = 1000./f
                spikeTimes.append(numpy.arange(1,tSim,isi))
            
            for target, st in zip(SNR, spikeTimes ) :
                    source = my_nest.Create('spike_generator', 
                                        params={'spike_times':st} )
                    my_nest.SetDefaults(syn, params={'delay':1.})
                    my_nest.Connect(source, [target], model=syn)
            
            my_nest.MySimulate(tSim)
                
            SNR.get_signal( 'g','g_GABAA_2', stop=tSim ) # retrieve signal
                  
            signal=SNR.signals['g_GABAA_2']
                
                
            for i, st in enumerate(spikeTimes, start=1):
                    
                    if SNR.mm_dt==0.1:  indecies=numpy.int64(numpy.ceil(st*10))+9
                    elif SNR.mm_dt==1.: indecies=numpy.int64(numpy.ceil(st))
                    
                    values=signal[i].signal[indecies]-signal[i].signal[indecies-1]
                    
                    ss=my_nest.GetDefaults(syn)       
                    synapticEficacy = ss['weight']*ss['U'] 
                    
                    if flag=='ss':  tmpSteadyState.append(values[-1]/synapticEficacy)
                    if flag=='max': tmpSteadyState.append(max(values)/synapticEficacy)
                    
            relativeFacilitation.append(tmpSteadyState)
            
        relativeFacilitation=numpy.array(relativeFacilitation)
        
    
        header=HEADER_SIMULATION_SETUP
        misc.text_save(header, save_header_at)
        misc.pickle_save([frequencies, relativeFacilitation], save_result_at)

        
    elif load: 
            frequencies, relativeFacilitation=misc.pickle_load(save_result_at)
        
    return frequencies, relativeFacilitation

def simulate_recovery(revoceryTimes, load=True):
    
    # Path were raw data is saved. For example the spike trains.
    save_result_at=OUTPUT_PATH+'/simulate_recovery.pkl'
    save_header_at=OUTPUT_PATH+'/simulate_recovery_header'   
    
    relativeRecovery=[]
    n=len(revoceryTimes)
    if not load:
        for syn in SYNAPSE_MODELS:
            my_nest.ResetKernel()  
            model_list, model_dict=models()     
            my_nest.MyLoadModels( model_list, NEURON_MODELS )
            my_nest.MyLoadModels( model_list, [syn])
            
            ss=my_nest.GetDefaults(syn)       
            synapticEficacy = ss['weight']*ss['U'] 
    
            SNR = MyGroup( NEURON_MODELS[0], n, mm=True, mm_dt = .1, 
                           params={'I_e':-150.}, record_from=['g_GABAA_2'])
            
            tSim=10000
            spikeTimes=[]
            for rt in revoceryTimes:
                #spikeTimes.append(numpy.array([1.,11.,21.,31.,41.,41+rt]))
                
                # Choosen so that it starts at a pairpulse ration of 0.2 
                spikeTimes.append(numpy.array([1.,11.,21.,31.,41.,
                                               51.,61.,71.,81.,91.,
                                               101.,111.,121.,131.,141.,
                                               151.,161.,171.,181.,191.,
                                               191+rt]))
     
            for target, st in zip(SNR, spikeTimes ) :
       
                source = my_nest.Create('spike_generator', 
                                    params={'spike_times':st} )
                my_nest.SetDefaults(syn, params={'delay':1.})
                my_nest.Connect(source, [target], model=syn)
        
            my_nest.MySimulate(tSim)
            SNR.get_signal( 'g','g_GABAA_2', stop=tSim ) # retrieve signal
            
            signal=SNR.signals['g_GABAA_2']
            
            tmpSteadyState=[]
            for i, st in enumerate(spikeTimes, start=1):
                
                if SNR.mm_dt==0.1:  indecies=numpy.int64(numpy.ceil(st*10))+9
                elif SNR.mm_dt==1.: indecies=numpy.int64(numpy.ceil(st))
                
                values=signal[i].signal[indecies]-signal[i].signal[indecies-1]
                
                tmpSteadyState.append(values[-1]/synapticEficacy)
                #tmpSteadyState.append(max(values)/synapticEficacy)
                
            relativeRecovery.append(tmpSteadyState)
            
        relativeRecovery=numpy.array(relativeRecovery)
        
        
        header=HEADER_SIMULATION_SETUP
        misc.text_save(header, save_header_at)    
        misc.pickle_save([revoceryTimes, relativeRecovery], save_result_at)
        

        
    elif load: 
            revoceryTimes, relativeRecovery=misc.pickle_load(save_result_at)
        
    return revoceryTimes, relativeRecovery

print 'Simulation'

# SIMULATION

# Steady state conductance all frequencies
frequenciesAll, steadyStateFacilitationAll = simulate_steady_state_freq(numpy.arange(1.,100.,1.), 
                                                                        load=True)
# Recovery spike
revoceryTimes, relativeRecovery = simulate_recovery(numpy.arange(100,5000,200), load=True)

# DISPLAY
plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0, fontsize=16)
font_size_text = 8
fig = pylab.figure( facecolor = 'w' )

ax_list = []
ax_list.append( MyAxes(fig, [ .075, .37, .135, .26 ] ) )    # text box
ax_list.append( MyAxes(fig, [ .26,  .6,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .53, .6,  .165, .34 ] ) )    #  
#ax_list.append( MyAxes(fig, [ .8,   .6,  .165-0.06, .34-0.12 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .26,  .1,  .165, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .53, .1,  .165, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .8,   .1,  .165, .34 ] ) )    # 

# Text
ax=ax_list[0]
plot_text(ax)

# Example steady state all
ax=ax_list[1]
plot_steady_state_freq(ax, frequenciesAll, steadyStateFacilitationAll)

# Example recovery spike
ax=ax_list[2]
plot_recovery(ax, revoceryTimes, relativeRecovery)


pylab.show()

fig.savefig( PICTURE_PATH + '/' + FILE_NAME  + '.svg', format = 'svg')
fig.savefig( PICTURE_PATH + '/' + FILE_NAME  + '.pdf', format = 'pdf')
