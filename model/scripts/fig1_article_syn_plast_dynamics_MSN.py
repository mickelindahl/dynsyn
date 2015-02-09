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
SYNAPSE_MODELS=['MSN_SNR_gaba_p1']

HEADER_SIMULATION_SETUP=( '**** BEGINNING GENERAL SCRIPT SETUP ****\n'+
                          'FILE_NAME:'+str(FILE_NAME)+'\n'+                         
                          'PICTURE_PATH:'+str(PICTURE_PATH)+'\n'+  
                          'OUTPUT_PATH:'+str(OUTPUT_PATH)+'\n\n'+    
                          'NEURON_MODEL:'+str(NEURON_MODELS)+'\n'+
                          'SYNAPSE_MODELS:'+str(SYNAPSE_MODELS)+'\n'+
                          '**** END GENERAL SCRIPT SETUP ****\n')

### END GLOBALS ###

def plot_steady_state_freq(ax, freq_, relative_fac ):
    colors = ['b','g','m']
    linestyles=['--','--','-']
    labels=[r'$ref_{init}^{MSN_{D1}}$', r'$ref_{max}^{MSN_{D1}}$' , r'$fac^{MSN_{D1}}$']
    coords=[[0.1, 0.22], [0.2, 0.74], [0.65, 0.45]]
    syn_static=['MSN_SNR_gaba_s_min', 'MSN_SNR_gaba_s_max']
    
    ax.plot(freq_,relative_fac[0,:],**{ 'color':colors[2]})
    ytext=ax.set_ylabel(r'$p_{ss}$/$p_1$') 
    #pylab.setp(ytext, fontsize=14.)
    
    ax.set_xlabel('Firing rate (spikes/s)')
    ax.my_set_no_ticks( yticks=5, xticks = 5 ) 
    
    
    my_nest.ResetKernel()
    model_list, model_dict=models() 
    syns=[SYNAPSE_MODELS[0]]
    syns.extend(syn_static)
    
    my_nest.MyLoadModels( model_list, syns)
    xdata=ax.lines[0].get_xdata()
    ss=my_nest.GetDefaults(SYNAPSE_MODELS[0])  
    synapticEficacy = ss['weight']*ss['U'] 
    
    for syn, color, ls in zip(syn_static, colors[0:2], linestyles[0:2]):
    
        sw=my_nest.GetDefaults(syn)  
        ax.plot([min(xdata), 48],[sw['weight']/synapticEficacy,
                                          sw['weight']/synapticEficacy],
                                          **{'color':color, 'linestyle':ls})

    for coord, label, color in zip(coords, labels, colors):
        ax.text( coord[0], coord[1], label , transform=ax.transAxes, 
                 fontsize=pylab.rcParams['font.size']+2,
                 **{ 'color' : color})  
    
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,48])      

    ax.set_xlim(misc.adjust_limit([0,50]))
    ax.set_ylim(misc.adjust_limit([0,5]))

    ax.my_set_no_ticks( yticks=6, xticks = 6 )

    pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend  
def plot_max_freq(ax, freq_, relative_fac ):
    colors = ['b','g','m']
    labels=[r'$ref_{init}^{MSN_{D1}}$', r'$ref_{max}^{MSN_{D1}}$' , r'$fac^{MSN_{D1}}$']
    coords=[[0.1, 0.02], [0.2, 0.55], [0.1, 0.8]]
    syn_static=['MSN_SNR_gaba_s_min', 'MSN_SNR_gaba_s_max']
    
    ax.plot(freq_,relative_fac[0,:],**{ 'color':colors[2]})
    
    ax.set_ylabel(r'$p_{max}$/$p_1$') 
    ax.set_xlabel('Firing rate (spikes/s)')
    ax.my_set_no_ticks( yticks=5, xticks = 5 ) 
    
    
    my_nest.ResetKernel()
    model_list, model_dict=models() 
    syns=[SYNAPSE_MODELS[0]]
    syns.extend(syn_static)
    
    my_nest.MyLoadModels( model_list, syns)
    xdata=ax.lines[0].get_xdata()
    ss=my_nest.GetDefaults(SYNAPSE_MODELS[0])  
    synapticEficacy = ss['weight']*ss['U'] 
    
    for syn, color in zip(syn_static, colors[0:2]):
    
        sw=my_nest.GetDefaults(syn)  
        ax.plot([min(xdata), max(xdata)],[sw['weight']/synapticEficacy,
                                          sw['weight']/synapticEficacy],
                                          **{'color':color})

    for coord, label, color in zip(coords, labels, colors):
        ax.text( coord[0], coord[1], label , transform=ax.transAxes, 
                 fontsize=pylab.rcParams['font.size']+2,
                 **{ 'color' : color})  
        
    ax.set_xlim(misc.adjust_limit([0,100]))
    ax.set_ylim(misc.adjust_limit([0,5]))

    ax.my_set_no_ticks( yticks=6, xticks = 6 )
    pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend
    
def plot_recovery(ax, revoceryTimes, relativeRecovery ):
    colors = ['m']
    ax.plot(revoceryTimes,relativeRecovery[0,:],**{'color':colors[0]})
    ax.text( 0.3, 0.3, r'$fac^{MSN_{D1}}$',
            transform=ax.transAxes, 
                 fontsize=pylab.rcParams['font.size']+2,
                 **{ 'color' : colors[0]})
    
    ytext=ax.set_ylabel(r'Pair pulse ratio') 
    #pylab.setp(ytext, fontsize=12.)
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
        
        if 'U' in ss.keys():                             
            tb = tb + ' %s %5s %3s \n' % ( 'U:', str ( ss['U'] ), '--' ) 
            tb = tb + ' %s %5s %3s \n' % ( 'tau_fac:', 
                                           str ( ss['tau_fac'] ), 'ms' )                                 
            tb = tb + ' %s %5s %3s \n' % ( 'tau_rec:', 
                                           str ( ss['tau_rec'] ), 'ms' )
            tb = tb + ' %s %5s %3s \n' % ( 'tau_psc:', 
                                           str ( ss['tau_psc'] ), 'ms' ) 
            synapticEficacy = ss['weight']*ss['U'] 
            tb = tb + ' %s %5s %3s \n' % ( 'P1:',synapticEficacy, 'nS' )

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
    
    relative_fac=[]
    n=len(frequencies)
    if not load: 
        for syn in SYNAPSE_MODELS:
            my_nest.ResetKernel()       
            model_list, model_dict=models()   
            my_nest.MyLoadModels( model_list, NEURON_MODELS )
            my_nest.MyLoadModels( model_list, [syn])
            
            ss=my_nest.GetDefaults(syn)       
            synapticEficacy = ss['weight']*ss['U'] 
                
            SNR = MyGroup( NEURON_MODELS[0], n, mm=True, mm_dt = .1, 
                           params={'I_e':-150.}, record_from=['g_GABAA_1'] )

            
            tSim=3*1000/frequencies[0]  
            spikeTimes=[]        
            for f in frequencies :
                isi  = 1000./f
                spikeTimes.append(numpy.arange(1,tSim,isi))
            
            

            for target, st in zip(SNR, spikeTimes ) :
                source = my_nest.Create('spike_generator', 
                                    params={'spike_times':st} )
                my_nest.SetDefaults(syn, params={'delay':1.})
                my_nest.Connect(source, [target], model=syn)
        
            my_nest.MySimulate(tSim)
            SNR.get_signal( 'g','g_GABAA_1', stop=tSim ) # retrieve signal
            SNR.save_signal( 'g','g_GABAA_1', stop=tSim )
            
            SNR.load_signal( 'g','g_GABAA_1')
    
            signal=SNR.signals['g_GABAA_1']
            
            tmpss=[]
            for i, st in enumerate(spikeTimes, start=1):
                
                if SNR.mm_dt==0.1:  indecies=numpy.int64(numpy.ceil(st*10))+9
                elif SNR.mm_dt==1.: indecies=numpy.int64(numpy.ceil(st))
                
                values=signal[i].signal[indecies]-signal[i].signal[indecies-1]
                
                if flag=='ss':  tmpss.append(values[-1]/synapticEficacy)
                if flag=='max': tmpss.append(max(values)/synapticEficacy)
                
            relative_fac.append(tmpss)
            
        relative_fac=numpy.array(relative_fac)
        max_rel_fac=str(max(relative_fac[0]))
        s='\n'
        s = s + ' %s %5s %3s \n' % ( flag+' Max:', max_rel_fac[0:6],  '--' )   
        info_string=s  
        
        header=HEADER_SIMULATION_SETUP+s
        misc.text_save(header, save_header_at)                        
        misc.pickle_save([frequencies, relative_fac, s], save_result_at)

        
    elif load: 
            revoceryTimes, relative_fac, info_string=misc.pickle_load(save_result_at)       
    
    return frequencies, relative_fac, info_string

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
            
    
    
            SNR = MyGroup( NEURON_MODELS[0], n, mm=True, mm_dt = .1, 
                           params={'I_e':-150.}, record_from=['g_GABAA_1'] )

            
            tSim=5000
            spikeTimes=[]
            for rt in revoceryTimes:
                spikeTimes.append(numpy.array([1.,11.,21.,31.,41.,41+rt]))
                
                
            for target, st in zip(SNR, spikeTimes ) :
       
                source = my_nest.Create('spike_generator', 
                                    params={'spike_times':st} )
                my_nest.SetDefaults(syn, params={'delay':1.})
                my_nest.Connect(source, [target], model=syn)
        
            my_nest.MySimulate(tSim)
            SNR.get_signal( 'g','g_GABAA_1', stop=tSim ) # retrieve signal
            
            signal=SNR.signals['g_GABAA_1']
            
            ss=my_nest.GetDefaults(syn)       
            synapticEficacy = ss['weight']*ss['U'] 
            
            tmpss=[]
            for i, st in enumerate(spikeTimes, start=1):
                
                if SNR.mm_dt==0.1:  indecies=numpy.int64(numpy.ceil(st*10))+9
                elif SNR.mm_dt==1.: indecies=numpy.int64(numpy.ceil(st))
                
                values=signal[i].signal[indecies]-signal[i].signal[indecies-1]
                
                tmpss.append(values[-1]/synapticEficacy)
                #tmpss.append(max(values)/synapticEficacy)
                
            relativeRecovery.append(tmpss)
            
        relativeRecovery=numpy.array(relativeRecovery) 
    
        header=HEADER_SIMULATION_SETUP
        misc.text_save(header, save_header_at)    
        misc.pickle_save([revoceryTimes, relativeRecovery], save_result_at)
   
    elif load: 
            revoceryTimes, relativeRecovery=misc.pickle_load(save_result_at)    

    return revoceryTimes, relativeRecovery

print 'Simulation'

# SIMULATION
info_string=''
# Steady state conductance all freq_
freq_all, ss_fac_all, s = simulate_steady_state_freq(numpy.arange(1.,100.,1.), load=True)
info_string=info_string+s

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
plot_text(ax, info_string)

# Example steady state all
ax=ax_list[1]
plot_steady_state_freq(ax, freq_all, ss_fac_all)

# Example steady state zoom
#ax=ax_list[2]
#plot_steady_state_freq_zoom(ax, freq_Zoom, relative_facZoom)

# Max
#ax=ax_list[2]
#plot_max_freq(ax, freq_all, max_fac_all)

# Example recovery spike
ax=ax_list[2]
plot_recovery(ax, revoceryTimes, relativeRecovery)
#ax.legend(numpoints=1, loc='best')

pylab.show()

fig.savefig( PICTURE_PATH + '/' + FILE_NAME  + '.svg', format = 'svg')
fig.savefig( PICTURE_PATH + '/' + FILE_NAME  + '.pdf', format = 'pdf')