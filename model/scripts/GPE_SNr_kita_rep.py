#! Imports
import math
import numpy
import pylab
import os
import sys
import time

if len(sys.argv) != 1: mpiRun = True
else:                  mpiRun = False
start = time.time() 

 
# Add directories to python path
sys.path.append(os.getcwd())                            
parent_dir='/'.join(os.getcwd().split('/')[0:-1])       
                   
model_dir=   '/'.join(os.getcwd().split('/')[0:-1])        
code_dir=  '/'.join(os.getcwd().split('/')[0:-2]) 
picture_dir=  '/'.join(os.getcwd().split('/')[0:-3]) + '/pictures'     
                
sys.path.append(model_dir) 
sys.path.append(code_dir+'/nest_toolbox') 
SPATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]
SNAME  = sys.argv[0].split('/')[-1].split('.')[0]

from model_params import models                               # Then import models     
from src import misc, my_nest, my_signals, plot_settings
from src.my_population import MyGroup, MyPoissonInput 
from src.my_axes import MyAxes 

SELECTION_THR=2.  # Hz
NEURON_MODELS=['SNR_aeif']
SYNAPSE_MODELS=['GPE_SNR_gaba_s_min', 'GPE_SNR_gaba_s_max', 'GPE_SNR_gaba_p']
SEL_ONSET = 2000.
ADJUST_XDATA_MS=1800.
I_E=-5.+70+120
    
def plot_example_SNR(ax, SNR_list):
    time_bin=5
    
    colors = misc.make_N_colors('Blues', 5)
    colors=['g','m']   
    labels=['Ref','Dyn GPe']
    
    SNR_list=[SNR_list[0], SNR_list[2]]
    
    
    for color, label, SNR in zip(colors, labels, SNR_list):
        signal=SNR.signals['spikes']
        signal.my_firing_rate(bin=time_bin, display=ax,
                          kwargs={'color':color})
    
    lines = ax.lines
    lines[0].set_xdata(lines[0].get_xdata()-ADJUST_XDATA_MS)
    lines[1].set_xdata(lines[1].get_xdata()-ADJUST_XDATA_MS)
    
    ax.set_title('bin=%i'%(time_bin),**{'fontsize':12})
    ax.set_ylabel('Frequency SNr (Hz)') 
    ax.set_xlabel('Time (ms)')
    #ax.set_ylim([0,55])
    ax.my_set_no_ticks( yticks=4, xticks = 4 ) 
    ax.set_xlim([0, 500])
    t_a = ax.transAxes
    ax.text( 0.1, 0.45, 'Ref' , transform=t_a, **{ 'color' : colors[0] })  
    #ax.text( 0.05, 0.25, 'Strong' , transform=t_a, **{ 'color' : colors[1] }) 
    ax.text( 0.1, 0.35, 'Dyn GPe' , transform=t_a, **{ 'color' : 'm' }) 
 
def plot_example_raster_GPE(ax, GPE_list):
    global ADJUST_XDATA_MS
    time_bin=100
    ax_twinx=ax.my_twinx()
    GPE=GPE_list[0]
    GPE.signals['spikes'].raster_plot(id_list=GPE, display=ax_twinx,
                                            kwargs={'color':'k', 'zorder':1})  
    
    line1 = ax_twinx.lines[0]
    #ax_twinx.set_ylim([0,500])
    ax_twinx.set_ylabel('Neuron id')

    GPE.signals['spikes'].my_firing_rate( bin=time_bin, display=ax,
                                          kwargs={'color':'k', 'linewidth':3,})
    line2 = ax.lines[0]
    
    line1.set_xdata(line1.get_xdata()-ADJUST_XDATA_MS)
    line2.set_xdata(line2.get_xdata()-ADJUST_XDATA_MS)


    ax.set_xlim([0,1000])
    ax_twinx.set_xlim([0,1000])
    ax.set_ylim([0,40])
    ax.my_set_no_ticks( yticks=4, xticks = 4 ) 
    ax_twinx.my_set_no_ticks( yticks=4, xticks = 4 ) 
    ax.set_title('bin=%i'%(time_bin),**{'fontsize':12})
    ax.set_ylabel('Frequency GPEs (Hz)')     
    #ax.text( 0.7, 0.75, 'All' , transform=ax.transAxes, **{ 'color' : 'k' })  
    #ax.text( 0.7, 0.85, 'Selected' , transform=ax.transAxes, **{ 'color' : 'grey' }) 
    
def plot_text(ax, info_string=''):
    
    my_nest.ResetKernel()
    MODEL_LIST=models()
    my_nest.MyLoadModels( MODEL_LIST, NEURON_MODELS )
    
    SNR = MyGroup( NEURON_MODELS[0], 1, mm_dt = 0.1)
    statusSNR = my_nest.GetStatus( SNR[:] )[0]
    
    tb = ''     
    tb = tb + info_string
    
    tb = tb + '\n'

    ax.text( 0.85, 0.5, tb , fontsize= font_size_text,
             horizontalalignment='right',
             verticalalignment='center',
             transform=ax.transAxes,     # to define coordinates in right scale
             **{ 'fontname' : 'monospace' })                           
    
    ax.my_remove_axis( xaxis=True, yaxis=True )
    ax.my_remove_spine(left=True,  bottom=True, right=True, top=True)

def simulate_example(hz=0, load=True):
    global NEURON_MODELS
    global SNAME
    global SPATH
    global SYNAPSE_MODELS
    global SEL_ONSET
    global I_E
    

    N_EXP = 200
    N_GPE = 50
    N_SEL = 30 # Number of selected GPE
    N_INH = 0 # Number of inhibited GPE
    RATE_BASE = 15 # Base rate
    RATE_SELE = hz # Selection rate 
    RATE_INHI = 0
    SAVE_AT = SPATH+'/'+NEURON_MODELS[0]+'-example.pkl'   
    SEL_TIME = 20.
    sim_time = SEL_TIME+SEL_ONSET+800.
    SNAME_NB = hz+1000    
    
    EXPERIMENTS=range(N_EXP)
   
    
    MODEL_LIST=models()
    my_nest.ResetKernel()       
    my_nest.MyLoadModels( MODEL_LIST, NEURON_MODELS )
    my_nest.MyLoadModels( MODEL_LIST, SYNAPSE_MODELS)       
 
    GPE_list=[] # GPE input for each experiment
    for i_exp in EXPERIMENTS:
        GPE = MyPoissonInput( n=N_GPE, sd=True, 
                              spath=SPATH, sname_nb=SNAME_NB+i_exp)
        GPE_list.append(GPE)
    
    SNR_list=[] # SNR groups for each synapse
    for i_syn, syn in enumerate(SYNAPSE_MODELS):
        SNR = MyGroup( NEURON_MODELS[0], n=N_EXP, params={'I_e':I_E}, 
                       sd=True, mm=False,
                       mm_dt=.1, record_from=[''], spath=SPATH, 
                       sname_nb=SNAME_NB+i_syn)
        SNR_list.append(SNR)
    
    
        
    if not load:
        for i_exp in EXPERIMENTS:    
            GPE = GPE_list[i_exp]
            
            # Set spike times
            # Base rate
            for id in GPE[:]:                 
                GPE.set_spike_times(id=id, rates=[RATE_BASE], times=[1], 
                                    t_stop=sim_time)               
      
            # Selection        
            for id in GPE[N_GPE-N_SEL:N_GPE+1]: 
                rates = [RATE_BASE, RATE_SELE, RATE_BASE]
                times = [1, SEL_ONSET, SEL_TIME + SEL_ONSET]
                t_stop = sim_time
                GPE.set_spike_times(id=id, rates=rates, times=times, 
                                    t_stop=t_stop)     
        
            # Inhibition        
            for id in GPE[N_GPE-N_SEL-N_INH:N_GPE+1-N_SEL]: 
                rates = [RATE_BASE, RATE_INHI, RATE_BASE]
                times = [1, SEL_ONSET, SEL_TIME + SEL_ONSET]
                t_stop = sim_time
                GPE.set_spike_times(id=id, rates=rates, times=times, 
                                    t_stop=t_stop) 
        
            # Connect         
            for i_syn, syn in enumerate(SYNAPSE_MODELS):       
                    target=SNR_list[i_syn][i_exp]
                    my_nest.ConvergentConnect(GPE[:], [target], model=syn)
                      
        my_nest.MySimulate( sim_time )

    
        for GPE in GPE_list: 
            GPE.get_signal( 's' )   
        for SNR in SNR_list: 
            SNR.get_signal( 's' ) 

        misc.pickle_save([GPE_list,SNR_list] , SAVE_AT)

    if load:
        GPE_list, SNR_list=misc.pickle_load(SAVE_AT)
        
    pre_ref=str(SNR_list[0].signals['spikes'].mean_rate(SEL_ONSET-500,SEL_ONSET)) 
    pre_dyn=str(SNR_list[2].signals['spikes'].mean_rate(SEL_ONSET-500,SEL_ONSET))   
      
    s='\n'
    s=s+'Example:\n'
    s = s + ' %s %5s %3s \n' % ( 'N experiments:', str ( len(EXPERIMENTS) ),  '#' )  
    s = s + ' %s %5s %3s \n' % ( 'Base rate:',   str ( RATE_BASE),'Hz' )     
    s = s + ' %s %5s %3s \n' % ( 'Selection rate:', str ( RATE_SELE ), 'Hz' )
    s = s + ' %s %5s %3s \n' % ( 'Selection time:', str ( SEL_TIME ), 'ms' )
    s = s + ' %s %5s %3s \n' % ( 'Pre sel rate Ref:', pre_ref[0:4], 'Hz' )
    s = s + ' %s %5s %3s \n' % ( 'Pre sel rate Dyn:', pre_dyn[0:4], 'Hz' )
      
    return GPE_list, SNR_list, s

stop = time.time()    
sec = stop - start
m = sec // 60
sec = sec - m*60
print 'Rank %i simulation time: %i minutes, %i seconds' % ( my_nest.Rank(), m, sec )

if not mpiRun:
    
    info_string=''
    GPE_list, SNR_list, s=simulate_example(hz=200, load=False)
    info_string=info_string+s  
    
        
    # DISPLAY
    plot_settings.set_mode(mode='by_fontsize', w = 1100.0, h = 450.0, fontsize=12)
    font_size_text = 8
    fig = pylab.figure( facecolor = 'w' )
    
    ax_list = []
    ax_list.append( MyAxes(fig, [ .075, .37, .135, .26 ] ) )    # text box
    ax_list.append( MyAxes(fig, [ .26,  .6,  .165, .34 ] ) )    # 
    ax_list.append( MyAxes(fig, [ .53,  .6,  .165, .34 ] ) )    # 
    #ax_list.append( MyAxes(fig, [ .8,   .6,  .165, .34 ] ) )    # 
    #ax_list.append( MyAxes(fig, [ .26,  .1,  .165, .34 ] ) )    # 
    #ax_list.append( MyAxes(fig, [ .53,  .1,  .165, .34 ] ) )    # 
    #ax_list.append( MyAxes(fig, [ .8,   .1,  .165, .34 ] ) )    # 
    
    # Text
    ax=ax_list[0]
    plot_text(ax, info_string)
    
    ax=ax_list[1]
    plot_example_raster_GPE(ax, GPE_list)
    
    ax=ax_list[2]
    plot_example_SNR(ax, SNR_list)


 
    pylab.show()
    
    fig.savefig( picture_dir + '/' + SNAME  + '.svg', dpi = 500, format = 'svg')