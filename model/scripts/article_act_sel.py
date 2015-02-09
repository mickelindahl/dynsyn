#! Imports
import copy
import math
import numpy
import pylab
import os
import sys
import time

start = time.time() 

numpy.random.seed(1) # set random seed
 
# Get directory where model and code resides 
model_dir=   '/'.join(os.getcwd().split('/')[0:-1])    
code_dir=  '/'.join(os.getcwd().split('/')[0:-2])     

# Add model, code and current directories to python path
sys.path.append(os.getcwd())  
sys.path.append(model_dir)
sys.path.append(code_dir+'/nest_toolbox') 
          
# Imports dependent on adding code model and nest_toolbox path
from model_params import models                                  
from src import misc, my_nest, my_signals, plot_settings
from src.my_population import MyGroup, MyPoissonInput 
from src.my_axes import MyAxes 
from simulation_utils import _simulate_model, simulate_model

### GLOBALS ###

# Paths and naming for saving data and picutes
FILE_NAME = sys.argv[0].split('/')[-1].split('.')[0]
model_name=os.getcwd().split('/')[-2]
PICTURE_PATH='/'.join(os.getcwd().split('/')[0:-3]) + '/pictures/'+model_name  
OUTPUT_PATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]

# Models
NEURON_MODELS=['SNR_aeif']
SYNAPSE_MODELS_TESTED=['MSN_SNR_gaba_s_min', 
                       'MSN_SNR_gaba_s_max', 
                       'MSN_SNR_gaba_p1']
SYNAPSE_MODELS_BACKGROUND=['GPE_SNR_gaba_p', 'STN_SNR_ampa_s']

# Neuron numbers
N_MSN=500
N_MSN_BURST=20
N_MSN_BASE=N_MSN-N_MSN_BURST
N_GPE = 50
N_STN = 50

IDS_MSN=numpy.arange(1,N_MSN+1)
IDS_MSN_BURST=numpy.arange(450,450+N_MSN_BURST)
IDS_MSN_NO_BURST=[id for id in IDS_MSN if id not in IDS_MSN_BURST]
IDS_GPE=numpy.arange(0,N_GPE)

# Neuron rate defaults
MSN_BASE_RATE=0.1
MSN_BURST_RATE=20.0
GPE_BASE_RATE=30
STN_BASE_RATE=10

# Misc
SELECTION_THR=5.  # spikes/s
SNR_INJECTED_CURRENT=280.0 # 280 pA
SEL_ONSET=2000
ADJUST_XDATA_MS=1500.
SEL_INTERVALS=[[0,100],
                [400,500],
                [0,500],
                [0,200],
                [300,500]]   

LIM_SYN_EVENTS=800

HEADER_SIMULATION_SETUP=( '**** BEGINNING GENERAL SCRIPT SETUP ****\n'+
                          'FILE_NAME:'+str(FILE_NAME)+'\n'+                         
                          'PICTURE_PATH:'+str(PICTURE_PATH)+'\n'+  
                          'OUTPUT_PATH:'+str(OUTPUT_PATH)+'\n\n'+                        
                          'NEURON_MODEL:'+str(NEURON_MODELS)+'\n'+
                          'SYNAPSE_MODELS_TESTED:'+str(SYNAPSE_MODELS_TESTED)+'\n'+
                          'SYNAPSE_MODELS_BACKGROUND:'+str(SYNAPSE_MODELS_BACKGROUND)+'\n\n'+
                          'N_MSN:'+str(N_MSN)+'\n'+
                          'N_MSN_BURST:'+str(N_MSN_BURST)+'\n'+
                          'N_GPE:'+str(N_GPE)+'\n\n'+
                          'MSN_BASE_RATE:'+str(MSN_BASE_RATE)+'\n'+
                          'MSN_BURST_RATE:'+str(MSN_BURST_RATE)+'\n'+
                          'GPE_BASE_RATE:'+str(GPE_BASE_RATE)+'\n\n'+
                           'STN_BASE_RATE:'+str(STN_BASE_RATE)+'\n\n'+

                          'SELECTION_THR:'+str(SELECTION_THR)+'\n'+
                          'SNR_INJECTED_CURRENT:'+str(SNR_INJECTED_CURRENT)+'\n'+
                          'SEL_ONSET:'+str(SEL_ONSET)+'\n' 
                          'ADJUST_XDATA_MS:'+str(ADJUST_XDATA_MS)+'\n\n'+
                          '**** END GENERAL SCRIPT SETUP ****\n')

### END GLOBALS ###


def plot_example_SNR(ax, SNR_list):
    time_bin=20
    
    colors=['b','g','m']   
    labels=[r'$ref_{init}^{MSN_{D1}}$' , r'$ref_{max}^{MSN_{D1}}$',  
            r'$fac^{MSN_{D1}}$']
    coords=[[0.4, 0.52], [ 0.05, 0.22], [0.1, 0.74]]
    
    for color, SNR in zip(colors, SNR_list):
        signal=SNR.signals['spikes']
        signal.my_firing_rate(bin=time_bin, display=ax,
                          kwargs={'color':color})
        
    lines = ax.lines
    for line in lines:
        line.set_xdata(line.get_xdata()-ADJUST_XDATA_MS)    
        misc.slice_line(line, xlim=[0,1490])
   
    ax.plot([0, 1490],[SELECTION_THR]*len([0, 1490]),**{'color':[0.6,0.6,0.6],
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
    ax.text( 0.9, 0.15,'Thr' , transform=ax.transAxes, **{ 'color' : [0.6,0.6,0.6]})   
        
    #ax.set_title('bin=%i'%(time_bin),**{'fontsize':12})
    ax.set_ylabel('Rate SNr (spikes/s)') 
    ax.set_xlabel('Time (ms)')
    ax.my_set_no_ticks( yticks=6, xticks = 8 ) 

    ax.set_xlim(misc.adjust_limit([0,1500]))
    ax.set_ylim(misc.adjust_limit([0,40]))

    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['text.fontsize'], 
                 **{'color': color})
 
def plot_example_raster_MSN(ax, MSN_spikes_and_ids):
    global ADJUST_XDATA_MS
    global N_MSN_BURST
    global N_MSN


    colors=['r','b']
    
    for c, data in zip(colors, MSN_spikes_and_ids):
        ax.plot(data[:,0], data[:,1], ',', color=c)
    
    lines1 = ax.lines
    #ax_twinx.set_ylim([0,500])
    ax.set_ylabel('MSN id')
    ax.set_xlabel('Time (ms)')

    for line in lines1:
        line.set_xdata(line.get_xdata()-ADJUST_XDATA_MS)
    
    ax.set_xlim(misc.adjust_limit([0,1500]))
    ax.set_ylim(misc.adjust_limit([0,N_MSN]))


    ax.my_set_no_ticks( yticks=6, xticks = 5 )   
    
    ax.text( 0.05, 0.05, 'Non-bursting MSNs' , backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : 'b' })  
    ax.text( 0.05, 0.15, 'Bursting MSNs' , backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : 'r' }) 
    
def plot_example_firing_frequency_MSN(ax, MSN_firing_rates):
    global ADJUST_XDATA_MS
    global N_GPE
    global N_MSN_BURST
    time_bin=20
    
    colors=['b','r','k']
    
    for c, data in zip(colors, MSN_firing_rates):
        ax.plot(data[0],data[1], color=c)
    
    ax.set_xlim(misc.adjust_limit([0,1500]))
    ax.set_ylim(misc.adjust_limit([0,35]))
    
    lines = ax.lines
    
    for line in lines:
        line.set_xdata(line.get_xdata()-ADJUST_XDATA_MS)
    
    misc.slice_line(lines[0], xlim=[0,1490])
    misc.slice_line(lines[1], xlim=[0,1490])
    misc.slice_line(lines[2], xlim=[0,1490])
    
    #ax.set_title('bin=%i'%(time_bin),**{'fontsize':12})
    ax.set_ylabel('Rate MSNs (spikes/s)')
    ax.set_xlabel('Time (ms)')     
    ax.my_set_no_ticks( yticks=6, xticks = 5 )
    
    ax.text( 0.05, 0.75, 'Non-bursting MSNs' , backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : 'b' })  
    ax.text( 0.05, 0.85, 'Bursting MSNs' , backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : 'r' }) 
    ax.text( 0.05, 0.65, 'Average over all' , backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : 'k' })  
    
def plot_selection_vs_neurons_full(ax, hzs, data):

    colors=['b','g','m']   
    labels=[r'$ref_{init}^{MSN_{D1}}$' , r'$ref_{max}^{MSN_{D1}}$',  
            r'$fac_{MSN_{D1}}$']  
    coords=[[0.03, 0.8], [ 0.05, 0.07], [0.03, 0.50]]   
    
    syn=SYNAPSE_MODELS_TESTED
    
    for id, label, color in zip([0,1,2],labels,colors):
        ax.plot(hzs,numpy.array(data[id]['msn_at_thr'][0])/float(N_MSN)*100.,**{'color':color}) 
    
    for id, label, color in zip([0,1,2],labels,colors):
        ax.plot(hzs,numpy.array(data[id]['msn_at_thr'][1])/float(N_MSN)*100.,**{'color':color,'linestyle':'--'})
    
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,48])    
        
    ax.set_xlabel(r'Rate $MSN_{D1}$ (spikes/s)') 
    ax.set_ylabel(r'Bursting $MSN_{D1}$ at thr (%)')
    ax.my_set_no_ticks( yticks=8, xticks = 6 ) 

    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['text.fontsize'], 
                 **{'color': color})
    
    line1=ax.plot(1,1,'-k')
    line2=ax.plot(1,1,'--k')
    leg=ax.legend([line1, line2],['First 100 ms of \nthe burst', 'Last 100 ms of \nthe burst'], loc='best')
    frame  = leg.get_frame() 
    #frame.set_visible(False) 
    pylab.setp(frame, edgecolor='w') 
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    pylab.setp(ltext, fontsize=font_size_text, backgroundcolor='w') 
    
  
    ax.set_xlim(misc.adjust_limit([6,48]))
    ax.set_ylim(misc.adjust_limit([0,30]))   
        
def plot_text(ax, info_string=''):
    
    my_nest.ResetKernel(threads=1)
    model_list, model_dict=models()
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    
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

def plot_SNr_rate_vs_syn_event1(ax, syn_events, meanRates, meanRatesStd):
    
    colors=['k']   
    labels=[ r'$fac^{MSN}$']
    coords=[[0.04, 0.17]]   
    
    for id, label, color in zip([2],labels,colors):
        ax.plot(syn_events, meanRates[id,:],**{'color':color}) 
        #ax.fill_between(syn_events, meanRates[id,:]+meanRatesStd[id,:], 
        #                meanRates[id,:]-meanRatesStd[id,:], 
         #               facecolor=color, alpha=0.1) 
 
    vec=[50,LIM_SYN_EVENTS]
    ax.plot(vec,[SELECTION_THR]*len(vec),**{'color':[0.6,0.6,0.6],
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
    
    
    
    ax.text( 0.8, 0.18,'Thr' , transform=ax.transAxes, **{ 'color' : [0.6,0.6,0.6]})
    ax.set_ylabel('Rate SNr (spikes/s)') 
    ax.set_xlabel('Input SNr (events/s)')
    ax.my_set_no_ticks( yticks=8, xticks = 5 ) 
    
    x_arrow=400./LIM_SYN_EVENTS
    y_arrow=2./50.
    
    ax.arrow(x_arrow, y_arrow, 0, 0.02, transform=ax.transAxes,
             width=0.01, head_width=0.03, head_starts_at_zero=True,
            head_length=0.02,**{ 'color':'k'})
    
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['text.fontsize'], 
                 **{'color': color})
    
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,LIM_SYN_EVENTS]) 
    
    ax.set_xlim(misc.adjust_limit([50,LIM_SYN_EVENTS]))
    ax.set_ylim(misc.adjust_limit([0,34]))

def plot_SNr_rate_vs_syn_event2(ax, nbNeurons, meanRates):
    colors=['k'] 
    ax.plot(nbNeurons,meanRates[2,:],**{'color':colors[0], 'linestyle':'-.'})  
 
    ax.set_ylabel('Rate SNr (spikes/s)') 
    ax.set_xlabel('Synaptic events (#/s)')
    ax.my_set_no_ticks( yticks=6, xticks = 5 ) 
    
    line1=ax.plot(1,1,'-k')
    line2=ax.plot(1,1,'-.k')
    leg=ax.legend([line1, line2],['Increased size of \nbursting subpopulation', 'Increased activity  in \nall MSNs'], 
                  loc='best')
    frame  = leg.get_frame() 
    frame.set_visible(False) 
    ltext  = leg.get_texts()  # all the text.Text instance in the legend

    
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,LIM_SYN_EVENTS])   
    
    pylab.setp(ltext, fontsize=font_size_text) 
    ax.set_xlim(misc.adjust_limit([50,LIM_SYN_EVENTS]))
    ax.set_ylim(misc.adjust_limit([0,34]))
  
  
def plot_SNr_rate_vs_syn_event_x3_early(ax, hzs):  
    
    colors=misc.make_N_colors('YlOrBr', 5)[2:]
    
    labels=['%i Hz'%hzs[0] , '%i Hz'%hzs[1], '%i Hz'%hzs[2]]
    coords=[[0.4, 0.52], [ 0.05, 0.22], [0.1, 0.74]]
    
    i_inverval=0 # 0-100 ms
    i_syn=2 # Plastic
    
    for hz, col in zip(hzs, colors):
        
        nb, rate_data, r_std_data, n_max_sel, s = simulate_selection_vs_neurons(SEL_INTERVALS,
                                                                       msn_burst_rate=hz,
                                                                       load=True)  
        syn_events=nb*hz+(N_MSN-nb)*0.1
        ax.plot(syn_events,rate_data[i_inverval][i_syn,:],**{'color':col}) 
    
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=font_size_text, 
                 **{'color': color})
    
    ax.set_ylabel('Rate SNr (spikes/s)') 
    ax.set_xlabel('Input SNr (events/s)')

    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,LIM_SYN_EVENTS])   
    
    ax.my_set_no_ticks( yticks=6, xticks = 5 ) 
    ax.set_xlim(misc.adjust_limit([50,LIM_SYN_EVENTS]))
    ax.set_ylim(misc.adjust_limit([0,34]))
    ax.text( 0.1, 0.85, 'First 100 ms of the burst', transform=ax.transAxes, 
                 fontsize=12, 
                 **{'color': 'k'})

def plot_SNr_rate_vs_syn_event_x3_late(ax, hzs):  
    
    colors=misc.make_N_colors('YlOrBr', 5)[2:]
    
    labels=['%i Hz'%hzs[0] , '%i Hz'%hzs[1], '%i Hz'%hzs[2]]
    coords=[[0.4, 0.52], [ 0.05, 0.22], [0.1, 0.74]]
    
    i_inverval=1 # 0-100 ms
    i_syn=2 # Plastic
    
    for hz, col in zip(hzs, colors):
        
        nb, rate_data, r_std_data, n_max_sel, s = simulate_selection_vs_neurons(SEL_INTERVALS,
                                                                       msn_burst_rate=hz,
                                                                       load=True)  
        syn_events=nb*hz+(N_MSN-nb)*0.1
        ax.plot(syn_events,rate_data[i_inverval][i_syn,:],**{'color':col}) 
    
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=font_size_text, 
                 **{'color': color})
    
    ax.set_ylabel('Rate SNr (spikes/s)') 
    ax.set_xlabel('Input SNr (events/s)')

    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,LIM_SYN_EVENTS])   
    
    ax.my_set_no_ticks( yticks=6, xticks = 5 ) 
    ax.set_xlim(misc.adjust_limit([50,LIM_SYN_EVENTS]))
    ax.set_ylim(misc.adjust_limit([0,34]))
    ax.text( 0.1, 0.85, 'Last 100 ms of the burst', transform=ax.transAxes, 
                 fontsize=font_size_text, 
                 **{'color': 'k'})
    
    
def plot_SNr_rate_std_vs_syn_event_x3_early(ax, hzs):  
    
    colors=misc.make_N_colors('YlOrBr', 5)[2:]
    
    labels=['%i Hz'%hzs[0] , '%i Hz'%hzs[1], '%i Hz'%hzs[2]]
    coords=[[0.4, 0.52], [ 0.05, 0.22], [0.1, 0.74]]
    
    i_inverval=0 # 0-100 ms
    i_syn=2 # Plastic
    
    for hz, col in zip(hzs, colors):
        
        nb, rate_data, r_std_data, n_max_sel, s = simulate_selection_vs_neurons(SEL_INTERVALS,
                                                                       msn_burst_rate=hz,
                                                                       load=True)  
        syn_events=nb*hz+(N_MSN-nb)*0.1
        ax.plot(syn_events,rate_data[i_inverval][i_syn,:]/r_std_data[i_inverval][i_syn,:],**{'color':col}) 
    
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=font_size_text, 
                 **{'color': color})
    
    ax.set_ylabel('Signal to noise ratio') 
    ax.set_xlabel('Input SNr (events/s)')

    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,LIM_SYN_EVENTS])   
    
    ax.my_set_no_ticks( yticks=6, xticks = 5 ) 
    ax.set_xlim(misc.adjust_limit([50,LIM_SYN_EVENTS]))
    ax.set_ylim(misc.adjust_limit([0,4]))
    ax.text( 0.1, 0.85, 'First 100 ms of the burst', transform=ax.transAxes, 
                 fontsize=font_size_text, 
                 **{'color': 'k'})
    
def plot_SNr_rate_std_vs_syn_event_x3_late(ax, hzs):  
    
    colors=misc.make_N_colors('YlOrBr', 5)[2:]
    
    labels=['%i Hz'%hzs[0] , '%i Hz'%hzs[1], '%i Hz'%hzs[2]]
    coords=[[0.4, 0.52], [ 0.05, 0.22], [0.1, 0.74]]
    
    i_inverval=1 # 0-100 ms
    i_syn=2 # Plastic
    
    for hz, col in zip(hzs, colors):
        
        nb, rate_data, r_std_data, n_max_sel, s = simulate_selection_vs_neurons(SEL_INTERVALS,
                                                                       msn_burst_rate=hz,
                                                                       load=True)  
        syn_events=nb*hz+(N_MSN-nb)*0.1
        ax.plot(syn_events,rate_data[i_inverval][i_syn,:]/r_std_data[i_inverval][i_syn,:],**{'color':col}) 
    
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=font_size_text, 
                 **{'color': color})
    
    ax.set_ylabel('Signal to noise ratio') 
    ax.set_xlabel('Input SNr (events/s)')

    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,LIM_SYN_EVENTS])   
    
    ax.my_set_no_ticks( yticks=6, xticks = 5 ) 
    ax.set_xlim(misc.adjust_limit([50,LIM_SYN_EVENTS]))
    ax.set_ylim(misc.adjust_limit([0,4]))
    ax.text( 0.1, 0.85, 'Last 100 ms of the burst', transform=ax.transAxes, 
                 fontsize=font_size_text, 
                 **{'color': 'k'})

def plot_thr_rate_vs_std(ax, data, hzs):
      
      colors=['b','g','m']  
      for syn, col in zip(data.keys(), colors):
          r=numpy.array(data[syn]['rates_thr'][0])
          r_std=numpy.array(data[syn]['rates_std_thr'][0])
          
          if syn==1:
              ax.plot(hzs, r, color=col)
              ax.fill_between(hzs, r-r_std, r+r_std, color=col,alpha=0.5)

      
def plot_example_3x_freq(ax, times, spk_mean, hzs=[10,20,40] ):
    spk_mean=spk_mean
    colors=misc.make_N_colors('YlOrBr', 5)[2:]
    
    labels=['%i Hz'%hzs[0] , '%i Hz'%hzs[1], '%i Hz'%hzs[2]]
    coords=[[0.4, 0.52], [ 0.05, 0.22], [0.1, 0.74]]
    for i, col in enumerate(colors):
        ax.plot(times, misc.convolve(spk_mean[i,:], 50, 'triangle',single=True)[0], color=col)
        #ax.plot(times, spk_mean[i,:], col)

    lines = ax.lines
    for line in lines:
        line.set_xdata(line.get_xdata()-ADJUST_XDATA_MS)    
        misc.slice_line(line, xlim=[0,1490])
    #ax.set_ylim([0,40])
    
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=font_size_text, 
                 **{'color': color})
    
    ax.set_ylabel('Rate SNr (spikes/s)') 
    ax.set_xlabel('Time (ms)')
    
    ax.my_set_no_ticks( yticks=6, xticks = 8 ) 
    ax.set_xlim(misc.adjust_limit([0,1500]))
    ax.set_ylim(misc.adjust_limit([0,40]))


def simulate_example(load=True):

    global GPE_BASE_RATE  
    global FILE_NAME
    global N_STN
    global N_GPE
    global N_MSN_BURST
    global N_MSN
    global NEURON_MODELS
    global OUTPUT_PATH
    global SEL_ONSET
    global SNR_INJECTED_CURRENT
    global SYNAPSE_MODELS_TESTED
    
    n_exp =200 # number of experiments   
    #n_exp =200 # number of experiments   
    
    # Path were raw data is saved. For example the spike trains.
    save_result_at=OUTPUT_PATH+'/simulate_example.pkl'
    save_header_at=OUTPUT_PATH+'/simulate_example_header'
    
    burst_time = 500.
    sim_time = burst_time+SEL_ONSET+1000.
      
    model_list, model_dict=models()
    my_nest.ResetKernel(threads=8)       
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    my_nest.MyLoadModels( model_list, SYNAPSE_MODELS_TESTED)       
    my_nest.MyLoadModels( model_list, SYNAPSE_MODELS_BACKGROUND)      
 
    SNR_list=[] # List with SNR groups for synapse. 
    if not load:    
        MSN_base=MyPoissonInput(n=N_MSN_BASE*n_exp, sd=True)
        MSN_burst=MyPoissonInput(n=N_MSN_BURST*n_exp, sd=True)     
        GPE=MyPoissonInput(n=N_GPE*n_exp, sd=True)
        STN=MyPoissonInput(n=N_STN*n_exp, sd=True)
                            
        # Set spike times MSN and GPe
        # Non bursting MSNs
        for id in MSN_base[:]:  
            seed=numpy.random.random_integers(0,1000000.0)        
            MSN_base.set_spike_times(id=id, rates=[MSN_BASE_RATE], times=[1], 
                        t_stop=sim_time, seed=seed)           
        
        # Background GPe    
        for id in GPE[:]:
            seed=numpy.random.random_integers(0,1000000.0)                  
            GPE.set_spike_times(id=id, rates=[GPE_BASE_RATE], times=[1], 
                                t_stop=sim_time, seed=seed)  
         # Background STN   
        for id in STN[:]:
            seed=numpy.random.random_integers(0,1000000.0)                  
            STN.set_spike_times(id=id, rates=[STN_BASE_RATE], times=[1], 
                                t_stop=sim_time, seed=seed)  
                                
        # Bursting MSNs        
        for id in MSN_burst[:]: 
            rates = [MSN_BASE_RATE, MSN_BURST_RATE, MSN_BASE_RATE]
            times = [1, SEL_ONSET, burst_time + SEL_ONSET]
            t_stop = sim_time
            seed=numpy.random.random_integers(0,1000000.0) 

            MSN_burst.set_spike_times(id=id, rates=rates, times=times, 
                                    t_stop=t_stop, seed=seed)       
        
        for i_syn in range(len(SYNAPSE_MODELS_TESTED)):
            
            params=[]
            I_e=my_nest.GetDefaults(NEURON_MODELS[0])['I_e']+SNR_INJECTED_CURRENT
            for i in range(n_exp):
                #params.append({'I_e':numpy.random.normal(I_e, 
                #                                         0.1*I_e)})
                params.append({'I_e':I_e})

            #{'I_e':SNR_INJECTED_CURRENT}
            SNR=MyGroup( NEURON_MODELS[0], n=n_exp, sd=True, 
               params=params, mm_dt=.1, 
               record_from=[''])
                    
            SNR_list.append(SNR)
    
        
        # Connect, experiment specific    
        sources_MSN_SNR_base=numpy.arange(0,n_exp*N_MSN_BASE)
        sources_MSN_SNR_burst=numpy.arange(0,n_exp*N_MSN_BURST)
        
        targets_MSN_SNR_base=numpy.mgrid[0:n_exp, 0:N_MSN_BASE][0].reshape(1,N_MSN_BASE*n_exp)[0]
        targets_MSN_SNR_burst=numpy.mgrid[0:n_exp, 0:N_MSN_BURST][0].reshape(1,N_MSN_BURST*n_exp)[0]
         
        sources_GPE_SNR=numpy.arange(0,n_exp*N_GPE)
        targets_GPE_SNR=numpy.mgrid[0:n_exp, 0:N_GPE][0].reshape(1,N_GPE*n_exp)[0]

        sources_STN_SNR=numpy.arange(0,n_exp*N_STN)
        targets_STN_SNR=numpy.mgrid[0:n_exp, 0:N_STN][0].reshape(1,N_STN*n_exp)[0]

        
        for i_syn, syn in enumerate(SYNAPSE_MODELS_TESTED):      
            syn= SYNAPSE_MODELS_TESTED[i_syn]
            SNR=SNR_list[i_syn]
            my_nest.Connect(MSN_base[sources_MSN_SNR_base], SNR[targets_MSN_SNR_base], model=syn)
            my_nest.Connect(MSN_burst[sources_MSN_SNR_burst], SNR[targets_MSN_SNR_burst], model=syn)
            my_nest.Connect(GPE[sources_GPE_SNR], SNR[targets_GPE_SNR], model=SYNAPSE_MODELS_BACKGROUND[0])
            my_nest.Connect(STN[sources_STN_SNR], SNR[targets_STN_SNR], model=SYNAPSE_MODELS_BACKGROUND[1])
                 
        my_nest.MySimulate( sim_time )
        
        MSN_base.get_signal( 's', start=0, stop=sim_time )   
        MSN_burst.get_signal( 's', start=0, stop=sim_time )   

        for SNR in SNR_list: 
            SNR.get_signal( 's', start=0, stop=sim_time ) 
        
        
        # Get firing rates of MSNs
        MSN_firing_rates=[]
        
        MSN_all=copy.deepcopy(MSN_base)
        MSN_all.merge(MSN_burst)
    
        time_bin=20.
        groups=[MSN_base, MSN_burst, MSN_all]
        for group in groups:
            timeAxis, firingRates=group.signals['spikes'].my_firing_rate(bin=time_bin, display=False)
            MSN_firing_rates.append([timeAxis, firingRates])
        
        # Pick out spikes for burst, base and all to use in scatter plot    
        MSN_spikes_and_ids=[]
        
        g1=MSN_burst.slice(MSN_burst[0:N_MSN_BURST])
        g2=MSN_base.slice(MSN_base[0:N_MSN_BASE])
            
        ids_MSN_burst=range(450,450+N_MSN_BURST)
        ids_MSN_base=[id for id in range(N_MSN) if id not in IDS_MSN_BURST]
        
        # Rename ids for plotting purpose        
        g1_dict=dict([[id1,id2] for id1 ,id2 in zip(g1.ids, ids_MSN_burst)  ])
        g2_dict=dict([[id1,id2] for id1 ,id2 in zip(g2.ids, ids_MSN_base)  ])
        
        groups=[g1,g2]
        dics=[g1_dict,g2_dict]
        for group, dic in zip(groups, dics):
            raw_data=group.signals['spikes'].raw_data()
            for i in range(raw_data.shape[0]):
                raw_data[i,1]=dic[raw_data[i,1]]
            MSN_spikes_and_ids.append(raw_data)
         
            
        pre_ref_1=str(SNR_list[0].signals['spikes'].mean_rate(SEL_ONSET-500,
                                                              SEL_ONSET)) 
        pre_ref_2=str(SNR_list[1].signals['spikes'].mean_rate(SEL_ONSET-500,
                                                              SEL_ONSET)) 
        pre_dyn=str(SNR_list[2].signals['spikes'].mean_rate(SEL_ONSET-500,
                                                            SEL_ONSET))   
          
        s='\n'
        s=s+'Simulate example:\n'
        s = s + '%s %5s %3s \n' % ( 'Simulation time', str ( sim_time ),  '#' ) 
        s = s + '%s %5s %3s \n' % ( 'N experiments:', str ( n_exp ),  '#' )  
        s = s + '%s %5s %3s \n' % ( 'Bin size MSN hz:',     str ( time_bin),'ms' )    
        s = s + '%s %5s %3s \n' % ( 'MSN base rate:',   str ( MSN_BASE_RATE),'spikes/s' )     
        s = s + '%s %5s %3s \n' % ( 'MSN burst rate:', str ( MSN_BURST_RATE ), 'spikes/s' )
        s = s + '%s %5s %3s \n' % ( 'GPe rate:', str ( GPE_BASE_RATE ), 'spikes/s' )
        s = s + '%s %5s %3s \n' % ( 'Burst time:', str ( burst_time ), 'ms' )
        s = s + '%s %5s %3s \n' % ( 'Pre sel rate Ref:', pre_ref_1[0:4], 'spikes/s' )
        s = s + '%s %5s %3s \n' % ( 'Pre sel rate Ref:', pre_ref_2[0:4], 'spikes/s' )
        s = s + '%s %5s %3s \n' % ( 'Pre sel rate Dyn:', pre_dyn[0:4], 'spikes/s' )
        
        header=HEADER_SIMULATION_SETUP+s
        misc.text_save(header, save_header_at)
        misc.pickle_save([MSN_firing_rates, 
                          MSN_spikes_and_ids, 
                          SNR_list, s], save_result_at)

           
    else:
        MSN_firing_rates, MSN_spikes_and_ids, SNR_list, s = misc.pickle_load(save_result_at)
      
    return MSN_firing_rates, MSN_spikes_and_ids, SNR_list, s
        
def simulate_selection_vs_neurons(selection_intervals=[0.0,500.0], msn_burst_rate=20, load=True, n_max_sel=100):    
    global SNR_INJECTED_CURRENT
    global NEURON_MODELS
    global N_GPE
    global N_STN
    global N_MSN_BURST
    global N_MSN
    global MSN_BASE_RATE
    global GPE_BASE_RATE  
    global FILE_NAME
    global OUTPUT_PATH
    global SYNAPSE_MODELS_TESTED
    global SEL_ONSET
          
    #n_exp=100
    n_exp=200
 
    save_result_at = (OUTPUT_PATH+'/simulate_selection_vs_neurons_'+str(msn_burst_rate)+'hz.pkl') 
    save_header_at = (OUTPUT_PATH+'/simulate_selection_vs_neurons_'+str(msn_burst_rate)+'hz_header') 

    burst_time = 500.
    sim_time = burst_time+SEL_ONSET+500.   
    
    EXPERIMENTS=range(n_exp)

    model_list, model_dict=models()
    my_nest.ResetKernel(threads=8)     
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    my_nest.MyLoadModels( model_list, SYNAPSE_MODELS_TESTED)   
    my_nest.MyLoadModels( model_list, SYNAPSE_MODELS_BACKGROUND)   
    
    if not load:
        MSN_list=[] # MSN input for each experiment
        for i_exp in EXPERIMENTS:
            MSN = MyPoissonInput( n=N_MSN+n_max_sel, sd=True)
            MSN_list.append(MSN)
        
        GPE_list=[] # GPE input for each experiment
        for i_exp in EXPERIMENTS:
            GPE = MyPoissonInput( n=N_GPE, sd=True)
            GPE_list.append(GPE)
        
        STN_list=[] # GPE input for each experiment
        for i_exp in EXPERIMENTS:
            STN = MyPoissonInput( n=N_STN, sd=True)
            STN_list.append(STN)
        
        
        SNR_list=[] # SNR groups for each synapse and number of selected MSN
        SNR_list_experiments=[]
        for i_syn, syn in enumerate(SYNAPSE_MODELS_TESTED):
            SNR = []
            for i_sel in range(n_max_sel+1): # Plus one to get no burst point
                
                I_e=my_nest.GetDefaults(NEURON_MODELS[0])['I_e']+SNR_INJECTED_CURRENT           
                SNR.append(MyGroup( NEURON_MODELS[0], n=n_exp, sd=True, params={'I_e':I_e}))
                
            SNR_list.append(SNR)    
           
        for i_exp in EXPERIMENTS:    
            MSN = MSN_list[i_exp]
            GPE = GPE_list[i_exp]
            STN = STN_list[i_exp]
            
            # Set spike times
            # Base rate
            for id in MSN[1:N_MSN]:                 
                MSN.set_spike_times(rates=[MSN_BASE_RATE], times=[1], 
                                    t_stop=sim_time,ids=[id], 
                                    seed=int(numpy.random.random()*10000.0))               
      
            # Selection        
            for id in MSN[N_MSN:N_MSN + n_max_sel]: 
                rates = [MSN_BASE_RATE, msn_burst_rate, MSN_BASE_RATE]
                times = [1, SEL_ONSET, burst_time + SEL_ONSET]
                t_stop = sim_time
                MSN.set_spike_times(rates=rates, times=times, 
                                    t_stop=t_stop, ids=[id], 
                                    seed=int(numpy.random.random()*10000.0))    
            
            # Base rate GPE
            for id in GPE[:]:                 
                GPE.set_spike_times(rates=[GPE_BASE_RATE], times=[1], 
                                    t_stop=sim_time, ids=[id], 
                                    seed=int(numpy.random.random()*10000.0))    

            # Base rate STN
            for id in STN[:]:                 
                STN.set_spike_times(rates=[STN_BASE_RATE], times=[1], 
                                    t_stop=sim_time, ids=[id], 
                                    seed=int(numpy.random.random()*10000.0))   
            
            # Connect, only for plastic     
            for i_syn, syn in enumerate(SYNAPSE_MODELS_TESTED):
                # i_sel goes over 0,..., n_max_sel
                for i_sel, n_sel in enumerate(range(0,n_max_sel+1)):       
                    target = SNR_list[i_syn][i_sel][i_exp]
                    
                    my_nest.ConvergentConnect(MSN[0:N_MSN-n_sel], [target], 
                                                  model=syn)
                    my_nest.ConvergentConnect(MSN[N_MSN:N_MSN+n_sel], 
                                                  [target], model=syn)
                    my_nest.ConvergentConnect(GPE[:], [target], 
                                          model=SYNAPSE_MODELS_BACKGROUND[0])
                    my_nest.ConvergentConnect(STN[:], [target], 
                                          model=SYNAPSE_MODELS_BACKGROUND[1])
                     
        my_nest.MySimulate( sim_time ) 
              
        for SNR_sel in SNR_list: 
            for SNR in SNR_sel: 
                SNR.get_signal( 's' ) 
        
           
        sel_interval_mean_rates=[]
        sel_interval_mean_rates_std=[]
        for i_interval, interval in enumerate(selection_intervals):        
            t1=selection_intervals[i_interval][0]
            t2=selection_intervals[i_interval][1]
            
            mean_rates=[]    
            mean_rates_std=[]
                
            # Time until arrival of spikes in SNr
            delay=my_nest.GetDefaults(SYNAPSE_MODELS_BACKGROUND[0])['delay']
            for SNR_sel in SNR_list: 
                m_r=[]
                m_r_std=[]
                for SNR in SNR_sel:
                                     
                     m_r.append(SNR.signals['spikes'].mean_rate(SEL_ONSET+t1+delay, SEL_ONSET+t2+delay))  
                     m_r_std.append(SNR.signals['spikes'].mean_rate_std(SEL_ONSET+t1+delay, SEL_ONSET+t2+delay)) 
                    
                mean_rates.append(m_r)
                mean_rates_std.append(m_r_std)
                
            mean_rates = numpy.array(mean_rates)
            mean_rates_std = numpy.array(mean_rates_std)
            
            sel_interval_mean_rates.append(mean_rates)
            sel_interval_mean_rates_std.append(mean_rates_std)

        nb_neurons = numpy.arange(0,n_max_sel+1,1)
            
        s='\n'
        s=s + 'simulate_selection_vs_neurons\n'
        s = s + '%s %5s %3s \n' % ( 'Simulation time', str ( sim_time ),  '#' ) 
        s = s + '%s %5s %3s \n' % ( 'N MSNs:', str ( N_MSN ),  '#' )     
        s = s + '%s %5s %3s \n' % ( 'N experiments:', str ( n_exp ),  '#' )    
        s = s + '%s %5s %3s \n' % ( 'MSN base rate:',   str ( MSN_BASE_RATE),'spikes/s' )     
        s = s + '%s %5s %3s \n' % ( 'MSN burst rate:', str ( MSN_BURST_RATE ), 'spikes/s' )
        s = s + '%s %5s %3s \n' % ( 'MSN burst time:', str ( burst_time ), 'ms' )
        s = s + '%s %5s %3s \n' % ( 'GPe base rate:', str ( GPE_BASE_RATE ), 'spikes/s' )
        s = s + '%s %5s %3s \n' % ( 'SNR injected current:', str ( SNR_INJECTED_CURRENT), 'pA' )
        for i_interval, interval in enumerate(selection_intervals):
            s = s + '%s %5s %3s \n' % ( 'Sel interval '+str(i_interval)+':', str ( selection_intervals), 'ms' )

        info_string=s
        
        header=HEADER_SIMULATION_SETUP+s
        misc.text_save(header, save_header_at)
        misc.pickle_save([nb_neurons, 
                          sel_interval_mean_rates, 
                          sel_interval_mean_rates_std,n_max_sel,
                           info_string], save_result_at)
        
    elif load:     
        nb_neurons, sel_interval_mean_rates, sel_interval_mean_rates_std, n_max_sel, info_string = misc.pickle_load(save_result_at)
    
    return nb_neurons, sel_interval_mean_rates, sel_interval_mean_rates_std, n_max_sel, info_string

def simulate_selection_vs_neurons_full(selRateInterval, load_pickle=True, load_raw=True):
    global OUTPUT_PATH

    save_result_at=OUTPUT_PATH+'/simulate_selection_vs_neurons_full.pkl'
    save_header_at=OUTPUT_PATH+'/simulate_selection_vs_neurons_full_header'
   

    # Range
    hzs=numpy.arange(7,49,1)
    
    #hzs=[8,20]
    if not load_pickle:
        data={}
        
        for syn in range(3):
            data[syn]={}
            data[syn]['rates_thr']=[ [] for k in range(len(SEL_INTERVALS))]
            data[syn]['rates_std_thr']=[ [] for k in range(len(SEL_INTERVALS))]
 
            data[syn]['msn_at_thr']=[ [] for k in range(len(SEL_INTERVALS))]
            data[syn]['n_max_sel']=[ [] for k in range(len(SEL_INTERVALS))]
        
        n_max_sel=218  
        progress=''         
        i_hz=0
        for hz in hzs:
                                  
            n, rate_data, r_std_data, n_max_sel, s = simulate_selection_vs_neurons(SEL_INTERVALS, hz, load_raw, n_max_sel=n_max_sel)
            n_sel_vec=numpy.arange(n_max_sel+1)
            
            # Clear n_max_cell
            n_max_sel=0
            for i_interval in range(len(rate_data)):
                for i_syn in range(len(rate_data[i_interval])):
                     
                    r_syn=rate_data[i_interval][i_syn]
                    r_std_syn=r_std_data[i_interval][i_syn]
                    
                    # Retrieve selection threshold passing
                    
                    r_std_syn_tmp=r_std_syn[r_syn<SELECTION_THR]
                    n_sel_vec_tmp=n_sel_vec[r_syn<SELECTION_THR]
                    r_syn_tmp=r_syn[r_syn<SELECTION_THR]
                    
                    data[i_syn]['rates_thr'][i_interval].append(r_syn_tmp[0])
                    data[i_syn]['rates_std_thr'][i_interval].append(r_std_syn_tmp[0])             
                    data[i_syn]['msn_at_thr'][i_interval].append(n_sel_vec_tmp[0])
                
                    # Find new n_max_sel
                    msn_at_thr=data[i_syn]['msn_at_thr'][i_interval][i_hz]
                    n_max_sel=int(numpy.ceil(max(msn_at_thr*2.0, n_max_sel)))
                    data[i_syn]['n_max_sel'][i_interval].append(n_max_sel)
            
            i_hz+=1 
            progress+=str(hz)+' hz finished, n_max_sel='+str(n_max_sel)+'\n' 
            print progress
        
        s='\n'
        s=s + 'simulate_selection_vs_neurons_full\n'
        s = s + ' %s %5s %s \n' % ( 'Range hz', str ( hzs[0])+'-'+str(hzs[-1]),  '#' )     


        header=HEADER_SIMULATION_SETUP+s
        misc.text_save(header, save_header_at)
        misc.pickle_save([data,s], save_result_at)
        info_string=s
    elif load_pickle:
        data, info_string=misc.pickle_load(save_result_at)
      
    
    return hzs, data, info_string

def simulate_example_3x_freq(hzs=[10,  20, 30], max_syn_event=400, res=1, load=True):

    n_msns=[max_syn_event/hzs[0], max_syn_event/hzs[1], max_syn_event/hzs[2]]



    save_result_at=OUTPUT_PATH+'/simulate_diff_3x_freq.pkl'
    n_exp=20
    if not load:  
        spk_mean=[]
        for hz, n_msn in zip(hzs,n_msns):
            sim_time= SEL_ONSET+1000
            params_msn={'base_rates':[MSN_BASE_RATE], 'base_times':[1], 'mod_rates': [MSN_BASE_RATE, hz, MSN_BASE_RATE],
                        'mod_times':[1,SEL_ONSET, SEL_ONSET+500], 'n_tot':N_MSN, 'n_mod':n_msn}
            params_gpe={'base_rates':[GPE_BASE_RATE], 'base_times':[1], 'n_tot':N_GPE, 'n_mod':0}
            params_stn={'base_rates':[STN_BASE_RATE], 'base_times':[1], 'n_tot':N_STN, 'n_mod':0}
            synapse_models={'MSN':'MSN_SNR_gaba_p1', 'GPE':'GPE_SNR_gaba_p',
                            'STN':'STN_SNR_ampa_s'}
            
            #times, spk_binned =_simulate_model([params_msn, params_gpe,params_gpe, 'SNR_izh', 
            #                                   synapse_models, sim_time, 0])
            
            t=time.time()
        
            times, spk_binned =simulate_model(params_msn, params_gpe, params_stn, 'SNR_aeif', 
                                              synapse_models, sim_time, res, n_exp=n_exp,threads=4)
            
            print 'Time:',time.time()-t
            spk_mean.append(numpy.mean(spk_binned,axis=0))
        spk_mean=numpy.array(spk_mean)*1000/res       
        misc.pickle_save([times, spk_mean], save_result_at)        
    else:        
        times, spk_mean = misc.pickle_load(save_result_at)   
    return times, spk_mean    

print 'Simulation'

# SIMULATION
info_string=''

sel_interval_4 = [0.0, 200.0]


# If simulating with MPI
# 1. Run with MPI load_pickle and load2 False
# 2. When finished run laod2=True and load_pickle=False to recover MPI saved data
# data and save it as pickled file.
# 3. Run with load_pickle and load2 both True and plot data

load_pickle = False # Load selection vs neuron full data
load_raw = False # Load selection vs neuron data used in full function, 
              # should be False when mpi run
hzs, data, s=simulate_selection_vs_neurons_full(sel_interval_4, load_pickle, load_raw)
info_string=info_string+s  

stop = time.time()    
sec = stop - start
m = sec // 60
sec = sec - m*60
print 'Rank %i simulation time: %i minutes, %i seconds' % ( my_nest.Rank(), m, sec )

   
MSN_firing_rates, MSN_spikes_and_ids, SNR_list, s=simulate_example( load=True )
info_string=info_string+s  

 

# DISPLAY
plot_settings.set_mode(pylab, mode='by_fontsize', 
                       w = 1100.0, h = 450.0, fontsize=16)
font_size_text = 12
fig = pylab.figure( facecolor = 'w' )

ax_list = []
ax_list.append( MyAxes(fig, [ .075, .37, .135, .26 ] ) )    # text box
ax_list.append( MyAxes(fig, [ .26,  .6,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .53,  .6,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .26,  .1,  .165*2.312, .34 ] ) )    #     
ax_list.append( MyAxes(fig, [ .8,   .6,  .165, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .53,  .1,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .1,  .165, .34 ] ) )    # 

# Text
ax=ax_list[0]
plot_text(ax, info_string)

ax=ax_list[1]
plot_example_raster_MSN(ax, MSN_spikes_and_ids)

ax=ax_list[2]
plot_example_firing_frequency_MSN(ax, MSN_firing_rates)
#plot_thr_rate_vs_std(ax, data, hzs)

ax=ax_list[3]
plot_example_SNR(ax, SNR_list)


ax=ax_list[4]  
plot_selection_vs_neurons_full(ax, hzs, data) 

# plot_selection_vs_neurons 

ax=ax_list[5]
# article_filtering.py has to be run before
data=misc.pickle_load(os.getcwd()+'/output/mean_rates_filtering'+
NEURON_MODELS[0])
MSNmeanRates=data['MSN_mean_rates']
SNRmeanRates=data['SNR_mean_rates']
print MSNmeanRates
print SNRmeanRates
syn_events=MSNmeanRates*N_MSN
SNRmeanRates[2,syn_events<LIM_SYN_EVENTS]
plot_SNr_rate_vs_syn_event2(ax, syn_events, SNRmeanRates)

hz=20
nb, rate_data, r_std_data, n_max_sel, s = simulate_selection_vs_neurons(SEL_INTERVALS,
                                                                       msn_burst_rate=hz,
                                                                       load=True)
i=2
syn_events=nb*hz+(N_MSN-nb)*0.1
plot_SNr_rate_vs_syn_event1(ax, syn_events, rate_data[i], r_std_data[i])




times, spk_mean=simulate_example_3x_freq(load=True, hzs=[10,20,40], max_syn_event=400)


# DISPLAY
plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0, fontsize=16)
font_size_text = 16
fig2 = pylab.figure( facecolor = 'w' )

ax_list = []
ax_list.append( MyAxes(fig2, [ .075, .37, .135, .26 ] ) )    # text box
ax_list.append( MyAxes(fig2, [ .26,  .6,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig2, [ .53,  .6,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig2, [ .8,   .6,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig2, [ .26,  .1,  .165*2.312, .34 ] ) )    #   
#ax_list.append( MyAxes(fig, [ .53,  .1,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig2, [ .8,   .1,  .165, .34 ] ) )    # 

ax=ax_list[1]
plot_SNr_rate_vs_syn_event_x3_early(ax, hzs=[10,20,40])

   
ax=ax_list[2]
plot_SNr_rate_vs_syn_event_x3_late(ax, hzs=[10,20,40])

ax=ax_list[3]
plot_SNr_rate_std_vs_syn_event_x3_early(ax, hzs=[10,20,40])

ax=ax_list[4]
plot_example_3x_freq(ax, times, spk_mean, hzs=[10,20,40] )

ax=ax_list[5]
plot_SNr_rate_std_vs_syn_event_x3_late(ax, hzs=[10,20,40])

pylab.show()

fig.savefig( PICTURE_PATH + '/' + FILE_NAME  + '.svg', format = 'svg')
fig.savefig( PICTURE_PATH + '/' + FILE_NAME  + '.pdf', format = 'pdf')
fig2.savefig( PICTURE_PATH + '/' + FILE_NAME  + '2.svg', format = 'svg')
fig2.savefig( PICTURE_PATH + '/' + FILE_NAME  + '2.pdf', format = 'pdf')

