#! Imports
import copy
import math
import numpy
import pylab
import os
import sys
import time

if len(sys.argv) != 1: mpiRun = True
else:                  mpiRun = False
start = time.time() 

numpy.random.seed(1) # set random seed
 
# Get directory where model and code resides 
model_dir=   '/'.join(os.getcwd().split('/')[0:-1])    
code_dir=  '/'.join(os.getcwd().split('/')[0:-2])     

# Add model, code and current directories to python path
sys.path.append(os.getcwd())  
sys.path.append(model_dir)
sys.path.append(code_dir+'/nest_toolbox') 

# Picture directory
model_name=os.getcwd().split('/')[-2]
picture_dir='/'.join(os.getcwd().split('/')[0:-3]) + '/pictures/'+model_name    
            
# Imports dependent on adding code model and nest_toolbox path
from model_params import models                                  
from src import misc, my_nest, my_signals, plot_settings
from src.my_population import MyGroup, MyPoissonInput 
from src.my_axes import MyAxes 

### GLOBALS ###
OUTPUT_PATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]
FILE_NAME  = sys.argv[0].split('/')[-1].split('.')[0]
SELECTION_THR=5.  # Hz
NEURON_MODELS=['SNR_izh']

# Neuron numbers
N_MSN=500
N_MSN_BURST=11
N_MSN_BASE=N_MSN-N_MSN_BURST
N_GPE = 22

IDS_MSN=numpy.arange(1,N_MSN+1)
IDS_MSN_BURST=numpy.arange(450,450+N_MSN_BURST)
IDS_MSN_NO_BURST=[id for id in IDS_MSN if id not in IDS_MSN_BURST]
IDS_GPE=numpy.arange(0,N_GPE)

# Neuron rate defaults
MSN_BASE_RATE=0.1
MSN_BURST_RATE=20.0
GPE_BASE_RATE=25

SNR_INJECTED_CURRENT=530.0
SEL_ONSET=2000
ADJUST_XDATA_MS=1500.

SYNAPSE_MODELS_TESTED=['MSN_SNR_gaba_s_min', 
                       'MSN_SNR_gaba_s_max', 
                       'MSN_SNR_gaba_p1']

SYNAPSE_MODELS_BACKGROUND=['GPE_SNR_gaba_p']

HEADER_SIMULATION_SETUP=( '**** BEGINNING GENERAL SCRIPT SETUP ****\n'+
                          'N_MSN:'+str(N_MSN)+'\n'+
                          'N_MSN_BURST:'+str(N_MSN_BURST)+'\n'+
                          'N_GPE:'+str(N_GPE)+'\n\n'+
                          'MSN_BASE_RATE:'+str(MSN_BASE_RATE)+'\n'+
                          'MSN_BURST_RATE:'+str(MSN_BURST_RATE)+'\n'+
                          'GPE_BASE_RATE:'+str(GPE_BASE_RATE)+'\n\n'+
                          'SNR_INJECTED_CURRENT:'+str(SNR_INJECTED_CURRENT)+'\n'+
                          'SEL_ONSET:'+str(SEL_ONSET)+'\n' 
                          'ADJUST_XDATA_MS:'+str(ADJUST_XDATA_MS)+'\n\n'+
                          'SYNAPSE_MODELS_TESTED:'+str(SYNAPSE_MODELS_TESTED)+'\n'+
                          'SYNAPSE_MODELS_BACKGROUND:'+str(SYNAPSE_MODELS_BACKGROUND)+'\n'+
                          '**** END GENERAL SCRIPT SETUP ****\n')


def plot_example_SNR(ax, SNR_list):
    time_bin=20
    
    colors=['b','g','m']   
    labels=[r'$\delta_{weak}^{MSN}$' , r'$\delta_{strong}^{MSN}$',  
            r'$\delta_{fac}^{MSN}$']
    coords=[[0.4, 0.49], [ 0.05, 0.27], [0.1, 0.75
                                         ]]
    
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
    ax.set_ylabel('Frequency SNr (Hz)') 
    ax.set_xlabel('Time (ms)')
    ax.my_set_no_ticks( yticks=6, xticks = 7 ) 

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
    ax.set_ylim(misc.adjust_limit([0,40]))
    
    lines = ax.lines
    
    for line in lines:
        line.set_xdata(line.get_xdata()-ADJUST_XDATA_MS)
    
    misc.slice_line(lines[0], xlim=[0,1490])
    misc.slice_line(lines[1], xlim=[0,1490])
    misc.slice_line(lines[2], xlim=[0,1490])
    
    #ax.set_title('bin=%i'%(time_bin),**{'fontsize':12})
    ax.set_ylabel('Frequency MSNs (Hz)')     
    ax.my_set_no_ticks( yticks=6, xticks = 5 )
    
    ax.text( 0.05, 0.75, 'Non-bursting MSNs' , backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : 'b' })  
    ax.text( 0.05, 0.85, 'Bursting MSNs' , backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : 'r' }) 
    ax.text( 0.05, 0.65, 'Average over all' , backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : 'k' })  
    
def plot_selection_vs_neurons_full(ax, hzs, data):

    colors=['b','g','m']   
    labels=[r'$\delta_{weak}^{MSN}$' , r'$\delta_{strong}^{MSN}$',  
            r'$\delta_{fac}^{MSN}$']  
    coords=[[0.03, 0.8], [ 0.05, 0.07], [0.03, 0.50]]   
    
    syn=SYNAPSE_MODELS_TESTED
    
    for id, label, color in zip([0,1,2],labels,colors):
        ax.plot(hzs,data[syn[id]]['thrVec'],**{'color':color}) 
    
    
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,48])    
        
    ax.set_xlabel('Frequency MSN (Hz)') 
    ax.set_ylabel('MSNs at thr (#)')
    ax.my_set_no_ticks( yticks=6, xticks = 6 ) 

    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['text.fontsize'], 
                 **{'color': color})
    
    ax.set_xlim(misc.adjust_limit([0,48]))
    ax.set_ylim(misc.adjust_limit([0,75]))    
    
def plot_selection_vs_neurons_full_extra(ax, hzs,data):   
    colors=['b','g','m']   
    labels=[r'$\delta_{weak}^{MSN}$' , r'$\delta_{strong}^{MSN}$',  
            r'$\delta_{fac}^{MSN}$']
    
    syn=SYNAPSE_MODELS_TESTED
    for id, label, color in zip([0,1,2],labels,colors):
        ax.plot(hzs,data[syn[id]]['thrVec'],**{'color':color,'linestyle':'--'})

    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,48])    


    line1=ax.plot(1,1,'-k')
    line2=ax.plot(1,1,'--k')
    leg=ax.legend([line1, line2],['First 200 ms of \nthe burst', '300-500 ms of \nthe burst'], loc='best')
    frame  = leg.get_frame() 
    #frame.set_visible(False) 
    pylab.setp(frame, edgecolor='w') 
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    pylab.setp(ltext, fontsize=10., backgroundcolor='w') 
    
  
    ax.set_xlim(misc.adjust_limit([-7,48]))
    ax.set_ylim(misc.adjust_limit([0,100]))  
    
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

def plot_SNr_rate_vs_syn_event1(ax, syn_events, meanRates, meanRatesStd):
    
    colors=['k']   
    labels=[ r'$\delta_{fac}^{MSN}$']
    coords=[[0.04, 0.16]]   
    
    for id, label, color in zip([2],labels,colors):
        ax.plot(syn_events,meanRates[id,:],**{'color':color}) 
#        ax.fill_between(syn_events, meanRates[id,:]+meanRatesStd[id,:], 
#                        meanRates[id,:]-meanRatesStd[id,:], 
#                        facecolor=color, alpha=0.5) 
 
    vec=[50,600]
    ax.plot(vec,[SELECTION_THR]*len(vec),**{'color':[0.6,0.6,0.6],
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
    
    
    
    ax.text( 0.8, 0.17,'Thr' , transform=ax.transAxes, **{ 'color' : [0.6,0.6,0.6]})
    ax.set_ylabel('Frequency SNr (Hz)') 
    ax.set_xlabel('Synaptic events (#/s)')
    ax.my_set_no_ticks( yticks=6, xticks = 5 ) 
    
    x_arrow=285./600.
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
        misc.slice_line(line, xlim=[0,600]) 
    
    ax.set_xlim(misc.adjust_limit([50,600]))
    ax.set_ylim(misc.adjust_limit([0,35]))

def plot_SNr_rate_vs_syn_event2(ax, nbNeurons, meanRates):
    colors=['k'] 
    ax.plot(nbNeurons,meanRates[2,:],**{'color':colors[0], 'linestyle':'-.'})  
 
    ax.set_ylabel('Frequency SNr (Hz)') 
    ax.set_xlabel('Synaptic events (#/s)')
    ax.my_set_no_ticks( yticks=5, xticks = 5 ) 
    
    line1=ax.plot(1,1,'-k')
    line2=ax.plot(1,1,'-.k')
    leg=ax.legend([line1, line2],['Bursting subpopulation', 'Increased activity  in \nall MSNs'], 
                  loc='best')
    frame  = leg.get_frame() 
    frame.set_visible(False) 
    ltext  = leg.get_texts()  # all the text.Text instance in the legend

    
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,600])   
    
    pylab.setp(ltext, fontsize=10.) 
    ax.set_xlim(misc.adjust_limit([50,600]))
    ax.set_ylim(misc.adjust_limit([0,35]))
    
def simulate_example(load=True):

    global SNR_INJECTED_CURRENT
    global NEURON_MODELS
    global N_GPE
    global N_MSN_BURST
    global N_MSN
    global GPE_BASE_RATE  
    global FILE_NAME
    global OUTPUT_PATH
    global SYNAPSE_MODELS_TESTED
    global SEL_ONSET
    
    #n_exp =200 # number of experiments   
    n_exp =20 # number of experiments   
    
    # Path were raw data is saved. For example the spike trains.
    save_result_at=OUTPUT_PATH+'/' + FILE_NAME + '-simulate_example.pkl'
    save_header_at=OUTPUT_PATH+'/' + FILE_NAME + '-simulate_example_header'


    
    burst_time = 500.
    sim_time = burst_time+SEL_ONSET+1000.
      
    MODEL_LIST=models()
    my_nest.ResetKernel()       
    my_nest.MyLoadModels( MODEL_LIST, NEURON_MODELS )
    my_nest.MyLoadModels( MODEL_LIST, SYNAPSE_MODELS_TESTED)       
    my_nest.MyLoadModels( MODEL_LIST, SYNAPSE_MODELS_BACKGROUND)      
 
    SNR_list=[] # List with SNR groups for synapse. 
    if not load:    
        MSN_base=MyPoissonInput(n=N_MSN_BASE*n_exp, sd=True)
        MSN_burst=MyPoissonInput(n=N_MSN_BURST*n_exp, sd=True)     
        GPE=MyPoissonInput(n=N_GPE*n_exp, sd=True)

                        
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
        
        for i_syn, syn in enumerate(SYNAPSE_MODELS_TESTED):      
            syn= SYNAPSE_MODELS_TESTED[i_syn]
            SNR=SNR_list[i_syn]
            my_nest.Connect(MSN_base[sources_MSN_SNR_base], SNR[targets_MSN_SNR_base], model=syn)
            my_nest.Connect(MSN_burst[sources_MSN_SNR_burst], SNR[targets_MSN_SNR_burst], model=syn)
            my_nest.Connect(GPE[sources_GPE_SNR], SNR[targets_GPE_SNR], 
                                          model=SYNAPSE_MODELS_BACKGROUND[0])
                 
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
         
        #times, binned_data=MSN_base.signals['spikes'].binned_raw_data(0, sim_time, res=1, clip=0)
        #filtered_binned_data=misc.time_resolved_rate(binned_data, 100, kernel_type='triangle', res=1)
            
        pre_ref_1=str(SNR_list[0].signals['spikes'].mean_rate(SEL_ONSET-500,
                                                              SEL_ONSET)) 
        pre_ref_2=str(SNR_list[1].signals['spikes'].mean_rate(SEL_ONSET-500,
                                                              SEL_ONSET)) 
        pre_dyn=str(SNR_list[2].signals['spikes'].mean_rate(SEL_ONSET-500,
                                                            SEL_ONSET))   
          
        s='\n'
        s=s+'Simulate example:\n'
        s = s + '%s %5s %3s \n' % ( 'N experiments:', str ( n_exp ),  '#' )  
        s = s + '%s %5s %3s \n' % ( 'Bin size MSN hz:',     str ( time_bin),'ms' )    
        s = s + '%s %5s %3s \n' % ( 'MSN base rate:',   str ( MSN_BASE_RATE),'Hz' )     
        s = s + '%s %5s %3s \n' % ( 'MSN burst rate:', str ( MSN_BURST_RATE ), 'Hz' )
        s = s + '%s %5s %3s \n' % ( 'GPe rate:', str ( GPE_BASE_RATE ), 'Hz' )
        s = s + '%s %5s %3s \n' % ( 'Burst time:', str ( burst_time ), 'ms' )
        s = s + '%s %5s %3s \n' % ( 'Pre sel rate Ref:', pre_ref_1[0:4], 'Hz' )
        s = s + '%s %5s %3s \n' % ( 'Pre sel rate Ref:', pre_ref_2[0:4], 'Hz' )
        s = s + '%s %5s %3s \n' % ( 'Pre sel rate Dyn:', pre_dyn[0:4], 'Hz' )
        
        header=HEADER_SIMULATION_SETUP+s
        
        misc.pickle_save([MSN_firing_rates, 
                          MSN_spikes_and_ids, 
                          SNR_list, s], save_result_at)
        misc.text_save(header, save_header_at)
           
    else:
        MSN_firing_rates, MSN_spikes_and_ids, SNR_list, s = misc.pickle_load(save_result_at)
      
    return MSN_firing_rates, MSN_spikes_and_ids, SNR_list, s
        
def simulate_selection_vs_neurons(selection_intervals=[0.0,500.0], hz=20, load=True):    
    global SNR_INJECTED_CURRENT
    global NEURON_MODELS
    global N_GPE
    global N_MSN_BURST
    global N_MSN
    global GPE_BASE_RATE  
    global FILE_NAME
    global OUTPUT_PATH
    global SYNAPSE_MODELS_TESTED
    global SEL_ONSET
          
    #n_exp=100
    n_exp=2
  
    if hz>7:
        n_max_sel=60
    if hz>20:
        n_max_sel=30
    else:
        n_max_sel=100
 
    RATE_BASE = 0.1
    RATE_SELE = hz
    save_result_at = (OUTPUT_PATH+'/'+FILE_NAME+'-simulate_selection_vs_neurons'+str(hz)+'-hz.pkl') 
    save_header_at = (OUTPUT_PATH+'/'+FILE_NAME+'-simulate_selection_vs_neurons'+str(hz)+'-hz_header') 

    burst_time = 500.
    sim_time = burst_time+SEL_ONSET+500.   
    
    EXPERIMENTS=range(n_exp)

    MODEL_LIST=models()
    my_nest.ResetKernel()     
    my_nest.MyLoadModels( MODEL_LIST, NEURON_MODELS )
    my_nest.MyLoadModels( MODEL_LIST, SYNAPSE_MODELS_TESTED)   
    my_nest.MyLoadModels( MODEL_LIST, SYNAPSE_MODELS_BACKGROUND)   
    
    MSN_list=[] # MSN input for each experiment
    for i_exp in EXPERIMENTS:
        MSN = MyPoissonInput( n=N_MSN+n_max_sel, sd=True)
        MSN_list.append(MSN)
    
    GPE_list=[] # GPE input for each experiment
    for i_exp in EXPERIMENTS:
        GPE = MyPoissonInput( n=N_GPE, sd=True)
        GPE_list.append(GPE)
    
    
    SNR_list=[] # SNR groups for each synapse and number of selected MSN
    SNR_list_experiments=[]
    for i_syn, syn in enumerate(SYNAPSE_MODELS_TESTED):
        SNR = []
        for i_sel in range(n_max_sel+1): # Plus one to get no burst point
            
            I_e=my_nest.GetDefaults(NEURON_MODELS[0])['I_e']+SNR_INJECTED_CURRENT           
            SNR.append(MyGroup( NEURON_MODELS[0], n=n_exp, sd=True, params={'I_e':I_e}))
            
        SNR_list.append(SNR)    
        
    if not load:
        for i_exp in EXPERIMENTS:    
            MSN = MSN_list[i_exp]
            GPE = GPE_list[i_exp]
        
            # Set spike times
            # Base rate
            for id in MSN[1:N_MSN]:                 
                MSN.set_spike_times(id=id, rates=[RATE_BASE], times=[1], 
                                    t_stop=sim_time,
                                    seed=int(numpy.random.random()*10000.0))               
      
            # Selection        
            for id in MSN[N_MSN:N_MSN + n_max_sel]: 
                rates = [RATE_BASE, RATE_SELE, RATE_BASE]
                times = [1, SEL_ONSET, burst_time + SEL_ONSET]
                t_stop = sim_time
                MSN.set_spike_times(id=id, rates=rates, times=times, 
                                    t_stop=t_stop,
                                    seed=int(numpy.random.random()*10000.0))    
            
            # Base rate GPE
            for id in GPE[:]:                 
                GPE.set_spike_times(id=id, rates=[GPE_BASE_RATE], times=[1], 
                                    t_stop=sim_time, 
                                    seed=int(numpy.random.random()*10000.0))    
            
            # Connect     
            for i_syn, syn in enumerate(SYNAPSE_MODELS_TESTED):
                # i_sel goes over 0,..., n_max_sel
                for i_sel in range(0,n_max_sel+1):       
                    target = SNR_list[i_syn][i_sel][i_exp]
                    
                    my_nest.ConvergentConnect(MSN[0:N_MSN-i_sel], [target], 
                                                  model=syn)
                    my_nest.ConvergentConnect(MSN[N_MSN:N_MSN+i_sel], 
                                                  [target], model=syn)
                    my_nest.ConvergentConnect(GPE[:], [target], 
                                          model=SYNAPSE_MODELS_BACKGROUND[0])
                     
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
        s = s + ' %s %5s %3s \n' % ( 'N MSNs:', str ( N_MSN ),  '#' )     
        s = s + ' %s %5s %3s \n' % ( 'N experiments:', str ( n_exp ),  '#' )    
        s = s + ' %s %5s %3s \n' % ( 'MSN base rate:',   str ( MSN_BASE_RATE),'Hz' )     
        s = s + ' %s %5s %3s \n' % ( 'MSN burst rate:', str ( MSN_BURST_RATE ), 'Hz' )
        s = s + ' %s %5s %3s \n' % ( 'GPe rate:', str ( GPE_BASE_RATE ), 'Hz' )
        s = s + ' %s %5s %3s \n' % ( 'Burst time:', str ( burst_time ), 'ms' )
        s = s + ' %s %5s %3s \n' % ( 'SNR_INJECTED_CURRENT:', str ( SNR_INJECTED_CURRENT), 'pA' )
        for i_interval, interval in enumerate(selection_intervals):
            s = s + ' %s %5s %3s \n' % ( 'Sel interval '+str(i_interval)+':', str ( selection_intervals), 'ms' )

        info_string=s
        
        header=HEADER_SIMULATION_SETUP+s
        misc.text_save(header, save_header_at)
        misc.pickle_save([nb_neurons, 
                          sel_interval_mean_rates, 
                          sel_interval_mean_rates_std, info_string], save_result_at)
        
    elif load:     
        nb_neurons, sel_interval_mean_rates, sel_interval_mean_rates_std, info_string = misc.pickle_load(save_result_at)
    
    return nb_neurons, sel_interval_mean_rates, sel_interval_mean_rates_std, info_string

def simulate_selection_vs_neurons_full(selRateInterval, load_pickle=True, load_raw=True):
    global N_MSN

    save_result_at=OUTPUT_PATH+'/'+FILE_NAME+'-simulate_selection_vs_neurons_full.pkl'
    save_result_at=OUTPUT_PATH+'/'+FILE_NAME+'-simulate_selection_vs_neurons_full_header'
   
   
   
    selection_intervals=[[0,200],
                     [300,500]]   
    #Range 1
    #hzs=range(5,8)

    # Range 2 
    #hzs=range(8,61,1) # MPI can cope with jump 7->8 when n max selected decreases
    
    # Range
    hzs=range(5,49,1)
    
    #hzs=[8,20]
    if not load_pickle:
        data={}
        for syn in SYNAPSE_MODELS_TESTED:
            data[syn]={}
            data[syn]['rates']=[ [] for k in range(len(selection_intervals))]
            data[syn]['selMat']=[ [] for k in range(len(selection_intervals))]
            data[syn]['thrVec']=[ [] for k in range(len(selection_intervals))]
            
            
        for hz in hzs:
            n, r, r_std, s = simulate_selection_vs_neurons(selection_intervals, hz, load_raw)

            print hz, 'hz finished'
            for i_syn, syn in enumerate(SYNAPSE_MODELS_TESTED):
                for k in range(len(selection_intervals)):  data[syn]['rates'][k].append(r[k][i_syn,:])

                
                # Create matricies   
                # Adjust rates if they are different length, add zeros at end of     
                # short vectors. OBS only true if last rate is zero in vectorn that
                # is elonged.
                  

                for k in range(len(selection_intervals)): 
                    rates=data[syn]['rates'][k]
                        
                    maxLen=0
                    for r in rates:
                            if len(r)>maxLen:
                                maxLen = len(r)
                    for i_r, r in  enumerate(rates):
                                rates[i_r]=numpy.append(r,numpy.zeros( (1,maxLen-len(r)) ))                  
                 
                    selMat=rates
                    thrVec=[]
                    
                    for i in range(rates.shape[0]):
                            p=True
                            for j in range(rates.shape[1]):
                                if SELECTION_THR<r[i,j]:
                                    selMat[i,j]=3
                                elif (SELECTION_THR>=rates[i,j]) and (SELECTION_THR<rates[i,j-1]) and p:
                                    selMat[i,j]=2
                                    thrVec.append(j+1)  # Neurons for threshold
                                    p=False
                                else:
                                    selMat[i,j]=1
                            if p: 
                                thrVec[j].append(100)
                
                        
                    data[syn]['selMat'][k]=selMat
                    data[syn]['thrVec'][k]=numpy.array(thrVec)
        
        if not mpiRun:
            header=HEADER_SIMULATION_SETUP+s
            misc.text_save(header, save_header_at)
            misc.pickle_save(data, save_result_at)
    
    elif load_pickle:
        data=misc.pickle_load(save_result_at)
      
    s='\n'

    info_string=s
    
    return hzs, data, info_string
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
hzs_1, data_1, s=simulate_selection_vs_neurons_full(sel_interval_4, load_pickle, load_raw)
info_string=info_string+s  

stop = time.time()    
sec = stop - start
m = sec // 60
sec = sec - m*60
print 'Rank %i simulation time: %i minutes, %i seconds' % ( my_nest.Rank(), m, sec )

if not mpiRun:
    
    sel_interval_5 = [200.0, 500.0]
    

    #hzs_2, data_2, s=simulate_selection_vs_neurons_full(sel_interval_5, load_pickle, True)
    
    MSN_firing_rates, MSN_spikes_and_ids, SNR_list, s=simulate_example(load=False )
    info_string=info_string+s  
    
    
    load_fast=True
    save_result_at=OUTPUT_PATH+'/'+NEURON_MODELS[0]+'simulate_selection_vs_neurons_full_load_fast.pkl'
    
    if not load_fast:

        sel_interval_1=[0.0, 500.0]
        nbNeurons1, meanRates1, meanRatesStd1, s = simulate_selection_vs_neurons(sel_interval_1,hz=20,
                                                                  load=True)      
        s = s + ' %s %5s %3s \n' % ( 'Sel thr:', str ( SELECTION_THR ), 'Hz' )
        s = s + ' \n%s %5s %3s \n' % ( 'Sel interval 1:', str ( sel_interval_1 ), 'pA' )
        info_string=info_string+s
        misc.pickle_save([nbNeurons1, meanRates1, meanRatesStd1,
                          info_string ], save_result_at)
    else:
        print 'hej'
        #nbNeurons1, meanRates1, meanRatesStd1, info_string = misc.pickle_load(save_result_at)
    
    
    # DISPLAY
    plot_settings.set_mode(pylab, mode='by_fontsize', 
                           w = 1100.0, h = 450.0, fontsize=12)
    font_size_text = 8
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
    
    ax=ax_list[3]
    plot_example_SNR(ax, SNR_list)

    '''
    ax=ax_list[4]  
    plot_selection_vs_neurons_full(ax, hzs_1, data_1)
    plot_selection_vs_neurons_full_extra(ax, hzs_2, data_2)  
    
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
    SNRmeanRates[2,syn_events<600]
    plot_SNr_rate_vs_syn_event2(ax, syn_events, SNRmeanRates)
       
    syn_events=nbNeurons1*20.0+(N_MSN-nbNeurons1)*0.1
    plot_SNr_rate_vs_syn_event1(ax, syn_events, meanRates1, meanRatesStd1)
    '''
    pylab.show()
    
    fig.savefig( picture_dir + '/' + FILE_NAME  + '.svg', format = 'svg')
    fig.savefig( picture_dir + '/' + FILE_NAME  + '.pdf', format = 'pdf')