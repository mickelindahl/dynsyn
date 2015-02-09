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
NEURON_MODELS=['SNR_izh']
SYNAPSE_MODELS_TESTED=['MSN_SNR_gaba_p1']
SYNAPSE_MODELS_BACKGROUND=['GPE_SNR_gaba_p', 'STN_SNR_ampa_s']

# Neuron numbers
N_MSN=500
N_MSN_BURST=20
N_MSN_BASE=N_MSN-N_MSN_BURST
N_GPE = 30
N_STN = 100

IDS_MSN=numpy.arange(1,N_MSN+1)
IDS_MSN_BURST=numpy.arange(450,450+N_MSN_BURST)
IDS_MSN_NO_BURST=[id for id in IDS_MSN if id not in IDS_MSN_BURST]
IDS_GPE=numpy.arange(0,N_GPE)

# Neuron rate defaults
MSN_BASE_RATE=0.1
MSN_BURST_RATE=20.0
GPE_BASE_RATE=25
STN_BASE_RATE=10

# Misc
SELECTION_THR=5.  # spikes/s
SNR_INJECTED_CURRENT=400.0
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


   
def plot_selection_vs_neurons_full(ax, hzs, data):

    colors=['m']   
    labels=[r'$\delta_{fac}^{MSN}$']# , r'$\delta_{strong}^{MSN}$',  
           # r'$\delta_{fac}^{MSN}$']  
    coords=[[0.03, 0.8]]#, [ 0.05, 0.07], [0.03, 0.50]]   
    
    syn=SYNAPSE_MODELS_TESTED
    
    for id, label, color in zip([0],labels,colors):
        ax.plot(hzs,data[id]['msn_at_thr'][0],**{'color':color}) 
    
    #for id, label, color in zip([0,1,2],labels,colors):
    #    ax.plot(hzs,data[id]['msn_at_thr'][1],**{'color':color,'linestyle':'--'})
    
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,48])    
        
    ax.set_xlabel('Firing rate MSN (spikes/s)') 
    ax.set_ylabel('Bursting MSNs at thr (#)')
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
    pylab.setp(ltext, fontsize=10., backgroundcolor='w') 
    
  
    ax.set_xlim(misc.adjust_limit([6,48]))
    ax.set_ylim(misc.adjust_limit([0,100]))   
        
def plot_text(ax, info_string=''):
    
    my_nest.ResetKernel(threads=1)
    model_list=models()
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
          
def simulate_selection_vs_neurons(selection_intervals=[0.0,500.0], 
                                  msn_burst_rate=20, load=True, 
                                  n_max_sel=100, rate_gpe=0, rate_stn=10):    
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
          
    n_exp=20
    #n_exp=200
 
    save_result_at = (OUTPUT_PATH+'/simulate_selection_vs_neurons_'+str(msn_burst_rate)+'hz.pkl') 
    save_header_at = (OUTPUT_PATH+'/simulate_selection_vs_neurons_'+str(msn_burst_rate)+'hz_header') 

    burst_time = 500.
    sim_time = burst_time+SEL_ONSET+500.   
    
    EXPERIMENTS=range(n_exp)

    model_list=models()
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
                MSN.set_spike_times(id=id, rates=[MSN_BASE_RATE], times=[1], 
                                    t_stop=sim_time,
                                    seed=int(numpy.random.random()*10000.0))               
      
            # Selection        
            for id in MSN[N_MSN:N_MSN + n_max_sel]: 
                rates = [MSN_BASE_RATE, msn_burst_rate, MSN_BASE_RATE]
                times = [1, SEL_ONSET, burst_time + SEL_ONSET]
                t_stop = sim_time
                MSN.set_spike_times(id=id, rates=rates, times=times, 
                                    t_stop=t_stop,
                                    seed=int(numpy.random.random()*10000.0))    
            
            # Base rate GPE
            for id in GPE[:]:                 
                GPE.set_spike_times(id=id, rates=[rate_gpe], times=[1], 
                                    t_stop=sim_time, 
                                    seed=int(numpy.random.random()*10000.0))    

            # Base rate STN
            for id in STN[:]:                 
                STN.set_spike_times(id=id, rates=[rate_stn], times=[1], 
                                    t_stop=sim_time, 
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

def simulate_selection_vs_neurons_full(selRateInterval, load_pickle=True, load_raw=True, rate_gpe=0, rate_stn=10, hz_start=7):
    global OUTPUT_PATH

    save_result_at=OUTPUT_PATH+'/simulate_selection_vs_neurons_full_'+str(1+rate_gpe)+str((1+rate_stn)*10)+'.pkl'
    save_header_at=OUTPUT_PATH+'/simulate_selection_vs_neurons_full_header'+str(1+rate_gpe)+str((1+rate_stn)*10)
   

    # Range
    hzs=numpy.arange(hz_start,49,1)
    
    #hzs=[8,20]
    if not load_pickle:
        data={}
        
        for syn in SYNAPSE_MODELS_TESTED:
            data[syn]={}
            data[syn]['rates_thr']=[ [] for k in range(len(SEL_INTERVALS))]
            data[syn]['rates_std_thr']=[ [] for k in range(len(SEL_INTERVALS))]
 
            data[syn]['msn_at_thr']=[ [] for k in range(len(SEL_INTERVALS))]
            data[syn]['n_max_sel']=[ [] for k in range(len(SEL_INTERVALS))]
        
        n_max_sel=218  
        progress=''         
        i_hz=0
        for hz in hzs:
                                  
            n, rate_data, r_std_data, n_max_sel, s = simulate_selection_vs_neurons(SEL_INTERVALS, hz, load_raw,  
                                                                                   n_max_sel, rate_gpe, rate_stn)
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
                    
                    try:
                        data[i_syn]['rates_thr'][i_interval].append(r_syn_tmp[0])
                        data[i_syn]['rates_std_thr'][i_interval].append(r_std_syn_tmp[0])             
                        data[i_syn]['msn_at_thr'][i_interval].append(n_sel_vec_tmp[0])
                    except:
                        data[i_syn]['rates_thr'][i_interval].append(0)
                        data[i_syn]['rates_std_thr'][i_interval].append(0)             
                        data[i_syn]['msn_at_thr'][i_interval].append(0)    
                    # Find new n_max_sel
                    msn_at_thr=data[i_syn]['msn_at_thr'][i_interval][i_hz]
                    n_max_sel=int(numpy.ceil(max(msn_at_thr*1.5, n_max_sel)))
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
hzs, data, s=simulate_selection_vs_neurons_full(sel_interval_4, load_pickle, 
                                                load_raw, rate_gpe=10, rate_stn=10, hz_start=20)
hzs2, data2, s=simulate_selection_vs_neurons_full(sel_interval_4, load_pickle, 
                                                  load_raw, rate_gpe=25, rate_stn=100, hz_start=20)


info_string=info_string+s  

stop = time.time()    
sec = stop - start
m = sec // 60
sec = sec - m*60
print 'Rank %i simulation time: %i minutes, %i seconds' % ( my_nest.Rank(), m, sec )

   

 

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


ax=ax_list[1]
plot_selection_vs_neurons_full(ax, hzs, data)

ax=ax_list[2]
plot_selection_vs_neurons_full(ax, hzs2, data2)
#plot_thr_rate_vs_std(ax, data, hzs)

ax=ax_list[3]



ax=ax_list[4]  


# plot_selection_vs_neurons 


pylab.show()

fig.savefig( PICTURE_PATH + '/' + FILE_NAME  + '.svg', format = 'svg')
fig.savefig( PICTURE_PATH + '/' + FILE_NAME  + '.pdf', format = 'pdf')

