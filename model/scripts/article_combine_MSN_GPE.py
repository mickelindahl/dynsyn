#! Imports
import math
import numpy
numpy.random.seed(1)
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

mode_name=os.getcwd().split('/')[-2]
picture_dir='/'.join(os.getcwd().split('/')[0:-3]) + '/pictures/'+mode_name 
                
sys.path.append(model_dir) 
sys.path.append(code_dir+'/nest_toolbox') 
SPATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]
SNAME  = sys.argv[0].split('/')[-1].split('.')[0]

from model_params import models                               # Then import models     
from src import misc, my_nest, my_signals, plot_settings
from src.my_population import MyGroup, MyPoissonInput 
from src.my_axes import MyAxes 

SELECTION_THR=5.  # spikes/s
#NEURON_MODELS=['SNR_aeif']
#SNR_INJECTED_CURRENT=[-5., -5.+80+120, -5.+95+120]

NEURON_MODELS=['SNR_izh']
c_add=0
#SNR_INJECTED_CURRENT=[120+20.+c_add, 120.+320.+c_add, 120.+20.+320.+c_add]
#SNR_INJECTED_CURRENT=[120.+20.+320.+c_add, 120.+20.+320.+c_add, 120.+20.+320.+c_add]
SNR_INJECTED_CURRENT=[530.0+c_add, 530.0+c_add, 530.0+c_add]

SYNAPSE_MODELS=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p']
SYNAPSE_MODELS_BACKGROUND=['STN_SNR_ampa_s']
SEL_ONSET = 2000.
ADJUST_XDATA_MS=1500.

N_GPE = 30
N_GPE_PAUSE=10
N_GPE_BURST=10

N_MSN = 500
N_MSN_BURST=20

N_STN=100

MSN_BASE_RATE=0.1
GPE_BASE_RATE=25
STN_BASE_RATE=10.

MSN_BURST_TIME=500
GPE_PAUSE_TIME=500
GPE_BURST_TIME=200

def plot_example_SNR(ax, SNR_list, flag):
    time_bin=20
    
    colors=misc.make_N_colors('cool',3)
    
    colors=['m','c',colors[1]]   
    #labels=[r'$\delta_{fac}^{MSN}$' , r'$\delta_{dep}^{GPe}$',  
    #        r'$\delta_{fac}^{MSN}$+$\delta_{dep}^{GPe}$']
    coords=[[0.05, 0.2], [ 0.05, 0.65], [0.05, 0.8]]


    
    SNR_list=[SNR_list[0], SNR_list[1], SNR_list[2]]
    
    
    for color, SNR in zip(colors, SNR_list):
        signal=SNR.signals['spikes']
        signal.my_firing_rate(bin=time_bin, display=ax, kwargs={'color':color})
    
    lines = ax.lines
    for line in lines:
        line.set_xdata(line.get_xdata()-ADJUST_XDATA_MS)
        misc.slice_line(line, xlim=[0,1490])
        
    if flag==0:
        leg=ax.legend(lines,['Burst in MSN subpopulation' , 'Burst in GPe subpopulation',  
            'Burst in both MSN- and GPe subpopulations'], loc=2)
    if flag==1:
        leg=ax.legend(lines,['Burst in MSN subpopulation' , 'Pause in GPe subpopulation',  
            'Burst in MSN- and pause in GPe subpopulations'], loc=2)
    
    frame  = leg.get_frame() 
    #frame.set_visible(False) 
    pylab.setp(frame, edgecolor='w') 
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    pylab.setp(ltext, fontsize=12., backgroundcolor='w') 
   
    ax.plot([0, 1490],[SELECTION_THR]*len([0, 1490]),**{'color':[0.6,0.6,0.6],
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
    ax.text( 0.9, 0.11,'Thr' , transform=ax.transAxes, **{ 'color' : [0.6,0.6,0.6]})   
       
    #ax.set_title('bin=%i'%(time_bin),**{'fontsize':12})
    ax.set_ylabel('Firing rate SNr (spikes/s)') 
    ax.set_xlabel('Time (ms)')
    ax.my_set_no_ticks( yticks=7, xticks=8 ) 

    ax.set_xlim(misc.adjust_limit([0,1500]))
    ax.set_ylim(misc.adjust_limit([0,115]))
    
    
    #for label, coord, color in zip(labels,coords,colors):
   #     ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
    #             fontsize=pylab.rcParams['text.fontsize']-4, 
    #             **{'color': color}) 

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

def simulate_example(MSN_hz=20, GPE_hz=0, load=True, n_gpe_sel=3, sel_time_GPE=500):
    global GPE_BASE_RATE
    global STN_BASE_RATE
    global MSN_BASE_RATE
    global MSN_BURST_TIME

    global NEURON_MODELS
    global N_GPE
    global N_STN
    global N_MSN
    global N_MSN_BURST

    global SNAME
    global SPATH
    global SYNAPSE_MODELS
    global SEL_ONSET
    global SNR_INJECTED_CURRENT
    
    n_exp = 200
   
    msn_rate_sel = MSN_hz # Selection rate     
    gpe_sel_rate = GPE_hz # Selection rate     

    sel_time_MSN = MSN_BURST_TIME
    sim_time = sel_time_MSN+SEL_ONSET+500.
    
    EXPERIMENTS=range(n_exp)
    
    MODEL_LIST=models()
    my_nest.ResetKernel()       
    my_nest.MyLoadModels( MODEL_LIST, NEURON_MODELS )
    my_nest.MyLoadModels( MODEL_LIST, SYNAPSE_MODELS)      
    my_nest.MyLoadModels( MODEL_LIST, SYNAPSE_MODELS_BACKGROUND)       
 
    
    MSN_list=[] # MSN input for each experiment
    for i_exp in EXPERIMENTS:
        MSN = MyPoissonInput( n=N_MSN+N_MSN_BURST, sd=True)
        MSN_list.append(MSN)
 
    GPE_list=[] # GPE input for each experiment
    for i_exp in EXPERIMENTS:
        GPE = MyPoissonInput( n=N_GPE+n_gpe_sel, sd=True)
        GPE_list.append(GPE)

    STN_list=[] # GPE input for each experiment
    for i_exp in EXPERIMENTS:
        STN = MyPoissonInput( n=N_STN, sd=True)
        STN_list.append(GPE)

    
    SNR_list=[] # SNR groups for each synapse
    
    
    for i, SNR_i_c in enumerate(SNR_INJECTED_CURRENT):
        I_e=my_nest.GetDefaults(NEURON_MODELS[0])['I_e']+SNR_i_c    
        SNR = MyGroup( NEURON_MODELS[0], n=n_exp, params={'I_e':I_e}, 
                       sd=True, mm=False,
                       mm_dt=.1, record_from=[''])
        SNR_list.append(SNR)

   
    if not load:
        for i_exp in EXPERIMENTS:    
            
            # MSN
            MSN = MSN_list[i_exp]
            
            # Set spike times
            # Base rate
            for id in MSN[0:N_MSN]:                 
                MSN.set_spike_times(id=id, rates=[MSN_BASE_RATE], times=[1], 
                                    t_stop=sim_time,
                                    seed=int(numpy.random.random()*10000.0))               
      
            # Selection MSN        
            for id in MSN[N_MSN:N_MSN+N_MSN_BURST]: 
                rates = [MSN_BASE_RATE, msn_rate_sel, MSN_BASE_RATE]
                times = [1, SEL_ONSET, sel_time_MSN + SEL_ONSET]
                t_stop = sim_time
                MSN.set_spike_times(id=id, rates=rates, times=times, 
                                    t_stop=t_stop, 
                                    seed=int(numpy.random.random()*10000.0))     
        
     
            # GPE
            GPE = GPE_list[i_exp]
            
            # Set spike times
            # Base rate
            for id in GPE[:]:                 
                GPE.set_spike_times(id=id, rates=[GPE_BASE_RATE], times=[1], 
                                    t_stop=sim_time,
                                    seed=int(numpy.random.random()*10000.0))               
      
            # Selection GPE        
            for id in GPE[N_GPE:N_GPE+n_gpe_sel]: 
                rates = [GPE_BASE_RATE, gpe_sel_rate, GPE_BASE_RATE]
                
                # If GPe excited smaller selection time
                times = [1, SEL_ONSET, sel_time_GPE + SEL_ONSET]
                t_stop = sim_time
                GPE.set_spike_times(id=id, rates=rates, times=times, 
                                    t_stop=t_stop, seed=int(numpy.random.random()*100000.0))     

            # Base rate STN
            for id in STN[:]:                 
                STN.set_spike_times(id=id, rates=[STN_BASE_RATE], times=[1], 
                                    t_stop=sim_time,
                                    seed=int(numpy.random.random()*10000.0))     
                
            idx_MSN_s=range(0,N_MSN-N_MSN_BURST)
            idx_MSN_s.extend(range(N_MSN,N_MSN+N_MSN_BURST))
            idx_GPE_s=range(0,N_GPE-n_gpe_sel)
            idx_GPE_s.extend(range(N_GPE,N_GPE+n_gpe_sel))
            
            # Connect with MSN burst         
            target=SNR_list[0][i_exp]
            my_nest.ConvergentConnect(MSN[idx_MSN_s], [target], model=SYNAPSE_MODELS[0])
            my_nest.ConvergentConnect(GPE[0:N_GPE], [target], model=SYNAPSE_MODELS[1])               
            my_nest.ConvergentConnect(STN[:], [target], model=SYNAPSE_MODELS_BACKGROUND[0]) 
      
            # With GPe pause
            target=SNR_list[1][i_exp]
            my_nest.ConvergentConnect(MSN[0:N_MSN], [target], model=SYNAPSE_MODELS[0])
            my_nest.ConvergentConnect(GPE[idx_GPE_s], [target], model=SYNAPSE_MODELS[1])                
            my_nest.ConvergentConnect(STN[:], [target], model=SYNAPSE_MODELS_BACKGROUND[0]) 
            
            # With MSN burst and GPe pause
            target=SNR_list[2][i_exp]
            my_nest.ConvergentConnect(MSN[idx_MSN_s], [target], model=SYNAPSE_MODELS[0])
            my_nest.ConvergentConnect(GPE[idx_GPE_s], [target], model=SYNAPSE_MODELS[1])         
            my_nest.ConvergentConnect(STN[:], [target], model=SYNAPSE_MODELS_BACKGROUND[0]) 
                      
        my_nest.MySimulate( sim_time )

        for MSN in MSN_list: 
            MSN.get_signal( 's' )      
        for GPE in GPE_list: 
            GPE.get_signal( 's' )   
        for SNR in SNR_list: 
            SNR.get_signal( 's' ) 

        misc.pickle_save([MSN_list, GPE_list,SNR_list] , save_at)

    if load:
        MSN_list, GPE_list, SNR_list=misc.pickle_load(save_at)
        
    pre_dyn_MSN=str(SNR_list[0].signals['spikes'].mean_rate(SEL_ONSET-500,
                                                            SEL_ONSET)) 
    pre_dyn_GPE=str(SNR_list[1].signals['spikes'].mean_rate(SEL_ONSET-500,
                                                            SEL_ONSET))   
      
    s='\n'
    s=s+'Example:\n'
    s = s + ' %s %5s %3s \n' % ( 'N experiments:', str ( len(EXPERIMENTS) ),  '#' )  
    s = s + ' %s %5s %3s \n' % ( 'N MSN:', str ( N_MSN ),  '#' )  
    s = s + ' %s %5s %3s \n' % ( 'N GPE:', str ( N_GPE ),  '#' )  
    s='\n'
    s = s + ' %s %5s %3s \n' % ( 'Base rate MSN:',   str ( MSN_BASE_RATE),'spikes/s' )     
    s = s + ' %s %5s %3s \n' % ( 'Sel rate MSN:', str ( msn_rate_sel ), 'spikes/s' )
    s = s + ' %s %5s %3s \n' % ( 'Sel time MSN:', str ( sel_time_MSN ), 'ms' )
    s='\n'
    s = s + ' %s %5s %3s \n' % ( 'Base rate GPe:',   str ( GPE_BASE_RATE),'spikes/s' )   
    s = s + ' %s %5s %3s \n' % ( 'Sel rate GPe:', str ( gpe_sel_rate ), 'spikes/s' )  
    s = s + ' %s %5s %3s \n' % ( 'Sel time GPe:', str ( sel_time_GPE ), 'ms' )
    s = s + ' %s %5s %3s \n' % ( 'Pre sel rate Dyn MSN:', pre_dyn_MSN[0:4], 'spikes/s' )
    s = s + ' %s %5s %3s \n' % ( 'Pre sel rate Dyn GPe:', pre_dyn_GPE[0:4], 'spikes/s' )
      
    return MSN_list, GPE_list, SNR_list, s

    info_string=s
    
    return MSN_hzs, GPE_hzs, data, info_string
print 'Simulation'

# SIMULATION
info_string=''
 
stop = time.time()    
sec = stop - start
m = sec // 60
sec = sec - m*60
print 'Rank %i simulation time: %i minutes, %i seconds' % ( my_nest.Rank(), m, sec )

if not mpiRun:
    
    load_fast=True
    save_at = SPATH+'/'+NEURON_MODELS[0]+'simulate_example'+'-load_fast.pkl' 
    if not load_fast:
        load=False        
        
        MSN_list_sel2, GPE_list_sel2, SNR_list_sel2, s2=simulate_example(MSN_hz=20, 
                                                                     GPE_hz=100, 
                                                                     load=load,
                                                                     n_gpe_sel=N_GPE_BURST,
                                                                     sel_time_GPE=GPE_BURST_TIME)
        
        
        MSN_list_inh, GPE_list_inh, SNR_list_inh, s3=simulate_example(MSN_hz=20, 
                                                                     GPE_hz=0, 
                                                                     load=load,
                                                                     n_gpe_sel=N_GPE_PAUSE,
                                                                      sel_time_GPE=GPE_PAUSE_TIME)
        info_string=info_string+s2+s3  
        misc.pickle_save([SNR_list_sel2,SNR_list_inh, info_string] , save_at)
    else:
        SNR_list_sel2,SNR_list_inh, info_string=misc.pickle_load(save_at) 

       
    # DISPLAY
    plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0, fontsize=12)
    font_size_text = 8
    fig = pylab.figure( facecolor = 'w' )
    
    ax_list = []
    ax_list.append( MyAxes(fig, [ .075, .37, .135, .26 ] ) )    # text box
    ax_list.append( MyAxes(fig, [ .26,  .6,  .165*2.312, .34 ] ) )    # 
    ax_list.append( MyAxes(fig, [ .26,  .1,  .165*2.312, .34 ] ) )    # 
    #ax_list.append( MyAxes(fig, [ .8,   .6,  .165, .34 ] ) )    # 
    #ax_list.append( MyAxes(fig, [ .26,  .1,  .165, .34 ] ) )    # 
    #ax_list.append( MyAxes(fig, [ .53,  .1,  .165, .34 ] ) )    # 
    #ax_list.append( MyAxes(fig, [ .8,   .1,  .165, .34 ] ) )    # 
    
    # Text
    ax=ax_list[0]
    plot_text(ax, info_string)
        
    ax=ax_list[1]
    plot_example_SNR(ax, SNR_list_sel2, flag=0)

    ax=ax_list[2]
    plot_example_SNR(ax, SNR_list_inh, flag=1)



    pylab.show()
    
    fig.savefig( picture_dir + '/' + SNAME  + '.svg', format = 'svg')
    fig.savefig( picture_dir + '/' + SNAME  + '.pdf', format = 'pdf')