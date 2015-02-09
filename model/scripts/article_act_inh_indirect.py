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

numpy.random.seed(1) # set random seed
 
# Add directories to python path
sys.path.append(os.getcwd())                            
parent_dir='/'.join(os.getcwd().split('/')[0:-1])                     
model_dir= '/'.join(os.getcwd().split('/')[0:-1])        
code_dir='/'.join(os.getcwd().split('/')[0:-2]) 

model_name=os.getcwd().split('/')[-2]
picture_dir='/'.join(os.getcwd().split('/')[0:-3]) + '/pictures/'+model_name 

sys.path.append(model_dir) 
sys.path.append(code_dir+'/nest_toolbox') 
SPATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]
SNAME  = sys.argv[0].split('/')[-1].split('.')[0]

from model_params import models                               # Then import models     
from src import misc, my_nest, my_signals, plot_settings
from src.my_population import MyGroup, MyPoissonInput 
from src.my_axes import MyAxes 


### GLOBALS ###
SELECTION_THR=5.  # spikes/s

#NEURON_MODELS=['SNR_aeif']
NEURON_MODELS=['SNR_izh']
SYNAPSE_MODELS_TESTED=['GPE_SNR_gaba_s_ref', 'GPE_SNR_gaba_p']
SYNAPSE_MODELS_BACKGROUND=['MSN_SNR_gaba_p1', 'STN_SNR_ampa_s']

SNR_INJECTED_CURRENT=400.0
SEL_ONSET=2000.
ADJUST_XDATA_MS=1500.
SEL_TIME = 500.
    
N_GPE = 30
N_SEL = 10
N_MSN = 500
N_STN = 100

MSN_RATE_BASE=0.1
GPE_BASE_RATE=25
STN_BASE_RATE=10

################
    
def plot_example_SNR(ax, SNR_list):
    global N_GPE
    time_bin=20
    
    colors=['r','c' ] 
    labels=[r'$\delta_{ref}^{GPe}$' , r'$\delta_{dep}^{GPe}$']
    coords=[[0.05, 0.7], [ 0.05, 0.2]]
    
    SNR_list=[SNR_list[0], SNR_list[1]]
    
    
    for color, SNR in zip(colors, SNR_list):
        signal=SNR.signals['spikes']
        signal.my_firing_rate(bin=time_bin, display=ax,
                          kwargs={'color':color})
    
    lines = ax.lines
    lines[0].set_xdata(lines[0].get_xdata()-ADJUST_XDATA_MS)
    lines[1].set_xdata(lines[1].get_xdata()-ADJUST_XDATA_MS)
    
    misc.slice_line(lines[0], xlim=[0,1490])
    misc.slice_line(lines[1], xlim=[0,1490])
    
    #ax.set_title('bin=%i'%(time_bin),**{'fontsize':12})
    ax.set_ylabel('Firing rate SNr (spikes/s)') 
    ax.set_xlabel('Time (ms)')
    
    ax.my_set_no_ticks( yticks=6, xticks = 8 ) 
    
    ax.set_xlim(misc.adjust_limit([0,1500]))
    ax.set_ylim(misc.adjust_limit([0,50]))

    
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['text.fontsize'], 
                 **{'color': color})   
 
def plot_example_raster_GPE(ax, GPE_list):
    global ADJUST_XDATA_MS
    global N_GPE
    global N_SEL
    
    GPE=GPE_list[0]
    
    GPE.signals['spikes'].raster_plot(id_list=GPE[N_GPE-N_SEL:], display=ax,
                                      kwargs={'color':'r', 'zorder':1})  
    GPE.signals['spikes'].raster_plot(id_list=GPE[0:N_GPE-N_SEL], display=ax,
                                      kwargs={'color':'b', 'zorder':1})  
    
    lines = ax.lines
    ax.set_ylabel('GPe id')


    for line in lines:
        line.set_xdata(line.get_xdata()-ADJUST_XDATA_MS)
    
    ax.set_xlim(misc.adjust_limit([0,1500]))
    ax.set_ylim(misc.adjust_limit([0,N_GPE]))
    ax.my_set_no_ticks( yticks=6, xticks = 5 )   
    
    ax.text( 0.05, 0.05, 'Non-pausing GPe neurons' , backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : 'b' })  
    ax.text( 0.05, 0.15, 'Pausing GPe neurons' , backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : 'r' }) 
    
def plot_example_firing_frequency_GPE(ax, GPE_list):
    global ADJUST_XDATA_MS
    global N_GPE
    global N_SEL
    time_bin=20
    
    GPE_a=GPE_list[0]
    GPE_s=GPE_a.slice(GPE_a.ids[N_GPE-N_SEL:])
    GPE_b=GPE_a.slice(GPE_a.ids[0:N_GPE-N_SEL])
    
    for g in GPE_list[1:]:
        GPE_a.merge(g)
        GPE_s.merge(g.slice(g.ids[N_GPE-N_SEL:]))
        GPE_b.merge(g.slice(g.ids[0:N_GPE-N_SEL]))
    
    GPE_a.signals['spikes'].my_firing_rate(bin=time_bin, display=ax,
                                          kwargs={'color':'k'})
    GPE_s.signals['spikes'].my_firing_rate( bin=time_bin, display=ax,
                                          kwargs={'color':'r'})
    GPE_b.signals['spikes'].my_firing_rate(bin=time_bin, display=ax,
                                          kwargs={'color':'b'})
    
    
    ax.set_xlim(misc.adjust_limit([0,1500]))
    ax.set_ylim(misc.adjust_limit([0,50]))
    
    lines = ax.lines
    
    for line in lines:
        line.set_xdata(line.get_xdata()-ADJUST_XDATA_MS)
    
    misc.slice_line(lines[0], xlim=[0,1490])
    misc.slice_line(lines[1], xlim=[0,1490])
    misc.slice_line(lines[2], xlim=[0,1490])
    
    #ax.set_title('bin=%i'%(time_bin),**{'fontsize':12})
    ax.set_ylabel('Firing rate GPe (spikes/s)')     
    ax.my_set_no_ticks( yticks=6, xticks = 5 )
    
    ax.text( 0.05, 0.75, 'Non-pausing GPe neurons' , backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : 'b' })  
    ax.text( 0.05, 0.85, 'Pausing GPe neurons' , backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : 'r' }) 
    ax.text( 0.05, 0.65, 'Average over all' , backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : 'k' })  
    
def plot_selection_vs_neurons(ax, nbNeurons, meanRates):
    global N_GPE    
    print nbNeurons
    print meanRates
    
    
    colors=['r','c' ] 
    labels=[r'$\delta_{ref}^{GPe}$' , r'$\delta_{dep}^{GPe}$']
    coords=[[0.1, 0.4], [ 0.1, 0.1]]
    
    line1 = ax.plot(nbNeurons,meanRates[0,:], **{'color':colors[0]})  
    line2 = ax.plot(nbNeurons,meanRates[1,:], **{ 'color':colors[1]})

    ax.set_ylabel('Firing rate SNr (spikes/s)') 
    ax.set_xlabel('Paused GPEs (#)')

    ax.my_set_no_ticks( yticks=8, xticks=5 ) 
    #ax.set_xlim(misc.adjust_limit([0, N_GPE]))
    #ax.set_ylim(misc.adjust_limit([0,110]))
    
 
    
    
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['text.fontsize'], 
                 **{'color': color})  

def plot_SNr_rate_vs_syn_event1(ax, syn_events, meanRates):
    
    colors=['k']   
    labels=[ r'$\delta_{dep}^{GPe}$']
    coords=[[0.1, 0.15]]   
    
    for id, label, color in zip([0],labels,colors):
        ax.plot(syn_events,meanRates[id,:],**{'color':color})  
 
    vec=[0,800]

    ax.set_ylabel('Firing rate SNr (spikes/s)') 
    ax.set_xlabel('Synaptic events (#/s)')
    ax.my_set_no_ticks( yticks=6, xticks = 5 ) 
    
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['text.fontsize'], 
                 **{'color': color})
    
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,800])  
    
    ax.set_xlim(misc.adjust_limit([0,800]))
    ax.set_ylim(misc.adjust_limit([20,115]))

def plot_SNr_rate_vs_syn_event2(ax, syn_events, meanRates):
    colors=['k'] 
    ax.plot(syn_events,meanRates[1,:],**{'color':colors[0], 'linestyle':'-.'})  
 
    ax.set_ylabel('Firing rate SNr (spikes/s)') 
    ax.set_xlabel('Synaptic events (#/s)')
    ax.my_set_no_ticks( yticks=5, xticks = 5 ) 
    
    line1=ax.plot(1,1,'-k')
    line2=ax.plot(1,1,'-.k')
    leg=ax.legend([line1, line2],['Pausing subpopulation', 'Increased activity in \nall GPe neurons'], loc='best')
    frame  = leg.get_frame() 
    frame.set_visible(False) 
    
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    pylab.setp(ltext, fontsize=10., backgroundcolor='w') 
    
    ax.set_xlim(misc.adjust_limit([0,800]))
    ax.set_ylim(misc.adjust_limit([20,110]))

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
    global SNR_INJECTED_CURRENT
    global NEURON_MODELS
    global N_GPE
    global N_SEL
    global N_MSN
    global N_STN
    global MSN_RATE_BASE  
    global STN_BASE_RATE      
    global SNAME
    global SPATH
    global SYNAPSE_MODELS
    global SEL_ONSET
    global GPE_BASE_RATE

    #n_exp = 20
    n_exp = 200

    RATE_SELE = hz # Selection rate 
    save_at = SPATH+'/'+NEURON_MODELS[0]+'-example.pkl'   

    sim_time = SEL_TIME+SEL_ONSET+500.
    SNAME_NB = hz+1000    
     
    experiments=range(n_exp)
       
    MODEL_LIST=models()
    my_nest.ResetKernel()       
    my_nest.MyLoadModels( MODEL_LIST, NEURON_MODELS)
    my_nest.MyLoadModels( MODEL_LIST, SYNAPSE_MODELS_TESTED)   
    my_nest.MyLoadModels( MODEL_LIST, SYNAPSE_MODELS_BACKGROUND)      
 
 
    GPE_list=[] # GPE input for each experiment
    for i_exp in experiments:
        GPE = MyPoissonInput( n=N_GPE, sd=True, spath=SPATH, 
                              sname_nb=SNAME_NB+i_exp)
        GPE_list.append(GPE)
    
    MSN_list=[] # MSN input for each experiment
    for i_exp in experiments:
        MSN = MyPoissonInput( n=N_MSN, sd=False )
        MSN_list.append(MSN)
    
    
    STN_list=[] # MSN input for each experiment
    for i_exp in experiments:
        STN = MyPoissonInput( n=N_STN, sd=False )
        STN_list.append(STN)
    
    SNR_list=[] # SNR groups for each synapse
    I_e=my_nest.GetDefaults(NEURON_MODELS[0])['I_e']+SNR_INJECTED_CURRENT
    for i_syn in range(len(SYNAPSE_MODELS_TESTED)):
        SNR = MyGroup( NEURON_MODELS[0], n=n_exp, params={'I_e':I_e}, sd=True)
        SNR_list.append(SNR)
        
    if not load:
        for i_exp in experiments:    
            GPE = GPE_list[i_exp]
            MSN = MSN_list[i_exp]
            STN = STN_list[i_exp]
            
            
            # Set spike times
            # Base rate MSN
            for id in MSN[:]:                 
                MSN.set_spike_times(id=id, rates=[MSN_RATE_BASE], times=[1], 
                                    t_stop=sim_time,
                                    seed=int(numpy.random.random()*10000.0))  
            # Base rate STN
            for id in STN[:]:                 
                STN.set_spike_times(id=id, rates=[STN_BASE_RATE], times=[1], 
                                    t_stop=sim_time,
                                    seed=int(numpy.random.random()*10000.0))  
                            
            # Set spike times
            # Base rate
            for id in GPE[0:N_GPE-N_SEL]:                 
                GPE.set_spike_times(id=id, rates=[GPE_BASE_RATE], times=[1], 
                                    t_stop=sim_time, 
                                    seed=int(numpy.random.random()*10000.0))               
      
            # Selection        
            for id in GPE[N_GPE-N_SEL:N_GPE]: 
                rates = [GPE_BASE_RATE, RATE_SELE, GPE_BASE_RATE]
                times = [1, SEL_ONSET, SEL_TIME + SEL_ONSET]
                t_stop = sim_time
                GPE.set_spike_times(id=id, rates=rates, times=times, 
                                    t_stop=t_stop, 
                                    seed=int(numpy.random.random()*10000.0))     
        
            # Connect         
            for i_syn, syn in enumerate(SYNAPSE_MODELS_TESTED):       
                target=SNR_list[i_syn][i_exp]
                my_nest.ConvergentConnect(GPE[:], [target], model=syn)
                my_nest.ConvergentConnect(MSN[:], [target], 
                                          model=SYNAPSE_MODELS_BACKGROUND[0])
                my_nest.ConvergentConnect(STN[:], [target], 
                                          model=SYNAPSE_MODELS_BACKGROUND[1])
                      
        my_nest.MySimulate( sim_time )

    
        for GPE in GPE_list: 
            GPE.get_signal( 's' )   
        for SNR in SNR_list: 
            SNR.get_signal( 's' ) 

        misc.pickle_save([GPE_list,SNR_list] , save_at)

    elif load:
        GPE_list, SNR_list=misc.pickle_load(save_at)
        
    pre_ref=str(SNR_list[0].signals['spikes'].mean_rate(SEL_ONSET-5000,SEL_ONSET)) 
    pre_dyn=str(SNR_list[1].signals['spikes'].mean_rate(SEL_ONSET-500,SEL_ONSET))    
    
    statusSynapes=[]
    for syn in SYNAPSE_MODELS_TESTED:
        statusSynapes.append( my_nest.GetDefaults(syn))
      
    s='\n'
    s=s+'Example:\n'
    s = s + ' %s %5s %3s \n' % ( 'N experiments:', str ( len(experiments) ),  '#' )  
    s = s + ' %s %5s %3s \n' % ( 'N GPEs:', str ( N_GPE ),  '#' )     
  
    s = s + ' %s %5s %3s \n' % ( 'Base rate:',   str ( GPE_BASE_RATE),'spikes/s' )     
    s = s + ' %s %5s %3s \n' % ( 'Selection rate:', str ( RATE_SELE ), 'spikes/s' )
    s = s + ' %s %5s %3s \n' % ( 'Selection time:', str ( SEL_TIME ), 'ms' )
    s = s + ' %s %5s %3s \n' % ( 'Pre sel rate Ref:', pre_ref[0:4], 'spikes/s' )
    s = s + ' %s %5s %3s \n' % ( 'Pre sel rate Dyn:', pre_dyn[0:4], 'spikes/s' )
    for ss in statusSynapes:
        s = s + '\n'  
        s = s + ' %s %10s\n' % ( 'Synapse', ss['synapsemodel'])   
        s = s + ' %s %5s %3s\n' % ( 'Weight', 
                                      str( round( ss['weight'], 1) ), 'nS')
    
      
    return GPE_list, SNR_list, s
    
def simulate_selection_vs_neurons(selRateInterval=[0.0,500.0], hz=0, load=True):    
    global SEL_ONSET
    global SNR_INJECTED_CURRENT
    global N_GPE
    global N_STN
    global MSN_RATE_BASE
    global GPE_BASE_RATE
    global STN_BASE_RATE
    SNAME_NB=hz

    #n_exp=200
    n_exp=20
    N_MAX_SEL=N_GPE+1 # Plus one to account for the case when all GPe have stopped 
      
    RATE_SELE = hz
    save_at = (SPATH+'/'+NEURON_MODELS[0]+'-' + str(SNAME_NB)+'-hz.pkl') 
    #SEL_TIME = 1000.

    sim_time = SEL_TIME+SEL_ONSET+500.   
    
    experiments=range(n_exp)

    MODEL_LIST=models()
    my_nest.ResetKernel()       
    my_nest.MyLoadModels( MODEL_LIST, NEURON_MODELS )
    my_nest.MyLoadModels( MODEL_LIST, SYNAPSE_MODELS_TESTED)   
    my_nest.MyLoadModels( MODEL_LIST, SYNAPSE_MODELS_BACKGROUND)   

    
    GPE_list=[] # GPE input for each experiment
    for i_exp in experiments:
        GPE = MyPoissonInput( n=N_GPE+N_MAX_SEL)
        GPE_list.append(GPE)
    
    MSN_list=[] # MSN input for each experiment
    for i_exp in experiments:
        MSN = MyPoissonInput( n=N_MSN, sd=True)
        MSN_list.append(MSN)

    STN_list=[] # STN input for each experiment
    for i_exp in experiments:
        STN = MyPoissonInput( n=N_STN, sd=True)
        STN_list.append(STN)
    
    SNR_list=[] # SNR groups for each synapse and number of selected GPE
    I_e=my_nest.GetDefaults(NEURON_MODELS[0])['I_e']+SNR_INJECTED_CURRENT
    for i, i_syn in enumerate(SYNAPSE_MODELS_TESTED):
        SNR = []
        for i_sel in range(N_MAX_SEL):
            SNR.append(MyGroup( NEURON_MODELS[0], n=n_exp, 
                       params={'I_e':I_e}, sd=True,
                              sd_params={'start':0., 'stop':sim_time}))

        SNR_list.append(SNR)    
        
    if not load:
        for i_exp in experiments:    
            GPE = GPE_list[i_exp]
            MSN = MSN_list[i_exp]
            STN = STN_list[i_exp]
        
            # Set spike times
            # Base rate
            for id in MSN[:]:                 
                MSN.set_spike_times(id=id, rates=[MSN_RATE_BASE], times=[1], 
                                    t_stop=sim_time,
                                    seed=int(numpy.random.random()*10000.0))  
            # Base rate STN
            for id in STN[:]:                 
                STN.set_spike_times(id=id, rates=[STN_BASE_RATE], times=[1], 
                                    t_stop=sim_time,
                                    seed=int(numpy.random.random()*10000.0))          
            # Set spike times
            # Base rate
            for id in GPE[0:N_GPE]:                 
                GPE.set_spike_times(id=id, rates=[GPE_BASE_RATE], times=[1], 
                                    t_stop=sim_time,
                                    seed=int(numpy.random.random()*10000.0))               
      
            # Selection        
            for id in GPE[N_GPE:N_GPE + N_MAX_SEL]: 
                rates = [GPE_BASE_RATE, RATE_SELE, GPE_BASE_RATE]
                times = [1, SEL_ONSET, SEL_TIME + SEL_ONSET]
                t_stop = sim_time
                GPE.set_spike_times(id=id, rates=rates, times=times, 
                                    t_stop=t_stop,
                                    seed=int(numpy.random.random()*10000.0))    
            
            # Connect     
            for i, syn in enumerate(SYNAPSE_MODELS_TESTED):  
                    for i_sel in range(N_MAX_SEL):       
                        target = SNR_list[i][i_sel][i_exp]
                        my_nest.ConvergentConnect(GPE[0:N_GPE-i_sel], [target], 
                                                  model=syn)
                        my_nest.ConvergentConnect(GPE[N_GPE:N_GPE+i_sel], 
                                                  [target], model=syn)
                        my_nest.ConvergentConnect(MSN[:], [target], 
                                          model=SYNAPSE_MODELS_BACKGROUND[0])
                        my_nest.ConvergentConnect(STN[:], [target], 
                                          model=SYNAPSE_MODELS_BACKGROUND[1])
                      
        my_nest.MySimulate( sim_time ) 
        
        
        for SNR_sel in SNR_list: 
            for SNR in SNR_sel: 
                SNR.get_signal( 's' ) 
        
        misc.pickle_save_groups(SNR_list, save_at)

    elif load:
        SNR_list=misc.pickle_load_groups(save_at)
     
    t1=selRateInterval[0]
    t2=selRateInterval[1]
    mean_rates=[]    
    
    delay=my_nest.GetDefaults(SYNAPSE_MODELS_TESTED[0])['delay'] 
    for SNR_sel in SNR_list: 
        m_r=[]
        for SNR in SNR_sel:  
            m_r_pre=SNR.signals['spikes'].mean_rate(SEL_ONSET-(t2-t1), 
                                                    SEL_ONSET)
            m_r_post=SNR.signals['spikes'].mean_rate(SEL_ONSET+t1+delay, 
                                                    SEL_ONSET+t2+delay)
            m_r.append(m_r_post)  
        mean_rates.append(m_r)
    mean_rates = numpy.array(mean_rates)
    nb_neurons = numpy.arange(0,N_MAX_SEL,1)
    
    s='\n'
    s = s + ' %s %5s %3s \n' % ( 'N GPEs:', str ( N_GPE ),  '#' )     
    s = s + ' %s %5s %3s \n' % ( 'N experiments:', str ( n_exp ),  '#' )    
    s = s + ' %s %5s %3s \n' % ( 'Base rate:',   str ( GPE_BASE_RATE),'spikes/s' )     
    s = s + ' %s %5s %3s \n' % ( 'Selection rate:', str ( RATE_SELE ), 'spikes/s' )
    s = s + ' %s %5s %3s \n' % ( 'Selection time:', str ( SEL_TIME ), 'ms' )
    s = s + ' %s %5s %3s \n' % ( 'SNR_INJECTED_CURRENT:', str ( SNR_INJECTED_CURRENT ), 'pA' )
    

    info_string=s
    
    return nb_neurons, mean_rates, info_string

print 'Simulation'

# SIMULATION
info_string=''

stop = time.time()    
sec = stop - start
m = sec // 60
sec = sec - m*60
print 'Rank %i simulation time: %i minutes, %i seconds' % ( my_nest.Rank(), m, sec )

if not mpiRun:
    
    load=True
    GPE_list, SNR_list, s=simulate_example(hz=0, load=True)
    info_string=info_string+s  
    

    sel_interval_2=[0.0, 500.0]
    nb, mr, s = simulate_selection_vs_neurons(sel_interval_2, load=True)
    
    info_string=info_string+s     
    
    # DISPLAY
    plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0, fontsize=12)
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
    plot_example_raster_GPE(ax, GPE_list)
    
    ax=ax_list[2]
    plot_example_firing_frequency_GPE(ax, GPE_list)
    
    ax=ax_list[3]
    plot_example_SNR(ax, SNR_list)

    ax=ax_list[4]   
    plot_selection_vs_neurons(ax, nb, mr)  

    ax=ax_list[5]   
    # Article filtering has to be run before
    data=misc.pickle_load(os.getcwd()+'/output/mean_rates_GPE_constant_supression'+
    NEURON_MODELS[0])
    GPEmeanRates=data['GPE_mean_rates']
    SNRmeanRates=data['SNR_mean_rates']
    print GPEmeanRates
    print SNRmeanRates
    syn_events=GPEmeanRates*N_GPE
    #SNRmeanRates[1,syn_events<600]
    plot_SNr_rate_vs_syn_event2(ax, syn_events, SNRmeanRates)
       
    #syn_events=nbNeurons1*20.0+(N_MSN-nbNeurons1)*0.1
    plot_SNr_rate_vs_syn_event1(ax, (max(nb)-nb)*GPE_BASE_RATE, mr)

    
    
    pylab.show()
    
    # dpi does not matter since svg and pdf are both vectorbased
    fig.savefig( picture_dir + '/' + SNAME  + '.svg', format = 'svg') 
    fig.savefig( picture_dir + '/' + SNAME  + '.pdf', format = 'pdf')