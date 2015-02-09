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
model_dir=   '/'.join(os.getcwd().split('/')[0:-1])        
code_dir=  '/'.join(os.getcwd().split('/')[0:-2]) 

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
NEURON_MODELS=['SNR_aeif']

SNR_INJECTED_CURRENT=400.0

SYNAPSE_MODELS_TESTED=['GPE_SNR_gaba_s_ref',  'GPE_SNR_gaba_p']
SYNAPSE_MODELS_BACKGROUND=['MSN_SNR_gaba_p1', 'STN_SNR_ampa_s']

SEL_ONSET = 2000.
ADJUST_XDATA_MS=1500.


N_GPE = 30
N_SEL = 10
N_MSN = 500
N_STN = 100
MSN_RATE_BASE=0.1
STN_RATE_BASE=10.

###############  
  
def plot_example_SNR(ax, SNR_list):
    time_bin=20
    
    colors=['r','c' ] 
    labels=[r'$\delta_{ref}^{GPe}$' , r'$\delta_{dep}^{GPe}$']
    coords=[[0.12, 0.60], [ 0.1, 0.2]]
    
    SNR_list=[SNR_list[0], SNR_list[1]]
    
    for color, label, SNR in zip(colors, labels, SNR_list):
        signal=SNR.signals['spikes']
        signal.my_firing_rate(bin=time_bin, display=ax,
                          kwargs={'color':color})
    
    lines = ax.lines
    lines[0].set_xdata(lines[0].get_xdata()-ADJUST_XDATA_MS)
    lines[1].set_xdata(lines[1].get_xdata()-ADJUST_XDATA_MS)
    
    ax.plot([0, 1190],[SELECTION_THR]*len([0, 2990]),
            **{'color':[0.6,0.6,0.6], 'label':'', 'linewidth':1,
               'linestyle':'--'})
    ax.text( 0.9, 0.13,'Thr' , transform=ax.transAxes, 
             **{ 'color' : [0.6,0.6,0.6]})   
    
    misc.slice_line(lines[0], xlim=[0,1190])
    misc.slice_line(lines[1], xlim=[0,1190])
    
    #ax.set_title('bin=%i'%(time_bin),**{'fontsize':12})
    ax.set_ylabel('Firing rate SNr (spikes/s)') 
    ax.set_xlabel('Time (ms)')
    ax.my_set_no_ticks( yticks=6, xticks = 8 ) 
    
    ax.set_xlim(misc.adjust_limit([0,1200]))
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
    #ax_twinx.set_ylim([0,500])
    ax.set_ylabel('GPe id')


    for line in lines:
        line.set_xdata(line.get_xdata()-ADJUST_XDATA_MS)
    
    ax.set_xlim(misc.adjust_limit([0,1200]))
    ax.set_ylim(misc.adjust_limit([0,N_GPE]))
    ax.my_set_no_ticks( yticks=6, xticks = 5 )   
    
    ax.text( 0.05, 0.05, 'Non-bursting GPe neurons' , backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : 'b' })  
    ax.text( 0.05, 0.15, 'Bursting GPe neurons' , backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : 'r' }) 
    
def plot_example_firing_frequency_GPE(ax, GPE_list):
    global ADJUST_XDATA_MS
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
    
    ax.set_xlim(misc.adjust_limit([0,1200]))
    ax.set_ylim(misc.adjust_limit([20,155]))
    
    lines = ax.lines
    
    for line in lines:
        line.set_xdata(line.get_xdata()-ADJUST_XDATA_MS)
    
    misc.slice_line(lines[0], xlim=[0,1190])
    misc.slice_line(lines[1], xlim=[0,1190])
    misc.slice_line(lines[2], xlim=[0,1190])
    
    #ax.set_title('bin=%i'%(time_bin),**{'fontsize':12})
    ax.set_ylabel('Firing rate GPe (spikes/s)')     
    ax.my_set_no_ticks( yticks=6, xticks = 5 )
    
    ax.text( 0.05, 0.75, 'Non-bursting GPe neurons' , backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : 'b' })  
    ax.text( 0.05, 0.85, 'Bursting GPe neurons' , backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : 'r' }) 
    ax.text( 0.05, 0.65, 'Average over all' , backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : 'k' })  
    
def plot_text(ax, info_string=''):
    
    my_nest.ResetKernel()
    MODEL_LIST,model_dic=models()
    my_nest.MyLoadModels( model_dic, NEURON_MODELS )
    
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

def simulate_example(hz=40, load=True):
    global SNR_INJECTED_CURRENT
    global NEURON_MODELS
    global N_GPE
    global N_STN
    global N_SEL
    global N_MSN
    global MSN_RATE_BASE  
    global SNAME
    global SPATH
    global SYNAPSE_MODELS
    global SEL_ONSET
    

    #N_EXP = 20
    N_EXP = 200

    RATE_BASE = 25 # Base rate
    RATE_SELE = hz # Selection rate 
    SAVE_AT = SPATH+'/'+NEURON_MODELS[0]+'-example.pkl'   
    SEL_TIME = 200.
    sim_time = SEL_TIME+SEL_ONSET+500.
    SNAME_NB = hz+1000    
    
    EXPERIMENTS=range(N_EXP)
   
    MODEL_LIST,model_dic=models()
    my_nest.ResetKernel()       
    my_nest.MyLoadModels( model_dic, NEURON_MODELS )
    my_nest.MyLoadModels( model_dic, SYNAPSE_MODELS_TESTED)
    my_nest.MyLoadModels( model_dic, SYNAPSE_MODELS_BACKGROUND)       
 
    GPE_list=[] # GPE input for each experiment
    for i_exp in EXPERIMENTS:
        GPE = MyPoissonInput( n=N_GPE, sd=True, 
                              spath=SPATH, sname_nb=SNAME_NB+i_exp)
        GPE_list.append(GPE)
    
    MSN_list=[] # MSN input for each experiment
    for i_exp in EXPERIMENTS:
        MSN = MyPoissonInput( n=N_MSN, sd=True )
        MSN_list.append(MSN)
        
    STN_list=[] # STN input for each experiment
    for i_exp in EXPERIMENTS:
        STN = MyPoissonInput( n=N_STN, sd=True )
        STN_list.append(STN)
    
    SNR_list=[] # SNR groups for each synapse
    for i_syn in SYNAPSE_MODELS_TESTED:
        I_e=my_nest.GetDefaults(NEURON_MODELS[0])['I_e']+SNR_INJECTED_CURRENT
            
        SNR = MyGroup( NEURON_MODELS[0], n=N_EXP, params={'I_e':I_e}, 
                       sd=True, mm=False,
                       mm_dt=.1, record_from=[''])
        SNR_list.append(SNR)
    
           
    if not load:
        for i_exp in EXPERIMENTS:    
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
                STN.set_spike_times(id=id, rates=[STN_RATE_BASE], times=[1], 
                                    t_stop=sim_time,
                                    seed=int(numpy.random.random()*10000.0)) 
            
            # Base rate
            for id in GPE[0:N_GPE-N_SEL]:                 
                GPE.set_spike_times(id=id, rates=[RATE_BASE], times=[1], 
                                    t_stop=sim_time, 
                                    seed=int(numpy.random.random()*10000.0))               
      
            # Selection        
            for id in GPE[N_GPE-N_SEL:N_GPE]: 
                rates = [RATE_BASE, RATE_SELE, RATE_BASE]
                times = [1, SEL_ONSET, SEL_TIME + SEL_ONSET]
                t_stop = sim_time
                GPE.set_spike_times(id=id, rates=rates, times=times, 
                                    t_stop=t_stop, 
                                    seed=int(numpy.random.random()*10000.0))    
        
            # Connect         
            for i, syn in enumerate(SYNAPSE_MODELS_TESTED):       
                target=SNR_list[i][i_exp]
                my_nest.ConvergentConnect(GPE[:], [target], 
                                          model=syn)
                my_nest.ConvergentConnect(MSN[:], [target], 
                                          model=SYNAPSE_MODELS_BACKGROUND[0])
                my_nest.ConvergentConnect(STN[:], [target], 
                                          model=SYNAPSE_MODELS_BACKGROUND[1])
                      
        my_nest.MySimulate( sim_time )

    
        for GPE in GPE_list: 
            GPE.get_signal( 's' )   
        for SNR in SNR_list: 
            SNR.get_signal( 's' ) 

        misc.pickle_save([GPE_list,SNR_list] , SAVE_AT)

    elif load:
        GPE_list, SNR_list=misc.pickle_load(SAVE_AT)
        
    pre_ref=str(SNR_list[0].signals['spikes'].mean_rate(SEL_ONSET-500,SEL_ONSET)) 
    pre_dyn=str(SNR_list[1].signals['spikes'].mean_rate(SEL_ONSET-500,SEL_ONSET))   
    
    statusSynapes=[]
    for syn in SYNAPSE_MODELS_TESTED:
        statusSynapes.append( my_nest.GetDefaults(syn ))
      
    s='\n'
    s=s+'Example:\n'
    s = s + ' %s %5s %3s \n' % ( 'N experiments:', str ( len(EXPERIMENTS) ),  '#' )  
    s = s + ' %s %5s %3s \n' % ( 'N GPEs:', str ( N_GPE ),  '#' )  
    
    s = s + ' %s %5s %3s \n' % ( 'Base rate:',   str ( RATE_BASE),'spikes/s' )     
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
    

print 'Simulation'

# SIMULATION
info_string='' 

stop = time.time()    
sec = stop - start
m = sec // 60
sec = sec - m*60
print 'Rank %i simulation time: %i minutes, %i seconds' % ( my_nest.Rank(), m, sec )

if not mpiRun:
    
    
    GPE_list, SNR_list, s=simulate_example(hz=100, load=True)
    info_string=info_string+s  
        
    # DISPLAY
    plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, 
                           h = 450.0, fontsize=16)
    font_size_text = 8
    fig = pylab.figure( facecolor = 'w' )
    
    ax_list = []
    ax_list.append( MyAxes(fig, [ .075, .37, .135, .26 ] ) )    # text box
    ax_list.append( MyAxes(fig, [ .26,  .6,  .165, .34 ] ) )    # 
    ax_list.append( MyAxes(fig, [ .53,  .6,  .165, .34 ] ) )    # 
    #ax_list.append( MyAxes(fig, [ .8,   .6,  .165, .34 ] ) )    # 
    ax_list.append( MyAxes(fig, [ .26,  .1,  .165*2.312, .34 ] ) )    # 
    #ax_list.append( MyAxes(fig, [ .53,  .1,  .165, .34 ] ) )    # 
    #ax_list.append( MyAxes(fig, [ .8,   .1,  .165, .34 ] ) )    # 
    
    # Text
    ax=ax_list[0]
    plot_text(ax, info_string)
    
    ax=ax_list[1]
    plot_example_raster_GPE(ax, GPE_list)
    
    ax=ax_list[2]
    plot_example_firing_frequency_GPE(ax, GPE_list)
    
    ax=ax_list[3]
    plot_example_SNR(ax, SNR_list)

    pylab.show()
    
    # dpi does not matter since svg and pdf are vetorbased
    fig.savefig( picture_dir + '/' + SNAME  + '.svg', format = 'svg') 
    fig.savefig( picture_dir + '/' + SNAME  + '.pdf', format = 'pdf')