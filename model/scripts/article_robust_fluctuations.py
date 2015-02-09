import numpy
import pylab
import os
import sys


if len(sys.argv) != 1: mpiRun = True
else:                  mpiRun = False

sys.path.append(os.getcwd())    # Add current directory to python path                                                
current_path=os.getcwd()

# First add parent directory to python path
model_dir=   '/'.join(os.getcwd().split('/')[0:-1])        
code_dir=  '/'.join(os.getcwd().split('/')[0:-2]) 

mode_name=os.getcwd().split('/')[-2]
picture_dir='/'.join(os.getcwd().split('/')[0:-3]) + '/pictures/'+mode_name   
                
sys.path.append(model_dir) 
sys.path.append(code_dir+'/nest_toolbox') 
spath  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]
sname_nb=0

# Then import model and network
from model_params import models
from src import misc, my_nest, my_signals, plot_settings
from src.my_population import MyGroup, MyPoissonInput  
from src.my_axes import MyAxes 


LOAD=True
ADJUST_XDATA_MS=300
SELECTION_THR=5.  # spikes/s

#NEURON_MODELS=['SNR_aeif']
#I_E_1=190
#I_E_2=0.

NEURON_MODELS=['SNR_izh']
I_E_1=120.+ 40+ 320  #
I_E_2=120.+40


synapseModels=['MSN_SNR_gaba_p0', 'MSN_SNR_gaba_p1','MSN_SNR_gaba_p2']

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
N_GPE=30
N_STN=100

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
SEL_INTERVAL_1 = [SEL_ONSET , SEL_ONSET+500.0]    # start, stop first selection
SEL_INTERVAL_2 = [SEL_ONSET+1000. , SEL_ONSET+ 1500.0]    # start, stop second selection


def plot_example_snr(ax, SNR):
    time_bin=20
    signal=SNR[0].signals['spikes']
    colors = ['k']
    
    sim_time=2500
    
    signal.my_firing_rate(bin=time_bin, 
                          display=ax, kwargs={'color':colors[0]})

    ax.my_set_no_ticks( yticks=6, xticks=7 ) 
    lines = ax.lines
    for line in lines:
        line.set_xdata(line.get_xdata()-ADJUST_XDATA_MS)
    
    
    ax.plot([0, sim_time],[SELECTION_THR]*len([0, sim_time]),**{'color':[0.6,0.6,0.6],
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
    ax.text( 0.9, 0.14,'Thr' , transform=ax.transAxes, **{ 'color' : [0.6,0.6,0.6]})   
    misc.slice_line(lines[0], xlim=[0,sim_time])
    
    ax.plot([1000,1500],[39,39],color='k', marker='|')
    ax.text( 0.47, 0.84,'stop' , transform=ax.transAxes, **{ 'color' : 'k' }) 

    ax.plot([500,700],[30,30],color='k', marker='|')
    ax.text( 0.23, 0.67,'1st' , transform=ax.transAxes, **{ 'color' : 'k' })
    
    ax.plot([1500,1700],[30,30],color='k', marker='|') 
    ax.text( 0.60, 0.67, '2nd' , transform=ax.transAxes, **{ 'color' : 'k' }) 

    ax.set_xlim(misc.adjust_limit([0,sim_time]))
    ax.set_ylim(misc.adjust_limit([0,45]))
    
    ax.set_ylabel('Firing rate SNr (spikes/s)') 
    
def plot_rate_first_and_second_bursts_full(ax, x, data):
    colors = ['k']
    max_delay=3100
    
    ax.text( 0.1, 0.4, r'$\delta_{fac}^{MSN}$' , 
             fontsize=pylab.rcParams['text.fontsize'], 
             transform=ax.transAxes, 
             **{ 'color' : colors[0]})
     
    line1=ax.plot(1,1,'-k')
    line2=ax.plot(1,1,'-.k')
    leg=ax.legend([line1, line2],['1st', '2nd'], loc='best')
    frame  = leg.get_frame() 
    frame.set_visible(False)
    
    ax.plot([0, max_delay],[SELECTION_THR]*len([0, max_delay]),**{'color':[0.6,0.6,0.6],
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
    ax.text( 0.8, 0.3,'Thr' , transform=ax.transAxes, **{ 'color' : [0.6,0.6,0.6]})   

    
    ax.plot(x,data['rates'][:,0],**{'label':'1st', 'color': colors[0]})
    ax.plot(x,data['rates'][:,1],**{'label':'2nd', 'color': colors[0], 
                                    'linestyle':'-.'})
    ax.set_ylabel('Firing rate SNr (spikes/s)') 
    ax.set_xlabel('Burst stop (ms)')
    ax.my_set_no_ticks( yticks=6, xticks = 5 )

    lines = ax.lines
    misc.slice_line(lines[2], xlim=[0,max_delay])
    misc.slice_line(lines[3], xlim=[0,max_delay])  

    ax.set_xlim(misc.adjust_limit([0,max_delay]))
    ax.set_ylim(misc.adjust_limit([3,10]))
    
def plot_text(ax, infoString=''):
    
    my_nest.ResetKernel()
    model_list=models()
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    
    SNR = MyGroup( NEURON_MODELS[0], 1, mm_dt = 0.1)
    statusSNR = my_nest.GetStatus( SNR[:] )[0]
    
    tb = ''
  
    tb = tb + ' %s %10s\n' % ( 'Neuron model', statusSNR['model'])   
    tb = tb + infoString
    
    tb = tb + '\n'

    ax.text( 0.85, 0.5, tb , fontsize= font_size_text,
             horizontalalignment='right',
             verticalalignment='center',
             transform=ax.transAxes,     # to define coordinates in right scale
             **{ 'fontname' : 'monospace' })                           
    
    ax.my_remove_axis( xaxis=True, yaxis=True )
    ax.my_remove_spine(left=True,  bottom=True, right=True, top=True)

def simulate_example(load=True):

    global GPE_BASE_RATE  
    global FILE_NAME
    global N_GPE
    global N_STN
    global N_MSN_BURST
    global N_MSN
    global NEURON_MODELS
    global OUTPUT_PATH
    global SEL_ONSET
    global SNR_INJECTED_CURRENT
    global SYNAPSE_MODELS_TESTED
    
    #n_exp =200 # number of experiments   
    n_exp =200 # number of experiments   
    
    # Path were raw data is saved. For example the spike trains.
    save_result_at=OUTPUT_PATH+'/simulate_example.pkl'
    save_header_at=OUTPUT_PATH+'/simulate_example_header'
    
    burst_time = 500.
    sim_time = SEL_INTERVAL_2[1]+500
      
    model_list=models()
    my_nest.ResetKernel(threads=8)       
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    my_nest.MyLoadModels( model_list, SYNAPSE_MODELS_TESTED)       
    my_nest.MyLoadModels( model_list, SYNAPSE_MODELS_BACKGROUND)      
 
    SNR_list=[] # List with SNR groups for synapse. 
    if not load:    
        MSN_base=MyPoissonInput(n=N_MSN_BASE*n_exp)
        MSN_burst=MyPoissonInput(n=N_MSN_BURST*n_exp)     
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
            rates = [MSN_BASE_RATE, MSN_BURST_RATE, MSN_BASE_RATE,
                     MSN_BURST_RATE, MSN_BASE_RATE]
            times = [1, SEL_INTERVAL_1[0], SEL_INTERVAL_1[1],
                     SEL_INTERVAL_2[0], SEL_INTERVAL_2[1]]
            t_stop = sim_time
            seed=numpy.random.random_integers(0,1000000.0) 

            MSN_burst.set_spike_times(id=id, rates=rates, times=times, 
                                    t_stop=t_stop, seed=seed)       
        
        for i_syn in range(len(SYNAPSE_MODELS_TESTED)):
            

            I_e=my_nest.GetDefaults(NEURON_MODELS[0])['I_e']+SNR_INJECTED_CURRENT
            SNR=MyGroup( NEURON_MODELS[0], n=n_exp, sd=True, 
               params={'I_e':I_e}, mm_dt=.1, 
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
            my_nest.Connect(GPE[sources_GPE_SNR], SNR[targets_GPE_SNR], 
                                          model=SYNAPSE_MODELS_BACKGROUND[0])
            my_nest.Connect(STN[sources_STN_SNR], SNR[targets_STN_SNR], 
                                          model=SYNAPSE_MODELS_BACKGROUND[1])
                 
        my_nest.MySimulate( sim_time )
        
        for SNR in SNR_list: 
            SNR.get_signal( 's', start=0, stop=sim_time )  
            
        pre_ref_1=str(SNR_list[0].signals['spikes'].mean_rate(SEL_ONSET-500,
                                                              SEL_ONSET)) 
        burst_1=str(SNR_list[0].signals['spikes'].mean_rate(SEL_ONSET,
                                                              SEL_ONSET+200)) 
        burst_2=str(SNR_list[0].signals['spikes'].mean_rate(SEL_ONSET+1000,
                                                              SEL_ONSET+1200))   
        s='\n'
        s=s+'Simulate example:\n'
        s = s + '%s %5s %3s \n' % ( 'Simulation time', str ( sim_time ),  '#' ) 
        s = s + '%s %5s %3s \n' % ( 'N experiments:', str ( n_exp ),  '#' )    
        s = s + '%s %5s %3s \n' % ( 'MSN base rate:',   str ( MSN_BASE_RATE),'spikes/s' )     
        s = s + '%s %5s %3s \n' % ( 'MSN burst rate:', str ( MSN_BURST_RATE ), 'spikes/s' )
        s = s + '%s %5s %3s \n' % ( 'GPe rate:', str ( GPE_BASE_RATE ), 'spikes/s' )
        s = s + '%s %5s %3s \n' % ( 'Burst time:', str ( burst_time ), 'ms' )
        s = s + '%s %5s %3s \n' % ( 'Pre sel rate Ref:', pre_ref_1[0:4], 'spikes/s' )
        s = s + '%s %5s %3s \n' % ( 'Burst 1:', burst_1[0:4], 'spikes/s' )
        s = s + '%s %5s %3s \n' % ( 'Burst 2:', burst_2[0:4], 'spikes/s' )
        header=s
        misc.text_save(header, save_header_at)
        misc.pickle_save([SNR_list, s], save_result_at)

           
    else:
        SNR_list, s = misc.pickle_load(save_result_at)
      
    return SNR_list, s

def simulate_rate_first_and_second_bursts(selection_intervals=[0.0,500.0,1000.,1500.], load=True):    
    global SNR_INJECTED_CURRENT
    global NEURON_MODELS
    global N_GPE
    global N_MSN_BURST
    global N_MSN
    global N_STN
    global MSN_BASE_RATE
    global GPE_BASE_RATE  
    global STN_BASE_RATE  
    global FILE_NAME
    global OUTPUT_PATH
    global SYNAPSE_MODELS_TESTED
    global SEL_ONSET
          
    #n_exp=20
    n_exp=200
    msn_burst_rate=20 
    n_msn_burst=N_MSN_BURST
    
    transient_stop=selection_intervals[2]-selection_intervals[1]
    save_result_at = (OUTPUT_PATH+'/simulate_rate_first_and_second_bursts_'+str(transient_stop)+'ms.pkl') 
    save_header_at = (OUTPUT_PATH+'/simulate_rate_first_and_second_bursts_'+str(transient_stop)+'ms_header') 

    burst_time = 500.
    sim_time = SEL_ONSET+selection_intervals[3]+500.   
    
    EXPERIMENTS=range(n_exp)

    model_list=models()
    my_nest.ResetKernel(threads=1)     
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    my_nest.MyLoadModels( model_list, SYNAPSE_MODELS_TESTED)   
    my_nest.MyLoadModels( model_list, SYNAPSE_MODELS_BACKGROUND)   
    if not load:
        MSN_list=[] # MSN input for each experiment
        for i_exp in EXPERIMENTS:
            MSN = MyPoissonInput( n=N_MSN+n_msn_burst, sd=True)
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
            for i_sel in range(1):
                
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
            for id in MSN[N_MSN:N_MSN + n_msn_burst]: 
                rates = [MSN_BASE_RATE, msn_burst_rate, MSN_BASE_RATE,
                         msn_burst_rate, MSN_BASE_RATE]
                
                t1=selection_intervals[0]
                t2=selection_intervals[1]
                t3=selection_intervals[2]
                t4=selection_intervals[3]
                times = [1, SEL_ONSET+t1, SEL_ONSET+t2, SEL_ONSET+t3,
                         SEL_ONSET+t4]
                t_stop = sim_time
                MSN.set_spike_times(id=id, rates=rates, times=times, 
                                    t_stop=t_stop,
                                    seed=int(numpy.random.random()*10000.0))    
            
            # Base rate GPE
            for id in GPE[:]:                 
                GPE.set_spike_times(id=id, rates=[GPE_BASE_RATE], times=[1], 
                                    t_stop=sim_time, 
                                    seed=int(numpy.random.random()*10000.0))    

            # Base rate GPE
            for id in STN[:]:                 
                STN.set_spike_times(id=id, rates=[STN_BASE_RATE], times=[1], 
                                    t_stop=sim_time, 
                                    seed=int(numpy.random.random()*10000.0))    
            
            # Connect     
            for i_syn, syn in enumerate(SYNAPSE_MODELS_TESTED):
                # i_sel goes over 0,..., n_max_sel
                for i_sel, n_sel in enumerate(range(n_msn_burst,n_msn_burst+1)):       
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
     
        t1=selection_intervals[0]
        t3=selection_intervals[2]
        
        mean_rates=[]    
        mean_rates_std=[]
        # Time until arrival of spikes in SNr
        delay=my_nest.GetDefaults(SYNAPSE_MODELS_BACKGROUND[0])['delay']
        for SNR_sel in SNR_list: 
            m_r=[]
            m_r_std=[]
            for SNR in SNR_sel:
                 
                 # Mean rate during first 200 ms                 
                 m_r.append(SNR.signals['spikes'].mean_rate(SEL_ONSET+t1+delay, SEL_ONSET+t1+200+delay))  
                 m_r.append(SNR.signals['spikes'].mean_rate(SEL_ONSET+t3+delay, SEL_ONSET+t3+200+delay))  

                 m_r_std.append(SNR.signals['spikes'].mean_rate_std(SEL_ONSET+t1+delay, SEL_ONSET+t1+200+delay))                 
                 m_r_std.append(SNR.signals['spikes'].mean_rate_std(SEL_ONSET+t3+delay, SEL_ONSET+t3+200+delay)) 

            mean_rates.append(m_r)
            mean_rates_std.append(m_r_std)
            
        mean_rates = numpy.array(mean_rates)
        mean_rates_std = numpy.array(mean_rates_std)
        
            
        s='\n'
        s=s + 'simulate_rate_first_and_second_bursts\n'
        s = s + '%s %5s %3s \n' % ( 'Simulation time', str ( sim_time ),  '#' ) 
        s = s + '%s %5s %3s \n' % ( 'N MSNs:', str ( N_MSN ),  '#' )     
        s = s + '%s %5s %3s \n' % ( 'N MSN_bursts:', str ( n_msn_burst ),  '#' ) 
        s = s + '%s %5s %3s \n' % ( 'N experiments:', str ( n_exp ),  '#' )    
        s = s + '%s %5s %3s \n' % ( 'MSN base rate:',   str ( MSN_BASE_RATE),'spikes/s' )     
        s = s + '%s %5s %3s \n' % ( 'MSN burst rate:', str ( MSN_BURST_RATE ), 'spikes/s' )
        s = s + '%s %5s %3s \n' % ( 'MSN burst time:', str ( burst_time ), 'ms' )
        s = s + '%s %5s %3s \n' % ( 'GPe base rate:', str ( GPE_BASE_RATE ), 'spikes/s' )
        s = s + '%s %5s %3s \n' % ( 'SNR injected current:', str ( SNR_INJECTED_CURRENT), 'pA' )
        for i_interval, interval in enumerate(selection_intervals):
            s = s + '%s %5s %3s \n' % ( 'Sel interval '+str(i_interval)+':', str ( selection_intervals), 'ms' )

        info_string=s
        
        header=s
        misc.text_save(header, save_header_at)
        misc.pickle_save([mean_rates, 
                          mean_rates_std, info_string], save_result_at)
        
    elif load:     
        mean_rates, mean_rates_std, info_string = misc.pickle_load(save_result_at)
    
    return mean_rates, mean_rates_std, info_string

def simulate_rate_first_and_second_bursts_full(load=True):
    global OUTPUT_PATH

    save_result_at=OUTPUT_PATH+'/simulate_rate_first_and_second_bursts_full.pkl'
    save_header_at=OUTPUT_PATH+'/simulate_rate_first_and_second_bursts_full_header'
   
    # Range
    transient_stops=numpy.arange(100,3200,500)
    
    #hzs=[8,20]
    if not load:
        data={}

        data['rates']=[]
        
        for stop in transient_stops:
                                  
            mean_rates, mean_rates_std, info_string=simulate_rate_first_and_second_bursts(selection_intervals=[0.0,500.0,500.+stop,1000.+stop], load=False)    
            data['rates'].append(mean_rates[0])

        
        s='\n'
        s=s + 'simulate_rate_first_and_second_bursts_full\n'
        s = s + ' %s %5s %s \n' % ( 'Transient stops', str ( transient_stops[0])+'-'+str(transient_stops[-1]),  'ms' )     

        header=s
        misc.text_save(header, save_header_at)
        misc.pickle_save([data,s], save_result_at)
        info_string=s
    elif load:
        data, info_string=misc.pickle_load(save_result_at)
      
    
    data['rates']=numpy.array(data['rates'])
    
    return transient_stops, data, info_string
print 'Simulation'


transient_stops, data, info_string=simulate_rate_first_and_second_bursts_full(load=False)

# SIMULATION
infoString=''

SNR, s =simulate_example(load=True)

infoString=infoString+s   

# DISPLAY

plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0, fontsize=12)
font_size_text = 8
fig = pylab.figure( facecolor = 'w' )

ax_list = []
ax_list.append( MyAxes(fig, [ .075, .37, .135, .26 ] ) )    # text box
#ax_list.append( MyAxes(fig, [ .26,  .6,  .165*2.312, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .53, .6,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .6,  .165, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .26,  .1,  .165, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .53, .1,  .165, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .8,   .1,  .165, .34 ] ) )    # 

import pprint
pprint.pprint(pylab.rcParams)

# Text
ax=ax_list[0]
plot_text(ax, infoString)

# Example msn
ax=ax_list[1]
plot_example_snr(ax, SNR)

# Example snr
ax=ax_list[2]
plot_rate_first_and_second_bursts_full(ax, transient_stops, data)


# plot_selection_vs_neurons
#ax=ax_list[3]
#plot_plastic_1_vs_2_rate(ax, firstMeanRates2, secondMeanRates2, shiftSecond2)
#ax.legend(numpoints=1, loc=[1.5,0.0])

# IF
#ax=ax_list[4]


pylab.show()

name = sys.argv[0].split('/')[-1].split('.')[0]
fig.savefig( picture_dir + '/' + name  + '.svg', format = 'svg')
fig.savefig( picture_dir + '/' + name  + '.pdf', format = 'pdf')