#! Imports
import math
import numpy
import pylab
import os
import sys
import time
import scipy.interpolate as sp

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

# Imports dependent on adding code model and nest_toolbox path
from model_params import models                                   
from src import misc, my_nest, my_signals, plot_settings
from src.my_population import MyGroup, MyPoissonInput 
from src.my_axes import MyAxes 

### GLOBALS ###

# Default params 
DP={}

# Paths and naming for saving data and picutes
DP['FILE_NAME'] = sys.argv[0].split('/')[-1].split('.')[0]
model_name=os.getcwd().split('/')[-2]
DP['PICTURE_PATH']=('/'.join(os.getcwd().split('/')[0:-3]) + '/pictures/'+
                    model_name)    
DP['OUTPUT_PATH']  = (os.getcwd()+'/output/' + 
                      sys.argv[0].split('/')[-1].split('.')[0])

# Models
DP['NEURON_MODELS']=['SNR_izh']
DP['SYNAPSE_MODELS_TESTED']=['MSN_SNR_gaba_s_min', 'MSN_SNR_gaba_s_max', 
                             'MSN_SNR_gaba_p1']
DP['SYNAPSE_MODELS_BACKGROUND']=['GPE_SNR_gaba_p']

# Neuron numbers
DP['N_MSN'] = 500
DP['N_MSN_BURST'] = 11    
DP['N_GPE'] = 22

# Rates
DP['MSN_BASE_RATE']=0.1 # Base rate
DP['GPE_BASE_RATE']=25 

# Misc
DP['LOWER_THR']=5.  # spikes/s
DP['UPPER_THR']=10.  # spikes/s
DP['SNR_INJECTED_CURRENT']=[530.0, 530.0, 530.0] 
DP['SELECTION_THR']=5.  # spikes/s
DP['SEL_ONSET']=2000
DP['SEL_TIME']=500
DP['ADJUST_XDATA_MS']=1500.
 

HEADER_SIMULATION_SETUP=( '**** BEGINNING GENERAL SCRIPT SETUP ****\n'+
                          str(DP)+                 
                          '**** END GENERAL SCRIPT SETUP ****\n')
### END GLOBALS ###

def plot_example_SNR(ax, SNR_list):
    time_bin=20
    
    colors=['b','g','m']   
    labels=[r'$\delta_{weak}^{MSN}$' , r'$\delta_{strong}^{MSN}$',  
            r'$\delta_{fac}^{MSN}$']
    coords=[[0.4, 0.45], [ 0.05, 0.12], [0.1, 0.6]]
    
    for color, SNR in zip(colors, SNR_list):
        signal=SNR.signals['spikes']
        signal.my_firing_rate(bin=time_bin, display=ax,
                          kwargs={'color':color})
        
    lines = ax.lines
    for line in lines:
        line.set_xdata(line.get_xdata()-DP['ADJUST_XDATA_MS'])    
        misc.slice_line(line, xlim=[0,1490])
   
    ax.plot([0, 1490],[DP['SELECTION_THR']]*len([0, 1490]),**{'color':[0.6,0.6,0.6],
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
    ax.text( 0.8, 0.11,'Thr' , transform=ax.transAxes, **{ 'color' : [0.6,0.6,0.6]})   
        
    ax.set_title('bin=%i'%(time_bin),**{'fontsize':12})
    ax.set_ylabel('Firing rate SNr (spikes/s)') 
    ax.set_xlabel('Time (ms)')
    ax.my_set_no_ticks( yticks=6, xticks = 5 ) 

    ax.set_xlim(misc.adjust_limit([0,1500]))
    ax.set_ylim(misc.adjust_limit([0,50]))

    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['text.fontsize'], 
                 **{'color': color})
        
def plot_signal_rates(ax, hzs, rates):

    colors=['b','m']   
    labels=[r'$\delta_{weak}^{MSN}$', r'$\delta_{fac}^{MSN}$']  
    coords=[[0.53, 0.49], [ 0.05, 0.21]]   
    
    syn=DP['SYNAPSE_MODELS_TESTED']
    
    for id, label, color in zip([0,2],labels,colors):
        ax.plot(hzs,rates[:,id],**{'color':color}) 
        
    ax.set_xlabel('Firing rate MSN burst (spikes/s)') 
    ax.set_ylabel('Firing rate SNr (spikes/s)')
    ax.my_set_no_ticks( yticks=8, xticks = 6 ) 

    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['text.fontsize'], 
                 **{'color': color})
    ax.set_yticks(hzs, hzs)
    
    ax.plot([min(hzs), max(hzs)],[DP['LOWER_THR']]*len([0, max(hzs)]),**{'color':[0.6,0.6,0.6],
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
    #ax.text( 0.8, 0.11,'Thr' , transform=ax.transAxes, **{ 'color' : [0.6,0.6,0.6]})   
    
    '''
    ax.semilogx([min(hzs), max(hzs)],[DP['UPPER_THR']]*len([0, max(hzs)]),**{'color':[0.6,0.6,0.6],
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
    #ax.text( 0.8, 0.11,'Thr' , transform=ax.transAxes, **{ 'color' : [0.6,0.6,0.6]})   
    '''
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0, 48])
        
    ax.set_xlim(misc.adjust_limit([0,48]))
    ax.set_ylim(misc.adjust_limit([0,30]))  

def plot_signal_performance(ax,  r_D, r_S, r_S_):

    colors=['b','m']   
    labels=[r'$\delta_{strong}^{MSN}$' , 
            r'$\delta_{fac}^{MSN}$']  
    coords=[[0.75, 0.51], [ 0.75, 0.36], [0.75, 0.21]]   
        
    
    syn=DP['SYNAPSE_MODELS_TESTED']
        
    #ax.set_xlabel('Firing rate MSN burst (spikes/s)') 
    ax.set_ylabel('Proportion (%)')
    
    ind=numpy.arange(3)
    width = 0.35       # the width of the bars
    
    
    rects1 = ax.bar(ind, [ r_D[1], r_S[1],r_S_[1]], width, color=colors[0])
    rects2 = ax.bar(ind+width, [ r_D[2], r_S[2],r_S_[2]], width, color=colors[1])
     
    ax.set_xticks(ind+width)
    ax.set_xticklabels( ['D','S','S*'] )
    
    
       
    #ax.text( 0.05, 0.85, labels[0] , backgroundcolor='w',
    #         transform=ax.transAxes, **{ 'color' : colors[0] })  
    #ax.text( 0.05, 0.75, labels[1] , backgroundcolor='w',
    #         transform=ax.transAxes, **{ 'color' : colors[1] })  

    ax.legend( (rects1[0], rects2[0]), (labels[0], labels[1]), loc='best' )
            
    #ax.set_xlim(misc.adjust_limit([0,3]))
    #ax.set_ylim(misc.adjust_limit([0,70]))     

def plot_text(ax, info_string=''):
    
    my_nest.ResetKernel()
    model_list=models()
    my_nest.MyLoadModels( model_list, DP['NEURON_MODELS'] )
    
    SNR = MyGroup( DP['NEURON_MODELS'][0], 1, mm_dt = 0.1)
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

def simulate_get_rates(msn_burst_rate=20, load=True, len_ms=500.):

    n_exp = 50
    sim_time = DP['SEL_TIME']+DP['SEL_ONSET']+1000.
   
    experiments=range(n_exp)
      
    model_list=models()
    my_nest.ResetKernel()       
    my_nest.MyLoadModels( model_list, DP['NEURON_MODELS'] )
    my_nest.MyLoadModels( model_list, DP['SYNAPSE_MODELS_TESTED'])       
    my_nest.MyLoadModels( model_list, DP['SYNAPSE_MODELS_BACKGROUND'])      
 
    MSN_list=[] # MSN input for each experiment
    for i_exp in experiments:
        MSN = MyPoissonInput( n=DP['N_MSN'], sd=True)
        MSN_list.append(MSN)
    
    GPE_list=[] # GPE input for each experiment
    for i_exp in experiments:
        GPE = MyPoissonInput( n=DP['N_GPE'], sd=True)
        GPE_list.append(GPE)
    
    SNR_list=[] # SNR groups for each synapse
    for i_syn, syn in enumerate(DP['SYNAPSE_MODELS_TESTED']):
        
        I_e=my_nest.GetDefaults(DP['NEURON_MODELS'][0])['I_e']+DP['SNR_INJECTED_CURRENT'][i_syn]          
        SNR = MyGroup( DP['NEURON_MODELS'][0], n=n_exp, sd=True, params={'I_e':I_e})
        SNR_list.append(SNR)
        
        
    for i_exp in experiments:    
        MSN = MSN_list[i_exp]
        GPE = GPE_list[i_exp]
        
        # Set spike times
        # Base rate
        for id in MSN[1:DP['N_MSN']]:                 
            MSN.set_spike_times(id=id, rates=[DP['MSN_BASE_RATE']], times=[1], 
                                t_stop=sim_time,
                                seed=int(numpy.random.random()*10000.0))           
  
  
        # Set spike times
        # Base rate
        for id in GPE[:]:                 
            GPE.set_spike_times(id=id, rates=[DP['GPE_BASE_RATE']], times=[1], 
                                t_stop=sim_time, 
                                seed=int(numpy.random.random()*10000.0))           
  
        # Selection        
        for id in MSN[DP['N_MSN']-DP['N_MSN_BURST']-50:DP['N_MSN']-50]: 
            rates = [DP['MSN_BASE_RATE'], msn_burst_rate, DP['MSN_BASE_RATE']]
            times = [1, DP['SEL_ONSET'], DP['SEL_TIME'] + DP['SEL_ONSET']]
            t_stop = sim_time
            MSN.set_spike_times(id=id, rates=rates, times=times, 
                                t_stop=t_stop,
                                seed=int(numpy.random.random()*10000.0))     
    
        # Connect  
        for i_syn, syn in enumerate(DP['SYNAPSE_MODELS_TESTED']):       
                target=SNR_list[i_syn][i_exp]
                my_nest.ConvergentConnect(MSN[:], [target], model=syn)
                my_nest.ConvergentConnect(GPE[:], [target], 
                                      model=DP['SYNAPSE_MODELS_BACKGROUND'][0])
                  
    my_nest.MySimulate( sim_time )

    for SNR in SNR_list: 
        SNR.get_signal( 's' ) 
        
    rate_ref_1=str(SNR_list[0].signals['spikes'].mean_rate(DP['SEL_ONSET'],
                                                          DP['SEL_ONSET']+len_ms)) 
    rate_ref_2=str(SNR_list[1].signals['spikes'].mean_rate(DP['SEL_ONSET'],
                                                          DP['SEL_ONSET']+len_ms)) 
    rate_dyn=str(SNR_list[2].signals['spikes'].mean_rate(DP['SEL_ONSET'],
                                                        DP['SEL_ONSET']+len_ms))   
      
    return [rate_ref_1,rate_ref_2, rate_dyn]

def simulate_signal_rates(load=True, hzs=[1,2]):
    
    # Path were raw data is saved. For example the spike trains.
    save_result_at=DP['OUTPUT_PATH']+'/simulate_signal_rates.pkl'
    save_header_at=DP['OUTPUT_PATH']+'/simulate_signal_rates_header' 
    
    rates=[]
    if not load:
        for hz in hzs:
            rates.append(simulate_get_rates(msn_burst_rate=hz,load=load))
        
        rates=numpy.array(rates)
        
        header=HEADER_SIMULATION_SETUP
        misc.text_save(header, save_header_at)
        misc.pickle_save(rates, save_result_at)
        
    else:
        rates=misc.pickle_load(save_result_at)
    
    return rates
        
def calulate_signal_preformance(hzs, rates):
    max_hz=max(hzs)
    r_D=[]
    r_S=[]
    r_S_=[]
    for i in [0,1,2]:
        r=rates[:,i]
        fl = sp.InterpolatedUnivariateSpline(hzs,r)
        #fl = sp.interp1d(hzs, r,kind='linear')
        x=numpy.arange(hzs[0],hzs[-1], 0.1)
        y=fl(x)
        x_cut=x[(y>DP['LOWER_THR']) * (y<DP['UPPER_THR'])]
        
        r_D.append((x_cut[-1]-x_cut[0])/(x[-1]-x[0])*100)
        r_S.append((x[-1]-x_cut[-1])/(x[-1]-x[0])*100)
        r_S_.append((x_cut[0]-x[0])/(x[-1]-x[0])*100)
    return r_D, r_S, r_S_
    
print 'Simulation'

# SIMULATION
info_string=''

sel_interval_4 = [0.0, 200.0]

# If simulating with MPI
# 1. Run with MPI load_pickle and load2 False
# 2. When finished run laod2=True and load_pickle=False to recover MPI saved data
# data and save it as pickled file.
# 3. Run with load_pickle and load2 both True and plot data


stop = time.time()    
sec = stop - start
m = sec // 60
sec = sec - m*60
print 'Rank %i simulation time: %i minutes, %i seconds' % ( my_nest.Rank(), m, sec )


#hzs=[0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 25.6, 51.2]

hzs=range(1,48,5)
signal_rates=simulate_signal_rates(load=True, hzs=hzs)


# DISPLAY
plot_settings.set_mode(pylab, mode='by_fontsize', 
                       w = 1100.0, h = 450.0, fontsize=12)
font_size_text = 8
fig = pylab.figure( facecolor = 'w' )

ax_list = []
ax_list.append( MyAxes(fig, [ .075, .37, .135, .26 ] ) )    # text box
ax_list.append( MyAxes(fig, [ .26,  .6,  .165, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .53,  .6,  .165, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .8,   .6,  .165, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .26,  .1,  .165, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .53,  .1,  .165, .34 ] ) )    # 
# ax_list.append( MyAxes(fig, [ .8,   .1,  .165, .34 ] ) )    # 

# Text
ax=ax_list[0]
plot_text(ax, info_string)
   
ax=ax_list[1]
plot_signal_rates(ax, hzs, signal_rates)


pylab.show()

fig.savefig( DP['PICTURE_PATH'] + '/' + DP['FILE_NAME']  + '.svg', format = 'svg')
fig.savefig( DP['PICTURE_PATH'] + '/' + DP['FILE_NAME']  + '.pdf', format = 'pdf')