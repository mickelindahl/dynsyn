
import numpy
import pylab
import os
import sys

numpy.random.seed(1) # set random seed

# First add parent directory to python path
model_dir=   '/'.join(os.getcwd().split('/')[0:-1])        
code_dir=  '/'.join(os.getcwd().split('/')[0:-2]) 

# Add model, code and current directories to python path
sys.path.append(os.getcwd())                                                                  
sys.path.append(model_dir) 
sys.path.append(code_dir+'/nest_toolbox') 

# Then import model and network
from model_params import models
from src import misc, my_nest, my_signals, plot_settings
from src.my_population import MyGroup, MyPoissonInput  
from src.my_axes import MyAxes 


### GLOBALS ###

# Paths and naming for saving data and picutes
FILE_NAME = sys.argv[0].split('/')[-1].split('.')[0]
model_name=os.getcwd().split('/')[-2]
PICTURE_PATH='/'.join(os.getcwd().split('/')[0:-3]) + '/pictures/'+model_name    
OUTPUT_PATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]

# Models
NEURON_MODELS=['SNR_aeif']
SYNAPSE_MODELS_TESTED=['MSN_SNR_gaba_s_min', 'MSN_SNR_gaba_s_max', 'MSN_SNR_gaba_p1']
SYNAPSE_MODELS_BACKGROUND=['GPE_SNR_gaba_p', 'STN_SNR_ampa_s']

# Neuron numbers
N_MSN=500
N_GPE=30
N_STN=100

# Neuron rate default
MSN_BURST_RATE=20.
GPE_BASE_RATE=25.
STN_BASE_RATE=10.
# Misc
SNR_INJECTED_CURRENT=400.0
SYN_EVENTS=400
N_MAX_BURSTING=20
SELECTION_THR=5.  # spikes/s
SEL_ONSET=2000
SEL_OFFSET=2500

HEADER_SIMULATION_SETUP=( '**** BEGINNING GENERAL SCRIPT SETUP ****\n'+
                          'FILE_NAME:'+str(FILE_NAME)+'\n'+                         
                          'PICTURE_PATH:'+str(PICTURE_PATH)+'\n'+  
                          'OUTPUT_PATH:'+str(OUTPUT_PATH)+'\n\n'+     
                          'N_MSN:'+str(N_MSN)+'\n'+
                          'N_STN:'+str(N_STN)+'\n'+
                          'N_GPE:'+str(N_GPE)+'\n\n'+
                          'GPE_BASE_RATE:'+str(GPE_BASE_RATE)+'\n'+
                          'STN_BASE_RATE:'+str(STN_BASE_RATE)+'\n\n'+
                          'SNR_INJECTED_CURRENT:'+str(SNR_INJECTED_CURRENT)+'\n'+
                          'SEL_ONSET:'+str(SEL_ONSET)+'\n' 
                          'SEL_OFFSET:'+str(SEL_OFFSET)+'\n' 
                          'SYN_EVENTS:'+str(SYN_EVENTS)+'\n\n'+
                          'SYNAPSE_MODELS_TESTED:'+str(SYNAPSE_MODELS_TESTED)+'\n'+
                          'SYNAPSE_MODELS_BACKGROUND:'+str(SYNAPSE_MODELS_BACKGROUND)+'\n'+
                          '**** END GENERAL SCRIPT SETUP ****\n')

### END GLOBALS ###


def plot_selection_vs_neurons(ax, MSNmeanRates, SNRmeanRates):
    colors=['b','g','m']   
    labels=[r'$\delta_{weak}^{MSN}$' , r'$\delta_{strong}^{MSN}$',  
            r'$\delta_{fac}^{MSN}$']
    coords=[[0.4, 0.42], [ 0.02, 0.001], [0.15, 0.78]]
    
    
    for i, color in enumerate(colors):
        ax.plot(MSNmeanRates,SNRmeanRates[i,:],**{'color':color})  

   
    ax.plot([0,2.5],[SELECTION_THR]*len([0,2]),**{'color':[0.6,0.6,0.6],
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
    ax.text( 0.8, 0.20,'Thr' , transform=ax.transAxes, **{ 'color' : [0.6,0.6,0.6]})
    
    ax.set_ylabel('Firing rate SNr (spikes/s)') 
    ax.set_xlabel('Firing rate MSN (spikes/s)')
    ax.my_set_no_ticks( yticks=8, xticks = 6 )
    
    ax.set_xlim(misc.adjust_limit([0, 2.5]))
    ax.set_ylim(misc.adjust_limit([0,28]))
    
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,2.5])   
                                  
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['text.fontsize'], 
                 **{'color': color})

def plot_MSN_vs_SNR_const_syn_events(ax, nMSN_range, SNRmeanRates):

    colors=['b','g','m' ] 
    labels=[r'$\delta_{weak}^{MSN}$' , r'$\delta_{strong}^{MSN}$', 
            r'$\delta_{fac}^{MSN}$']
    coords=[[0.1, 0.68], [ 0.05, 0.12], [ 0.35, 0.45]]
    
    for i, color in enumerate(colors):
        ax.plot(nMSN_range,SNRmeanRates[i,:],**{'color':color})  
    
    ax.set_ylabel('Firing rate SNr (spikes/s)') 
    ax.set_xlabel('Bursting MSN (#)')
    ax.my_set_no_ticks( yticks=7, xticks=6 )
    
    
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,N_MAX_BURSTING])   
    
    ax.plot([0,N_MAX_BURSTING],[SELECTION_THR]*len([0,2]),**{'color':[0.6,0.6,0.6],
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
        
    ax.text( 0.1, 0.33,'Thr' , transform=ax.transAxes, **{ 'color' : [0.6,0.6,0.6]})
    
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['text.fontsize'], 
                 **{'color': color})

    ax.set_xlim(misc.adjust_limit([0, N_MAX_BURSTING]))      
    ax.set_ylim(misc.adjust_limit([0,18]))

def plot_text(ax, infoString=''):
    
    my_nest.ResetKernel()
    model_list, model_dict=models()
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

def simulate_MSN_vs_SNR_rate(load=True):
    global SNR_INJECTED_CURRENT
    global N_MSN
    global N_GPE
    global N_STN
    global GPE_BASE_RATE
    
    # Path were raw data is saved. For example the spike trains.
    save_result_at=OUTPUT_PATH+'/simulate_MSN_vs_SNR_rate.pkl'
    save_header_at=OUTPUT_PATH+'/simulate_MSN_vs_SNR_rate_header'
    
    MSNmeanRates=numpy.arange(0.1,3.1,0.1)
    SNRmeanRates=[]

    sim_time=100000.
   
    if not load:
        for r in MSNmeanRates:
            my_nest.ResetKernel(threads=3)
            model_list, model_dict=models()
            my_nest.MyLoadModels( model_list, NEURON_MODELS )
            my_nest.MyLoadModels( model_list, SYNAPSE_MODELS_TESTED )
            my_nest.MyLoadModels( model_list, SYNAPSE_MODELS_BACKGROUND )
            
            MSN = MyPoissonInput( n=N_MSN )           
            GPE = MyPoissonInput( n=N_GPE )
            STN = MyPoissonInput( n=N_STN )
            
            I_e=my_nest.GetDefaults(NEURON_MODELS[0])['I_e']+SNR_INJECTED_CURRENT             
            SNR = MyGroup( NEURON_MODELS[0], n=len(SYNAPSE_MODELS_TESTED), params={'I_e':I_e},
                           sd=True)
   

            for id in MSN[:]:    
                MSN.set_spike_times(id=id, rates=numpy.array([r]), 
                                        times=numpy.array([1]), 
                                        t_stop=sim_time,
                                        seed=int(numpy.random.random()*10000.0)) 
      
                  
            # Base rate GPE
            for id in GPE[:]:                 
                    GPE.set_spike_times(id=id, rates=[GPE_BASE_RATE], times=[1], 
                                    t_stop=sim_time, 
                                    seed=int(numpy.random.random()*10000.0))
                    
            # Base rate STN
            for id in STN[:]:                 
                    STN.set_spike_times(id=id, rates=[STN_BASE_RATE], times=[1], 
                                    t_stop=sim_time, 
                                    seed=int(numpy.random.random()*10000.0))
                    
            for i, syn in enumerate(SYNAPSE_MODELS_TESTED):
                    my_nest.ConvergentConnect(MSN[:],[SNR[i]], model=syn)
                    my_nest.ConvergentConnect(GPE[:],[SNR[i]], model=SYNAPSE_MODELS_BACKGROUND[0])
                    my_nest.ConvergentConnect(STN[:],[SNR[i]], model=SYNAPSE_MODELS_BACKGROUND[1])
            
            my_nest.MySimulate( sim_time )

            SNR.get_signal( 's') # retrieve signal
                  
            SNRmeanRates.append(SNR.signals['spikes'].mean_rates(1000.0,sim_time))   
        
        SNRmeanRates=numpy.array(SNRmeanRates).transpose()
        MSNmeanRates=numpy.array(MSNmeanRates)
        
        rateAtThr=''
        for SNRr in SNRmeanRates:
            tmp=str(MSNmeanRates[SNRr>=SELECTION_THR][-1])
            rateAtThr+=' '+tmp[0:4]
        
            
        s='\n'
        s =s + 'simulate_MSN_vs_SNR_rate:\n'   
        s = s + ' %s %5s %3s \n' % ( 'N MSNs:', str ( N_MSN ),  '#' )  
        s = s + ' \n%s \n%5s %3s \n' % ( 'MSN rates:', str ( MSNmeanRates[0] ) + '-'+ 
                                         str ( MSNmeanRates[-1] ),  'spikes/s' ) 
        s = s + ' %s %5s %3s \n' % ( 'N GPes:', str ( N_GPE ),  '#' )   
        s = s + ' %s %5s %3s \n' % ( 'Threshold SNr:', str ( SELECTION_THR ),  'spikes/s' )
        s = s + ' \n%s \n%5s %3s \n' % ( 'MSN rate right before threshold SNr:', str ( rateAtThr ),  'spikes/s' )   
        s = s + ' \n%s %5s %3s \n' % ( 'Simulation time:', str ( sim_time), 'ms' )
        s = s + ' %s %5s %3s \n' % ( 'Injected current:', str ( SNR_INJECTED_CURRENT ), 'pA' )
        infoString=s

        header=HEADER_SIMULATION_SETUP+s
        misc.text_save(header, save_header_at)    
        misc.pickle_save([MSNmeanRates, SNRmeanRates, infoString], 
                                save_result_at)
    elif load:
        MSNmeanRates, SNRmeanRates, infoString = misc.pickle_load(save_result_at)       
        
    return MSNmeanRates, SNRmeanRates, infoString
 
def simulate_MSN_vs_SNR_const_syn_events(load=True):
    global SNR_INJECTED_CURRENT
    global N_MSN
    global N_GPE
    global MSN_BURST_RATE
    global GPE_BASE_RATE
    
    # Path were raw data is saved. For example the spike trains.
    save_result_at=OUTPUT_PATH+'/simulate_MSN_vs_SNR_const_syn_events.pkl'
    save_header_at=OUTPUT_PATH+'/simulate_MSN_vs_SNR_const_syn_events_header'
    
    # REMARK can not be more than rate=const_syn_events/burst_rate
    n_MSN_bursting=numpy.arange(0,N_MAX_BURSTING+1) 

    n_exp=200
    #n_exp=20

    # Solve (500-n)*x + 20*n=600, where 500 is total number of MSNs, 20 is burst
    # activation, x is MSN mean rate and n is number of bursters. 
    # Then x=(600-20*n)/(500-n)
    MSNmeanRates=(SYN_EVENTS-MSN_BURST_RATE*n_MSN_bursting)/(N_MSN-n_MSN_bursting)
    
    SNRmeanRates=[]

    sim_time=3000.

    if not load:   
        for r, n_MSN_b in zip(MSNmeanRates, n_MSN_bursting):
            my_nest.ResetKernel(threads=4)
            model_list, model_dict=models()
            my_nest.MyLoadModels( model_list, NEURON_MODELS )
            my_nest.MyLoadModels( model_list, SYNAPSE_MODELS_TESTED )
            my_nest.MyLoadModels( model_list, SYNAPSE_MODELS_BACKGROUND )
            
            MSN=[]
            SNR=[]
            GPE=[]
            STN=[]
            
            for i in range(n_exp):
                MSN.append(MyPoissonInput( n=N_MSN, sd=True))
                GPE.append(MyPoissonInput( n=N_GPE, sd=True))
                STN.append(MyPoissonInput( n=N_STN, sd=True))
                
                I_e=my_nest.GetDefaults(NEURON_MODELS[0])['I_e']+SNR_INJECTED_CURRENT                   
                SNR.append(MyGroup( NEURON_MODELS[0], 
                                    n=len(SYNAPSE_MODELS_TESTED), 
                                    params={'I_e':I_e}, sd=True))
                 
            for i_exp in range(n_exp):
                for id in MSN[i_exp][:N_MSN-n_MSN_b]:    
                    MSN[i_exp].set_spike_times(id=id, rates=numpy.array([r]), 
                                        times=numpy.array([1]), 
                                        t_stop=sim_time,
                                        seed=int(numpy.random.random()*10000.0)) 
                 
                for id in MSN[i_exp][N_MSN-n_MSN_b:]:    
                    MSN[i_exp].set_spike_times(id=id, rates=numpy.array([r,MSN_BURST_RATE,r]), 
                                        times=numpy.array([1,SEL_ONSET,SEL_OFFSET]), 
                                        t_stop=sim_time,
                                        seed=int(numpy.random.random()*10000.0))     
                  
                # Base rate GPE
                for id in GPE[i_exp][:]:                 
                    GPE[i_exp].set_spike_times(id=id, rates=[GPE_BASE_RATE], times=[1], 
                                    t_stop=sim_time, 
                                    seed=int(numpy.random.random()*10000.0))
                # Base rate STN
                for id in STN[i_exp][:]:                 
                    STN[i_exp].set_spike_times(id=id, rates=[STN_BASE_RATE], times=[1], 
                                    t_stop=sim_time, 
                                    seed=int(numpy.random.random()*10000.0))
                    
                for j, syn in enumerate(SYNAPSE_MODELS_TESTED):
                    my_nest.ConvergentConnect(MSN[i_exp][:],[SNR[i_exp][j]], model=syn)
                    my_nest.ConvergentConnect(GPE[i_exp][:],[SNR[i_exp][j]], model=SYNAPSE_MODELS_BACKGROUND[0])                
                    my_nest.ConvergentConnect(STN[i_exp][:],[SNR[i_exp][j]], model=SYNAPSE_MODELS_BACKGROUND[1])
                
            my_nest.MySimulate( sim_time )
            
            delay=my_nest.GetDefaults(SYNAPSE_MODELS_BACKGROUND[0])['delay']
            SNRmeanRates_tmp=[]    
            for i in range(n_exp): 
                SNR[i].get_signal( 's') # retrieve signal
      
                SNRmeanRates_tmp.append(SNR[i].signals['spikes'].mean_rates(SEL_ONSET+delay,SEL_OFFSET+delay))   
            
            SNRmeanRates.append(numpy.mean(SNRmeanRates_tmp,axis=0))

        SNRmeanRates=numpy.array(SNRmeanRates).transpose()

        s='\n'
        s =s + 'simulate_MSN_vs_SNR_const_syn_events:\n'   
        s = s + '%s %5s %3s \n' % ( 'Syn events:', str ( SYN_EVENTS ),  '#' )  
        s = s + '%s %5s %3s \n' % ( 'n_exp:', str ( n_exp ),  '#' )  
        infoString=s
        
        header=HEADER_SIMULATION_SETUP+s
        misc.text_save(header, save_header_at)    
        misc.pickle_save([SNRmeanRates, infoString], save_result_at)
    
    elif load:
        SNRmeanRates, infoString = misc.pickle_load(save_result_at)       
    
    
    return n_MSN_bursting, MSNmeanRates, SNRmeanRates, infoString    
print 'Simulation'

# SIMULATION
infoString=''

# simulate_MSN_vs_SNR_rate
MSNmeanRates, SNRmeanRates, s = simulate_MSN_vs_SNR_rate(load=True)
infoString=infoString+s    


n_MSN_bursting, MSNmeanRates_cse, SNRmeanRates_cse, s = simulate_MSN_vs_SNR_const_syn_events(load=True)
infoString=infoString+s   

# DISPLAY
plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0, fontsize=12)
font_size_text = 8
fig = pylab.figure( facecolor = 'w' )

ax_list = []
ax_list.append( MyAxes(fig, [ .075, .37, .135, .26 ] ) )    # text box
ax_list.append( MyAxes(fig, [ .26,  .6,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .53, .6,  .165, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .8,   .6,  .165, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .26,  .1,  .165, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .53, .1,  .165, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .8,   .1,  .165, .34 ] ) )    # 

# Text
ax=ax_list[0]
plot_text(ax, infoString)

# plot_selection_vs_neurons
ax=ax_list[1]
plot_selection_vs_neurons(ax, MSNmeanRates, SNRmeanRates)

data={}
data['MSN_mean_rates']=MSNmeanRates
data['SNR_mean_rates']=SNRmeanRates
misc.pickle_save(data, os.getcwd()+'/output/mean_rates_filtering'+
            NEURON_MODELS[0])

# IF
ax=ax_list[2]
plot_MSN_vs_SNR_const_syn_events(ax, n_MSN_bursting, SNRmeanRates_cse)

pylab.show()

fig.savefig( PICTURE_PATH + '/' + FILE_NAME  + '.svg', format = 'svg')
fig.savefig( PICTURE_PATH + '/' + FILE_NAME  + '.pdf', format = 'pdf')

