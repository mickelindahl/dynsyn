
import numpy
import pylab
import os
import sys


if len(sys.argv) != 1: mpiRun = True
else:                  mpiRun = False

numpy.random.seed(1) # set random seed

sys.path.append(os.getcwd())    # Add current directory to python path                                                
current_path=os.getcwd()

# First add parent directory to python path
model_dir=   '/'.join(os.getcwd().split('/')[0:-1])        
code_dir=  '/'.join(os.getcwd().split('/')[0:-2]) 

model_name=os.getcwd().split('/')[-2]
picture_dir='/'.join(os.getcwd().split('/')[0:-3]) + '/pictures/'+model_name      
                
sys.path.append(model_dir) 
sys.path.append(code_dir+'/nest_toolbox') 
SPATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]

# Then import model and network
from model_params import models
from src import misc, my_nest, my_signals, plot_settings
from src.my_population import MyGroup, MyPoissonInput  
from src.my_axes import MyAxes 


### GLOBALS ###
SELECTION_THR=5.  # spikes/s

N_GPE =30
N_MSN=500
N_STN=100

MSN_BASE_RATE=0.1
GPE_BASE_RATE=25
STN_BASE_RATE=10

# Fix number of synaptic events at the backround rate times number
# of presynaptic GPe neurons
CONSTANT_SYN_EVENTS=N_GPE*GPE_BASE_RATE 
SNR_INJECTED_CURRENT= 400.0   

NEURON_MODELS=['SNR_izh']
SYNAPSE_MODELS_TESTED=['GPE_SNR_gaba_s_ref', 'GPE_SNR_gaba_p']
SYNAPSE_MODELS_BACKGROUND=['MSN_SNR_gaba_p1', 'STN_SNR_ampa_s']
### END GLOBALS ###

def plot_selection_vs_neurons(ax, GPEmeanRates, SNRmeanRates):
    colors=['r','c' ] 
    labels=[r'$\delta_{ref}^{GPe}$' , r'$\delta_{dep}^{GPe}$']
    coords=[[0.15, 0.6], [ 0.4, 0.17]]
    
    for i, color in enumerate(colors):
        ax.plot(GPEmeanRates,SNRmeanRates[i,:],**{'color':color})  
    
    
    vec=[0,100]
    ax.plot(vec,[SELECTION_THR]*len(vec),**{'color':[0.6,0.6,0.6],
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
    ax.text( 0.05, 0.08,'Thr' , transform=ax.transAxes, **{ 'color' : [0.6,0.6,0.6]})   
    ax.set_ylabel('Firing rate SNr (spikes/s)') 
    ax.set_xlabel('Firing rate GPe (spikes/s)')
    ax.my_set_no_ticks( yticks=8, xticks=8 )
    
    
        
    SNRmr=SNRmeanRates[1,SNRmeanRates[1,:]<25]
    GPEmr=GPEmeanRates[SNRmeanRates[1,:]<25]
    GPEmr=GPEmr[SNRmr>=6.5]
    SNRmr=SNRmr[SNRmr>=6.5]+7
    
    #x_arrow=GPEmr[-1]
    #y_arrow=SNRmr[-1]
    
    #ax.arrow(x_arrow, y_arrow, 0, -2,
    #width=1, head_width=3, head_starts_at_zero=True,
    #head_length=2,**{ 'color':'k'})
    

    ax.set_xlim(misc.adjust_limit([0,100]))
    ax.set_ylim(misc.adjust_limit([0,110]))
    
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,100]) 
    
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['text.fontsize'], 
                 **{'color': color})
              
def plot_selection_vs_neurons_percent(ax, GPEmeanRates, SNRmeanRates):
    colors=['k' ] 
    labels=[r'$\delta_{dep}^{GPe}$']
    coords=[[ 0.1, 0.78]]
    max_SNR=max(SNRmeanRates[1,:])
    
    SNRmeanRates[1,:]=SNRmeanRates[1,:]/max_SNR*100.0
    
    
    print SNRmeanRates[1,:]
    print GPEmeanRates[:]
    
    
    SNRmr=SNRmeanRates[1,SNRmeanRates[1,:]>=50]
    GPEmr=GPEmeanRates[SNRmeanRates[1,:]>=50]
    lines=ax.plot(GPEmr,SNRmr,**{'color':[1,0.5,0.5], 'linewidth':5})  
    ax.text( GPEmr[-1]+4, SNRmr[-1], str(int(GPEmr[-1]))+' spikes/s',
                 fontsize=pylab.rcParams['text.fontsize']-6, 
                 **{'color': 'k'})
    
    SNRmr=SNRmeanRates[1,SNRmeanRates[1,:]<50]
    GPEmr=GPEmeanRates[SNRmeanRates[1,:]<50]
    GPEmr=GPEmr[SNRmr>=25]
    SNRmr=SNRmr[SNRmr>=25]
    lines.append(ax.plot(GPEmr,SNRmr,**{'color':[0.1,1,0.1], 'linewidth':5}))
    ax.text( GPEmr[-1]+2, SNRmr[-1], str(int(GPEmr[-1]))+' spikes/s',
                 fontsize=pylab.rcParams['text.fontsize']-6, 
                 **{'color': 'k'})
    
    
    SNRmr=SNRmeanRates[1,SNRmeanRates[1,:]<25]
    GPEmr=GPEmeanRates[SNRmeanRates[1,:]<25]
    GPEmr=GPEmr[SNRmr>=12.5]
    SNRmr=SNRmr[SNRmr>=12.5]
    lines.append(ax.plot(GPEmr,SNRmr,**{'color':[0.5,0.5,1], 'linewidth':5}))
    ax.text( GPEmr[-1]+2, SNRmr[-1], str(int(GPEmr[-1]))+' spikes/s',
                 fontsize=pylab.rcParams['text.fontsize']-6, 
                 **{'color': 'k'})

    
    ax.plot(GPEmeanRates,SNRmeanRates[1,:],**{'color':colors[0]})  
    
    ax.legend(linesnumpoints=1, loc='best')

    leg=ax.legend(lines,['100-50%', '50-25%','25-12.5%'], loc='best')
    frame  = leg.get_frame() 
    frame.set_visible(False) 

    ax.set_ylabel('Max firing rate SNr (%)') 
    ax.set_xlabel('Firing rate GPe (spikes/s)')
    ax.my_set_no_ticks( yticks=8, xticks=8 )
    ax.set_xlim(misc.adjust_limit([0,100]))
    ax.set_ylim(misc.adjust_limit([0,110]))
    
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,100]) 
    
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['text.fontsize'], 
                 **{'color': color})
        
def plot_GPE_vs_SNR_const_syn_events(ax, nGPE_range, SNRmeanRates):
    global N_GPE
    colors=['r','c' ] 
    labels=[r'$\delta_{ref}^{GPe}$' , r'$\delta_{dep}^{GPe}$']
    coords=[[0.35, 0.18], [ 0.2, 0.35]]
    
    for i, color in enumerate(colors):
        ax.plot(abs(nGPE_range-N_GPE),SNRmeanRates[i,:],**{'color':color})  
    
    ax.set_ylabel('Firing rate SNr (spikes/s)') 
    ax.set_xlabel('Pausing GPe (#)')
    ax.my_set_no_ticks( yticks=8, xticks=7 )

    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,25])
    ax.set_xlim(misc.adjust_limit([0,25]))
    ax.set_ylim(misc.adjust_limit([20,80]))
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['text.fontsize'], 
                 **{'color': color})

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

def simulate_GPE_vs_SNR_rate(load=True):
    global N_GPE
    global N_MSN
    global N_STN
    global MSN_BASE_RATE
   
    # Path were data is saved. For example the spike trains.
    save_at = (SPATH+'/'+NEURON_MODELS[0]+'-' + '-GPE_vs_SNR_rate.pkl') 
    
    GPEmeanRates=numpy.arange(0,150,1)
    SNRmeanRates=[]

    sim_time=50000.
    I_e=0.
    
    if not load:
        for r in GPEmeanRates:
            my_nest.ResetKernel()
            model_list=models()
            my_nest.MyLoadModels( model_list, NEURON_MODELS )
            my_nest.MyLoadModels( model_list, SYNAPSE_MODELS_TESTED )
            my_nest.MyLoadModels( model_list, SYNAPSE_MODELS_BACKGROUND )
            
            MSN = MyPoissonInput( n=N_MSN)           
            GPE = MyPoissonInput( n=N_GPE)          
            STN = MyPoissonInput( n=N_STN) 
            
            I_e=my_nest.GetDefaults(NEURON_MODELS[0])['I_e']+SNR_INJECTED_CURRENT
            SNR = MyGroup( NEURON_MODELS[0], n=len(SYNAPSE_MODELS_TESTED), 
                           sd=True,params={'I_e':I_e})
            
            for id in GPE[:]:
                GPE.set_spike_times(id=id, rates=[r], times=[1], 
                                    t_stop=sim_time, 
                                    seed=int(numpy.random.random()*10000.0))  
            for id in MSN[:]:
                MSN.set_spike_times(id=id, rates=[MSN_BASE_RATE], times=[1], 
                                    t_stop=sim_time, 
                                    seed=int(numpy.random.random()*10000.0))              

            for id in STN[:]:
                STN.set_spike_times(id=id, rates=[STN_BASE_RATE], times=[1], 
                                    t_stop=sim_time, 
                                    seed=int(numpy.random.random()*10000.0))
        
                
            for i, syn in enumerate(SYNAPSE_MODELS_TESTED):
                my_nest.ConvergentConnect(GPE[:],[SNR[i]], model=syn)
                my_nest.ConvergentConnect(MSN[:],[SNR[i]], model=SYNAPSE_MODELS_BACKGROUND[0])
                my_nest.ConvergentConnect(STN[:],[SNR[i]], model=SYNAPSE_MODELS_BACKGROUND[1])
                
            my_nest.MySimulate( sim_time )
            SNR.get_signal( 's') # retrieve signal
  
                
            SNRmeanRates.append(SNR.signals['spikes'].mean_rates( 5000, sim_time))   
        
        SNRmeanRates=numpy.array(SNRmeanRates).transpose()
        GPEmeanRates=numpy.array(GPEmeanRates)
        

        rateAtThr=''
        for SNRr in SNRmeanRates:
            tmp=str(GPEmeanRates[SNRr>=SELECTION_THR][-1])
            rateAtThr+=' '+tmp[0:4]
        
            
        
        
        s='\n'
        s =s + 'GPE vs SNr rate:\n'   
        s = s + ' %s %5s %3s \n' % ( 'N GPEs:', str ( N_GPE ),  '#' )  
        s = s + ' %s %5s %3s \n' % ( 'Max SNr rate:', str ( SNRmeanRates[0] ),  '#' )  
        s = s + ' \n%s \n%5s %3s \n' % ( 'GPE rates:', str ( GPEmeanRates[0] ) + '-'+ 
                                         str ( GPEmeanRates[-1] ),  'spikes/s' ) 
        s = s + ' %s %5s %3s \n' % ( 'Threshold SNr:', str ( SELECTION_THR ),  'spikes/s' )
        s = s + ' \n%s \n%5s %3s \n' % ( 'GPE rate at threshold SNr:', str ( rateAtThr ),  'spikes/s' )   
        s = s + ' \n%s %5s %3s \n' % ( 'Simulation time:', str ( sim_time), 'ms' )
        s = s + ' %s %5s %3s \n' % ( 'I_e:', str ( I_e ), 'pA' )
        infoString=s
        
        
        misc.pickle_save([GPEmeanRates, SNRmeanRates, infoString], 
                                save_at)
    
    elif load:
        GPEmeanRates, SNRmeanRates, infoString = misc.pickle_load(save_at)
    return GPEmeanRates, SNRmeanRates, infoString

def simulate_GPE_vs_SNR_const_syn_events(load=True):
    global N_GPE
    global N_MSN
    global MSN_BASE_RATE
    global SNR_INJECTED_CURRENT
    
    save_at = (SPATH+'/'+NEURON_MODELS[0]+'-' + '-GPE_vs_SNR_const_syn_events.pkl') 
    
    nGPE_range=numpy.arange(N_GPE,4,-1)
    
    # To maintain CONSTANT_SYN_EVENTS in to SNr while changing number of pausing 
    # GPe we have to increase the mean rate of the non-pausing GPe's
    
    GPEmeanRates=CONSTANT_SYN_EVENTS/nGPE_range
    SNRmeanRates=[]

    sim_time=10000.
    I_e=0.
    
    if not load:
        for r, n_gpe in zip(GPEmeanRates,nGPE_range):
            my_nest.ResetKernel()
            model_list=models()
            my_nest.MyLoadModels( model_list, NEURON_MODELS )
            my_nest.MyLoadModels( model_list, SYNAPSE_MODELS_TESTED )      
            my_nest.MyLoadModels( model_list, SYNAPSE_MODELS_BACKGROUND )
            
            GPE = MyPoissonInput( n=n_gpe)          
            MSN = MyPoissonInput( n=N_MSN)          
            STN = MyPoissonInput( n=N_STN)    
            
            I_e=my_nest.GetDefaults(NEURON_MODELS[0])['I_e']+SNR_INJECTED_CURRENT
            SNR = MyGroup( NEURON_MODELS[0], n=len(SYNAPSE_MODELS_TESTED), 
                           sd=True,params={'I_e':I_e})

            for id in GPE[:]:
                GPE.set_spike_times(id=id, rates=[r], times=[1], 
                                    t_stop=sim_time, 
                                    seed=int(numpy.random.random()*10000.0))  
            for id in MSN[:]:
                MSN.set_spike_times(id=id, rates=[MSN_BASE_RATE], times=[1], 
                                    t_stop=sim_time, 
                                    seed=int(numpy.random.random()*10000.0))              
            for id in STN[:]:
                STN.set_spike_times(id=id, rates=[STN_BASE_RATE], times=[1], 
                                    t_stop=sim_time, 
                                    seed=int(numpy.random.random()*10000.0))          
                
            for i, syn in enumerate(SYNAPSE_MODELS_TESTED):
                my_nest.ConvergentConnect(GPE[:],[SNR[i]], model=syn)
                my_nest.ConvergentConnect(MSN[:],[SNR[i]], model=SYNAPSE_MODELS_BACKGROUND[0])
                my_nest.ConvergentConnect(STN[:],[SNR[i]], model=SYNAPSE_MODELS_BACKGROUND[1])
    
            
            my_nest.MySimulate( sim_time )

            SNR.get_signal( 's') # retrieve signal
                
            SNRmeanRates.append(SNR.signals['spikes'].mean_rates(5000,sim_time))   
        
        SNRmeanRates=numpy.array(SNRmeanRates).transpose()
        GPEmeanRates=numpy.array(GPEmeanRates)

        rateAtThr=''
        for SNRr in SNRmeanRates:
            tmp=str(GPEmeanRates[SNRr>=SELECTION_THR][-1])
            rateAtThr+=' '+tmp[0:4]
        
            
        
        
        s='\n'
        s =s + 'GPE vs SNr rate:\n'   
        s = s + ' %s %5s %3s \n' % ( 'N GPEs:', str ( N_GPE ),  '#' )  
        s = s + ' \n%s %5s %3s \n' % ( 'GPE rates:', str ( GPEmeanRates[0] ) + '-'+ 
                                         str ( GPEmeanRates[-1] ),  'spikes/s' ) 
        s = s + ' %s %5s %3s \n' % ( 'Threshold SNr:', str ( SELECTION_THR ),  'spikes/s' )
        s = s + ' \n%s %5s %3s \n' % ( 'GPE rate at threshold SNr:', str ( rateAtThr ),  'spikes/s' )   
        s = s + ' \n%s %5s %3s \n' % ( 'Simulation time:', str ( sim_time), 'ms' )
        s = s + ' %s %5s %3s \n' % ( 'I_e:', str ( I_e ), 'pA' )
        s = s + ' %s %5s %3s \n' % ( 'Steady state rate ref:', str ( round(SNRmeanRates[0][0],1) ), 'pA' )
        s = s + ' %s %5s %3s \n' % ( 'Steady state rate dyn:', str ( round(SNRmeanRates[1][0],1) ), 'pA' )
        statusSynapse=[]
        for syn in SYNAPSE_MODELS_TESTED:
            statusSynapse.append( my_nest.GetDefaults(syn) )

            for ss in statusSynapse:
                s = s + '\n'  
                s = s + ' %s %10s\n' % ( 'Synapse', ss['synapsemodel'])   
                s = s + ' %s %5s %3s\n' % ( 'Weight', 
                                      str( round( ss['weight'], 1) ), 'nS')
        
        infoString=s
        
        
        misc.pickle_save([GPEmeanRates, SNRmeanRates, infoString], 
                                save_at)
    
    elif load:
        GPEmeanRates, SNRmeanRates, infoString = misc.pickle_load(save_at)
    return nGPE_range, GPEmeanRates, SNRmeanRates, infoString
  
print 'Simulation'

# SIMULATION
infoString=''

# simulate_GPE_vs_SNR_rate
GPEmeanRates, SNRmeanRates, s = simulate_GPE_vs_SNR_rate(load=True)
infoString=infoString+s    

# simulate_GPE_vs_SNR_rate
nGPE_range, GPEmeanRates_cse, SNRmeanRates_cse, s = simulate_GPE_vs_SNR_const_syn_events(load=True)
infoString=infoString+s  

# DISPLAY
plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0, fontsize=12)
font_size_text = 8
fig = pylab.figure( facecolor = 'w' )

ax_list = []
ax_list.append( MyAxes(fig, [ .075, .37, .135, .26 ] ) )    # text box
ax_list.append( MyAxes(fig, [ .26,  .6,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .53, .6,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .6,  .165, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .26,  .1,  .165, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .53, .1,  .165, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .8,   .1,  .165, .34 ] ) )    # 

# Text
ax=ax_list[0]
plot_text(ax, infoString)


# plot_selection_vs_neurons
ax=ax_list[1]
plot_selection_vs_neurons(ax, GPEmeanRates, SNRmeanRates)

data={}
data['GPE_mean_rates']=GPEmeanRates
data['SNR_mean_rates']=SNRmeanRates
misc.pickle_save(data, os.getcwd()+'/output/mean_rates_GPE_constant_supression'+
            NEURON_MODELS[0])

# plot_selection_vs_neurons
ax=ax_list[2]
plot_selection_vs_neurons_percent(ax, GPEmeanRates, SNRmeanRates)
# IF
ax=ax_list[3]
plot_GPE_vs_SNR_const_syn_events(ax, nGPE_range, SNRmeanRates_cse)

pylab.show()

name = sys.argv[0].split('/')[-1].split('.')[0]
fig.savefig( picture_dir + '/' + name  + '.svg', format = 'svg')
fig.savefig( picture_dir + '/' + name  + '.pdf', format = 'pdf')
