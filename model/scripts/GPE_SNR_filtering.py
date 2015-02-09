
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
picture_dir=  '/'.join(os.getcwd().split('/')[0:-3]) + '/pictures'     
                
sys.path.append(model_dir) 
sys.path.append(code_dir+'/nest_toolbox') 
spath  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]

# Then import model and network
from model_params import models
from src import misc, my_nest, my_signals, plot_settings
from src.my_population import MyGroup  
from src.my_axes import MyAxes 

# Simulate or use stored data
LOAD=False

neuronModels=['SNR_aeif']
synapseModels=['GPE_SNR_gaba_s_min','GPE_SNR_gaba_s_max', 'GPE_SNR_gaba_p']

def plot_example_GPE(ax, GPE):
    ax_twinx=ax.my_twinx()
    GPE.signals['spikes'].raster_plot(display=ax_twinx,kwargs={'color':'k',
                                                               'zorder':1})  
    ax_twinx.set_ylabel('Neuron id')
    time_bin=500
   # GPE.signals['spikes'].my_firing_rate( bin=time_bin, display=ax,
   #                                       kwargs={'color':'r',
   #                                               'linewidth':9,
   #                                               'zorder':20})
    GPE.signals['spikes'].my_firing_rate( bin=time_bin, display=ax,
                                          kwargs={'color':'k',
                                                  'linewidth':3,})
    ax.my_set_no_ticks( yticks=4, xticks = 4 ) 
    ax_twinx.my_set_no_ticks( yticks=4, xticks = 4 ) 
    ax.set_title('bin=%i'%(time_bin),**{'fontsize':12})
    ax.set_ylabel('Frequency GPEs (Hz)')
    
def plot_example_snr(ax, SNR):
    time_bin=500
    
    colors = misc.make_N_colors('Blues', 5)
    colors=['g','r', colors[1], colors[2], colors[3]]   
    
    signal=SNR.signals['spikes']
    signal.my_firing_rate(id_list=[SNR[0]], bin=time_bin, display=ax,
                          kwargs={'color':colors[0]})
    signal.my_firing_rate(id_list=[SNR[1]], bin=time_bin, display=ax,
                          kwargs={'color':colors[1]})
    signal.my_firing_rate(id_list=[SNR[2]], bin=time_bin, display=ax,
                          kwargs={'color':colors[2]})
    #signal.my_firing_rate(id_list=[SNR[3]], bin=time_bin, display=ax,
    #                      kwargs={'color':colors[3]})
    #signal.my_firing_rate(id_list=[SNR[4]], bin=time_bin, display=ax,
    #                      kwargs={'color':colors[4]})
    
    ax.my_set_no_ticks( yticks=4, xticks = 4 ) 
    ax.set_ylim([0,40])
    ax.text( 7000, 33, 'Weak' , **{ 'color' : colors[0] })  
    ax.text( 4000, 7, 'Strong' , **{ 'color' : colors[1] }) 
    ax.text( 7000, 24, 'Set 1' , **{ 'color' : colors[2] }) 
    #ax.text( 7000, 20, 'Set 1+2' , **{ 'color' : colors[3] }) 
    #ax.text( 7000, 16, 'Set 2' , **{ 'color' : colors[4] }) 
    
    #ax.legend(('Weak','Strong','Data 1','Data 1+2','Data 2'), loc='best')
    ax.set_title('bin=%i'%(time_bin),**{'fontsize':12})
    ax.set_ylabel('Frequency SNr (Hz)') 

def plot_selection_vs_neurons(ax, GPEmeanRates, SNRmeanRates):
    colors = misc.make_N_colors('Blues', 5)
    colors=['g','r', colors[1], colors[2], colors[3]] 
    
    ax.plot(GPEmeanRates,SNRmeanRates[0,:],**{'label':'Weak','color':colors[0]})  
    ax.plot(GPEmeanRates,SNRmeanRates[1,:],**{'label':'Strong','color':colors[1]})  
    ax.plot(GPEmeanRates,SNRmeanRates[2,:],**{'label':'Data 1','color':colors[2]})  
    #ax.plot(GPEmeanRates,SNRmeanRates[3,:],**{'label':'Data 1+2','color':colors[3]})
    #ax.plot(GPEmeanRates,SNRmeanRates[4,:],**{'label':'Data 3','color':colors[4]})
    ax.set_ylabel('Frequency SNr (Hz)') 
    ax.set_xlabel('Frequency GPEs (Hz)')
    ax.my_set_no_ticks( yticks=4, xticks = 4 )
    ax.set_ylim([0,40])
    #ax.set_xlim([0,1.4])
    ax.text( 0.7, 0.85, 'Weak' , transform=ax.transAxes, **{ 'color' : colors[0] })  
    ax.text( 0.7, 0.75, 'Strong' , transform=ax.transAxes, **{ 'color' : colors[1] }) 
    ax.text( 0.7, 0.65,'Set 1' , transform=ax.transAxes, **{ 'color' : colors[2] }) 
    #ax.text( 0.7, 0.55, 'Set 1+2' , transform=ax.transAxes, **{ 'color' : colors[3] }) 
    #ax.text( 0.7, 0.45, 'Set 2' , transform=ax.transAxes, **{ 'color' : colors[4] })  

def plot_text(ax, infoString=''):
    
    my_nest.ResetKernel()
    model_list=models()
    my_nest.MyLoadModels( model_list, neuronModels )
    
    SNR = MyGroup( neuronModels[0], 1, mm_dt = 0.1)
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
    
def simulate_example_GPE_snr():  
    nFun=0  # Function number
    nSim=0  # Simulation number within function
    
    rates=numpy.array([20,30])
    times=numpy.array([0.,5000.])
    nGPE =10
    simTime=10000.
    I_e=0.
    
    my_nest.ResetKernel()
    model_list=models()
    my_nest.MyLoadModels( model_list, neuronModels )
    my_nest.MyLoadModels( model_list, synapseModels )

    GPE = MyGroup( 'spike_generator', nGPE, mm_dt=1.0, mm=False, sd=False,
                   spath=spath, 
                   siter=str(nFun)+str(nSim))  
    SNR = MyGroup( neuronModels[0], n=len(synapseModels), params={'I_e':I_e},
                   mm_dt = .1, mm=False, spath=spath, 
                   siter=str(nFun)+str(nSim) )
    nSim+=1
    if not LOAD:  
        spikeTimes=[]
        for i in range(nGPE):
            spikes=misc.inh_poisson_spikes( rates, times,                        
                                        t_stop=simTime, 
                                        n_rep=1, seed=i )
            my_nest.SetStatus([GPE[i]], params={ 'spike_times':spikes } ) 
            for spk in spikes: spikeTimes.append((i,spk))   
        # add spike list for GPE to GPE spike list
        GPE.signals['spikes'] = my_signals.MySpikeList(spikeTimes, GPE.ids)     
        GPE.save_signal( 's') 
       
    
        
        for i, syn in enumerate(synapseModels):
            my_nest.ConvergentConnect(GPE[:],[SNR[i]], model=syn)

        my_nest.MySimulate( simTime )
        SNR.save_signal( 's') 
        SNR.get_signal( 's' ) # retrieve signal
    elif LOAD: 
        GPE.load_signal( 's')
        SNR.load_signal( 's')
   
    
    SNR_rates=[SNR.signals['spikes'].mean_rates(0,5000), 
               SNR.signals['spikes'].mean_rates(5000, 10000)]     
    for i in range(0, len(SNR_rates)):      
        for j in range(0, len(SNR_rates[0])):
            SNR_rates[i][j]=int(SNR_rates[i][j])
    s='\n'
    s =s + 'Example plot GPE and SNr:\n' 
    s =s + 'Synapse models:\n'
    for syn in synapseModels:
        s = s + ' %s\n' % (syn )    
    s = s + ' %s %5s %3s \n' % ( 'N GPE:', str ( nGPE ),  '#' )    
    s = s + ' %s %5s %3s \n' % ( 'GPE Rates:',   str ( [str(round(r,1)) 
                                                        for r in rates]),'Hz' )     
    s = s + ' %s %5s %3s \n' % ( '\nSNR Rates 0-5000:\n',   
                                 str ( SNR_rates [0]) ,'Hz' )   
    s = s + ' %s %5s %3s \n' % ( '\nSNR Rates 10000-5000:\n',  
                                  str ( SNR_rates [1]) ,'Hz' )   
    s = s + ' %s %5s %3s \n' % ( '\nTimes:', str ( times), 'ms' )
    s = s + ' %s %5s %3s \n' % ( 'I_e:', str ( I_e ), 'pA' )
    
    infoString=s
    
    
    return GPE, SNR, infoString

def simulate_GPE_vs_SNR_rate():
    nFun=1  # Function number
    nSim=0  # Simulation number within function
    
    GPEmeanRates=numpy.arange(1,50,2)
    SNRmeanRates=[]
    nGPE =10
    simTime=10000.
    I_e=0.
    
    
    for r in GPEmeanRates:
        my_nest.ResetKernel()
        model_list=models()
        my_nest.MyLoadModels( model_list, neuronModels )
        my_nest.MyLoadModels( model_list, synapseModels )
    
        GPE = MyGroup( 'spike_generator', nGPE, mm_dt=1.0, mm=False, sd=False,
                       spath=spath, siter=str(nFun)+str(nSim))  
        SNR = MyGroup( neuronModels[0], n=len(synapseModels), params={'I_e':I_e},
                       mm_dt=.1, mm=False, spath=spath, 
                   siter=str(nFun)+str(nSim))
        nSim+=1
        
        if not LOAD:
            spikeTimes=[]
            for i in range(nGPE):
                spikes=misc.inh_poisson_spikes( numpy.array([r]), 
                                                numpy.array([0]),                        
                                            t_stop=simTime, 
                                            n_rep=1, seed=i )
                my_nest.SetStatus([GPE[i]], params={ 'spike_times':spikes } ) 
                for spk in spikes: spikeTimes.append((i,spk))   
            # add spike list for GPE to GPE spike list
            GPE.signals['spikes'] = my_signals.MySpikeList(spikeTimes, GPE.ids)     
            GPE.save_signal( 's') 
           
    
            
            for i, syn in enumerate(synapseModels):
                my_nest.ConvergentConnect(GPE[:],[SNR[i]], model=syn)
        
            my_nest.MySimulate( simTime )
            SNR.save_signal( 's') 
            SNR.get_signal( 's') # retrieve signal
        elif LOAD: 
            SNR.load_signal( 's')     
            
        SNRmeanRates.append(SNR.signals['spikes'].mean_rates(0,simTime))   
    
    SNRmeanRates=numpy.array(SNRmeanRates).transpose()
    GPEmeanRates=numpy.array(GPEmeanRates)
    
    THR=2.
    rateAtThr=''
    for SNRr in SNRmeanRates:
        tmp=str(GPEmeanRates[SNRr>=THR][-1])
        rateAtThr+=' '+tmp[0:4]
    
        
    s='\n'
    s =s + 'GPE vs SNr rate:\n'   
    s = s + ' %s %5s %3s \n' % ( 'N GPEs:', str ( nGPE ),  '#' )  
    s = s + ' \n%s \n%5s %3s \n' % ( 'GPE rates:', str ( GPEmeanRates[0] ) + '-'+ 
                                     str ( GPEmeanRates[-1] ),  'Hz' ) 
    s = s + ' %s %5s %3s \n' % ( 'Threshold SNr:', str ( THR ),  'Hz' )
    s = s + ' \n%s \n%5s %3s \n' % ( 'GPE rate at threshold SNr:', str ( rateAtThr ),  'Hz' )   
    s = s + ' \n%s %5s %3s \n' % ( 'Simulation time:', str ( simTime), 'ms' )
    s = s + ' %s %5s %3s \n' % ( 'I_e:', str ( I_e ), 'pA' )
    infoString=s
    
    return GPEmeanRates, SNRmeanRates, infoString
    
print 'Simulation'

# SIMULATION
infoString=''

# simulate_example_GPE_snr
GPE, SNR, s = simulate_example_GPE_snr()
infoString=infoString+s    

# simulate_GPE_vs_SNR_rate
GPEmeanRates, SNRmeanRates, s = simulate_GPE_vs_SNR_rate()
infoString=infoString+s    


# DISPLAY
plot_settings.set_mode(mode='by_fontsize', w = 750.0, h = 400.0, fontsize=12)
font_size_text = 9
fig = pylab.figure( facecolor = 'w' )

ax_list = []
ax_list.append( MyAxes(fig, [ .1, .37, .18,  .26 ] ) )    # text box
ax_list.append( MyAxes(fig, [ .35, .6, .24, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .73, .6, .24, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .35, .1, .24, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .73, .15, .24, .34 ] ) )    # 

# Text
ax=ax_list[0]
plot_text(ax, infoString)

# Example GPE
ax=ax_list[1]
plot_example_GPE(ax, GPE)

# Example snr
ax=ax_list[2]
plot_example_snr(ax, SNR)


# plot_selection_vs_neurons
ax=ax_list[3]
plot_selection_vs_neurons(ax, GPEmeanRates, SNRmeanRates)



# IF
#ax=ax_list[4]


pylab.show()

name = sys.argv[0].split('/')[-1].split('.')[0]
fig.savefig( picture_dir + '/' + name  + '.svg', dpi = 500, format = 'svg')
