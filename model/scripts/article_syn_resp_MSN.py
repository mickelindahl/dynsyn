#! Imports
import numpy
import pylab
import os
import sys


# Add directories to python path
sys.path.append(os.getcwd())                            
parent_dir='/'.join(os.getcwd().split('/')[0:-1])       
                   
model_dir=   '/'.join(os.getcwd().split('/')[0:-1])        
code_dir=  '/'.join(os.getcwd().split('/')[0:-2]) 
picture_dir=  '/'.join(os.getcwd().split('/')[0:-3]) + '/pictures'     
                
sys.path.append(model_dir) 
sys.path.append(code_dir+'/nest_toolbox') 
spath  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]

from model_params import models                               # Then import models     
from src import misc, my_nest, my_signals, plot_settings
from src.my_population import MyGroup 
from src.my_axes import MyAxes 


NEURON_MODELS=['SNR_aeif']
SYNAPES_MODELS=['MSN_SNR_gaba_s_min', 'MSN_SNR_gaba_s_max', 'GPE_SNR_gaba_s_min']

def plot_response_example_clamped_silent(ax, SNR):
    SNR.signals['V_m'].plot(id_list=[SNR[0]], display=ax,kwargs={'color':'b'})
    SNR.signals['V_m'].plot(id_list=[SNR[1]], display=ax,kwargs={'color':'g'})
    SNR.signals['V_m'].plot(id_list=[SNR[1]], display=ax,kwargs={'color':'m'})
    
    ax.text( 0.7, 0.4,r'$\delta_r^{GPe}$''Ref 1' , 
             transform=ax.transAxes, **{ 'color' : 'b' }) 
    ax.text( 0.7, 0.3, r'$\delta_r^{GPe}$''Ref 1', 
             transform=ax.transAxes, **{ 'color' : 'g' }) 
    ax.text( 0.7, 0.2, r'$\delta_r^{GPe}$''Ref 1' , 
              transform=ax.transAxes, **{ 'color' : 'm' })  
    ax.my_set_no_ticks( yticks=4, xticks = 4 )   
    
def plot_response_example_clamped_spiking(ax, SNR):
    SNR.signals['V_m'].plot(id_list=[SNR[0]], display=ax,kwargs={'color':'g'})
    SNR.signals['V_m'].plot(id_list=[SNR[1]], display=ax,kwargs={'color':'r'})
    SNR.signals['V_m'].plot(id_list=[SNR[2]], display=ax,kwargs={'color':'k'})
    ax.my_set_no_ticks( yticks=4, xticks = 4 )   
    ax.text( 0.1, 0.6,'Ref 1' , transform=ax.transAxes, **{ 'color' : 'g' }) 
    ax.text( 0.1, 0.5, 'Ref 2' , transform=ax.transAxes, **{ 'color' : 'r' }) 

def plot_voltage_ipsc(ax, voltage, ipsc):
    ax.plot(voltage,ipsc[0,:],**{'color':'g'})
    ax.plot(voltage,ipsc[1,:],**{'color':'r'})
    ax.my_set_no_ticks( yticks=4, xticks = 4 ) 
    ax.set_ylabel('IPSC amplitude (pA)') 
    ax.set_xlabel('Holding potential (mV)')
    ax.text( 0.1, 0.4,'Ref 1' , transform=ax.transAxes, **{ 'color' : 'g' }) 
    ax.text( 0.2, 0.8, 'Ref 2' , transform=ax.transAxes, **{ 'color' : 'r' }) 

def plot_voltage_ipsp(ax, voltage, ipsp):
    colors = ['g', 'r']
    ax.plot(voltage,ipsp[0,:],**{'color':colors[0]})
    ax.plot(voltage,ipsp[1,:],**{'color':colors[1]})
    ax.my_set_no_ticks( yticks=4, xticks = 4 ) 
    ax.set_ylabel('IPSP amplitude (mV)')
    ax.set_xlabel('Holding potential (mV)')
    ax.set_ylim()
    ax.text( 0.55, 0.8, 'Ref 1' , transform=ax.transAxes, **{ 'color' : colors[0] }) 
    ax.text( 0.3, 0.45, 'Ref 2' , transform=ax.transAxes, **{ 'color' : colors[1] }) 

def plot_text(ax, infoString):
    
    my_nest.ResetKernel()
    model_list=models()
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    my_nest.MyLoadModels( model_list, SYNAPES_MODELS )
    
    SNR = MyGroup( NEURON_MODELS[0], 2, mm_dt = 0.1)
    statusSNR = my_nest.GetStatus( SNR[:] )[0]
    
    statusSynapes=[]
    for s in SYNAPES_MODELS:
        statusSynapes.append( my_nest.GetDefaults(s) )
    tb = ''     
    tb = tb + infoString
         
    tb = tb + '\n'
    tb = tb + 'Neuron models:\n '                                           
    tb = tb + ' %s \n' % ( NEURON_MODELS[0] )
    
    tb = tb + '\n'
    tb = tb + ' %s %5s %3s \n' % ( 'E rev:', 
                                   str( statusSNR['GABAA_1_E_rev' ] ), 
                                   'mV' )
    tb = tb + ' %s %5s %3s \n' % ( 'Tau_decay:', 
                                   str( statusSNR['GABAA_1_Tau_decay' ] ), 
                                   'mV' )  
    
    for ss in statusSynapes:
        tb = tb + '\n'  
        tb = tb + ' %s %10s\n' % ( 'Synapse', ss['synapsemodel'])   
        tb = tb + ' %s %5s %3s\n' % ( 'Weight', 
                                      str( round( ss['weight'], 1) ), 'nS')
        
    ax.text( 0.85, 0.5, tb , fontsize= font_size_text,
             horizontalalignment='right',
             verticalalignment='center',
             transform=ax.transAxes,     # to define coordinates in right scale
             **{ 'fontname' : 'monospace' })                           
    
    ax.my_remove_axis( xaxis=True, yaxis=True )
    ax.my_remove_spine(left=True,  bottom=True, right=True, top=True)
     
def simulate_response_example_clamped_silent(I_e):
    '''
    Response when SNR is clamped to a voltage for Ref 2 and Ref 1 synapse
    '''
    spike_at = 500.  # ms
    simTime  = 700.  # ms

    my_nest.ResetKernel()
    model_list=models()
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    my_nest.MyLoadModels( model_list, SYNAPES_MODELS )
    
    SNR = MyGroup( NEURON_MODELS[0], len(SYNAPES_MODELS),
                   mm=True, mm_dt = 0.1, params={'I_e':I_e} )

    SG = my_nest.Create('spike_generator', params={'spike_times':[spike_at]} )

    for i in range(len(SYNAPES_MODELS)):
        my_nest.Connect(SG, [SNR[i]], model=SYNAPES_MODELS[i])

    my_nest.MySimulate(simTime)
    
    SNR.get_signal( 'v','V_m', stop=simTime ) # retrieve signal
    SNR.signals['V_m']=SNR.signals['V_m'].my_time_slice(400, 700)
    
    clamped_at  = SNR.signals['V_m'][1].signal[-1]
    size_weak   = min(SNR.signals['V_m'][1].signal)-clamped_at
    size_strong = min(SNR.signals['V_m'][2].signal)-clamped_at
    s=''
    s = s + ' %s %5s %3s \n' % ( 'Clamped at:', str ( round(clamped_at,1) ),  'mV' )    
    s = s + ' %s %5s %3s \n' % ( 'Size Ref 1:',   str ( round(size_weak,1) ),   'mV' )     
    s = s + ' %s %5s %3s \n' % ( 'Size Ref 2:', str ( round(size_strong,1) ), 'mV' )

    infoString=s
    
    return SNR, infoString

def simulate_response_example_clamped_spiking(I_e):
    '''
    Response when SNR is clamped to a voltage for Ref 2 and Ref 1 synapse
    '''
    spike_at = 500.  # ms
    simTime  = 700.  # ms

    my_nest.ResetKernel()
    model_list=models()
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    my_nest.MyLoadModels( model_list, SYNAPES_MODELS )
    
    SNR = MyGroup( NEURON_MODELS[0], len(SYNAPES_MODELS)+1, 
                   sd=True, mm=True,  mm_dt = 1., params={'I_e':I_e} )

    SG = my_nest.Create('spike_generator', params={'spike_times':[spike_at]} )

    for i in range(len(SYNAPES_MODELS)):
        my_nest.Connect(SG, [SNR[i]], model=SYNAPES_MODELS[i])

    my_nest.MySimulate(simTime)
    
    SNR.get_signal( 'v','V_m', stop=simTime ) # retrieve signal
    SNR.get_signal( 's') # retrieve signal
    SNR.signals['V_m'].my_set_spike_peak( 21, spkSignal= SNR.signals['spikes'] ) 
    SNR.signals['V_m']=SNR.signals['V_m'].my_time_slice(500, 560)
   
    return SNR

def simulate_voltage_ipsc(I_vec):

    simTime  = 700.  # ms
    spikes_at = numpy.arange(500., len(I_vec)*simTime,simTime)  # ms
    

    voltage     = []    # mV
    ipsc_weak   = []    # mV
    ipsc_strong = []    # mV

    my_nest.ResetKernel()
    model_list=models()
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    my_nest.MyLoadModels( model_list, SYNAPES_MODELS )
    
    SNR = MyGroup( NEURON_MODELS[0], len(SYNAPES_MODELS), mm=True,
                    mm_dt = 0.1 )

    SG = my_nest.Create('spike_generator', params={'spike_times':spikes_at} )

    for i in range(len(SYNAPES_MODELS)):
        my_nest.Connect(SG, [SNR[i]], model=SYNAPES_MODELS[i])
    
    simTimeTot=0
    for I_e in I_vec:
        
        my_nest.SetStatus(SNR[:], params={'I_e':float(I_e)})
        my_nest.MySimulate(simTime)
        simTimeTot+=simTime

    SNR.get_signal( 'c','I_GABAA_1', stop=simTimeTot ) # retrieve signal
    simTimeAcum=0
    
    
    for I_e in I_vec:
        
        signal=SNR.signals['I_GABAA_1'].my_time_slice(400+simTimeAcum, 700+simTimeAcum)
        simTimeAcum+=simTime
    
        clamped_at  = signal[1].signal[-1]
        minV=min(signal[1].signal)
        maxV=max(signal[1].signal)
        if abs(minV-clamped_at)<abs(maxV-clamped_at):
        
            size_weak   = max(signal[1].signal)-clamped_at
            size_strong = max(signal[2].signal)-clamped_at
        else:
            size_weak   = min(signal[1].signal)-clamped_at
            size_strong = min(signal[2].signal)-clamped_at
        
        voltage.append(clamped_at)
        ipsc_weak.append(size_weak)
        ipsc_strong.append(size_strong)
        
    ipsc=numpy.array([ipsc_weak, ipsc_strong])
    return voltage, ipsc

def simulate_voltage_ipsp(I_vec):

    simTime  = 700.  # ms
    spikes_at = numpy.arange(500., len(I_vec)*simTime,simTime)  # ms
    

    voltage     = []    # mV
    ipsp_weak   = []    # mV
    ipsp_strong = []    # mV

    my_nest.ResetKernel()
    model_list=models()
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    my_nest.MyLoadModels( model_list, SYNAPES_MODELS )
    
    SNR = MyGroup( NEURON_MODELS[0], len(SYNAPES_MODELS), mm=True,  
                   mm_dt = 0.1 )

    SG = my_nest.Create('spike_generator', params={'spike_times':spikes_at} )

    for i in range(len(SYNAPES_MODELS)):
        my_nest.Connect(SG, [SNR[i]], model=SYNAPES_MODELS[i])

    
    simTimeTot=0
    for I_e in I_vec:
        
        my_nest.SetStatus(SNR[:], params={'I_e':float(I_e)})
        my_nest.MySimulate(simTime)
        simTimeTot+=simTime

    SNR.get_signal( 'v','V_m', stop=simTimeTot ) # retrieve signal
    simTimeAcum=0
    
    
    for I_e in I_vec:
        
        signal=SNR.signals['V_m'].my_time_slice(400+simTimeAcum, 700+simTimeAcum)
        simTimeAcum+=simTime
    
        clamped_at  = signal[1].signal[-1]
        minV=min(signal[1].signal)
        maxV=max(signal[1].signal)
        if abs(minV-clamped_at)<abs(maxV-clamped_at):
        
            size_weak   = max(signal[1].signal)-clamped_at
            size_strong = max(signal[2].signal)-clamped_at
        else:
            size_weak   = min(signal[1].signal)-clamped_at
            size_strong = min(signal[2].signal)-clamped_at
        
        voltage.append(clamped_at)
        ipsp_weak.append(size_weak)
        ipsp_strong.append(size_strong)
        
    ipsp=numpy.array([ipsp_weak, ipsp_strong])
    return voltage, ipsp
        
        

print 'Simulation'   
    
# SIMULATION

# Example synaptic response silent
SNR_plot1, infoStringFig1=simulate_response_example_clamped_silent(-145.)


# Example synaptic response spiking
SNR_plot2 = simulate_response_example_clamped_spiking(0.)

# Voltage against ipsc size
voltage, ipsc=simulate_voltage_ipsc(range(-200,-130,10))

# Voltage against ipsp size
voltage, ipsp=simulate_voltage_ipsp(range(-200,-130,10))

# DISPLAY
plot_settings.set_mode(mode='by_fontsize', w = 1100.0, h = 450.0, fontsize=12)
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
plot_text(ax, infoStringFig1)


# Example synaptic response 
ax=ax_list[1]
plot_response_example_clamped_silent(ax, SNR_plot1)


# Voltage vs ipsp
ax=ax_list[2]
plot_voltage_ipsp(ax, voltage, ipsp)

# Example synaptic response spiking
#ax=ax_list[3]
#plot_response_example_clamped_spiking(ax, SNR_plot2)

# Current vs ipsc 
#ax=ax_list[3]
#plot_voltage_ipsc(ax, voltage, ipsc)




pylab.show()

name = sys.argv[0].split('/')[-1].split('.')[0]
fig.savefig( picture_dir + '/' + name  + '.svg' , format = 'svg')