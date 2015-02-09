#! MSN to SNR
#! ===============================
#! Investigate how static vs plastic synapse work between 
#! striatal population of 50 neurons and one SNR neuron. n 
#! experiments is run and average activity is displayed. 

#! Experiment 1:
#! --------------
#! Simulate 50 MSN firing at 1 Hz and one SNR firng at 
#! 10-20 Hz. Use static and plastic synapses. Tune static 
#! synapse such that SNR is not inhibit when all 50 neurons 
#! fire at 1 Hz. Let one MSN fire at 25 Hz. Is it enough to
#! inhibit  SNR or not. If not fire up another MSN. I can 
#! basically create several 50 - 1 MSN - SNR networks,

#! Experiment 2:
#! --------------
#! Simulate 50 MSN with facilitating synapse tuned after 
#! Sims 2008. What happens at 5, 10, 15, 20, 25 Hz?

#! Imports

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


def create_models(addModels):
    ''' copy model'''
    model_list = models()                               # Get model list
    for model in model_list:                            # Create models
        if model[1] in addModels : 
            my_nest.CopyModel(model[0], model[1], model[2] )

def create_input_population( nInput, nRep=1, simTime=20000 ):
    ''' Define input as inhomogenous poisson processes '''
    
    spike_times = []
    inputName  = 'spike_generator'
    Input = MyGroup( inputName, nInput, mm_dt=1.0, 
                   spath=spath, sname='MSN', mm=False, sd=False)  
    for i in range(nInput): 
        rates=meanInputRates[:]
        times=meanInputTimings[:]
        
        if Input.ids[i] in selectedInputIds:
            rates.extend(selectedInputRates)
            times.extend(selectedInputTimings)
            
        rates.append(inbetweenInputRate)
        times.append(inbetweenInputTime)
            
        spikes=misc.inh_poisson_spikes( rates, times,                        
                                       t_stop=simTime, 
                                       n_rep=nRep, seed=i )
        
        # create spike list for input
        for spk in spikes: spike_times.append((i,spk))                          
        my_nest.SetStatus([Input.ids[i]], params={ 'spike_times':spikes } ) 
        
 
    # add spike list for input to input spike list
    Input.signals['spikes'] = my_signals.MySpikeList(spike_times, Input.ids) 
    
    return Input

def create_output_population( nOutput, outputAddCurrent, outputName, sname, spath):
    
    create_models( outputName )
    
    Output = MyGroup( outputName, nOutput, mm_dt=1.0, 
                      sname=sname, spath=spath )          
           
    I_e    = my_nest.GetStatus(Output.local_ids,'I_e')[0]
    
    # add  output current
    my_nest.SetStatus( Output[:], { 'I_e':I_e + outputAddCurrent } )
    
    return Output
 
def plot_input_raster(ax, Input, nYticks, simTime):
    ''' Raster plot of input'''
    
    Input.signals['spikes'].raster_plot( Input[:], t_start=0, t_stop=simTime, 
                        display=ax, 
                        kwargs={ 'marker':'o', 'markersize':1., 
                                 'markeredgecolor':'b' } )

    ax.remove_axis( xaxis=True ) 
    ax.set_no_ticks(yticks=nYticks )
    
def plot_input_mean_rate(ax, binSize, Input, nRep, nYticks, simTime):
    Input.signals['spikes'].firing_rate( id_list = Input[:], bin = binSize,  
                              n_rep = 1, display = ax )
    ax.set_xlabel( 'Time (ms)' )
    ax.set_ylabel( 'Frequency (spikes/s)' )
    ax.set_title( '')
    ax.set_xlim([0, simTime/nRep])

    ax.set_no_ticks( yticks = nYticks, xticks = 3 )     

def plot_output_example(axes, colors, Output,  synapses ):    
    
    kwargs={}
    
    for ax, color, neuron, synapse in zip(axes, colors, Output[:], synapses ):
        
        kwargs['color']=color
        kwargs['linewidth']=1
        
        # Data visualization                                                                                                
        Output.signals['V_m'].plot( id_list = [ neuron ], display = ax, kwargs = kwargs  )
    
        ax.set_xlim(tLim)
        ax.remove_axis( xaxis=True, yaxis=True ) 
        ax.set_no_ticks(yticks=3 )

  
def plot_output_mean_rate( ax, binSize, colors, labels, nRep, Output, nOutputPerSynapse, simTime, synapses ):
    
    nSynapses = len(synapses)
    
    for i in range(nSynapses):   
      
      
        iId   = nOutputPerSynapse*i  
        jId   = nOutputPerSynapse*( i + 1 )

        
        # Display     
        Output.signals['spikes'].firing_rate( id_list=Output.ids[ iId:jId ], 
                                   bin=binSize, 
                                   n_rep=nRep, display=ax, 
                                   kwargs={ 'label':labels[i], 
                                         'color':colors[i] } )
        
    ax.set_xlabel('Time (ms)'); ax.set_ylabel('Frequency (spikes/s)')
    ax.set_xlim([0, simTime/nRep]); 
    ax.set_no_ticks( yticks=4, xticks = 3 )    

def multiply_weight_factor( synapses, weightFactor ):
    
    for i in range( len(synapses) ):   
      
      # Multiply with weight factor
      defaults=my_nest.GetDefaults(synapses[i])
      my_nest.SetDefaults(synapses[i], 
                       { 'weight':defaults['weight']*weightFactor } )

def simulate_output_example( Input, Output, simTime, synapses ):
    
    for synapse, target in zip(synapses, Output):
        my_nest.CC( Input[:], [target], model=synapse)                                   
                           
    Output.gaussian_model_par(sead=1)   # model param gaussian distributed
    Output.gaussian_conn_par(sead=1)    # connection param gaussian distributed
    
    my_nest.MySimulate( simTime )
   
    return Output
  
def simulate_output_mean_rate( Input, Output, nOutputPerSynapse, simTime, synapses ):
    
    nSynapses = len(synapses)               
    
    for i in range(nSynapses):   
      
      # connect input and output   
      my_nest.CC( Input,                                                        
                 Output[ i*nOutputPerSynapse:( i + 1 )*nOutputPerSynapse ], 
                 model=synapses[i])                                   
                           
    Output.gaussian_model_par(sead=1)   # model param gaussian distributed
    Output.gaussian_conn_par(sead=1)    # connection param gaussian distributed
    
    my_nest.MySimulate( simTime )
    
    return Output

def text(ax, binSize,  nRep, Output, synapses):
    
    statusOutput = my_nest.GetStatus( Output[:] )[0]
    
    statusSynapes=[]
    for s in synapses:
        statusSynapes.append( my_nest.GetDefaults(s) )
    
    
    tb      = ''   
    tb = tb + ' %s %10s\n' % ( 'Model', statusOutput['model'])                   
    tb = tb + '\n'                                             
    tb = tb + ' %s %5s %3s\n' % ( 'Bin size', binSize, 'ms')      
    tb = tb + ' %s %5s %3s\n' % ( 'Experiments', nRep, '# ')  
    tb = tb + ' %s %5s %3s\n' % ( 'Mean rates', meanInputRates, 'Hz ')   
    tb = tb + ' %s %5s %3s\n' % ( 'selected number', nSelectedInput, '# ')  
    tb = tb + ' %s %5s %3s\n' % ( 'selected rate', selectedInputRates, 'Hz ')
    tb = tb + ' %s %5s %3s\n' % ( 'Inbetween time', inbetweenInputTime, 'ms ')
    tb = tb + ' %s %5s %3s\n' % ( 'Inbetween rate', inbetweenInputRate, 'Hz ')
    tb = tb + '\n'
    tb = tb + ' %s %5s %3s \n' % ( 'E rev:', 
                                   str( statusOutput['GABAA_1_E_rev' ] ), 
                                   'mV' )
    tb = tb + ' %s %5s %3s \n' % ( 'Tau_decay:', 
                                   str( statusOutput['GABAA_1_Tau_decay' ] ), 
                                   'mV' )  

    for ss in statusSynapes:
        tb = tb + '\n'  
        tb = tb + ' %s %10s\n' % ( 'Synapse', ss['synapsemodel'])   
        tb = tb + ' %s %5s %3s\n' % ( 'Weight', 
                                      str( round( ss['weight'], 1) ), 'nS')
        
        if 'U' in ss.keys():
            tb = tb + '\n '                             
            tb = tb + ' %s %5s %3s \n' % ( 'U:', str ( ss['U'] ), '  ' ) 
            tb = tb + ' %s %5s %3s \n' % ( 'tau_fac:', 
                                           str ( ss['tau_fac'] ), 'ms' )                                 
            tb = tb + ' %s %5s %3s \n' % ( 'tau_rec:', 
                                           str ( ss['tau_rec'] ), 'ms' )
            tb = tb + ' %s %5s %3s \n' % ( 'tau_psc:', 
                                           str ( ss['tau_psc'] ), 'ms' ) 
     
    
    ax.text( 0.85, 0.5, tb , fontsize= font_size_text,
             horizontalalignment='right',
             verticalalignment='center',
             transform=ax.transAxes,     # to define coordinates in right scale
             **{ 'fontname' : 'monospace' })                           
    
    ax.remove_axis( xaxis=True, yaxis=True )
    ax.remove_spine(left=True,  bottom=True, right=True, top=True)

synapses=['MSN_SNR_gaba_s_min', 'MSN_SNR_gaba_s_max', 'MSN_SNR_gaba_p0',
          'MSN_SNR_gaba_p1','MSN_SNR_gaba_p2']
weightFactor=0.33
nSynapse=len(synapses)  


switch=3
if switch==1:                           	    # Figure 1 mean rate increase
    mode = 'mean-increase'
    outputName               = 'SNR_aeif'
    meanInputRates           = [0.1,0.3]
    meanInputTimings         = range(0, 1000, 500)
    selectedInputIds         = [] # range(20,28)# Specific input to MSNs   
    selectedInputRates       = [] # [15, 0.1, 15 ]
    selectedInputTimings     = [] #range( 1500, 3000, 500 )
    inbetweenInputRate       = 0.1
    inbetweenInputTime       = 1000  
    nSelectedInput           = len(selectedInputIds)
    experimentSimulationTime = 2000.0 
    nRep    = 50                                # number of experiments 
    tLim    =[0, inbetweenInputTime]
    simTime =  experimentSimulationTime*nRep     
    nInput = 500                             # Number of MSN
    
    nOutputPerSynapse = 1     
    nOutput = nSynapse*nOutputPerSynapse           
    outputAddCurrent = 100.          
elif switch==2:                                 # Figure 2 selection x1
    mode                     = 'selection-x1'
    outputName               = 'SNR_aeif'
    meanInputRates           = [0.1]
    meanInputTimings         = [0]
    selectedInputIds         = range(20,28)# Specific input to MSNs   
    selectedInputRates       = [15]
    selectedInputTimings     = [500]
    inbetweenInputRate       = 0.1
    inbetweenInputTime       = 1000  
    nSelectedInput           = len(selectedInputIds)
    experimentSimulationTime = 10000.0 
    inbetweenSimulationTime  = 0
    nRep    = 50                                # number of experiments 
    tLim    =[0, inbetweenInputTime]
    simTime =  experimentSimulationTime*nRep     
    nInput = 500                             # Number of MSN
    
    nOutputPerSynapse = 1     
    nOutput = nSynapse*nOutputPerSynapse           
    outputAddCurrent = 100.    
elif switch == 3:                               # Figure 3 selection x2
    mode                     = 'selection-x2'  
    outputName               = 'SNR_aeif'
    meanInputRates           = [0.1]
    meanInputTimings         = [0]
    selectedInputIds         = range(20,28)# Specific input to MSNs   
    selectedInputRates       = [15, 0.1, 15 ]
    selectedInputTimings     = range(500,2000,500)
    inbetweenInputRate       = 0.1
    inbetweenInputTime       = 2000  
    nSelectedInput           = len(selectedInputIds)
    experimentSimulationTime = 10000.0 
    inbetweenSimulationTime  = 0
    nRep    = 10                                # number of experiments 
    tLim    =[0, inbetweenInputTime]
    simTime =  experimentSimulationTime*nRep     
    nInput = 500                             # Number of MSN
    
    nOutputPerSynapse = 1     
    nOutput = nSynapse*nOutputPerSynapse           
    outputAddCurrent = 100.    
elif switch == 4:                               # Figure 4 selection x2-long
    mode                     = 'selection-x2-long'  
    outputName               = 'SNR_aeif'
    meanInputRates           = [0.1]
    meanInputTimings         = [0]
    selectedInputIds         = range(20,28)# Specific input to MSNs   
    selectedInputRates       = [15, 0.1, 15 ]
    selectedInputTimings     = [500, 1000, 2000]
    inbetweenInputRate       = 0.1
    inbetweenInputTime       = 2500  
    nSelectedInput           = len(selectedInputIds)
    experimentSimulationTime = 10000.0 
    inbetweenSimulationTime  = 0
    nRep    = 50                                # number of experiments 
    tLim    =[0, inbetweenInputTime]
    simTime =  experimentSimulationTime*nRep     
    nInput = 500                             # Number of MSN
    
    nOutputPerSynapse = 1     
    nOutput = nSynapse*nOutputPerSynapse           
    outputAddCurrent = 100.  
elif switch == 5:                               # Figure 5 selection x2 long-long
    mode                     = 'selection-x2-long-long'  
    outputName               = 'SNR_aeif'
    meanInputRates           = [0.1]
    meanInputTimings         = [0]
    selectedInputIds         = range(20,28)# Specific input to MSNs   
    selectedInputRates       = [15, 0.1, 15 ]
    selectedInputTimings     = [500, 1000, 7000]
    inbetweenInputRate       = 0.1
    inbetweenInputTime       = 7500  
    nSelectedInput           = len(selectedInputIds)
    experimentSimulationTime = 20000.0 
    inbetweenSimulationTime  = 0
    nRep    = 50                                # number of experiments 
    tLim    =[0, inbetweenInputTime]
    simTime =  experimentSimulationTime*nRep     
    nInput = 500                                # Number of MSN
    
    nOutputPerSynapse = 1     
    nOutput = nSynapse*nOutputPerSynapse           
    outputAddCurrent = 100.   
                                        


# Simulation
# ----------

# Example plot
useSynapses=synapses[0:5]
MSN=create_input_population( nInput, nRep=1, 
                         simTime=experimentSimulationTime  )                         
SNR_example_plot=create_output_population( len(useSynapses), outputAddCurrent, 
                                           outputName, 
                          sname=outputName, spath=spath )
create_models(synapses)                                 # create synapse models
multiply_weight_factor(synapses, weightFactor)          # change weights
if not mpiRun or sys.argv[1]=='simulate':
    SNR_example_plot = simulate_output_example( MSN, SNR_example_plot, 
                                                simTime, synapses )
# SNR mean rate
my_nest.ResetKernel()
MSN=create_input_population( nInput, nRep=nRep, simTime=simTime  )       
SNR_mean_rate=create_output_population( nOutput, outputAddCurrent, outputName, 
                                        sname=outputName, spath=spath )
create_models(synapses)                             # create synapse models
multiply_weight_factor(synapses, weightFactor)      # change weights

if not mpiRun or sys.argv[1]=='simulate':
    SNR_mean_rate = simulate_output_mean_rate( MSN, SNR_mean_rate, 
                                               nOutputPerSynapse, 
                                               simTime, synapses )
    
if mpiRun and sys.argv[1]=='simulate':    
    SNR_example_plot.save_signal( 'v','V_m', stop=simTime ) 
    SNR_mean_rate.save_signal( 's', start=1, stop=simTime )  
    sys.exit(0)      # stop script


# Display
# -------
print mpiRun and sys.argv[1]=='plot'
if mpiRun and sys.argv[1]=='plot':
        SNR_example_plot.load_signal( 'v','V_m') 
        SNR_mean_rate.load_signal( 's' ) 
else:
        SNR_example_plot.get_signal( 'v','V_m', stop=simTime ) 
        SNR_mean_rate.get_signal( 's', start=1, stop=simTime ) 
        


plot_settings.set_mode(mode='dynamic', w = 1100.0, h = 400.0)
font_size_text = 7
fig = pylab.figure( facecolor = 'w' )
#fig = pylab.figure( figsize = ( 34, 12 ), dpi = 50, facecolor='w' )

nb_yticks_rp = 4     # number of ticks raster plots
nb_yticks_sr = 3     # number of ticks spike rate plots
binSize      = 50    # spike rate bin    


ax_list = []
ax_list.append( MyAxes(fig, [ 0.05, .37, .09, .26 ] ) )   # text box
ax_list.append( MyAxes(fig, [ 0.22, .6,  .2, .3   ] ) )   # MSNs raster plot
ax_list.append( MyAxes(fig, [ 0.22, .12,  .2, .3  ] ) )   # MSNs rates
ax_list.append( MyAxes(fig, [ 0.57, .87, .2, .05  ] ) )   # SNR trace plastic
ax_list.append( MyAxes(fig, [ 0.57, .81, .2, .05  ] ) )   # SNR trace static w
ax_list.append( MyAxes(fig, [ 0.57, .74, .2, .05  ] ) )   # SNR trace static s
ax_list.append( MyAxes(fig, [ 0.57, .67, .2, .05  ] ) )   # SNR trace static s
ax_list.append( MyAxes(fig, [ 0.57, .6,  .2, .05  ] ) )   # SNR trace static s
ax_list.append( MyAxes(fig, [ 0.57, .12,  .2, .3 ] ) )    # SNR rates
#for ax in ax_list: fig.add_axes(ax)

# Text 
ax=ax_list[0]
my_nest.ResetKernel()
SNR=create_output_population( nOutput, outputAddCurrent,  # create output
               outputName, sname=outputName, 
               spath=spath )
create_models(synapses)                           # create synapse models
multiply_weight_factor(synapses, weightFactor)    # change weights
text(ax, binSize,  nRep, SNR, synapses)           # plot text

# MSN raster plot
ax=ax_list[1]
my_nest.ResetKernel()
MSN=create_input_population( nInput, nRep=1,              # create input
           simTime=experimentSimulationTime )                      
plot_input_raster( ax, MSN, nb_yticks_rp, 
                   experimentSimulationTime )
ax.set_xlim(tLim)

# MSN mean rate
ax=ax_list[2]
my_nest.ResetKernel()
MSN=create_input_population( nInput, nRep=nRep, 
           simTime=simTime )
plot_input_mean_rate(ax, binSize=binSize, Input=MSN, nRep=nRep,
                nYticks=nb_yticks_sr, 
                simTime=simTime)
ax.set_xlim(tLim)


# Example plots
axes=[ax_list[3], ax_list[4], ax_list[5], ax_list[6], ax_list[7]]
my_nest.ResetKernel()
colors = misc.make_N_colors('Blues', 5)
colors=['g','r', colors[1], colors[2], colors[3]]
plot_output_example(axes, colors, SNR_example_plot,  useSynapses )
  


# SNR mean rate
ax=ax_list[8]
labels=['Weak static', 'Strong static', 'Plastic data I', 'Plastic data I+II',
        'Plastic data II']

plot_output_mean_rate( ax, binSize, colors, labels, nRep, SNR_mean_rate, 
                       nOutputPerSynapse, simTime, synapses )
ax.legend( numpoints=1 , loc=[ 1.1, 1.6 ])
ax=ax_list[8]
ax.set_xlim(tLim)


#fig.savefig( path + '/GPE.eps' , dpi=300,papertype='a4',format='eps',transparent=True)
name = sys.argv[0].split('/')[-1].split('.')[0]
fig.savefig( picture_dir + '/' + name + '-'+ outputName +'-' +mode +'.svg' , 
             dpi = 500, format = 'svg')

pylab.show()
