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
                
sys.path.append(model_dir) 
sys.path.append(code_dir+'/nest_toolbox') 
spath  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]

from model_params import models                               # Then import models     
from src import misc, my_nest, my_signals, plot_settings
from src.my_population import MyGroup 
from src.my_axes import MyAxes 

### GLOBALS ###

# Paths and naming for saving data and picutes
FILE_NAME = sys.argv[0].split('/')[-1].split('.')[0]
model_name=os.getcwd().split('/')[-2]
PICTURE_PATH='/'.join(os.getcwd().split('/')[0:-3]) + '/pictures/'+model_name    
OUTPUT_PATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]

# Models
NEURON_MODELS=['GPE_aeif']
SYNAPSE_MODELS=['MSN_GPE_gaba_s_min', 'MSN_GPE_gaba_s_max', 
                'GPE_GPE_gaba_s','STN_GPE_ampa_s']

# Misc
SNR_INJECTED_CURRENT=-13.5
I_RANGE=range(-200, -20,1)

HEADER_SIMULATION_SETUP=( '**** BEGINNING GENERAL SCRIPT SETUP ****\n'+
                          'FILE_NAME:'+str(FILE_NAME)+'\n'+                         
                          'PICTURE_PATH:'+str(PICTURE_PATH)+'\n'+  
                          'OUTPUT_PATH:'+str(OUTPUT_PATH)+'\n\n'+    
                          'NEURON_MODEL:'+str(NEURON_MODELS)+'\n'+
                          'SYNAPSE_MODELS:'+str(SYNAPSE_MODELS)+'\n'+
                          'SNR_INJECTED_CURRENT:'+str(SNR_INJECTED_CURRENT)+'\n'+
                          'I_RANGE:'+str(I_RANGE)+'\n'+
                          '**** END GENERAL SCRIPT SETUP ****\n')
### END GLOBALS ###
    
def plot_response_example_clamped_silent(ax, SNR):
    
    colors = ['b', 'g','r','k']
    labels=[r'$ref_{init}^{MSN}$' , r'$ref_{max}^{MSN}$',  
            r'$ref_{25 Hz}^{GPe}$', r'$ref^{STN}$']
    coords=[[0.05, 0.4], [0.05, 0.14], [0.05, 0.64],[0.05, 0.84]]
    
    for i, color in enumerate(colors):   
        SNR.signals['V_m'].plot(id_list=[SNR[i]], display=ax,
                                kwargs={'color':color})
    
    
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['text.fontsize']+2, 
                 **{'color': color}) 
    hold_at=-68 
    lines = ax.lines    
    for line in lines:        
        line.set_xdata(line.get_xdata()-400)    
        #line.set_ydata(line.get_ydata()+68.) 
        
    #ax.text( 0.3, 0.05, 'Holding potential '+ str(hold_at), transform=ax.transAxes, 
    #             fontsize=pylab.rcParams['text.fontsize']-4, 
    #             **{'color': 'k'})
        
    ax.set_xlim(misc.adjust_limit([0,300]))
    #ax.set_ylim(misc.adjust_limit([-73,-67]))
    
    ax.my_set_no_ticks( yticks=8, xticks =7 ) 
    ax.set_ylabel('IPSP amplitude (mV)')

def plot_voltage_ipsp(ax, v_mat, ipsp_mat):
    colors = ['b', 'g','r']
    labels=[r'$ref_{init}^{MSN}$' , r'$ref_{max}^{MSN}$',  r'$ref_{25 Hz}^{GPe}$']
    
    coords=[[0.64, 0.3], [0.1, 0.29], [0.55, 0.7]]
    for i, color in enumerate(colors):
        ax.plot(v_mat[i,:],ipsp_mat[i,:],**{'color':color})

    ax.set_ylabel('IPSP amplitude (mV)')
    ax.set_xlabel('Holding potential (mV)')
    
    ax.set_xlim(misc.adjust_limit([-80,-65]))
    #ax.set_ylim(misc.adjust_limit([-4,2]))
    
    lines = ax.lines    
    for line in lines:
        misc.slice_line(line, xlim=[-80,-65])   
    
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['text.fontsize']+2, **{'color': color}) 
    ax.my_set_no_ticks( yticks=8, xticks = 6 )
    
 
    
def plot_text(ax, infoString):
    
    my_nest.ResetKernel()
    model_list, model_dict=models()
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    my_nest.MyLoadModels( model_list, SYNAPSE_MODELS )
    
    SNR = MyGroup( NEURON_MODELS[0], 2, mm_dt = 0.1)
    statusSNR = my_nest.GetStatus( SNR[:] )[0]
    
    statusSynapes=[]
    for s in SYNAPSE_MODELS:
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
    tb = tb + ' %s %5s %3s \n' % ( 'E rev:', 
                                   str( statusSNR['AMPA_E_rev' ] ), 
                                   'mV' )
    tb = tb + ' %s %5s %3s \n' % ( 'Tau_decay:', 
                                   str( statusSNR['AMPA_Tau_decay' ] ), 
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
    model_list, model_dict=models()
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    my_nest.MyLoadModels( model_list, SYNAPSE_MODELS )
    
    SNR = MyGroup( NEURON_MODELS[0], len(SYNAPSE_MODELS),
                   mm=True, mm_dt = 0.1, params={'I_e':I_e} )

    SG = my_nest.Create('spike_generator', params={'spike_times':[spike_at]} )

    for i in range(len(SYNAPSE_MODELS)):
        my_nest.Connect(SG, [SNR[i]], model=SYNAPSE_MODELS[i])

    my_nest.MySimulate(simTime)
    
    SNR.get_signal( 'v','V_m', stop=simTime ) # retrieve signal
    SNR.signals['V_m']=SNR.signals['V_m'].my_time_slice(400, 700)
    
    clamped_at  = SNR.signals['V_m'][1].signal[-1]
    size_MSN_weak   = min(SNR.signals['V_m'][1].signal)-clamped_at
    size_MSN_strong = min(SNR.signals['V_m'][2].signal)-clamped_at
    size_GPE_ref    = min(SNR.signals['V_m'][3].signal)-clamped_at
    s=''
    s = s + ' %s %5s %3s \n' % ( 'Clamped at:', str ( round(clamped_at,1) ),  'mV' )    
    s = s + ' %s %5s %3s \n' % ( r'$\delta_w^{MSN}$', 
                                 str ( round(size_MSN_weak,1) ),   'mV' )     
    s = s + ' %s %5s %3s \n' % ( r'$\delta_s^{MSN}$', 
                                 str ( round(size_MSN_strong,1) ), 'mV' )
    s = s + ' %s %5s %3s \n' % ( r'$\delta_s^{MSN}$', 
                                 str ( round(size_GPE_ref,1) ), 'mV' )

    infoString=s
    
    return SNR, infoString

def simulate_voltage_ipsp(I_vec):

    v_mat=[]
    ipsp_mat=[]
            
    for receptor, syn_model in zip( ['V_m', 'V_m','V_m'],
                                    SYNAPSE_MODELS):
        
        my_nest.ResetKernel()
        model_list, model_dict=models()
        my_nest.MyLoadModels( model_list, NEURON_MODELS )
        my_nest.MyLoadModels( model_list, [syn_model] )
        
        SNR = MyGroup( NEURON_MODELS[0], 1, mm=True,
                        mm_dt = 0.1 )
        
        voltage, ipsp=SNR.I_PSE(I_vec, synapse_model=syn_model, id=0, 
                                receptor=receptor)
        v_mat.append(voltage)
        ipsp_mat.append(ipsp)
        
    v_mat=numpy.array(v_mat)
    ipsp_mat=numpy.array(ipsp_mat) 
    
    return v_mat, ipsp_mat
        

print 'Simulation'   
    
# SIMULATION

# Example synaptic response silent
SNR_plot1, infoStringFig1=simulate_response_example_clamped_silent(SNR_INJECTED_CURRENT)

# Voltage against ipsp size
v_mat, ipsp_mat=simulate_voltage_ipsp(I_RANGE)

# DISPLAY
plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0, fontsize=16)
font_size_text = 8
fig = pylab.figure( facecolor = 'w' )

ax_list = []
ax_list.append( MyAxes(fig, [ .075, .37, .135, .26 ] ) )    # text box
ax_list.append( MyAxes(fig, [ .26,  .6,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .53, .6,  .165, .34 ] ) )    #  
#ax_list.append( MyAxes(fig, [ .8,   .6,  .165-0.06, .34-0.12 ] ) )    # 
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
plot_voltage_ipsp(ax, v_mat, ipsp_mat)

# Example synaptic response spiking
#ax=ax_list[3]
#plot_response_example_clamped_spiking(ax, SNR_plot2)

# Current vs ipsc 
#ax=ax_list[3]
#plot_voltage_ipsc(ax, voltage, ipsc)

pylab.show()

fig.savefig( PICTURE_PATH + '/' + FILE_NAME  + '.svg' , format = 'svg')
fig.savefig( PICTURE_PATH + '/' + FILE_NAME  + '.pdf' , format = 'pdf')