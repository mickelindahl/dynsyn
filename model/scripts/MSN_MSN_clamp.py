#! ===========
#! MSN to MSN
#! ===========


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

my_nest.ResetKernel()
  
model_list = models()                                         # Get model list
for model in model_list: 
    my_nest.CopyModel( model[ 0 ], model[ 1 ], model[ 2 ] )   # Create models
neuron_model='MSN_izh'
MSN = MyGroup( neuron_model, 3, mm_dt = 0.1 )

#! Spike train experiment 1. Train of 8 spikes at 20 Hz and the recovery spike 
#! at 550 ms as in Planert 2009
spike_times = range( 10, 430, 50 )                                                  
spike_times.extend( [ 430 + 550 ] )

# input
sgs = my_nest.Create('spike_generator',                               
                  params={'spike_times':[float(sp) for sp in spike_times]})

syn_model = 'MSN_MSN_gaba_s'
my_nest.Connect( sgs, [MSN[0]], model = syn_model )          # connect MSNs
T = 2000                                                     # simulation time
my_nest.Simulate( T )                                        # simulate
    
MSN.get_signal('v', 'V_m')  


pylab.close('all')                                          # display
# Create figure where figsize(width,height) and figure dimenstions window 
# width = figsize(width) x dpi and window hight = figsize(hight) x dpi
plot_settings.set_mode(mode='dynamic', w = 700.0, h = 400.0)
font_size_text = 10
fig = pylab.figure( facecolor = 'w' )
pylab.suptitle('MSN to MSN')

ax_list = []
ax_list.append( MyAxes(fig, [ .1, .37, .2,  .26 ] ) )     # text box
ax_list.append( MyAxes(fig, [ .50, .15, .40, .65 ] ) )       # voltage trace

ds=my_nest.GetDefaults( syn_model )
sn=my_nest.GetStatus( MSN[:] )[ 0 ]

ax         = ax_list[ 0 ]                                 # Text box
tb         = ''                                                
tb = tb + ' %6s %7s %3s \n' % ( 'Synapse model:',  ds[ 'synapsemodel'      ], ' ' )
tb = tb + ' \n '
tb = tb + ' %6s %7s %3s \n' % ( 'Delay:',          ds[ 'delay'             ], 'ms' )
tb = tb + ' %6s %7s %3s \n' % ( 'Weight:',         ds[ 'weight'            ], 'nS' )
tb = tb + ' \n '
tb = tb + ' %6s %7s %3s \n' % ( 'Decay:',          sn[ 'GABAA_2_Tau_decay' ], 'ms' )
tb = tb + ' %6s %7s %3s \n' % ( 'E_rev:',          sn[ 'GABAA_2_E_rev'     ], 'mV' )


ax.text( 0.9, 0.5, tb ,
         fontsize            = font_size_text,
         horizontalalignment = 'right',
         verticalalignment   = 'center',
         transform           = ax.transAxes,     # to define coordinates in scale
         **{ 'fontname' : 'monospace' })                           
ax.remove_axis( xaxis = True, yaxis = True )
ax.remove_spine(left=True,  bottom=True, right=True, top=True)

ax = ax_list[1]
MSN.signals['V_m'].plot( id_list = [ MSN[ 0 ] ], display = ax, 
                         kwargs = { 'color' : 'b' } )
ax.legend( numpoints = 1 )#, 'Theta'])
ax.set_ylabel( 'Membrane potential [mV]' )
ax.set_no_ticks(xticks=3, yticks=3)
pylab.show()


name = sys.argv[0].split('/')[-1].split('.')[0]
fig.savefig( picture_dir + '/' + name + '-'+ neuron_model + '.svg' , 
             dpi = 500, format = 'svg')