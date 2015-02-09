#! MSN to SNR static and plastic synapses
#! =======================================
#! In paper by Connelly 2010 short term plasticity of SNR->SNR and MSN->SNR are 
#! investigated. Inspired by this we here set out to create facilitating 
#! MSN->SNR synapses. Since Connelly does not provide parameters for Tsodyks 
#! short-term-plasticity synapse we decided to use facilitating parameters from
#! Markram 1998. It turns out that the parameters produce similare facilitating
#! response as shown in Connelly 2010. Here Connelly 2010 data are reproduced 
#! using Markram 1998 parameters. Synaptic strength will be tuned to network, 
#! or to 2-5 mV when hyperpolarised 10 mV below reversal potential.

#! Data Connelly 2010
#! * Synsptic time constant 5.2
#! * Max curr approximatly 4-5x200 pA = 800-1000 pA. Enough to silance the cell they conclude
#! * Facilitation recovery time constant 2000 ms
#! * Pn/P1 10 hz, 50 and 100hz ~4. 
#! * Pair pulse ration 4 when fastes
#! * Held at -70
#! * Reversal potential 5+-4, that is -65+-4
#!
#!
#! Connelly simulated at 
#! 
#! 1. Stimulation at 10 Hz
#!
#! #. Stimulation at 50 Hz
#!
#! #. Stimulation at 100 Hz


                                  
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

 
# Add directories to python path
sys.path.append(os.getcwd())                            
parent_dir='/'.join(os.getcwd().split('/')[0:-1])       
                   
model_dir=   '/'.join(os.getcwd().split('/')[0:-1])        
code_dir=  '/'.join(os.getcwd().split('/')[0:-2]) 
picture_dir=  '/'.join(os.getcwd().split('/')[0:-3]) + '/pictures'     
                
sys.path.append(model_dir) 
sys.path.append(code_dir+'/nest_toolbox') 
SPATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]
SNAME  = sys.argv[0].split('/')[-1].split('.')[0]

from model_params import models                               # Then import models     
from src import misc, my_nest, my_signals, plot_settings
from src.my_population import MyGroup 
from src.my_axes import MyAxes 
#! Neuron
#! =======
neuron = 'SNR_aeif'
letter = 'A'

model_list = models()                                                           # get model list
for model in model_list: my_nest.CopyModel( model[ 0 ], model[ 1 ], model[ 2 ] )   # create models
I_e=-43.5
SNR = MyGroup( neuron, 3 , mm_dt = 0.1, params={'I_e':I_e})                                         # create neurons, recording interval dt = 0.1 to get exact current readings



syn='MSN_SNR_gaba_p1'
#! Configurables
#! =============
# SNR_SNR_g_gaba = [1.4, 10.1]  (Sims et al. 2008)
   
plastic_w = round( my_nest.GetDefaults( syn )[ 'weight' ] , 1)       


#! Input
#! =====                  
#! Train at 3, 10, 50 and 100 Hz  and then a recovery spike at 

keys       = ['10 Hz', '50 hz', '100 hz' ]                                      # stimulation frequencies
nb_spikes  = [   50.,     50.,     50. ]                                        # number of spikes for each stimulation
shift      = 50.                                                                # first spike          
spikes     = {}                                                                 # store spike time lists
sg         = {}                                                                 # store spike generators
rec1, rec2 = 500, 1500
for key, nb in zip( keys, nb_spikes ):
    hz             = float( key.split()[ 0 ] )
    spikes[ key ]  = numpy.linspace( shift, shift + nb*1000/hz, nb+1 )            # create frequency train
    rec_spk=numpy.array([rec1, rec2])+max(spikes[ key ])
    
    spikes[ key ]=numpy.append(spikes[ key ],rec_spk) 
    sg[ key ]      = my_nest.Create( 'spike_generator', params = { 'spike_times' : spikes[ key ] } )

#! Connections
#! ===========
for i, key in enumerate(keys):
    my_nest.Connect( sg[ key ],    [SNR.ids[ i ] ], params = {'weight' : 1. }, model =syn)

#! Simulate
#! ========
T = int( spikes[ keys[ 0 ] ][ -1 ] + 2000 )                                     # simulation time
my_nest.Simulate( T )                                                              # simulate


#! Plot
#! ====
SNR.get_signal( 'v','g_GABAA_1', stop=T ) # retrieve signal

delay      = my_nest.GetDefaults(syn)['delay']
cond_peaks = {}
for id, analog_signal in SNR.signals['g_GABAA_1'].analog_signals.iteritems(): # find conductance peaks

    cond_peaks[ keys[ id-1 ] ] = []
    for i in spikes[ keys[ id-1 ] ]:
        i=int(10*i)+delay*10
        valley = numpy.min( analog_signal.signal[ i - 20 : i + 10 ] )
        cond_peaks[ keys[ id-1 ] ].append( numpy.max( analog_signal.signal[ i - 20 : i + 10 ]-valley ) )


# Create figure where figsize(width,height) and figure dimenstions window 
# width = figsize(width) x dpi and window hight = figsize(hight) x dpi
plot_settings.set_mode(mode='by_fontsize', w = 750.0, h = 400.0, fontsize=12)
font_size_text = 4
fig = pylab.figure( facecolor = 'w' )

#fig = pylab.figure( figsize = ( 20, 7 ), dpi = 50, facecolor = 'w' )

ax_list = []
ax_list.append( pylab.axes( [ .05, .37, .08, .43 ] ) )                          # text box
ax_list.append( pylab.axes( [ .27, .72, .28, .08 ] ) )                          # conductance
ax_list.append( pylab.axes( [ .27, .47, .28, .08 ] ) )                          # conductance
ax_list.append( pylab.axes( [ .27, .2,  .28, .08 ] ) )                          # conductance
ax_list.append( pylab.axes( [ .69, .2,  .28, .6  ] ) )                          # p1/p#


ax = ax_list[ 0 ]                                                               # Text box
tb = ''                                                
tb = tb + ' %s %5s %3s \n' % ( 'plastic w:',         str ( plastic_w ),                                           'pS' )
tb = tb + '\n '                             
tb = tb + ' %s %5s %3s \n' % ( 'E rev:',             str ( my_nest.GetStatus(SNR.ids,'GABAA_1_E_rev')[0] ),          'mV' )
tb = tb + ' %s %5s %3s \n' % ( 'Tau_decay:',         str ( my_nest.GetStatus(SNR.ids,'GABAA_1_Tau_decay')[0] ),      'mV' )  
tb = tb + '\n '                             
tb = tb + ' %s %5s %3s \n' % ( 'U:',                 str ( my_nest.GetDefaults(syn)['U'] ),             '  ' ) 
tb = tb + ' %s %5s %3s \n' % ( 'tau_fac:',           str ( my_nest.GetDefaults(syn)['tau_fac'] ),       'ms' )                                 
tb = tb + ' %s %5s %3s \n' % ( 'tau_rec:',           str ( my_nest.GetDefaults(syn)['tau_rec'] ),       'ms' )
tb = tb + ' %s %5s %3s \n' % ( 'tau_psc:',           str ( my_nest.GetDefaults(syn)['tau_psc'] ),       'ms' )     
tb = tb + '\n '                             
tb = tb + ' %s %5s %3s \n' % ( 'recovery spike 1:',  str ( rec1 ),                                                'ms' )
tb = tb + ' %s %5s %3s \n' % ( 'recovery spike 2:',  str ( rec2 ),                                                'ms' )
tb = tb + '\n '  
tb = tb + ' %s %5s %3s \n' % ( keys[ 0 ] + ' first:', str ( round( cond_peaks[ keys[ 0 ] ][ 0 ],1) ),             'nS' )
tb = tb + ' %s %5s %3s \n' % ( keys[ 0 ] + ' last:',  str ( round( cond_peaks[ keys[ 0 ] ][ 9 ],1) ),             'nS' )
tb = tb + ' %s %5s %3s \n' % ( keys[ 1 ] + ' first:', str ( round( cond_peaks[ keys[ 1 ] ][ 0 ],1) ),             'nS' )
tb = tb + ' %s %5s %3s \n' % ( keys[ 1 ] + ' last:',  str ( round( cond_peaks[ keys[ 1 ] ][ 9 ],1) ),             'nS' )

ax.text( 0.8, 0.5, tb ,
         fontsize            = font_size_text,
         horizontalalignment = 'right',
         verticalalignment   = 'center',
         transform           = ax.transAxes,                                    # to define coordinates in scale
         **{ 'fontname' : 'monospace' })                           
#ps.push_remove_axis( ax, xaxis = 'remove', yaxis = 'remove' )
ax.set_title('MSN-SNR')

ax   = ax_list[ 1 ]                                                             # Conductances
SNR.signals[ 'g_GABAA_1' ].plot( id_list=[ SNR.ids[ 0 ] ], display=ax ) 
#ax.set_xlim([0,1250])
ax.text( -0.3, 1.75, letter+'i',fontsize = 20, style = 'oblique', transform = ax.transAxes) # graph numbering A, B,...
#ps.push_remove_axis( ax, xaxis = 'push', yaxis = 'push', nb_yticks = 2 )
ax.set_ylabel('')

ax   = ax_list[ 2 ]                                                             # signalsuctances
SNR.signals[ 'g_GABAA_1' ].plot( id_list=[ SNR.ids[ 1 ] ], display=ax,
                                   kwargs = {'color':'g'}  ) 
#ax.set_xlim([0,350])
#ps.push_remove_axis( ax, xaxis = 'push', yaxis = 'push', nb_yticks = 2 )

ax   = ax_list[ 3 ]                                                             # Conductances
SNR.signals[ 'g_GABAA_1' ].plot( id_list=[ SNR.ids[ 2 ] ], display=ax,
                                  kwargs = {'color':'r'} ) 
#ax.set_xlim([0,200])
#ps.push_remove_axis( ax, xaxis = 'push', yaxis = 'push', nb_yticks = 2 )
ax.set_ylabel('')

ax   = ax_list[ 4 ]                                                             # Plot p1/p#
for key in keys:
    val = cond_peaks[key]
    ax.plot( [v/val[0] for v in val] )                            
ax.set_xlabel('Pulse (#)') 
ax.set_ylabel('Pn/P1')   
ax.legend(keys, loc='best')
#ps.push_remove_axis( ax, xaxis = 'push', yaxis = 'push' )
ax.text( -0.3, 1.1, letter+'ii' ,fontsize = 20, style = 'oblique', transform = ax.transAxes) # graph numbering A, B,...


path = '/home/mikael/activity-phd/projects/igor10/dev/lines10/pictures'
#fig.savefig( path + '/SNR.eps' , dpi=300,papertype='a4',format='eps',transparent=True)
#name = sys.argv[0].split('/')[-1].split('.')[0]
#fig.savefig( path + '/'+name +'.png' , dpi = 500, format = 'png')


pylab.show()
