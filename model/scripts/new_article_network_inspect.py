import numpy
import pylab
import os
import sys
import time as ttime
import pprint
import nest
# Get directory where model and code resides 
model_dir=   '/'.join(os.getcwd().split('/')[0:-1])    
code_dir=  '/'.join(os.getcwd().split('/')[0:-2])  

# Add model, code and current directories to python path
sys.path.append(os.getcwd())  
sys.path.append(model_dir)
sys.path.append(code_dir+'/nest_toolbox') 
from src import my_nest, misc, my_topology, plot_settings
from src.my_axes import MyAxes 
import nest.topology as tp
from simulation_utils import inspect_network



def plot_layer(source, target,type, layer_dic):

    layer_dic[target].plot(ax=ax, nodesize=10,nodecolor='b')
    n=layer_dic[target].plot_targets(source_layer=layer_dic[source], type=type,
                                  ax=ax, src_color='red', src_size=50, 
                                  tgt_color='red', tgt_size=5)
    m, std=layer_dic[target].stat_connections(source_layer=layer_dic[source], type=type)
    ax.set_title(source+'-'+target+': ' +str(n))

    
    n_source=len(layer_dic[source].ids)
    n_target=len(layer_dic[target].ids)
    
    
    s='\n'
    s=s+target+' layer:\n'
    s=s+source+':\n'
    s = s + ' %s %5s %3s \n' % ( 'Conn', str ( round(m,1) ) + r' $\pm$'+ str(round(std,1)),  '#' )  
    s = s + ' %s %5s %3s \n' % ( 'source per target', str ( round(n_source*m/n_target,1) )+ r'+ $\pm$'+ str(round(n_source*std/n_target,1)),  '#' )  
    return s

def plot_text(ax, info_string=''):
    

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
    
layer_dic=inspect_network()    
#Inspect network
 


plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 650.0, fontsize=12)
font_size_text = 8
fig = pylab.figure( facecolor = 'w' )
ax_list = []
ax_list.append( MyAxes(fig, [ .075, .37, .135, .26 ] ) )    # text box
ax_list.append( MyAxes(fig, [ .26,  .7,  .165, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .53,  .7,  .165, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .7,  .165, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .26,  .4,  .165, .2 ] ) )    #     
ax_list.append( MyAxes(fig, [ .53,  .4,  .165, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .4,  .165, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .26,  .1,  .165, .2 ] ) )    #     
ax_list.append( MyAxes(fig, [ .53,  .1,  .165, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .1,  .165, .2 ] ) )    # 


ax=ax_list[1]
s=''
s+=plot_layer('MSN_D1', 'SNR','gaba',layer_dic)
ax=ax_list[2]
s+=plot_layer('STN', 'SNR','gaba',layer_dic)
ax=ax_list[3]
s+=plot_layer('GPE', 'STN','gaba',layer_dic)
ax=ax_list[4]
#s+=plot_layer('MSN_D2', 'GPE','gaba',layer_dic)
ax=ax_list[5]
s+=plot_layer('STN', 'GPE','ampa',layer_dic)
ax=ax_list[6]
s+=plot_layer('GPE', 'GPE','gaba',layer_dic)
ax=ax_list[7]
s+=plot_layer('GPE', 'SNR','gaba',layer_dic)

ax=ax_list[0]
plot_text(ax, info_string=s)
pylab.show()

