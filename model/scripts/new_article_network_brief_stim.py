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
from simulation_utils import simulate_network

OUTPUT_PATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]

def plot_example_firing_rate(ax, layer, name,  color='b', ylim=[]):
    time_bin=20


    signal=layer.signals['spikes']
    #signal.my_firing_rate(bin=time_bin, display=ax,
    #                      kwargs={'color':color})

    print name, 'CV:', numpy.mean(signal.cv_isi())
    print name, 'mean:', signal.mean_rate()
    print name, 'std:', signal.mean_rate_std()
    hist=signal.spike_histogram(time_bin=1, normalized=True)
    spk_mean=numpy.mean(hist, axis=0)
    spk_mean=misc.convolve(spk_mean, 20, 'triangle',single=True)[0]
    print spk_mean.shape
    time=numpy.arange(1,len(spk_mean)+1)
    ax.plot(time,spk_mean)
    
    
    if layer.id_mod: 
        print numpy.array(layer.id_mod)-layer.ids[0]-1
        spk_mean=numpy.mean(hist[numpy.array(layer.id_mod)-layer.ids[0],:], axis=0)
        spk_mean=misc.convolve(spk_mean, 20, 'triangle',single=True)[0]
        print spk_mean.shape
        time=numpy.arange(1,len(spk_mean)+1)
        ax.plot(time,spk_mean,'r')
        
    
    #spk_rate=misc.convolve(spk_mean[i,:], 50, 'triangle',single=True)[0]
    #plot()
    #meanRate=round(layer.signals['spikes'].mean_rate(1000,sim_time),1)
    #spk=SNR.signals['spikes'].time_slice(1000,sim_time).raw_data()
    #CV=numpy.std(numpy.diff(spk[:,0],axis=0))/numpy.mean(numpy.diff(spk[:,0],axis=0))

    ax.set_ylabel('Firing rate '+name+' (spikes/s)') 
    ax.set_xlabel('Time (ms)')    
    ax.my_set_no_ticks( yticks=6, xticks = 8 ) 
    ax.set_ylim(ylim)
    ax.set_xlim([450,750])
    
def plot_example_raster(ax, layer, name):
    global ADJUST_XDATA_MS
    
    
    layer.signals['spikes'].raster_plot(display=ax,
                                      kwargs={'color':'b', 'zorder':1})  

    if len(layer.id_mod):layer.signals['spikes'].raster_plot(display=ax, id_list=layer.id_mod,
                                      kwargs={'color':'r', 'zorder':1})  
 
    lines = ax.lines
    ax.set_ylabel(name+' id')
    ax.my_set_no_ticks( yticks=6, xticks = 5 )
    ax.set_ylim([layer.ids[0],layer.ids[-1]])
    
    '''
    for line in lines:
        line.set_xdata(line.get_xdata()-ADJUST_XDATA_MS)
    
    ax.set_xlim(misc.adjust_limit([0,1500]))
    ax.set_ylim(misc.adjust_limit([0,N_GPE]))
       
    
    ax.text( 0.05, 0.05, 'Non-pausing GPe neurons' , backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : 'b' })  
    ax.text( 0.05, 0.15, 'Pausing GPe neurons' , backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : 'r' }) 
    
    '''
        
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


params_msn_d1={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, 100, 0.1],
            'mod_times':[1,1000+15, 1000+25],  'n_mod':300}    
params_msn_d2={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, 100, 0.1],
            'mod_times':[1.,1000.+15, 1000.+25.],  'n_mod':300, 'focus':False, 
            'skip':1} 
params_stn={'rate':350., 'mod':True,'mod_rate':5000., 'mod_times':[1000., 1000.+10.]} 
sim_time=1200.
N_MSN=15000
#synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p_stoc']
synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_s_ref']
model_params={'neurons':{'MSN_D1':{'n':N_MSN},
                         'MSN_D2':{'n':N_MSN},
                         'GPE': {'paused':0.}}}

save_result_at=OUTPUT_PATH+'/simulate_network_brief_stim.plk'
if 1:
    layer_dic=simulate_network(params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, sim_time=sim_time, seed=1,
                           I_e_add={'SNR':300, 'STN':0,'GPE':30}, threads=4, 
                           start_rec=500.,model_params=model_params)    
    misc.pickle_save(layer_dic, save_result_at)  
else:
    layer_dic=misc.pickle_load(save_result_at)  
 
 #Inspect results
plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0+275.0, fontsize=12)
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
plot_example_raster(ax, layer_dic['MSN_D1'], 'MSN D1',)
ax=ax_list[2]
plot_example_raster(ax, layer_dic['MSN_D2'], 'MSN D2',)
ax=ax_list[3]
plot_example_firing_rate(ax, layer_dic['MSN_D1'], 'MSN D1',ylim=[0,35])
#plot_example_firing_rate(ax, layer_dic['MSN_D2'], 'MSN D2', color='r', ylim=[0,1])
ax=ax_list[4]
plot_example_raster(ax, layer_dic['SNR'], 'SNr')

ax=ax_list[5]
plot_example_raster(ax, layer_dic['GPE'], 'GPe')
ax=ax_list[6]
plot_example_raster(ax, layer_dic['STN'], 'STN')

ax=ax_list[7]
plot_example_firing_rate(ax, layer_dic['SNR'], 'SNr',ylim=[0,70])
ax=ax_list[8]
plot_example_firing_rate(ax, layer_dic['GPE'], 'GPe',ylim=[0,70])
ax=ax_list[9]
plot_example_firing_rate(ax, layer_dic['STN'], 'STN', ylim=[0,35])


ax=ax_list[0]
#plot_text(ax, info_string=s)
pylab.show()

