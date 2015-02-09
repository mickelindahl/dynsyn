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
from simulation_utils import simulate_network_direct_sep_freq


model_name=os.getcwd().split('/')[-2]
picture_dir='/'.join(os.getcwd().split('/')[0:-3]) + '/pictures/'+model_name 
SNAME  = sys.argv[0].split('/')[-1].split('.')[0]

OUTPUT_PATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]

def plot_example_firing_rate(ax, spk, name,  color=['b','g','r'], ylim=[]):
    
    
    for col, sp in zip(color, spk):
        time=numpy.arange(1,len(sp)+1)
        sp=misc.convolve(sp, 100, 'triangle',single=True)[0]
        ax.plot(time,sp, **{'color':col})

    #spk_rate=misc.convolve(spk_mean[i,:], 50, 'triangle',single=True)[0]
    #plot()
    #meanRate=round(layer.signals['spikes'].mean_rate(1000,sim_time),1)
    #spk=SNR.signals['spikes'].time_slice(1000,sim_time).raw_data()
    #CV=numpy.std(numpy.diff(spk[:,0],axis=0))/numpy.mean(numpy.diff(spk[:,0],axis=0))

    ax.set_ylabel('Firing rate '+name+' (spikes/s)') 
    ax.set_xlabel('Time (ms)')    
    ax.my_set_no_ticks( yticks=6, xticks = 8 ) 
    ax.set_ylim(ylim)
    
def plot_example_raster(ax, layer, name):
    global ADJUST_XDATA_MS
    global N_GPE
    global N_SEL
    
    GPE=layer
    
    GPE.signals['spikes'].raster_plot(display=ax,
                                      kwargs={'color':'b', 'zorder':1})  

    
    lines = ax.lines
    ax.set_ylabel(name+' id')
    ax.my_set_no_ticks( yticks=6, xticks = 5 )
    
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




params_msn_d1={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, 20, 0.1],
            'mod_times':[1,1000, 1000+500],  'n_mod':0}    
params_msn_d2={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, 20, 0.1],
            'mod_times':[1,1000, 1000+500],  'n_mod':0} 
params_stn={'rate':350., 'mod':False,'mod_rate':0., 'mod_times':[1000., 1000.+500.]} 
synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p']


seed=[1,1,1]
sim_time=2000.
freq=[15., 30, 45.]
mod=[150,150,150]
#synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p_stoc']



save_result_at=OUTPUT_PATH+'/simulate_network_direct_sep_freq.plk'
if 1:
    spk=simulate_network_direct_sep_freq(mod, freq, params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, sim_time=sim_time, seed=seed,
                           I_e_add={'SNR':200, 'STN':0,'GPE':30}, threads=4, start_rec=500.)
    misc.pickle_save(spk, save_result_at)
else: 
    spk=misc.pickle_load(save_result_at)
 
 
mod=[150,int(150*15./30.),int(150*15./45.)]
save_result_at=OUTPUT_PATH+'/simulate_network_direct_sep_freq_const_syn_event.plk'
if 1:
    spk2=simulate_network_direct_sep_freq(mod, freq, params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, sim_time=sim_time, seed=seed,
                           I_e_add={'SNR':200, 'STN':0,'GPE':30}, threads=4, start_rec=500.)
    misc.pickle_save(spk2, save_result_at)
else: 
    spk2=misc.pickle_load(save_result_at) 
 


synapse_models_msn=['MSN_SNR_gaba_s_min','MSN_SNR_gaba_s_max','MSN_SNR_gaba_p1']
mods=numpy.arange(50,1000,50 )
freq=numpy.arange(1,47,1)
save_result_at=OUTPUT_PATH+'/simulate_network_direct_combo.plk'
if 1:
    spk3={}
    
    
    for syn in synapse_models_msn:
        spk3[syn]={}
        for n_mod in mods: 
            synapse_models=[syn, 'GPE_SNR_gaba_p']
            
            mod=numpy.ones(len(freq))*n_mod
            seed=range(len(freq))
        
            spk3[syn][n_mod]=simulate_network_direct_sep_freq(mod, freq, params_msn_d1, params_msn_d2, params_stn,
                                   synapse_models, sim_time=sim_time, seed=seed,
                                   I_e_add={'SNR':280, 'STN':1,'GPE':20}, threads=4, start_rec=500.)
    misc.pickle_save(spk3, save_result_at)

else: 
    spk3=misc.pickle_load(save_result_at) 
 
spikes={}
for key1 in spk3.keys():
 
    spikes=spk3[key1].values()
    mods=spk3[key1].keys()
    

 #Inspect results
plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0+225.0, fontsize=16)
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
plot_example_firing_rate(ax, spk, 'SNr',ylim=[0,40])

ax=ax_list[2]
plot_example_firing_rate(ax, spk2, 'SNr',ylim=[0,40])

ax=ax_list[0]
#plot_text(ax, info_string=s)
pylab.show()

# dpi does not matter since svg and pdf are both vectorbased
fig.savefig( picture_dir + '/' + SNAME  + '.svg', format = 'svg') 
fig.savefig( picture_dir + '/' + SNAME  + '.pdf', format = 'pdf')