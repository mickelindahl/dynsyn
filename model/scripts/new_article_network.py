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
from simulation_utils import simulate_network, simulate_network_poisson

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
    spk_mean=misc.convolve(spk_mean, 100, 'triangle',single=True)[0]
    print spk_mean.shape
    time=numpy.arange(1,len(spk_mean)+1)
    ax.plot(time,spk_mean)
    
    
    if layer.id_mod: 
        print numpy.array(layer.id_mod)-layer.ids[0]-1
        spk_mean=numpy.mean(hist[numpy.array(layer.id_mod)-layer.ids[0],:], axis=0)
        spk_mean=misc.convolve(spk_mean, 100, 'triangle',single=True)[0]
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
    
def plot_example_raster(ax, layer, name):
    global ADJUST_XDATA_MS
    
    
    layer.signals['spikes'].raster_plot(display=ax,
                                      kwargs={'color':'b', 'zorder':1})  

    #if layer.id_mod:layer.signals['spikes'].raster_plot(display=ax, id_list=layer.id_mod,
    #                                  kwargs={'color':'r', 'zorder':1})  
 
    lines = ax.lines
    ax.set_ylabel(name+' id')
    ax.my_set_no_ticks( yticks=6, xticks = 5 )
    ax.set_ylim([layer.ids[0],layer.ids[-1]])
       
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

def simulate_poisson_15000(params_msn_d1, params_msn_d2, params_stn,
                                       synapse_models, sim_time, seed,
                                       I_e_add, threads, start_rec):
    N_MSN=15000
    model_params={'misc':{'N_MSN':N_MSN},
              'conns':{'MSN_D1_SNR':{'p':500./float(N_MSN)},
                       'MSN_D2_GPE':{'p':500./float(N_MSN),
                                     'lines':False}},
              'neurons':{'MSN_D1':{'n':0},
                         'MSN_D2':{'n':0},
                         'MSN_D1_bg':{'n':300, 'lesion':False},
                         'MSN_D2_bg':{'n':300, 'lesion':False},  
                         'GPE': {'paused':False}}}
    layer_dic=simulate_network_poisson(params_msn_d1, params_msn_d2, params_stn,
                                       synapse_models, sim_time, seed,
                                       I_e_add, threads, start_rec, 
                                       model_params)    

    mr=[]
    signal=layer_dic['SNR'].signals['spikes']
    m_hist=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
    mr.append(numpy.mean(m_hist, axis=0))
    signal=layer_dic['GPE'].signals['spikes']
    m_hist=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
    mr.append(numpy.mean(m_hist, axis=0))
    signal=layer_dic['STN'].signals['spikes']
    m_hist=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
    mr.append(numpy.mean(m_hist, axis=0))
    
    return layer_dic, mr

def simulate_15000(params_msn_d1, params_msn_d2, params_stn, synapse_models, sim_time, seed,
                                       I_e_add, threads, start_rec):
    N_MSN=15000
    model_params={'misc':{'N_MSN':N_MSN},
              'conns':{'MSN_D2_GPE':{'lines':False}},
             'neurons':{'MSN_D1':{'n':N_MSN},
                         'MSN_D2':{'n':N_MSN},
                         'MSN_D1_bg':{'n':0},
                         'MSN_D2_bg':{'n':0},  
                         'GPE': {'paused':False}}}
    layer_dic=simulate_network(params_msn_d1, params_msn_d2, params_stn,
                                       synapse_models, sim_time, seed,
                                       I_e_add, threads, start_rec, 
                                       model_params)    

    mr=[]
    signal=layer_dic['SNR'].signals['spikes']
    m_hist=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
    mr.append(numpy.mean(m_hist, axis=0))
    signal=layer_dic['GPE'].signals['spikes']
    m_hist=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
    mr.append(numpy.mean(m_hist, axis=0))
    signal=layer_dic['STN'].signals['spikes']
    m_hist=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
    mr.append(numpy.mean(m_hist, axis=0))
    
    return layer_dic, mr

def simulate_1500(params_msn_d1, params_msn_d2, params_stn, synapse_models, sim_time, seed,
                                       I_e_add, threads, start_rec):
    N_MSN=1500
    model_params={'misc':{'N_MSN':N_MSN},
              'conns':{'MSN_D2_GPE':{'lines':False},
                       'GPE_STN':{'lesion':False},
                       'GPE_GPE':{'lesion':False},
                       'GPE_SNR':{'lesion':False},
                       'STN_GPE':{'lesion':False},
                       'MSN_D1_SNR':{'lesion':False},
                       'MSN_D2_GPE':{'lines':False, 
                                     'lesion':False},
                       'STN_SNR':{'lesion':False}},
              'neurons':{'MSN_D1':{'n':N_MSN},
                         'MSN_D2':{'n':N_MSN},
                         'MSN_D1_bg':{'n':0},
                         'MSN_D2_bg':{'n':0},  
                         'GPE': {'paused':False}}}

    layer_dic=simulate_network(params_msn_d1, params_msn_d2, params_stn,
                                       synapse_models, sim_time, seed,
                                       I_e_add, threads, start_rec, 
                                       model_params)    

    mr=[]
    signal=layer_dic['SNR'].signals['spikes']
    m_hist=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
    mr.append(numpy.mean(m_hist, axis=0))
    signal=layer_dic['GPE'].signals['spikes']
    m_hist=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
    mr.append(numpy.mean(m_hist, axis=0))
    signal=layer_dic['STN'].signals['spikes']
    m_hist=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
    mr.append(numpy.mean(m_hist, axis=0))
    
    return layer_dic, mr

def simulate_poisson_15000_eval(load, save_at, n_exp, params_msn_d1, params_msn_d2, params_stn,
                                       synapse_models, sim_time, seed,
                                       I_e_add, threads, start_rec):

    rates=[]
    if not load:
        for i in range(n_exp):
            seed=i
            layer_dic, r=simulate_poisson_15000(params_msn_d1, params_msn_d2, params_stn,
                                       synapse_models, sim_time, seed,
                                       I_e_add, threads, start_rec)
            rates.append(r)
        
        rates=numpy.array(rates)    


        misc.pickle_save(rates,save_at)   
    else:
        rates=misc.pickle_load(save_at)   
    
    mr=numpy.mean(rates, axis=0)
    stdr=numpy.std(rates, axis=0)
    
    return mr, stdr

def simulate_15000_eval(load, save_at, n_exp, params_msn_d1, params_msn_d2, params_stn,
                                       synapse_models, sim_time, seed,
                                       I_e_add, threads, start_rec):

    rates=[]
    if not load:
        for i in range(n_exp):
            seed=i
            layer_dic, r=simulate_15000(params_msn_d1, params_msn_d2, params_stn,
                                       synapse_models, sim_time, seed,
                                       I_e_add, threads, start_rec)
            rates.append(r)
        
        rates=numpy.array(rates)    


        misc.pickle_save(rates,save_at)   
    else:
        rates=misc.pickle_load(save_at)   
    
    mr=numpy.mean(rates, axis=0)
    stdr=numpy.std(rates, axis=0)
    
    return mr, stdr

def simulate_1500_eval(load, save_at, n_exp, params_msn_d1, params_msn_d2, params_stn,
                                       synapse_models, sim_time, seed,
                                       I_e_add, threads, start_rec):

    rates=[]
    if not load:
        for i in range(n_exp):
            seed=i
            layer_dic, r=simulate_1500(params_msn_d1, params_msn_d2, params_stn,
                                       synapse_models, sim_time, seed,
                                       I_e_add, threads, start_rec)
            rates.append(r)
        
        rates=numpy.array(rates)    


        misc.pickle_save(rates,save_at)   
    else:
        rates=misc.pickle_load(save_at)   
    
    mr=numpy.mean(rates, axis=0)
    stdr=numpy.std(rates, axis=0)
    
    return mr, stdr
br=0.1
bg_rate=br*500.0
params_msn_d1={'base_rates':[br], 'base_times':[1.], 'mod_rates': [],
            'mod_times':[],  'n_mod':0, 
            'bg_rate':[bg_rate]}    
params_msn_d2={'base_rates':[br], 'base_times':[1.], 'mod_rates': [],
            'mod_times':[],  'n_mod':0, 'focus':False, 
            'skip':1, 'bg_rate':[bg_rate]} 
params_stn={'rate':219., 'mod':False,'mod_rate':0., 'mod_times':[]} 

sim_time=4500.
seed=1
start_rec=500.0
#synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p_stoc']
#synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_s_ref']
synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p']
n_exp=5
threads=4


save_at=OUTPUT_PATH+'/simulate_poisson_15000_eval.plk'
mr, std=simulate_poisson_15000_eval(1, save_at, n_exp, params_msn_d1, params_msn_d2, params_stn,
                                       synapse_models, sim_time, seed,
                                       {}, threads, start_rec)
save_at=OUTPUT_PATH+'/simulate_15000_eval.plk'
mr2, std2=simulate_15000_eval(1, save_at, n_exp, params_msn_d1, params_msn_d2, params_stn,
                                       synapse_models, sim_time, seed,
                                       {}, threads, start_rec)
save_at=OUTPUT_PATH+'/simulate_1500_eval.plk'
mr3, std3=simulate_1500_eval(1, save_at, n_exp, params_msn_d1, params_msn_d2, params_stn,
                                       synapse_models, sim_time, seed,
                                       {}, threads, start_rec)


save_result_at=OUTPUT_PATH+'/simulate_poisson_15000.plk2'
if 1:
    layer_dic, mr=simulate_poisson_15000(params_msn_d1, params_msn_d2, params_stn,
                                       synapse_models, sim_time, seed,
                                       {}, threads, start_rec)
    misc.pickle_save([layer_dic, mr], save_result_at)  
else:
    layer_dic, mr=misc.pickle_load(save_result_at)  



save_result_at=OUTPUT_PATH+'/simulate_15000.plk'
if 1:
    layer_dic2, mr=simulate_15000(params_msn_d1, params_msn_d2, params_stn, synapse_models, sim_time, seed,
                                       {}, threads, start_rec)
    misc.pickle_save([layer_dic2, mr], save_result_at)  
else:
    layer_dic2, mr=misc.pickle_load(save_result_at)  


save_result_at=OUTPUT_PATH+'/simulate_1500.plk2'+nest.version()
if 1:
    layer_dic3, mr=simulate_1500(params_msn_d1, params_msn_d2, params_stn, synapse_models, sim_time, seed,
                                       {}, threads, start_rec)
    misc.pickle_save([layer_dic3, mr], save_result_at)  
else:
    layer_dic3, mr=misc.pickle_load(save_result_at)  

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
plot_example_firing_rate(ax, layer_dic['SNR'], 'SNr',ylim=[0,90])
ax=ax_list[8]
plot_example_firing_rate(ax, layer_dic['GPE'], 'GPe',ylim=[0,90])
ax=ax_list[9]
plot_example_firing_rate(ax, layer_dic['STN'], 'STN', ylim=[0,35])


ax=ax_list[0]
#plot_text(ax, info_string=s)

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

layer_dic=layer_dic2
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
plot_example_firing_rate(ax, layer_dic['SNR'], 'SNr',ylim=[0,90])
ax=ax_list[8]
plot_example_firing_rate(ax, layer_dic['GPE'], 'GPe',ylim=[0,90])
ax=ax_list[9]
plot_example_firing_rate(ax, layer_dic['STN'], 'STN', ylim=[0,35])
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

layer_dic=layer_dic3
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
plot_example_firing_rate(ax, layer_dic['SNR'], 'SNr',ylim=[0,90])
ax=ax_list[8]
plot_example_firing_rate(ax, layer_dic['GPE'], 'GPe',ylim=[0,90])
ax=ax_list[9]
plot_example_firing_rate(ax, layer_dic['STN'], 'STN', ylim=[0,35])

print 'poisson 15000',mr, std
print '15000',mr2, std2
print '1500',mr3, std3

pylab.show()

