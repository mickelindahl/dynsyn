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
import scipy.optimize as opt

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
        
    ax.set_ylabel('Firing rate '+name+' (spikes/s)') 
    ax.set_xlabel('Time (ms)')    
    ax.my_set_no_ticks( yticks=6, xticks = 8 ) 
    ax.set_ylim(ylim)
    
def plot_example_raster(ax, layer, name):
    global ADJUST_XDATA_MS
    
    
    layer.signals['spikes'].raster_plot(display=ax,
                                      kwargs={'color':'b', 'zorder':1})  

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

def simulate(params_msn_d1, params_msn_d2, params_stn, synapse_models, sim_time, seed,
                               I_e_add, threads,  start_rec, model_params, p_weights):
    
    layer_dic=simulate_network_poisson(params_msn_d1, params_msn_d2, params_stn,
                               synapse_models, sim_time, seed, I_e_add, threads, 
                               start_rec, model_params, p_weights=p_weights)    

        
    signal=layer_dic['SNR'].signals['spikes']
    m_hist=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
    mr=numpy.mean(m_hist[1000:], axis=0)
    
    signal=layer_dic['GPE'].signals['spikes']
    m_hist=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
    r_GPE=numpy.mean(m_hist[1000:], axis=0)
    
    signal=layer_dic['STN'].signals['spikes']
    m_hist=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
    r_STN=numpy.mean(m_hist[1000:], axis=0)

    
    return mr,  r_GPE, r_STN

def error_fun(x, n_exp, r_target, params_msn_d1, params_msn_d2, params_stn, sim_time, 
                               I_e_add, threads,  start_rec, model_params, p_weights):

    synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_s_ref']
    
    p_weights=p_weights=numpy.ones(17)
    p_GPE_SNR=16
    p_weights[p_GPE_SNR]=x
    mrs=[]
    for i in range(n_exp):
        seed=i
        mr, r_GPE, r_STN=simulate(params_msn_d1, params_msn_d2, params_stn, synapse_models, sim_time, seed,
                               I_e_add, threads,  start_rec, model_params, p_weights)
        mrs.append(mr)
    
    e=mr-r_target     

    return e**2

def fmin(load, save_at, x0, n_exp,  r_target, params_msn_d1, params_msn_d2, params_stn, sim_time,
                               I_e_add, threads,  start_rec, model_params, p_weights):
    

     #[current, w_GPE_STN]
    args=(n_exp,  r_target, params_msn_d1, params_msn_d2, params_stn, sim_time,
                               I_e_add, threads,  start_rec, model_params, p_weights )
    if not load:
        [xopt,fopt, iter, funcalls , warnflag, allvecs] = opt.fmin(error_fun, x0, args=args, maxiter=20, maxfun=20,full_output=1, retall=1)

        misc.pickle_save([xopt,fopt, iter, funcalls , warnflag, allvecs], save_at)
    else:
        [xopt,fopt, iter, funcalls , warnflag, allvecs]=misc.pickle_load(save_at)        
    return xopt,  fopt   
    
N_MSN=15000.0
bg_rate=0.1*500.0
params_msn_d1={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [],
            'mod_times':[],  'n_mod':0, 
            'bg_rate':[bg_rate]}    
params_msn_d2={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [],
            'mod_times':[],  'n_mod':0, 'focus':False, 
            'skip':1, 'bg_rate':[bg_rate]} 
params_stn={'rate':219., 'mod':False,'mod_rate':0., 'mod_times':[]} 
sim_time=4000.
start_rec=1000.
threads=4
model_params={'misc':{'N_MSN':N_MSN},
              'conns':{'MSN_D1_SNR':{'p':500./float(N_MSN)},
                       'MSN_D2_GPE':{'p':500./float(N_MSN),
                                     'lines':False}},
              'neurons':{'MSN_D1':{'n':0},
                         'MSN_D2':{'n':0},
                         'MSN_D1_bg':{'n':300,
                                      'lesion':False},
                         'MSN_D2_bg':{'n':300,
                                      'lesion':False},  
                         'GPE': {'paused':False}}}

synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p']
save_result_at=OUTPUT_PATH+'/simulate_target.plk'
p_weights=numpy.ones(17)
n_exp=5 

if 0:
    r_target=[]
    for i in range(n_exp):
        seed=i
        r_SNR, r_GPE, r_STN=simulate(params_msn_d1, params_msn_d2, params_stn, synapse_models, sim_time, seed,
                               {}, threads,  start_rec, model_params, p_weights)    
    
        r_target.append(r_SNR)
    misc.pickle_save([r_target, r_GPE, r_STN], save_result_at)  
else:
    r_target, r_GPE, r_STN=misc.pickle_load(save_result_at)  

r_target=numpy.array(r_target)
mr_target=numpy.mean(r_target, axis=0)
# Will find the value that GPE_SNR_ref synapse should be multiplied with 
# so that it have the same effect in the network as the plastic at base-level
# network activty    
x0=1.0
save_at=OUTPUT_PATH+'/simulate_network_fmin'+str(n_exp)+'.plk' 
x, e=fmin(0, save_at, x0, n_exp,  mr_target, params_msn_d1, params_msn_d2, params_stn, sim_time,
                               {}, threads,  start_rec, model_params, p_weights)



print x, e, r_target, r_GPE, r_STN
