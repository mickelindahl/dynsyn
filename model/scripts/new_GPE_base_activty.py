import numpy
import pylab
import os
import sys
import time as ttime

# Get directory where model and code resides 
model_dir=   '/'.join(os.getcwd().split('/')[0:-1])    
code_dir=  '/'.join(os.getcwd().split('/')[0:-2])  

# Add model, code and current directories to python path
sys.path.append(os.getcwd())  
sys.path.append(model_dir)
sys.path.append(code_dir+'/nest_toolbox') 
from src import my_nest, misc, my_topology, plot_settings
from src.my_axes import MyAxes 

from simulation_utils import simulate_basa_line_GPe
from src import misc
import scipy.optimize as opt
OUTPUT_PATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]


# Assume n_stn=30, n_gpe=30, n_msn=500 and weight_msn
# Find: gpe_curr, STN-GPe weights and GPe-GPE.
#
# 1. 15 Hz (50 % decrease) without STN ( must put GPe collateral at 15 Hz then 
# 2. 46 ~Hz (55 % increase) without GABA (from MSN and GPe). (And STN=4 Hz since
#    Found that STN fired at 4 Hz when GPe increased with 55 in model simulation.
# 3. 30 Hz with all connections
# 
#
# Sol:x0=[54, 0.9, 0.16] 

neuron_model=['GPE_aeif']
syn_models=['MSN_GPE_gaba_p','STN_GPE_ampa_s','GPE_GPE_gaba_s' ]
msn_rate, stn_rate, gpe_rate=.1, 10.0, 30.0        
n_msn, n_stn, n_gpe =500, 30, 30

def plot_example(ax, GPE_target, sim_time, x, type):
    meanRate=round(GPE_target.signals['spikes'].mean_rate(1000,sim_time),1)
    spk=GPE_target.signals['spikes'].time_slice(1000,sim_time).raw_data()
    CV=numpy.std(numpy.diff(spk[:,0],axis=0))/numpy.mean(numpy.diff(spk[:,0],axis=0))

    
    GPE_target.signals['V_m'].my_set_spike_peak( 15, spkSignal= GPE_target.signals['spikes'] )          
    pylab.rcParams.update( {'path.simplify':False}    )
    
    GPE_target.signals['V_m'].plot(display=ax)
    ax.set_title( type+'\n'+str(meanRate)+ 'Hz, CV='+str(round(CV,4))+
                  ', curr='+str(round(x[0],4))+', w_GPE_GPE='+str(round(x[1],4))+ 
                  ' w_STN_GPE='+str(round(x[2],4)) )


def restriction_1(msn_rate, stn_rate, gpe_rate, n_msn, n_stn, n_gpe,x, neuron_model, syn_models):
    # 1. 15 Hz (50 % decrease) without STN (n_stn=0), OBS! this means that gpe_rate=15 Hz    
    c,w1,w2=x
    target_rate1=15.
    n_stn_ch=0
    gpe_rate_ch=15.0
    GPE_target=simulate_basa_line_GPe(msn_rate, stn_rate, gpe_rate_ch, n_msn, n_stn_ch, n_gpe, 
                               neuron_model, syn_models, c, # This is the current to add to I_e base
                               sim_time, 8, w_GPE_GPE=w1, w_STN_GPE=0) 

    e=round(GPE_target.signals['spikes'].mean_rate(1000,sim_time),1)-target_rate1
    return GPE_target, e

def restriction_2(msn_rate, stn_rate, gpe_rate, n_msn, n_stn, n_gpe,x, neuron_model, syn_models):
        # 2. 46 ~Hz (55 % increase) without GABA (from MSN and GPe). (And STN=4 Hz since
        #    Found that STN fired at 4 Hz when GPe increased with 55 in model simulation.
        c,w1,w2=x
        target_rate2=30*1.55
        n_gpe_ch=0
        n_msn_ch=0
        stn_rate_ch=4.
        GPE_target=simulate_basa_line_GPe(msn_rate, stn_rate_ch, gpe_rate, n_msn_ch, n_stn, n_gpe_ch, 
                               neuron_model, syn_models, c, # This is the current to add to I_e base
                               sim_time, 8, w_GPE_GPE=0, w_STN_GPE=w2 )
 
        e=round(GPE_target.signals['spikes'].mean_rate(1000,sim_time),1)-target_rate2
        
        return GPE_target, e

def restriction_3(msn_rate, stn_rate, gpe_rate, n_msn, n_stn, n_gpe,x, neuron_model, syn_models):
        # 3. 30 Hz with all connections 
        c,w1,w2=x
        target_rate3=30.0
        GPE_target=simulate_basa_line_GPe(msn_rate, stn_rate, gpe_rate, n_msn, n_stn, n_gpe, 
                               neuron_model, syn_models, c, # This is the current to add to I_e base
                               sim_time, 8, w_GPE_GPE=w1, w_STN_GPE=w2 )
   
        e=round(GPE_target.signals['spikes'].mean_rate(1000,sim_time),1)-target_rate3
        return GPE_target, e



def error_fun(x, sim_time):
        
        GPE_target, e1=restriction_1(msn_rate, stn_rate, gpe_rate, n_msn, n_stn, n_gpe,x, neuron_model, syn_models)
        GPE_target, e2=restriction_2(msn_rate, stn_rate, gpe_rate, n_msn, n_stn, n_gpe,x, neuron_model, syn_models)
        GPE_target, e3=restriction_3(msn_rate, stn_rate, gpe_rate, n_msn, n_stn, n_gpe,x, neuron_model, syn_models)

        print e1**2+e2**2+e3**2
        return e1**2+e2**2+e3**2
        
def fmin(load, save_at):
    sim_time=10000
  
    x0=[42, 1.3, 0.35] #[50, 0.9, 0.2]  #[current, w_GPE_GPE, w_STN_GPE]
    if not load:
        [xopt,fopt, iter, funcalls , warnflag, allvecs] = opt.fmin(error_fun, x0, args=([sim_time]), maxiter=50, maxfun=50,full_output=1, retall=1)

        misc.pickle_save([xopt,fopt, iter, funcalls , warnflag, allvecs], save_at)
    else:
        [xopt,fopt, iter, funcalls , warnflag, allvecs]=misc.pickle_load(save_at)        
    return xopt  



sim_time=10000
  

save_at=OUTPUT_PATH+'/simulate_network_fmin.plk' 
x=fmin(1,save_at)
#x=[53, 0.8855, 0.15]
#x=[42, 1.3, 0.35] 
GPE_target1, e1=restriction_1(msn_rate, stn_rate, gpe_rate, n_msn, n_stn, n_gpe,x, neuron_model, syn_models)
GPE_target2, e2=restriction_2(msn_rate, stn_rate, gpe_rate, n_msn, n_stn, n_gpe,x, neuron_model, syn_models)
GPE_target3, e3=restriction_3(msn_rate, stn_rate, gpe_rate, n_msn, n_stn, n_gpe,x, neuron_model, syn_models)

plot_settings.set_mode(pylab, mode='by_fontsize', w = 1000.0, h = 1000.0, fontsize=16)
font_size_text = 8
fig = pylab.figure( facecolor = 'w')

ax_list=[]
ax_list.append( MyAxes(fig, [ .1, .7, .8, .2 ] ) )    # text box
ax_list.append( MyAxes(fig, [ .1, .4, .8, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .1, .1, .8, .2 ] ) )    # 
 
ax=ax_list[0]
plot_example(ax, GPE_target1,sim_time, x, type='No STN')  

ax=ax_list[1]
plot_example(ax, GPE_target2, sim_time, x, type='No GABA and STN rate=4 Hz')  

ax=ax_list[2]
plot_example(ax, GPE_target3, sim_time, x, type='Normal')  
        
pylab.show()
# Knas fixa
