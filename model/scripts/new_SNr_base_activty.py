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

from simulation_utils import simulate_basa_line_SNr, simulate_basa_line_SNr_multiple
from src import misc
import scipy.optimize as opt
OUTPUT_PATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]

# Assume: n_msn, n_stn, weight_msn_snr, weight_gpe_snr
# Find: n_gpe, weight_stn_snr and snr_current
#
# 1. With no GPe input (still with MSN input) and STN at 10 Hz (not assuming that
#    GPe are silent) the fire at 312 % if baseline (30*3=90 Hz)s. 
#    REMARK: If STN at 20 Hz (assuming GPe to STN is silent) then to be able to
#    to come close to restriction STN needs to be depressing onto SNr  
# 2. 15 Hz without STN OBS Since stn is of the GPE should fire at 15 Hz
# 3. 30 Hz normal, GPe at 30 Hz, STN at 10 Hz
# 
#
# At GPe = base*0.5 Hz and no stn then SNr =base*0.5  0
#simulate_basa_line_SNr(msn_rate=0.1, gpe_rate=30., stn_rate=10., n_msn=500,n_gpe=32, n_stn=30, 
#                       neuron_model=['SNR_aeif'],
#                       snr_current=280,# This is the current to add to I_e base
#                       sim_time=10000, threads=8, stn_syn='STN_SNR_ampa_p3')

neuron_model=['SNR_aeif']
syn_models=[]

msn_rate, gpe_rate, stn_rate=0.1, 30., 10.
n_msn, n_stn=500, 30 
syn_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p','STN_SNR_ampa_s']

def plot_example(ax, SNR_target, sim_time, x,  data, type):
    meanRate=round(SNR_target.signals['spikes'].mean_rate(2000,sim_time),1)
    spk=SNR_target.signals['spikes'].time_slice(2000,sim_time).raw_data()
    CV=numpy.std(numpy.diff(spk[:,0],axis=0))/numpy.mean(numpy.diff(spk[:,0],axis=0))

    
    SNR_target.signals['V_m'].my_set_spike_peak( 15, spkSignal= SNR_target.signals['spikes'] )          
    pylab.rcParams.update( {'path.simplify':False}    )
    
    SNR_target.signals['V_m'].plot(display=ax)
    ax.set_title( type+' mean='+str(round(data[0],3))+'mean='+str(meanRate)+ 'Hz, std='+str(round(data[1],3))+
                  'Hz, CV='+str(round(CV,3))+
                 ', snr_curr='+str(round(x[0],3))+ 
                  ' w_STN_SNR='+str(round(x[1],3)) )
        
def restriction_1(msn_rate, gpe_rate, stn_rate_ch, n_msn, n_gpe, x,  neuron_model, syn_models, n_exp,  record_vm=True):
# 1. With no GPe input (still with MSN input) and STN at 10 Hz (since in the
#    experiment the firing rate of STN were not recorded we can not know if it
#    was increased or not/and to what spatial extend, could be that SNr recieves 
#    increased input due to increased STN. However fitting with constatn 
#    synapses STN to SNr makes it difficult the hit the target rates) to  are silent) the fire at
#    312 % if baseline (30*3=90 Hz)s. REMARK: If STN at 20 then to be able to
#    to come close to restriction STN needs to be depressing onto SNr  
    c,w=x
    target_rate1=30*3.12
    n_gpe_ch=0
    #stn_rate_ch=20.0
    SNR_list=simulate_basa_line_SNr_multiple(msn_rate, gpe_rate, stn_rate_ch, 
                                             n_msn, n_gpe_ch, n_stn, 
                                             neuron_model,syn_models, c,# This is the current to add to I_e base
                                             sim_time, 8,  w_STN_SNR= w,seed=n_exp,
                                             record_vm=record_vm, n_neurons=n_exp) 
 
    r_SNR=[]
    for SNR in SNR_list:
        r_SNR.append(SNR.signals['spikes'].mean_rate(1000,sim_time))
    
    r=numpy.array(r_SNR)
    mr=numpy.mean(r,axis=0)
    std=numpy.std(r,axis=0)
    e=mr-target_rate1
    
    #Normalize error
    e=e/target_rate1
        
    return SNR_list[0], mr, std, e

def restriction_2(msn_rate, gpe_rate, stn_rate, n_msn, n_gpe, x, neuron_model, syn_models, n_exp, record_vm=True):
        # 2. 15 Hz without STN. OBS Since stn is of the GPE should fire at 15 Hz
        c, w=x
        target_rate2=15.0
        stn_rate_ch=0.
        gpe_rate_ch=15.0 # Since stn is of the GPE should fire at 15 Hz
        r_SNR=[]
        
        SNR_list=simulate_basa_line_SNr_multiple(msn_rate, gpe_rate_ch, stn_rate_ch, 
                                              n_msn, n_gpe, n_stn, 
                                              neuron_model, syn_models, c,# This is the current to add to I_e base
                                              sim_time, 8, w_STN_SNR= w, seed=n_exp,
                                              record_vm=record_vm, n_neurons=n_exp) 
        r_SNR=[]
        for SNR in SNR_list:
            r_SNR.append(SNR.signals['spikes'].mean_rate(1000,sim_time))
        
        r=numpy.array(r_SNR)
        mr=numpy.mean(r,axis=0)
        std=numpy.std(r,axis=0)
        e=mr-target_rate2
        
        #Normalize error
        e=e/target_rate2
        
        return SNR_list[0], mr, std, e
    
def restriction_3(msn_rate, gpe_rate, stn_rate, n_msn, n_gpe, x, neuron_model, syn_models, n_exp, record_vm=True):
        # 3. 30 Hz normal, GPe at 30 Hz, STN at 10 Hz
        c, w=x
        target_rate3=30.0
        
        SNR_list=simulate_basa_line_SNr_multiple(msn_rate, gpe_rate, stn_rate, n_msn, n_gpe, n_stn, 
                                     neuron_model, syn_models, 
                                      c,# This is the current to add to I_e base
                                      sim_time, 8, w_STN_SNR= w, seed=n_exp,
                                      record_vm=record_vm, n_neurons=n_exp) 
        r_SNR=[]
        for SNR in SNR_list:
            r_SNR.append(SNR.signals['spikes'].mean_rate(1000,sim_time))
        
        r=numpy.array(r_SNR)
        mr=numpy.mean(r,axis=0)
        std=numpy.std(r,axis=0)
        e=mr-target_rate3
        
        #Normalize error
        e=e/target_rate3
        
        return SNR_list[0], mr, std, e


def GPE_SRN_ref(msn_rate, gpe_rate, stn_rate, n_msn, x,   neuron_model, syn_models):
        # 32 Hz with GPE to SNr ref
        w=x
        target_rate2=32.0
        SNR_target=simulate_basa_line_SNr(msn_rate, gpe_rate, stn_rate, n_msn, int(n), n_stn, 
                                      neuron_model,syn_models, 
                                      c,# This is the current to add to I_e base
                                      sim_time, 8, w_STN_SNR= w) 
 
        e=round(SNR_target.signals['spikes'].mean_rate(1000,sim_time),1)-target_rate2
        
        return SNR_target, e

def error_fun(x, sim_time,n_exp, n_gpe, res1_stn_rate_ch):
        
        record_vm=False
        SNR_target, mr, std, e1=restriction_1(msn_rate, gpe_rate, res1_stn_rate_ch, n_msn, n_gpe, x, neuron_model, syn_models, n_exp,record_vm)
        SNR_target, mr, std, e2=restriction_2(msn_rate, gpe_rate, stn_rate, n_msn, n_gpe, x, neuron_model, syn_models, n_exp, record_vm)
        SNR_target, mr, std, e3=restriction_3(msn_rate, gpe_rate, stn_rate, n_msn, n_gpe, x, neuron_model, syn_models, n_exp, record_vm)
             
        print e1**2 + e2**2 + e3**2
        return e1**2 + e2**2 + e3**2
        
def fmin(load, save_at,  x0, n_exp, n_gpe, res1_stn_rate_ch):
    sim_time=10000
  

    if not load:
        [xopt,fopt, iter, funcalls , warnflag, allvecs] = opt.fmin(error_fun, x0, 
                                                                   args=([sim_time, n_exp, n_gpe,res1_stn_rate_ch]), maxiter=10, maxfun=10,full_output=1, retall=1)

        misc.pickle_save([xopt,fopt, iter, funcalls , warnflag, allvecs], save_at)
    else:
        [xopt,fopt, iter, funcalls , warnflag, allvecs]=misc.pickle_load(save_at)        
    return xopt, fopt  



sim_time=20000
n_exp=4.


# Solution when STN is at 10 while GPe is silent
  
res1_stn_rate_ch=10.0
save_at=OUTPUT_PATH+'/simulate_network_fmin.plk' 
n_conn=[31,33]
load=1
s1=''
s2=''
for n in range(n_conn[0], n_conn[1]):
    save_at=OUTPUT_PATH+'/simulate_network_fmin'+str(n)+str(res1_stn_rate_ch)+str(n_exp)+'.plk'     
    x0=[220, 0.75]
    x, e=fmin(load,save_at, x0, n_exp, n,  res1_stn_rate_ch)
    s1=s1+'\ne={0:.5f}, n_gpe={1}, c={2:.2f}, w_STN_GPE={3:.3f}, '.format(e, n, x[0], x[1])
    if load:
        SNR_target1, mr1, std1,  e1=restriction_1(msn_rate, gpe_rate, res1_stn_rate_ch, n_msn, n, x, neuron_model, syn_models, n_exp)
        SNR_target2, mr2, std2,  e2=restriction_2(msn_rate, gpe_rate, stn_rate, n_msn, n, x, neuron_model, syn_models, n_exp)
        SNR_target3, mr3, std3,  e3=restriction_3(msn_rate, gpe_rate, stn_rate, n_msn, n, x, neuron_model, syn_models, n_exp)
        s2=s2+'\nm1-3={0:.2f}, {1:.2f}, {2:.2f} std1-3={3:.2f}, {4:.2f}, {5:.2f} e1-3={6:.4f}, {7:.4f}, {8:.4f}'.format(mr1,mr2,mr3,std1,std2,std3,e1,e2,e3)
print s1
print s2
        
    
#res1_stn_rate_ch=20
#save_at=OUTPUT_PATH+'/simulate_network_fmin.plk' 
#n_conn=[23,27]
#for n in range(n_conn[0], n_conn[1]):
#    save_at=OUTPUT_PATH+'/simulate_network_fmin'+str(n)+str(res1_stn_rate_ch)+str(n_exp)+'.plk'     
#    x0=[190, 0.75]
#    x=fmin(1,save_at, x0, n_exp, n,  res1_stn_rate_ch)    
    
#x=[25, 190, 0.75]
#x=[33, 228, 1.083]
#x=[228, 1.083]
n_gpe=33
x=[239, 0.91]
n_gpe=32
stn_rate=10.0
n_exp=4.
SNR_target1, mr1, std1,  e1=restriction_1(msn_rate, gpe_rate, stn_rate, n_msn, n_gpe, x, neuron_model, syn_models, n_exp)
SNR_target2, mr2, std2,  e2=restriction_2(msn_rate, gpe_rate, stn_rate, n_msn, n_gpe, x, neuron_model, syn_models, n_exp)
SNR_target3, mr3, std3,  e3=restriction_3(msn_rate, gpe_rate, stn_rate, n_msn, n_gpe, x, neuron_model, syn_models, n_exp)


plot_settings.set_mode(pylab, mode='by_fontsize', w = 1000.0, h = 1000.0, fontsize=16)
font_size_text = 8
fig = pylab.figure( facecolor = 'w')

ax_list=[]
ax_list.append( MyAxes(fig, [ .1, .7, .8, .2 ] ) )    # text box
ax_list.append( MyAxes(fig, [ .1, .4, .8, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .1, .1, .8, .2 ] ) )    # 
 
ax=ax_list[0]
plot_example(ax, SNR_target1,sim_time, x, [mr1, std1], type='No GABA')  

ax=ax_list[1]
plot_example(ax, SNR_target2, sim_time, x, [mr2,std2], type='No STN')  

ax=ax_list[2]
plot_example(ax, SNR_target3, sim_time, x, [mr3, std3], type='Normal')  


        
pylab.show()