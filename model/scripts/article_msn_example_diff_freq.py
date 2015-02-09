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

from simulation_utils import _simulate_model, simulate_model
from src import misc

OUTPUT_PATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]
ADJUST_XDATA_MS=1500.
SEL_ONSET=2000


N_MSN=500
N_MSN_BURST=17
N_MSN_BASE=N_MSN-N_MSN_BURST
N_GPE = 30
N_STN = 100

MSN_BASE_RATE=0.1
GPE_BASE_RATE=25
STN_BASE_RATE=10

def plot_diff_3x_freq(ax, times, spk_mean ):
    spk_mean=spk_mean
    colors=['b','g','r']
    for i, col in enumerate(colors):
        #ax.plot(times, misc.convolve(spk_mean, 50, 'triangle',single=True)[0]*1000,col)
        ax.plot(times, spk_mean[i,:], col)

    lines = ax.lines
    for line in lines:
        line.set_xdata(line.get_xdata()-ADJUST_XDATA_MS)    
        misc.slice_line(line, xlim=[0,1490])
    #ax.set_ylim([0,40])
    
    ax.set_xlim(misc.adjust_limit([0,1500]))
    ax.set_ylim(misc.adjust_limit([0,40]))
    pylab.show()

def simulate_diff_3x_freq(hzs=[10,20, 30],n_msns=[40/2, 20/2, 13/2], res=20, load=True):

    save_result_at=OUTPUT_PATH+'/simulate_diff_3x_freq.pkl'
    n_exp=100
    if not load:  
        spk_mean=[]
        for hz, n_msn in zip(hzs,n_msns):
            sim_time= SEL_ONSET+1000
            params_msn={'base_rates':[MSN_BASE_RATE], 'base_times':[1], 'mod_rates': [MSN_BASE_RATE, hz, MSN_BASE_RATE],
                        'mod_times':[1,SEL_ONSET, SEL_ONSET+500], 'n_tot':N_MSN, 'n_mod':n_msn}
            params_gpe={'base_rates':[GPE_BASE_RATE], 'base_times':[1], 'n_tot':N_GPE, 'n_mod':0}
            params_stn={'base_rates':[STN_BASE_RATE], 'base_times':[1], 'n_tot':N_STN, 'n_mod':0}
            synapse_models={'MSN':'MSN_SNR_gaba_p1', 'GPE':'GPE_SNR_gaba_p',
                            'STN':'STN_SNR_ampa_s'}
            
            #times, spk_binned =_simulate_model([params_msn, params_gpe,params_gpe, 'SNR_izh', 
            #                                   synapse_models, sim_time, 0])
            
            t=ttime.time()
        
            times, spk_binned =simulate_model(params_msn, params_gpe, params_stn, 'SNR_izh', 
                                              synapse_models, sim_time, res, n_exp=n_exp,threads=4)
            
            print 'Time:',ttime.time()-t
            spk_mean.append(numpy.mean(spk_binned,axis=0))
        spk_mean=numpy.array(spk_mean)*1000/res       
        misc.pickle_save([times, spk_mean], save_result_at)        
    else:        
        times, spk_mean = misc.pickle_load(save_result_at)   
    return times, spk_mean
            
times, spk_mean=simulate_diff_3x_freq(load=False)
ax=pylab.subplot(111)
plot_diff_3x_freq(ax, times, spk_mean )
pylab.show()




