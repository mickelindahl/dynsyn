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



hzs=numpy.arange(7,49,1)

hzs=[10,20, 30]
n_msns=[40/2, 20/2, 13/2]
colors=['b','g','r']
n_exp=200

for hz,n_msn, col in zip(hzs,n_msns,colors ):
    sim_time=2500
    params_msn={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, hz, 0.1],
                'mod_times':[1,1000, 2000], 'n_tot':500, 'n_mod':n_msn}
    params_gpe={'base_rates':[25], 'base_times':[1], 'mod_rates': [25, 75, 25],
                'mod_times':[1,1900, 2000], 'n_tot':30, 'n_mod':15}
    params_stn={'base_rates':[10], 'base_times':[1], 'n_tot':100, 'n_mod':0}
    synapse_models={'MSN':'MSN_SNR_gaba_p1', 'GPE':'GPE_SNR_gaba_p',
                    'STN':'STN_SNR_ampa_s'}
    
    #times, spk_binned =_simulate_model([params_msn, params_gpe,params_gpe, 'SNR_izh', 
    #                                   synapse_models, sim_time, 0])
    
    t=ttime.time()

    times, spk_binned =simulate_model(params_msn, params_gpe, params_stn, 'SNR_izh', 
                                      synapse_models, sim_time, n_exp=n_exp,threads=4)
    
    print 'Time:',ttime.time()-t
    spk_mean=numpy.sum(spk_binned,axis=0)
    #spk_std=numpy.std(spk_binned,axis=0)
    #pylab.plot(times, spk_mean*1000/n_exp)
    pylab.plot(times, misc.convolve(spk_mean, 50, 'triangle',single=True)[0]*1000/n_exp,col)
pylab.ylim([0,40])
pylab.show()



