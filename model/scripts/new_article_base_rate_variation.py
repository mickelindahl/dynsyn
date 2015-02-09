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
from model_params import models, network  
from src.my_population import MyGroup

OUTPUT_PATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]

my_nest.ResetKernel(threads=4)       
model_list, model_dict=models({})
      



sim_time=10000.0
n=100
neuron_models=['SNR_aeif','GPE_aeif','STN_75_aeif']
norm_std=[0.1/0.3, 0.1/0.22, 0.1/0.61] # This shoud give 1.0 Hz in std for each 
                                       # of them
f=[0.122,0.108,0.19]
mrs=[13.96*f[0], 15.33*f[1], 9.63*f[2]]

save_result_at=OUTPUT_PATH+'/simulate.plk'
if 0:
    neuron_list=[]
    for i, model in enumerate(neuron_models):
        my_nest.MyLoadModels( model_dict, [model] )
        I_in_vitro=my_nest.GetDefaults(model)['I_e']
        neuron=MyGroup( model, n=n, sd=True, mm_dt=.1, mm=False)  
        for id in neuron.ids:
            I=numpy.random.normal(I_in_vitro, I_in_vitro*norm_std[i]*mrs[i])  
            my_nest.SetStatus([id],{'I_e':I})
        neuron_list.append(neuron)  
        noise=my_nest.Create('noise_generator', params={'mean':0.,'std':1.})
        rec=my_nest.GetStatus(neuron[:])[0]['receptor_types']
        
        for id in neuron.ids:
            my_nest.Connect( noise, [id], params = { 'receptor_type' : rec['CURR'] } )    
    
    my_nest.MySimulate(sim_time)         
            
    mr_list=[]
    for neuron in neuron_list:
        neuron.get_signal( 's', start=0, stop=sim_time )   
        signal=neuron.signals['spikes']
        mr=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=1)
        mr_list.append(mr) 
    
    misc.pickle_save(mr_list, save_result_at)  
else:
    mr_list=misc.pickle_load(save_result_at)      

for i in range(len(neuron_models)):
    try: my_nest.MyLoadModels( model_dict, [neuron_models[i]] )
    except: print 'Already loaded'
    I_in_vitro=my_nest.GetDefaults(neuron_models[i])['I_e']
    ax=pylab.subplot(3,1,i+1)   
    ax.hist(mr_list[i], 20)   
    m=numpy.mean(mr_list[i], axis=0)
    std=numpy.std(mr_list[i], axis=0)
    ax.set_title('m={0:.2f} std={1:.2f}, norm std by mean={2:.2f}, curr std={3:.2f}'.format(m, std, std/m, I_in_vitro*norm_std[i]*mrs[i])) 
        
pylab.show()  
        
        
        
        
        
        