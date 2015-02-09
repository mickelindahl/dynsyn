
import numpy
import pylab
import os
import sys
import random

# Get directory where model and code resides 
model_dir=   '/'.join(os.getcwd().split('/')[0:-1])    
code_dir=  '/'.join(os.getcwd().split('/')[0:-2])  

# Add model, code and current directories to python path
sys.path.append(os.getcwd())  
sys.path.append(model_dir)
sys.path.append(code_dir+'/nest_toolbox') 

          
# Imports dependent on adding code model and nest_toolbox path
from model_params import models, network                                  
from src import my_nest, misc
from src.my_population import MyGroup, MyPoissonInput, MyLayerGroup,MyLayerPoissonInput 
from multiprocessing import Pool
from numpy.random import random_integers
import nest.topology as tp
import nest
numpy.random.seed(0)
N_GPE=30
N_MSN=500
N_STN=100

BASE_RATE_GPE=25
BASE_RATE_MSN=0.1
BASE_RATE_STN=10.0

SYNAPSE_MODEL_BACKROUND_GPE=['GPE_SNR_gaba_p']
SYNAPSE_MODEL_BACKROUND_MSN=['MSN_SNR_gaba_p1']
SYNAPSE_MODEL_BACKROUND_STN=['STN_SNR_ampa_s']


I_IN_VIVO={'SNR':239, 'STN':0,'GPE':42}
BASE_RATE_CTX_STN=189.
I_E_VARIATION={'SNR':8.5,'GPE':3.8, 'STN':1.8} # see script new_article_base_rates, generate gauss dist 0.2 std out of in_vivo mean rates
#I_E_VARIATION={'SNR':0,'GPE':0, 'STN':0}

def simulate_basa_line_SNr(msn_rate, gpe_rate, stn_rate, n_msn, n_gpe, n_stn, 
                           neuron_model, syn_models, snr_current, sim_time=1000, 
                           threads=8,  w_STN_SNR=0, seed=0, record_vm=True, multiple=False):
    
    if not multiple:
        Params_in={}
        if w_STN_SNR: Params_in['STN_SNR_ampa_s']=w_STN_SNR
        
        model_list, model_dict=models(Params_in)
        
        my_nest.ResetKernel(threads=threads)       
        my_nest.MyLoadModels( model_dict, neuron_model )
        my_nest.MyLoadModels( model_dict, syn_models)       
       
    SNR_list=[] # List with SNR groups for synapse. 
    
    if n_msn>0: MSN=MyPoissonInput(n=n_msn, sd=True)
    if n_gpe>0: GPE=MyPoissonInput(n=n_gpe, sd=False)
    if n_stn>0: STN=MyPoissonInput(n=n_stn, sd=False)
    
    if n_msn>0: MSN.set_spike_times(rates=[ msn_rate], times=[1], t_stop=sim_time, seed=seed)    
    if n_gpe>0: GPE.set_spike_times(rates=[gpe_rate], times=[1], t_stop=sim_time, seed=seed)     
    if n_stn>0: STN.set_spike_times(rates=[stn_rate], times=[1], t_stop=sim_time, seed=seed)     
           
    I_e=my_nest.GetDefaults(neuron_model[0])['I_e']+snr_current

    SNR=MyGroup( neuron_model[0], n=1, sd=True, params={'I_e':I_e}, mm_dt=.1, mm=record_vm)    
       
    if n_msn>0: my_nest.Connect(MSN[:],SNR[:]*len(MSN[:]), model=syn_models[0])
    if n_gpe>0: my_nest.Connect(GPE[:],SNR[:]*len(GPE[:]), model=syn_models[1])   
    if n_stn>0: my_nest.Connect(STN[:],SNR[:]*len(STN[:]), model=syn_models[2])
                    
   
    return SNR
def simulate_basa_line_SNr_multiple(msn_rate, gpe_rate, stn_rate, n_msn, n_gpe, n_stn, 
                           neuron_model, syn_models, snr_current, sim_time=1000, 
                           threads=8,  w_STN_SNR=0, seed=0, record_vm=True, n_neurons=1):
    
    Params_in={}
    if w_STN_SNR: Params_in['STN_SNR_ampa_s']=w_STN_SNR
    
    model_list, model_dict=models(Params_in)
    
    my_nest.ResetKernel(threads=threads)       
    my_nest.MyLoadModels( model_dict, neuron_model )
    my_nest.MyLoadModels( model_dict, syn_models)     
    SNR_list=[]
    for i in range(int(n_neurons)):
        seed=seed+i
        SNR_list.append(simulate_basa_line_SNr(msn_rate, gpe_rate, stn_rate, n_msn, n_gpe, n_stn, 
                           neuron_model, syn_models, snr_current, sim_time, 
                           threads,  w_STN_SNR, seed, record_vm, True)) 
    
    
    my_nest.MySimulate(sim_time)    
    
    for i in range(int(n_neurons)):
        SNR_list[i].get_signal( 's', start=0, stop=sim_time )
        if record_vm: SNR_list[i].get_signal( 'v',recordable='V_m', start=0, stop=sim_time ) 
        
    return SNR_list

def simulate_basa_line_STN(ctx_rate, gpe_rate, n_ctx, n_gpe, neuron_model, syn_models, stn_current, sim_time=1000, threads=8, w_GPE_STN=0):
    
    Params_in={}
    if w_GPE_STN: Params_in['GPE_STN_gaba_s']=w_GPE_STN
    
    model_list, model_dict=models(Params_in)
    my_nest.ResetKernel(threads=threads)       
    my_nest.MyLoadModels( model_dict, neuron_model )
    my_nest.MyLoadModels( model_dict, syn_models)       
    
    SNR_list=[] # List with SNR groups for synapse. 
    
    if n_ctx>0: CTX=MyPoissonInput(n=n_ctx, sd=True)
    if n_gpe>0: GPE=MyPoissonInput(n=n_gpe, sd=False)
    
    if n_ctx>0: CTX.set_spike_times(rates=[ ctx_rate], times=[1], t_stop=sim_time, seed=0)    
    if n_gpe>0: GPE.set_spike_times(rates=[gpe_rate], times=[1], t_stop=sim_time, seed=0)     
           
    I_e=my_nest.GetDefaults(neuron_model[0])['I_e']+stn_current

    STN=MyGroup( neuron_model[0], n=1, sd=True, params={'I_e':I_e}, 
                 mm_dt=.1, mm=True)    
        
    if n_ctx>0: my_nest.Connect(CTX[:], STN[:]*len(CTX[:]), model=syn_models[0])
    if n_gpe>0: my_nest.Connect(GPE[:],STN[:]*len(GPE[:]), model=syn_models[1])   

                    
    my_nest.MySimulate(sim_time)    
    
    STN.get_signal( 's', start=0, stop=sim_time )
    STN.get_signal( 'v',recordable='V_m', start=0, stop=sim_time )  
     
    return  STN   
 
def simulate_basa_line_GPe(msn_rate, stn_rate, gpe_rate,  n_msn, n_stn, n_gpe, neuron_model, syn_models, gpe_current, sim_time=1000, threads=8, w_GPE_GPE=False, w_STN_GPE=False):
    
    Params_in={}
    if w_GPE_GPE: Params_in['GPE_GPE_gaba_s']=w_GPE_GPE
    if w_STN_GPE: Params_in['STN_GPE_ampa_s']=w_STN_GPE
    
    model_list, model_dict=models(Params_in)
    my_nest.ResetKernel(threads=threads)       
    my_nest.MyLoadModels( model_dict, neuron_model )
    my_nest.MyLoadModels( model_dict, syn_models)       
    
    SNR_list=[] # List with SNR groups for synapse. 
    
    if n_msn>0: MSN=MyPoissonInput(n=n_msn, sd=False)
    if n_stn>0: STN=MyPoissonInput(n=n_stn, sd=False)
    if n_gpe>0: GPE=MyPoissonInput(n=n_gpe, sd=False)
    
    if n_msn>0: MSN.set_spike_times(rates=[msn_rate], times=[1], t_stop=sim_time, seed=0)    
    if n_stn>0: STN.set_spike_times(rates=[stn_rate], times=[1], t_stop=sim_time, seed=0)     
    if n_gpe>0: GPE.set_spike_times(rates=[gpe_rate], times=[1], t_stop=sim_time, seed=0)  
               
    I_e=my_nest.GetDefaults(neuron_model[0])['I_e']+gpe_current

    GPE_target=MyGroup( neuron_model[0], n=1, sd=True, params={'I_e':I_e}, 
                 mm_dt=.1, mm=True)    
        
    if n_msn>0: my_nest.Connect(MSN[:], GPE_target[:]*len(MSN[:]), model=syn_models[0])
    if n_stn>0: my_nest.Connect(STN[:],GPE_target[:]*len(STN[:]), model=syn_models[1])   
    if n_gpe>0: my_nest.Connect(GPE[:],GPE_target[:]*len(GPE[:]), model=syn_models[2])   
                    
    my_nest.MySimulate(sim_time)    
    
    GPE_target.get_signal( 's', start=0, stop=sim_time )
    GPE_target.get_signal( 'v',recordable='V_m', start=0, stop=sim_time )     
    
    return GPE_target       
    
def simulate_network(params_msn_d1, params_msn_d2, params_stn,
                     synapse_models, sim_time, seed, I_e_add, threads=1, 
                     start_rec=0, model_params={}, params_in={}, dis_conn_GPE_STN=False):    
    '''
        params_msn_d1 - dictionary with timing and burst freq setup for msn
                     {'base_rates':0.1, 
                      'base_times':[1], 
                      'mod_rates': 20,
                      'mod_times':[1,200], 
                      'mod_units':list()
                      'n_tot':500, 
                       n_mod=20}
        params_msn_d2 - dictionary with timing and burst freq setup for gpe
        params_stn    - dictionary {'rate':50}
                     same as params_msn
        neuron_model - string, the neuron model to use 
        synapse_models - dict, {'MSN':'...', 'GPE':,'...', 'STN':'...'}
        sim_time - simulation time
        seed - seed for random generator
        I_e_add - diabled
        start_rec - start recording from
        model_params - general model paramters
    '''

  
    my_nest.ResetKernel(threads=8) 
    numpy.random.seed(seed)
    
    params = {'conns':{'MSN_D1_SNR':{'syn':synapse_models[0]},   
                       'GPE_SNR':{'syn':synapse_models[1]}}}   
    
    params=misc.dict_merge(model_params, params)
    params=misc.dict_merge({'neurons':{'GPE':{'paused':0}}}, params)
               
    model_list, model_dict = models(params_in)
    layer_list, connect_list = network(model_dict, params)
    
    
    # Create neurons and synapses
    layer_dic={}  
    for name, model, props  in layer_list:
        
        # Update input current
        my_nest.MyLoadModels( model_dict, [model[1]] )
        if name in I_IN_VIVO.keys():
            I_e=my_nest.GetDefaults(model[1])['I_e']+I_IN_VIVO[name]
            my_nest.SetDefaults(model[1],{'I_e':I_e})  
        
        
        #! Create layer, retrieve neurons ids per elements and p
        if model[0]=='spike_generator':
            layer=MyLayerPoissonInput(layer_props=props, sd=True, sd_params={'start':start_rec, 'stop':sim_time})
 
        else:
            layer=MyLayerGroup(layer_props=props, sd=True, mm=False, mm_dt = 0.1,
                               sd_params={'start':start_rec, 'stop':sim_time} )
               
            Im=0
            for iter, id in enumerate(layer[:]):

                if name=='GPE' and params_msn_d2['n_mod'] and iter<params['neurons']['GPE']['paused']:
                        scg = my_nest.Create( 'step_current_generator',n=1)  
                        rec=my_nest.GetStatus([id])[0]['receptor_types']
                        my_nest.SetStatus(scg, {'amplitude_times':params_msn_d2['mod_times'],
                                                'amplitude_values':[0.,-300.,0.]})
                        my_nest.Connect( scg, [id],  params = { 'receptor_type' : rec['CURR'] } )
                    
                
                I_e=my_nest.GetDefaults(model[1])['I_e']  
                if  I_E_VARIATION[name]:I=numpy.random.normal(I_e, I_E_VARIATION[name])
                else:I=I_e
                #I=I_e
                
                my_nest.SetStatus([id],{'I_e':I})
                Im+=my_nest.GetStatus([id],'I_e')[0]
            Im/=len(layer[:])
            print name, Im
        layer_dic[name]=layer
    
    
    # Connect populations
    for conn in connect_list:
        name=conn[0] + '_' + conn[1]
        my_nest.MyLoadModels( model_dict, [conn[2]['synapse_model']] )
        
        
        if dis_conn_GPE_STN=='GPE' and (name in ['GPE_SNR']):
            r,syn=32*30.0, 'GPE_SNR_gaba_s_ref'
            if not syn in my_nest.Models():    
                my_nest.MyLoadModels( model_dict, [syn] )
            pg=my_nest.Create('poisson_generator',1,{'rate':r,'start':1.})
            my_nest.DivergentConnect(pg, layer_dic[conn[1]].ids, model=syn)
        elif dis_conn_GPE_STN=='STN' and (name in ['STN_SNR']):
            r, syn=30*10.0, 'STN_SNR_ampa_s'
            if not syn in my_nest.Models():               
                my_nest.MyLoadModels( model_dict, [syn] )
            pg=my_nest.Create('poisson_generator',1,{'rate':r,'start':1.})
            my_nest.DivergentConnect(pg, layer_dic[conn[1]].ids, model=syn)
            
            
        else:
            name=name+'_'+conn[3]    
            tp.ConnectLayers(layer_dic[conn[0]].layer_id, layer_dic[conn[1]].layer_id, conn[2])
            layer_dic[conn[1]].add_connection(source=layer_dic[conn[0]], type=conn[3], props=conn[2])
    
    # Sort MSN D2 such that the closest to center is first in ids list.
    # Do this to we can get focused inhibition in GPe    
    if params_msn_d2['focus']:
        MSN_D2_idx=layer_dic['MSN_D2'].sort_ids()
    else:
        MSN_D2_idx=range(len(numpy.array(layer_dic['MSN_D2'].ids)))
        
    n_mod_msn_d1=params_msn_d1['n_mod']
    n_mod_msn_d2=params_msn_d2['n_mod']
 
    MSN_D1_ids=layer_dic['MSN_D1'].ids
    MSN_D2_ids=layer_dic['MSN_D2'].ids
    
    MSN_D1_mod,MSN_D2_mod=[],[]
    if params_msn_d1['n_mod']:MSN_D1_mod=MSN_D1_ids[0:n_mod_msn_d1]
    if params_msn_d2['n_mod']:MSN_D2_mod=MSN_D2_ids[0:n_mod_msn_d2*params_msn_d2['skip']:params_msn_d2['skip']]
    
       
    MSN_D1_base=list(set(MSN_D1_ids).difference(MSN_D1_mod))
    MSN_D2_base=list(set(MSN_D2_ids).difference(MSN_D2_mod))
    
    #layer_dic['MSN_D1'].ids[0:n_base_msn_d1]
    
    #MSN_D2_ids=numpy.array(layer_dic['MSN_D2'].ids)
    #MSN_D2_base=MSN_D2_ids#[MSN_D2_idx[0:n_base_msn_d1]]
    
    
    #set().difference(t) 
    
    layer_dic['MSN_D1'].set_spike_times(params_msn_d1['base_rates'], params_msn_d1['base_times'], sim_time, ids=MSN_D1_base)
    layer_dic['MSN_D2'].set_spike_times(params_msn_d2['base_rates'], params_msn_d2['base_times'], sim_time, ids=MSN_D2_base)
    
    if params_msn_d1['n_mod']: layer_dic['MSN_D1'].set_spike_times(params_msn_d1['mod_rates'], params_msn_d1['mod_times'], sim_time, ids=MSN_D1_mod)        
    if params_msn_d2['n_mod']: layer_dic['MSN_D2'].set_spike_times(params_msn_d2['mod_rates'], params_msn_d2['mod_times'], sim_time, ids=MSN_D2_mod) 
     
    STN_CTX_input_base=my_nest.Create('poisson_generator',params={'rate':BASE_RATE_CTX_STN, 'start':0., 'stop':sim_time})
    my_nest.MyLoadModels( model_dict, ['CTX_STN_ampa_s'] )
    my_nest.DivergentConnect(STN_CTX_input_base, layer_dic['STN'].ids, model='CTX_STN_ampa_s')
    
    
    if params_stn['mod']:
        STN_CTX_input_mod=my_nest.Create('poisson_generator',params={'rate':params_stn['mod_rate'], 
                                                                     'start':params_stn['mod_times'][0],
                                                                     'stop':params_stn['mod_times'][1]})
        my_nest.DivergentConnect(STN_CTX_input_mod, layer_dic['STN'].ids, model='CTX_STN_ampa_s')
    
    #tar=[]
    #for id in layer_dic['MSN_D1'].ids:
    #    tar.extend(sorted(nest.GetStatus(my_nest.FindConnections([id]),'target'))[:-1])
        
    #pylab.subplot(211).hist(tar, 1500) 
#    
#    tar=[]
#    for id in layer_dic['MSN_D2'].ids:
#        tar.extend(sorted(nest.GetStatus(my_nest.FindConnections([id]),'target'))[1:])
#        
#    pylab.subplot(212).hist(tar, 1500) 
    #pylab.show()   
#    
#    
    my_nest.MySimulate(sim_time)    
       
    if params_msn_d1['n_mod']:layer_dic['MSN_D1'].id_mod=MSN_D1_mod
    if params_msn_d2['n_mod']:layer_dic['MSN_D2'].id_mod=MSN_D2_mod
        
    layer_dic['MSN_D1'].get_signal( 's', start=start_rec, stop=sim_time )
    layer_dic['MSN_D2'].get_signal( 's', start=start_rec, stop=sim_time )
    layer_dic['GPE'].get_signal( 's', start=start_rec, stop=sim_time )    
    layer_dic['SNR'].get_signal( 's', start=start_rec, stop=sim_time )    
    layer_dic['STN'].get_signal( 's', start=start_rec, stop=sim_time )    
    
    return layer_dic

def simulate_network_test(params_msn_d1, params_msn_d2, params_stn,
                     synapse_models, sim_time, seed, I_e_add, threads=1, 
                     start_rec=0, model_params={}, params_in={}, dis_conn_GPE_STN=False):    
    '''
        params_msn_d1 - dictionary with timing and burst freq setup for msn
                     {'base_rates':0.1, 
                      'base_times':[1], 
                      'mod_rates': 20,
                      'mod_times':[1,200], 
                      'mod_units':list()
                      'n_tot':500, 
                       n_mod=20}
        params_msn_d2 - dictionary with timing and burst freq setup for gpe
        params_stn    - dictionary {'rate':50}
                     same as params_msn
        neuron_model - string, the neuron model to use 
        synapse_models - dict, {'MSN':'...', 'GPE':,'...', 'STN':'...'}
        sim_time - simulation time
        seed - seed for random generator
        I_e_add - diabled
        start_rec - start recording from
        model_params - general model paramters
    '''

  
    my_nest.ResetKernel(threads=8) 
    numpy.random.seed(seed)
    
    params = {'conns':{'MSN_D1_SNR':{'syn':synapse_models[0]},   
                       'GPE_SNR':{'syn':synapse_models[1]}}}   
    
    params=misc.dict_merge(model_params, params)
    params=misc.dict_merge({'neurons':{'GPE':{'paused':0}}}, params)
               
    model_list, model_dict = models(params_in)
    layer_list, connect_list = network(model_dict, params)
    
    
    # Create neurons and synapses
    layer_dic={}  
    for name, model, props  in layer_list:
        
        # Update input current
        my_nest.MyLoadModels( model_dict, [model[1]] )
        if name in I_IN_VIVO.keys():
            I_e=my_nest.GetDefaults(model[1])['I_e']+I_IN_VIVO[name]
            my_nest.SetDefaults(model[1],{'I_e':I_e})  
        
        
        #! Create layer, retrieve neurons ids per elements and p
        if model[0]=='spike_generator':
            layer=MyLayerPoissonInput(layer_props=props, sd=True, sd_params={'start':start_rec, 'stop':sim_time})
 
        else:
            layer=MyLayerGroup(layer_props=props, sd=True, mm=False, mm_dt = 0.1,
                               sd_params={'start':start_rec, 'stop':sim_time} )
               
            
            for iter, id in enumerate(layer[:]):

                if name=='GPE' and params_msn_d2['n_mod'] and iter<params['neurons']['GPE']['paused']:
                        scg = my_nest.Create( 'step_current_generator',n=1)  
                        rec=my_nest.GetStatus([id])[0]['receptor_types']
                        my_nest.SetStatus(scg, {'amplitude_times':params_msn_d2['mod_times'],
                                                'amplitude_values':[0.,-300.,0.]})
                        my_nest.Connect( scg, [id],  params = { 'receptor_type' : rec['CURR'] } )
                    
                
                I_e=my_nest.GetDefaults(model[1])['I_e']  
                if  I_E_VARIATION[name]:I=numpy.random.normal(I_e, I_E_VARIATION[name])
                else:I=I_e
                #I=I_e
                my_nest.SetStatus([id],{'I_e':I})
        layer_dic[name]=layer
    
    mm=nest.Create('multimeter', 1)
    recodables=['V_m', 'I', 'g_AMPA', 'g_NMDA', 'g_GABAA_1', 'g_GABAA_2']
    
    my_nest.SetStatus(mm, {'interval': 0.1, 'record_from': recodables})
    my_nest.Connect(mm, [layer_dic['STN'].ids[0]])
    
    # Connect populations
    for conn in connect_list:
        name=conn[0] + '_' + conn[1]
        my_nest.MyLoadModels( model_dict, [conn[2]['synapse_model']] )
        
        
        if dis_conn_GPE_STN=='GPE' and (name in ['GPE_SNR']):
            r,syn=32*30.0, 'GPE_SNR_gaba_s_ref'
            if not syn in my_nest.Models():    
                my_nest.MyLoadModels( model_dict, [syn] )
            pg=my_nest.Create('poisson_generator',1,{'rate':r,'start':1.})
            my_nest.DivergentConnect(pg, layer_dic[conn[1]].ids, model=syn)
        elif dis_conn_GPE_STN=='STN' and (name in ['STN_SNR']):
            r, syn=30*10.0, 'STN_SNR_ampa_s'
            if not syn in my_nest.Models():               
                my_nest.MyLoadModels( model_dict, [syn] )
            pg=my_nest.Create('poisson_generator',1,{'rate':r,'start':1.})
            my_nest.DivergentConnect(pg, layer_dic[conn[1]].ids, model=syn)
            
            
        else:
            name=name+'_'+conn[3]    
            tp.ConnectLayers(layer_dic[conn[0]].layer_id, layer_dic[conn[1]].layer_id, conn[2])
            layer_dic[conn[1]].add_connection(source=layer_dic[conn[0]], type=conn[3], props=conn[2])
    
    # Sort MSN D2 such that the closest to center is first in ids list.
    # Do this to we can get focused inhibition in GPe    
    if params_msn_d2['focus']:
        MSN_D2_idx=layer_dic['MSN_D2'].sort_ids()
    else:
        MSN_D2_idx=range(len(numpy.array(layer_dic['MSN_D2'].ids)))
        
    n_mod_msn_d1=params_msn_d1['n_mod']
    n_mod_msn_d2=params_msn_d2['n_mod']
 
    MSN_D1_ids=layer_dic['MSN_D1'].ids
    MSN_D2_ids=layer_dic['MSN_D2'].ids
    
    MSN_D1_mod,MSN_D2_mod=[],[]
    if params_msn_d1['n_mod']:MSN_D1_mod=MSN_D1_ids[0:n_mod_msn_d1]
    if params_msn_d2['n_mod']:MSN_D2_mod=MSN_D2_ids[0:n_mod_msn_d2*params_msn_d2['skip']:params_msn_d2['skip']]
    
       
    MSN_D1_base=list(set(MSN_D1_ids).difference(MSN_D1_mod))
    MSN_D2_base=list(set(MSN_D2_ids).difference(MSN_D2_mod))
    
    #layer_dic['MSN_D1'].ids[0:n_base_msn_d1]
    
    #MSN_D2_ids=numpy.array(layer_dic['MSN_D2'].ids)
    #MSN_D2_base=MSN_D2_ids#[MSN_D2_idx[0:n_base_msn_d1]]
    
    
    #set().difference(t) 
    
    layer_dic['MSN_D1'].set_spike_times(params_msn_d1['base_rates'], params_msn_d1['base_times'], sim_time, ids=MSN_D1_base)
    layer_dic['MSN_D2'].set_spike_times(params_msn_d2['base_rates'], params_msn_d2['base_times'], sim_time, ids=MSN_D2_base)
    
    if params_msn_d1['n_mod']: layer_dic['MSN_D1'].set_spike_times(params_msn_d1['mod_rates'], params_msn_d1['mod_times'], sim_time, ids=MSN_D1_mod)        
    if params_msn_d2['n_mod']: layer_dic['MSN_D2'].set_spike_times(params_msn_d2['mod_rates'], params_msn_d2['mod_times'], sim_time, ids=MSN_D2_mod) 
     
    STN_CTX_input_base=my_nest.Create('poisson_generator',params={'rate':BASE_RATE_CTX_STN, 'start':0., 'stop':sim_time})
    my_nest.MyLoadModels( model_dict, ['CTX_STN_ampa_s'] )
    my_nest.DivergentConnect(STN_CTX_input_base, layer_dic['STN'].ids, model='CTX_STN_ampa_s')
    
    
    if params_stn['mod']:
        STN_CTX_input_mod=my_nest.Create('poisson_generator',params={'rate':params_stn['mod_rate'], 
                                                                     'start':params_stn['mod_times'][0],
                                                                     'stop':params_stn['mod_times'][1]})
        my_nest.DivergentConnect(STN_CTX_input_mod, layer_dic['STN'].ids, model='CTX_STN_ampa_s')
    
    #tar=[]
    #for id in layer_dic['MSN_D1'].ids:
    #    tar.extend(sorted(nest.GetStatus(my_nest.FindConnections([id]),'target'))[:-1])
        
    #pylab.subplot(211).hist(tar, 1500) 
#    
#    tar=[]
#    for id in layer_dic['MSN_D2'].ids:
#        tar.extend(sorted(nest.GetStatus(my_nest.FindConnections([id]),'target'))[1:])
#        
#    pylab.subplot(212).hist(tar, 1500) 
    #pylab.show()   
#    
#    
    my_nest.MySimulate(sim_time)    
       
    if params_msn_d1['n_mod']:layer_dic['MSN_D1'].id_mod=MSN_D1_mod
    if params_msn_d2['n_mod']:layer_dic['MSN_D2'].id_mod=MSN_D2_mod
        
    #layer_dic['MSN_D1'].get_signal( 's', start=start_rec, stop=sim_time )
    #layer_dic['MSN_D2'].get_signal( 's', start=start_rec, stop=sim_time )
    #layer_dic['GPE'].get_signal( 's', start=start_rec, stop=sim_time )    
    #layer_dic['SNR'].get_signal( 's', start=start_rec, stop=sim_time )    
    #layer_dic['STN'].get_signal( 's', start=start_rec, stop=sim_time )    
    
    st_mm=my_nest.GetStatus(mm)[0]
    pylab.plot(st_mm['events']['g_AMPA'])
    pylab.plot(st_mm['events']['g_GABAA_1'])
    pylab.plot(st_mm['events']['g_NMDA'])
    pylab.plot(st_mm['events']['g_GABAA_2'])
    m_ampa=numpy.mean(st_mm['events']['g_AMPA'])
    m_gaba=numpy.mean(st_mm['events']['g_GABAA_1'])
    pylab.title("{0} m_ampa:{1:2.1f} m_gaba:{2:2.1f}".format(my_nest.version(),m_ampa,m_gaba ))
    pylab.show()
    return layer_dic


def simulate_network_poisson(params_msn_d1, params_msn_d2, params_stn,
                     synapse_models, sim_time, seed, I_e_add, threads=1, 
                     start_rec=0, model_params={}, params_in={}, p_weights=False, 
                     p_conn=False, p_I_e=False):    
    '''
    
    Assume that the background MSN are static weak, then can use poisson process 
    for them,
        params_msn_d1 - dictionary with timing and burst freq setup for msn
                     {'base_rates':0.1, 
                      'base_times':[1], 
                      'mod_rates': 20,
                      'mod_times':[1,200], 
                      'mod_units':list()
                      'n_tot':500, 
                       n_mod=20}
        params_msn_d2 - dictionary with timing and burst freq setup for gpe
        params_stn    - dictionary {'rate':50}
                     same as params_msn
        neuron_model - string, the neuron model to use 
        synapse_models - dict, {'MSN':'...', 'GPE':,'...', 'STN':'...'}
        sim_time - simulation time
        seed - seed for random generator
        I_e_add - diabled
        start_rec - start recording from
        model_params - general model paramters
    '''
    
    params = {'conns':{'MSN_D1_SNR':{'syn':synapse_models[0]},   
                       'GPE_SNR':{'syn':synapse_models[1]}}}  
    
    
    my_nest.ResetKernel(threads=8) 
    numpy.random.seed(seed)
    
 
    
    params=misc.dict_merge(model_params, params)
    params=misc.dict_merge({'neurons':{'GPE':{'paused':0}}}, params)
               
    model_list, model_dict = models({}, p_weights)
    layer_list, connect_list = network(model_dict, params, p_conn)
    
    dic_p_I_e={'SNR':1.,'GPE':1.,'STN':1.}
    if p_I_e is not False:
        dic_p_I_e['SNR']*=p_I_e[0]
        dic_p_I_e['GPE']*=p_I_e[1]
        dic_p_I_e['STN']*=p_I_e[2]
            
    # Create neurons and synapses
    layer_dic={}  
    for name, model, props  in layer_list:
        
        # Update input current
        my_nest.MyLoadModels( model_dict, [model[1]] )
        if name in I_IN_VIVO.keys():
            I_in_vitro=my_nest.GetDefaults(model[1])['I_e']
            I_e=I_in_vitro+I_IN_VIVO[name]
            my_nest.SetDefaults(model[1],{'I_e':I_e*dic_p_I_e[name]})  
                
        #! Create layer, retrieve neurons ids per elements and p
        if model[0]=='spike_generator':
            layer=MyLayerPoissonInput(layer_props=props, sd=True, 
                                      sd_params={'start':start_rec, 
                                                 'stop':sim_time})
        elif model[0]=='poisson_generator':
            layer=MyPoissonInput(model[0],props['columns'],
                                 sd=True, sd_params={'start':start_rec, 'stop':sim_time})
        
        
        else:
            layer=MyLayerGroup(layer_props=props, sd=True, mm=False, mm_dt = 0.1,
                               sd_params={'start':start_rec, 'stop':sim_time} )
               
            
            for iter, id in enumerate(layer[:]):

                if name=='GPE' and params_msn_d2['n_mod'] and iter<params['neurons']['GPE']['paused']:
                        scg = my_nest.Create( 'step_current_generator',n=1)  
                        rec=my_nest.GetStatus([id])[0]['receptor_types']
                        my_nest.SetStatus(scg, {'amplitude_times':params_msn_d2['mod_times'],
                                                'amplitude_values':[0.,-300.,0.]})
                        my_nest.Connect( scg, [id],  params = { 'receptor_type' : rec['CURR'] } )
                    
                
                I_e=my_nest.GetDefaults(model[1])['I_e']              
                if  I_E_VARIATION[name]:I=numpy.random.normal(I_e, I_E_VARIATION[name]) #I_E_VARIATION[name])
                else:I=I_e
                my_nest.SetStatus([id],{'I_e':I})
                
        layer_dic[name]=layer
    
    
    # Connect populations
    for conn in connect_list:
        print [conn[2]['synapse_model']]
        if not conn[2]['synapse_model'] in nest.Models():
            my_nest.MyLoadModels( model_dict, [conn[2]['synapse_model']] )
        
        if layer_dic[conn[0]].model == 'poisson_generator':
            my_nest.Connect(layer_dic[conn[0]].ids, layer_dic[conn[1]].ids,
                                model=conn[2]['synapse_model'])
        else:
           
            name=conn[0] + '_' + conn[1]+'_'+conn[3]    
            tp.ConnectLayers(layer_dic[conn[0]].layer_id, layer_dic[conn[1]].layer_id, conn[2])
            layer_dic[conn[1]].add_connection(source=layer_dic[conn[0]], type=conn[3], props=conn[2])
    
    # Sort MSN D2 such that the closest to center is first in ids list.
    # Do this to we can get focused inhibition in GPe
    
    if params_msn_d2['focus']:
        MSN_D2_idx=layer_dic['MSN_D2'].sort_ids()
    else:
        MSN_D2_idx=range(len(numpy.array(layer_dic['MSN_D2'].ids)))
        
    n_mod_msn_d1=params_msn_d1['n_mod']
    n_mod_msn_d2=params_msn_d2['n_mod']
 
    MSN_D1_ids=layer_dic['MSN_D1'].ids
    MSN_D2_ids=layer_dic['MSN_D2'].ids
    
    MSN_D1_mod,MSN_D2_mod=[],[]
    if params_msn_d1['n_mod']:MSN_D1_mod=MSN_D1_ids[0:n_mod_msn_d1]
    if params_msn_d2['n_mod']:MSN_D2_mod=MSN_D2_ids[0:n_mod_msn_d2*params_msn_d2['skip']:params_msn_d2['skip']]  
    
    MSN_D1_base=list(set(MSN_D1_ids).difference(MSN_D1_mod))
    MSN_D2_base=list(set(MSN_D2_ids).difference(MSN_D2_mod))
    
    
    layer_dic['MSN_D1'].set_spike_times(params_msn_d1['base_rates'], params_msn_d1['base_times'], sim_time, ids=MSN_D1_base)
    layer_dic['MSN_D2'].set_spike_times(params_msn_d2['base_rates'], params_msn_d2['base_times'], sim_time, ids=MSN_D2_base)
    
    if params_msn_d1['n_mod']: layer_dic['MSN_D1'].set_spike_times(params_msn_d1['mod_rates'], params_msn_d1['mod_times'], sim_time)        
    if params_msn_d2['n_mod']: layer_dic['MSN_D2'].set_spike_times(params_msn_d2['mod_rates'], params_msn_d2['mod_times'], sim_time, ids=MSN_D2_mod) 
    
    # If background poisson are use
    if params_msn_d1['bg_rate']:layer_dic['MSN_D1_bg'].set_spike_times(params_msn_d1['bg_rate'], [1.], sim_time)
    if params_msn_d2['bg_rate']:layer_dic['MSN_D2_bg'].set_spike_times(params_msn_d2['bg_rate'], [1.], sim_time)
    
    STN_CTX_input_base=my_nest.Create('poisson_generator',params={'rate':BASE_RATE_CTX_STN, 'start':0., 'stop':sim_time})
    my_nest.MyLoadModels( model_dict, ['CTX_STN_ampa_s'] )
    
    if 'STN' in layer_dic.keys(): my_nest.DivergentConnect(STN_CTX_input_base, layer_dic['STN'].ids, model='CTX_STN_ampa_s')
    
    
    if params_stn['mod'] and 'STN' in layer_dic.keys():
        STN_CTX_input_mod=my_nest.Create('poisson_generator',params={'rate':params_stn['mod_rate'], 
                                                                     'start':params_stn['mod_times'][0],
                                                                     'stop':params_stn['mod_times'][1]})
        my_nest.DivergentConnect(STN_CTX_input_mod, layer_dic['STN'].ids, model='CTX_STN_ampa_s')
    
    my_nest.MySimulate(sim_time)    
        
    if params_msn_d1['n_mod']:layer_dic['MSN_D1'].id_mod=MSN_D1_mod
    if params_msn_d2['n_mod']:layer_dic['MSN_D2'].id_mod=MSN_D2_mod
        
    if 'MSN_D1' in layer_dic.keys():
        layer_dic['MSN_D1'].get_signal( 's', start=start_rec, stop=sim_time )
    if 'MSN_D2' in layer_dic.keys():
        layer_dic['MSN_D2'].get_signal( 's', start=start_rec, stop=sim_time )
    if 'GPE' in layer_dic.keys():
        layer_dic['GPE'].get_signal( 's', start=start_rec, stop=sim_time )    
    if 'SNR' in layer_dic.keys():
        layer_dic['SNR'].get_signal( 's', start=start_rec, stop=sim_time )    
    if 'STN' in layer_dic.keys():
        layer_dic['STN'].get_signal( 's', start=start_rec, stop=sim_time )    
    
    return layer_dic

     
def simulate_network_indirect_sep_freq(skip,focus, mod, freq, params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, sim_time=2000., seed=1,
                           I_e_add={'SNR':300, 'STN':0,'GPE':20}, threads=4, start_rec=0):
    
    
    
    spk=[]
    spk2=[]
    spk3=[]
    for f, m, se, fo, sk in zip(freq, mod, seed, focus, skip):
        
        params_msn_d2.update({ 'mod_rates': [0.1, f, 0.1],
                               'mod_times':[1,1000, 1000+500],  
                               'n_mod':int(m),
                               'focus':fo,
                               'skip':sk})
        layer_dic=simulate_network(params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, sim_time=sim_time, seed=se,
                           I_e_add=I_e_add, threads=4,start_rec=start_rec,
                           model_params={'misc':{'N_MSN':N_MSN}}) 
    
        layer_dic['SNR'].get_signal( 's', start=start_rec, stop=sim_time )
        signal=layer_dic['SNR'].signals['spikes']
        spk_mean=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
        spk.append(spk_mean)
        
        layer_dic['GPE'].get_signal( 's', start=start_rec, stop=sim_time )
        signal=layer_dic['GPE'].signals['spikes']
        spk_mean=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
        spk2.append(spk_mean)
        
        layer_dic['STN'].get_signal( 's', start=start_rec, stop=sim_time )
        signal=layer_dic['STN'].signals['spikes']
        spk_mean=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
        spk3.append(spk_mean)
        
        
    return spk, spk2, spk3    
 

def simulate_network_direct_indirect_onoff_vs_rate(m_d1, m_d2, params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, I_e_add,sim_time=2000., seed=1,
                            threads=4, start_rec=0,
                           model_params={}, flag_bg=False):
    
    
    if flag_bg:
        N_MSN=model_params['neurons']['MSN_D1']['n']
        
        params_msn_d1.update({'bg_rate':[0.1*(500-500*m_d1/float(N_MSN))]})
        params_msn_d2.update({'bg_rate':[0.1*(500-500*m_d2/float(N_MSN))]})
        
        # Change paramters to use poisson background
        params={'conns':{'MSN_D1_SNR':{'p':500./float(N_MSN)},
                         'MSN_D2_GPE':{'p':500./float(N_MSN),
                              'lines':False}},
                      'neurons':{'MSN_D1':{'n':m_d1},
                                 'MSN_D2':{'n':m_d2},
                                 'MSN_D1_bg':{'n':300,
                                              'lesion':False},
                                 'MSN_D2_bg':{'n':300,
                                              'lesion':False},  
                                 'GPE': {'paused':False}}}
        model_params=misc.dict_merge(model_params, params)
    
    
    params_msn_d1.update({'n_mod':int(m_d1)})
    params_msn_d2.update({'n_mod':int(m_d2)})
    
    
    layer_dic=simulate_network_poisson(params_msn_d1, params_msn_d2, params_stn,
                       synapse_models, sim_time=sim_time, seed=seed,
                       I_e_add=I_e_add, threads=4,start_rec=start_rec,
                       model_params=model_params) 

    layer_dic['SNR'].get_signal( 's', start=start_rec, stop=sim_time )
    signal=layer_dic['SNR'].signals['spikes']
    spk_mean=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
    spk=spk_mean

    layer_dic['GPE'].get_signal( 's', start=start_rec, stop=sim_time )
    signal=layer_dic['GPE'].signals['spikes']
    spk_mean=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
    spk2=spk_mean
        
    return spk, spk2 

def simulate_network_hyper_direct_onoff_vs_rate(freq, params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, sim_time=2000., seed=1,
                           I_e_add={'SNR':300, 'STN':0,'GPE':20}, threads=4, start_rec=0,
                           model_params={}):
    
    
    
    spk, spk2, spk3=[],[],[]
    for f, se in zip(freq, seed):
        
        params_stn.update({'rate':350., 'mod':True,'mod_rate':f})
        
        layer_dic=simulate_network(params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, sim_time=sim_time, seed=se,
                           I_e_add=I_e_add, threads=4,start_rec=start_rec,
                           model_params=model_params) 
    
        layer_dic['SNR'].get_signal( 's', start=start_rec, stop=sim_time )
        signal=layer_dic['SNR'].signals['spikes']
        spk_mean=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
        spk.append(spk_mean)

        layer_dic['GPE'].get_signal( 's', start=start_rec, stop=sim_time )
        signal=layer_dic['GPE'].signals['spikes']
        spk_mean=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
        spk2.append(spk_mean)
        
        layer_dic['STN'].get_signal( 's', start=start_rec, stop=sim_time )
        signal=layer_dic['STN'].signals['spikes']
        spk_mean=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
        spk3.append(spk_mean)
        
    return spk, spk2, spk3 



def inspect_network():    
    model_list, model_dict = models()
    layer_list, connect_list = network(model_dict, {'misc':{'N_MSN':15000},
                                                    'conns':{'MSN_D2_GPE':{'lines':False}}})
    
    
    # Create neurons and synapses
    layer_dic={}  
    for name, model, props  in layer_list:
        my_nest.MyLoadModels( model_dict, [model[1]] )

    
        #! Create layer, retrieve neurons ids per elements and p
        if model[0]=='spike_generator':
            layer=MyLayerPoissonInput(layer_props=props, sd=False)
        else:  
            layer=MyLayerGroup(layer_props=props, sd=True, mm=True, mm_dt = 0.1 )
        layer_dic[name]=layer
    
    
    # Connect populations
    for conn in connect_list:
        my_nest.MyLoadModels( model_dict, [conn[2]['synapse_model']] )
        name=conn[0] + '_' + conn[1]+'_'+conn[3]    
        tp.ConnectLayers(layer_dic[conn[0]].layer_id, layer_dic[conn[1]].layer_id, conn[2])
        layer_dic[conn[1]].add_connection(source=layer_dic[conn[0]], type=conn[3], props=conn[2])
        
    return layer_dic
    
    
