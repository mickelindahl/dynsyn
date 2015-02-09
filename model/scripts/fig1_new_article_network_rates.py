import copy
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
model_name=os.getcwd().split('/')[-2]
picture_dir='/'.join(os.getcwd().split('/')[0:-3]) + '/pictures/'+model_name 
SNAME  = sys.argv[0].split('/')[-1].split('.')[0]

OUTPUT_PATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]

def plot_control(ax, data, data2):
    y=data[0]
    y_std=data[1]
    y_control=copy.deepcopy(y[0,:])

    #N=len(xticklabels)
    ind = 0
    width=0.4
    alpha=0.5
    
    
    x=range(len(y))
    rects1=ax.bar(0, y[:1,0], width, yerr=y_std[:1,0], color='b', alpha=alpha )
    rects2=ax.bar(1, y[:1,1], width, yerr=y_std[:1,1], color='g',alpha=alpha)
    rects3=ax.bar(2, y[:1,2], width,yerr=y_std[:1,2],  color='r',alpha=alpha)

    alpha=1.
    y=data2[0]
    y_std=data2[1]
    y_control=copy.deepcopy(y[0,:])

    
    x=range(len(y))
    rects1=ax.bar(0+width, y[:1,0], width, yerr=y_std[:1,0], color='b', alpha=alpha )
    rects2=ax.bar(1+width, y[:1,1], width, yerr=y_std[:1,1], color='g',alpha=alpha)
    rects3=ax.bar(2+width, y[:1,2], width,yerr=y_std[:1,2],  color='r',alpha=alpha)
    
    
    ax.set_ylabel('Rate (spikes/s')
    #ax.set_xlabel('Nucleus')
    #ax.set_title('Control')
    
    labels=['Shaded: static STN-SNr','Solid: depressing STN-SNr']
    colors=['k','k']
    coords=[[0.02, 0.7], [ 0.02, 0.8]]
    for label, coord, color in zip(labels,coords,colors):
            ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['font.size']-7, 
                 **{'color': color})
    
    ax.set_xticks( numpy.arange(0.4, 3 + 0.4, 1) )
    ax.set_xticklabels( ['SNr','GPe', 'STN'], rotation=55, ha='right')
    ax.set_ylim(misc.adjust_limit([0,70]))
              
def plot_control_only_one(ax, data, data2):
    y=data[0]
    y_std=data[1]
    y_control=copy.deepcopy(y[0,:])

    #N=len(xticklabels)
    ind = 0
    width=0.8
    alpha=1.
    
    
    x=range(len(y))
    rects1=ax.bar(0, y[:1,0], width, color='b', alpha=alpha )
    rects2=ax.bar(1, y[:1,1], width, yerr=y_std[:1,1], color='g',alpha=alpha)
    rects3=ax.bar(2, y[:1,2], width,yerr=y_std[:1,2],  color='r',alpha=alpha)
    (_, caplines, _) =ax.errorbar(0.4, y[:1,0], yerr=y_std[:1,0], fmt='o', color='k', markersize=5, linewidth=1.0,  capsize=5.0,markeredgewidth=1.0 )
    (_, caplines, _) =ax.errorbar(1.4, y[:1,1], yerr=y_std[:1,1], fmt='o', color='k', markersize=5, linewidth=1.0,  capsize=5.0,markeredgewidth=1.0 )
    (_, caplines, _) =ax.errorbar(2.4, y[:1,2], yerr=y_std[:1,2], fmt='o', color='k', markersize=5, linewidth=1.0,  capsize=5.0,markeredgewidth=1.0 )
    
    #for cap in caplines:
    #    cap.set_linewidth(20)
    #    cap.set_markeredgewidth(2)
    pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend
    
    
    alpha=1.
    y=data2[0]
    y_std=data2[1]
    y_control=copy.deepcopy(y[0,:])

    '''
    x=range(len(y))
    rects1=ax.bar(0+width, y[:1,0], width, yerr=y_std[:1,0], color='b', alpha=alpha )
    rects2=ax.bar(1+width, y[:1,1], width, yerr=y_std[:1,1], color='g',alpha=alpha)
    rects3=ax.bar(2+width, y[:1,2], width,yerr=y_std[:1,2],  color='r',alpha=alpha)
    '''
    
    ax.set_ylabel('Rate spike/s')
    #ax.set_xlabel('Nucleus')
    #ax.set_title('Control')
    '''
    labels=['Shaded: static STN-SNr','Solid: depressing STN-SNr']
    colors=['k','k']
    coords=[[0.02, 0.7], [ 0.02, 0.8]]
    for label, coord, color in zip(labels,coords,colors):
            ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['font.size']-7, 
                 **{'color': color})
    '''
    ax.set_xticks( numpy.arange(0.4, 3 + 0.4, 1) )
    ax.set_xticklabels( ['SNr','GPe', 'STN'], rotation=0, ha='center')
    ax.set_ylim(misc.adjust_limit([0,40]))
                             
def plot_lesion(ax, data, data2, xticklabels):
    
    y=copy.deepcopy(data[0])
    y_std=copy.deepcopy(data[1])
    y_control=copy.deepcopy(y[0,:])
    
    for i in range(y.shape[0]):
        y[i,:]/=y_control
        #y[i,:]/=y_control
        y[i,:]*=100.0
    for i in range(y.shape[0]):
        y_std[i,:]/=y_control
        #y[i,:]/=y_control
        y_std[i,:]*=100.0   
        
        
    y-=100    
    N=len(xticklabels)
    ind = numpy.arange(N)
    width=0.15
    alpha=0.5
    
    
    x=range(len(y))
    rects1=ax.bar(ind, y[1:,0], width,  color='b', alpha=alpha )
    rects2=ax.bar(ind+width*2, y[1:,1], width,  color='g',alpha=alpha)
    rects3=ax.bar(ind+width*4, y[1:,2], width,color='r',alpha=alpha)

    alpha=1.
    y=copy.deepcopy(data2[0])
    y_std=copy.deepcopy(data2[1])
    y_control=copy.deepcopy(y[0,:])
    
    for i in range(y.shape[0]):
        y[i,:]/=y_control
        #y[i,:]/=y_control
        y[i,:]*=100.0
    y-=100    
    
    x=range(len(y))
    rects1=ax.bar(ind+width, y[1:,0], width, color='b', alpha=alpha )
    rects2=ax.bar(ind+width*3, y[1:,1], width, color='g',alpha=alpha)
    rects3=ax.bar(ind+width*5, y[1:,2], width,  color='r',alpha=alpha)
    
    
    ax.set_ylabel('Rate change (%)')
    ax.set_xlabel('Lesioned nucleus')
    labels=['$SNr$','$GPe$','$STN$']
    colors=['b','g','r']
    coords=[[0.05, 0.5], [ 0.05, 0.6], [ 0.05, 0.7]]
    for label, coord, color in zip(labels,coords,colors):
            ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['font.size']-2, 
                 **{'color': color})
   
   
    labels=['Shaded: static \nSTN-SNr','Solid: depressing \nSTN-SNr']
    colors=['k','k']
    coords=[[0.4, 0.6], [ 0.4, 0.8]]
    for label, coord, color in zip(labels,coords,colors):
            ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['font.size']-7, 
                 **{'color': color})
    ax.set_xticks( numpy.arange(0.4,len(y)-1+0.4,1) )
    ax.set_xticklabels( xticklabels, rotation=0, ha='center')
    ax.set_ylim(misc.adjust_limit([-100,400]))


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
                                       I_e_add, threads, start_rec,
                                       model_params):
    N_MSN=15000
    
    params={'misc':{'N_MSN':N_MSN},
              'conns':{'MSN_D1_SNR':{'p':500./float(N_MSN)},
                       'MSN_D2_GPE':{'p':500./float(N_MSN),
                                     'lines':False}},
              'neurons':{'MSN_D1':{'n':0},
                         'MSN_D2':{'n':0},
                         'MSN_D1_bg':{'n':300, 'lesion':False},
                         'MSN_D2_bg':{'n':300, 'lesion':False},  
                         'GPE': {'paused':False}}}
    
    model_params=misc.dict_merge(params, model_params)
    
    layer_dic=simulate_network_poisson(params_msn_d1, params_msn_d2, params_stn,
                                       synapse_models, sim_time, seed,
                                       I_e_add, threads, start_rec, 
                                       model_params)    

    mr=[]
    std=[]
    if 'SNR' in layer_dic.keys():
        signal=layer_dic['SNR'].signals['spikes']
        m_hist=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=1)
        mr.append(numpy.mean(m_hist, axis=0))
        std.append(numpy.std(m_hist, axis=0))
    else: 
        mr.append(0)
        std.append(0)
    
    if 'GPE' in layer_dic.keys():
        signal=layer_dic['GPE'].signals['spikes']
        m_hist=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=1)
        mr.append(numpy.mean(m_hist, axis=0))
        std.append(numpy.std(m_hist, axis=0))
    else:
        mr.append(0)
        std.append(0)
        
    if 'STN' in layer_dic.keys():
        signal=layer_dic['STN'].signals['spikes']
        m_hist=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=1)
        mr.append(numpy.mean(m_hist, axis=0))
        std.append(numpy.std(m_hist, axis=0))
    else:        
        mr.append(0)
        std.append(0)
    return layer_dic, mr, std

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
    if 'SNR' in layer_dic.keys():
        signal=layer_dic['SNR'].signals['spikes']
        m_hist=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
        mr.append(numpy.mean(m_hist, axis=0))
    else: mr.append(0)
    
    if 'GPE' in layer_dic.keys():
        signal=layer_dic['GPE'].signals['spikes']
        m_hist=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
        mr.append(numpy.mean(m_hist, axis=0))
    else:mr.append(0)
    if 'STN' in layer_dic.keys():
        signal=layer_dic['STN'].signals['spikes']
        m_hist=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
        mr.append(numpy.mean(m_hist, axis=0))
    else:mr.append(0)
    
    return layer_dic, mr

def simulate_1500(params_msn_d1, params_msn_d2, params_stn, synapse_models, sim_time, seed,
                                       I_e_add, threads, start_rec):
    N_MSN=1500
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
    if 'SNR' in layer_dic.keys():
        signal=layer_dic['SNR'].signals['spikes']
        m_hist=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
        mr.append(numpy.mean(m_hist, axis=0))
    else: mr.append(0)
    
    if 'GPE' in layer_dic.keys():
        signal=layer_dic['GPE'].signals['spikes']
        m_hist=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
        mr.append(numpy.mean(m_hist, axis=0))
    else:mr.append(0)
    if 'STN' in layer_dic.keys():
        signal=layer_dic['STN'].signals['spikes']
        m_hist=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
        mr.append(numpy.mean(m_hist, axis=0))
    else:mr.append(0)
    
    return layer_dic, mr

def simulate_poisson_15000_eval(load, save_at, n_exp, params_msn_d1, params_msn_d2, params_stn,
                                       synapse_models, sim_time, seed,
                                       I_e_add, threads, start_rec, model_params={}):

    rates=[]
    stds=[]
    if not load:
        for i in range(n_exp):
            seed=i
            layer_dic, r, std=simulate_poisson_15000(params_msn_d1, params_msn_d2, params_stn,
                                       synapse_models, sim_time, seed,
                                       I_e_add, threads, start_rec, model_params)
            rates.append(r)
            stds.append(std)
        
        rates=numpy.array(rates)    
        stds=numpy.array(stds)    

        misc.pickle_save([rates, stds],save_at)   
    else:
        
        rates, stds=misc.pickle_load(save_at)   
    
    mr=numpy.mean(rates, axis=0)
    stdr=numpy.mean(stds, axis=0)
    
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
                                       I_e_add, threads, start_rec, model_param)
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

def simulate_lesion(load, save_at, lesions, threads, stn_syn):
    
    n_exp=5
    
    data=[]
    br=0.1
    bg_rate=br*500.0
    
    params_msn_d1={'base_rates':[br], 'base_times':[1], 'mod_rates': [],
                'mod_times':[],  'n_mod':0, 
                'bg_rate':[bg_rate]}    
    params_msn_d2={'base_rates':[br], 'base_times':[1], 'mod_rates': [],
                'mod_times':[],  'n_mod':0, 'focus':False, 
                'skip':1, 'bg_rate':[bg_rate]} 
    
    params_stn={'rate':210., 'mod':False,'mod_rate':0., 'mod_times':[]} 
    synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p']
    
    sim_time=2500.
    seed=1
    start_rec=500.0
    
    
    for l in lesions:
        model_params={}
        model_params['conns']={'STN_SNR':{'syn':stn_syn}}
        
        if l=='Control':
            pass
        if l=='STN lesion':
            model_params['neurons']= {'STN':{'lesion':True}}

        if l=='GPe lesion':
            model_params['neurons']= {'GPE':{'lesion':True}}
        
        # If lesion MSN neurons, then this must be handeled in 
        # simulate_network, which is currently not the case
        if l=='MSN D1 lesion': 
            model_params['conns']= {'MSN_D1_SNR':{'lesion':True},
                                    'MSN_D1_bg_SNR':{'lesion':True}}

        if l=='MSN D2 lesion': 
            model_params['conns']= {'MSN_D2_GPE':{'lesion':True},
                                    'MSN_D2_bg_GPE':{'lesion':True}}        
        
        save_at_tmp=save_at+l
        r_m, r_std=simulate_poisson_15000_eval(load, save_at_tmp, n_exp, params_msn_d1, 
                                               params_msn_d2, params_stn,
                                               synapse_models, sim_time, seed,
                                               {}, threads, start_rec, model_params)
        data.append([r_m, r_std])
        
    data=zip(*data)
    
    for i in range(len(data)):
             data[i]=numpy.array(data[i] )     
        
    return data
        


#synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p_stoc']
#synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_s_ref']


threads=4
stn_syn='STN_SNR_ampa_s'
lesions=['Control','STN lesion','GPe lesion','MSN D1 lesion', 'MSN D2 lesion']
save_at=OUTPUT_PATH+'/simulate_lesion'+ stn_syn +'.plk'
data=simulate_lesion(1, save_at, lesions, threads, stn_syn)
#ax=pylab.subplot(211)



threads=4
stn_syn='STN_SNR_ampa_p3'
lesions=['Control','STN lesion','GPe lesion','MSN D1 lesion', 'MSN D2 lesion']
save_at=OUTPUT_PATH+'/simulate_lesion'+ stn_syn +'.plk'
data2=simulate_lesion(1, save_at, lesions, threads, stn_syn)
#ax=pylab.subplot(212)



#pylab.show()



 #Inspect results
plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0+275.0, fontsize=16)
font_size_text = 10
fig = pylab.figure( facecolor = 'w' )
ax_list = []
ax_list.append( MyAxes(fig, [ .075, .37, .135, .26 ] ) )    # text box
ax_list.append( MyAxes(fig, [ .26,  .7,  .165, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .53,  .7,  .165, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .7,  .165, .2 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .26,  .4,  .165, .2 ] ) )    #     
#ax_list.append( MyAxes(fig, [ .53,  .4,  .165, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .4,  .165, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .26,  .1,  .165, .2 ] ) )    #     
ax_list.append( MyAxes(fig, [ .53,  .1,  .165, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .1,  .165, .2 ] ) )    # 
ax=ax_list[1]
plot_control_only_one(ax, data, data2)
ax=ax_list[2]
lesions=['Control','STN','GPe','$MSN_{D1}$', '$MSN_{D2}$']
plot_lesion(ax, data, data2, lesions[1:])
ax=ax_list[3]
plot_control(ax, data, data2)
pylab.show()

# dpi does not matter since svg and pdf are both vectorbased
fig.savefig( picture_dir + '/' + SNAME  + '.svg', format = 'svg') 
fig.savefig( picture_dir + '/' + SNAME  + '.pdf', format = 'pdf')