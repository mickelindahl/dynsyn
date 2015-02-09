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
from simulation_utils import  simulate_network_poisson

model_name=os.getcwd().split('/')[-2]
picture_dir='/'.join(os.getcwd().split('/')[0:-3]) + '/pictures/'+model_name 
SNAME  = sys.argv[0].split('/')[-1].split('.')[0]

OUTPUT_PATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]

def plot_data(ax, data, title, xticklabels, alpha=1, show_label=True):
    
    N=len(xticklabels)
    ind = numpy.arange(N)
    width=0.25
    
    y=data[0]
    y_std=data[1]
    
    x=range(len(y))
    #rects1=ax.bar(ind, y[:,0], width, yerr=y_std[:,0], color='b', alpha=alpha )
    #rects2=ax.bar(ind+width, y[:,1], width, yerr=y_std[:,1], color='g',alpha=alpha)
    #rects3=ax.bar(ind+width*2, y[:,2], width, yerr=y_std[:,2], color='r',alpha=alpha)
    rects1=ax.bar(ind, y[:,0], width, color='b', alpha=alpha )
    rects2=ax.bar(ind+width, y[:,1], width,  color='g',alpha=alpha)
    rects3=ax.bar(ind+width*2, y[:,2], width,  color='r',alpha=alpha)

    
    if alpha==1:
        ax.set_xticks( numpy.arange(0.4,len(y)+0.4,1) )
        if show_label:
            ax.set_xticklabels( xticklabels, rotation=55, ha='right')
        else:
            ax.set_xticklabels( '' ) 
        ax.set_title(title, fontsize=16)
        
        ax.set_ylabel('Rate change (%)')
        ax.set_xlabel('Parameter')
        
        
        labels=['$SNr$','$GPe$','$STN$']
        colors=['b','g','r']
        coords=[[0.1, 0.1], [ 0.1, 0.2], [ 0.1, 0.3]]
        for label, coord, color in zip(labels,coords,colors):
            ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['text.fontsize'], 
                 **{'color': color})
        labels=['Solid: Effect of increasing parameter with 20 %','Shaded: Effect of decreasing parameter with 20 %']
        colors=['k','k']
        coords=[[0.5, 0.6], [ 0.5, 0.8]]
        for label, coord, color in zip(labels,coords,colors):
            ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['text.fontsize']-4, 
                 **{'color': color})
    #if alpha!=1:
       
    ax.set_xlim([0,15])
    ax.set_ylim([-30,35])
    
        
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

def simulate_sensitivity_fun(n_exp, params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, sim_time=1000., seed=2,
                           I_e_add={}, threads=4, 
                           start_rec=500.,model_params={}, p_weights=False, 
                           p_conn=False, p_I_e=False):  
    r=[]
    for i in range(n_exp): 
        seed=i
        layer_dic=simulate_network_poisson(params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, sim_time, seed, I_e_add, threads, 
                           start_rec, model_params, {}, p_weights, p_conn, p_I_e)    
    
    
        layer_dic['SNR'].get_signal( 's', start=start_rec, stop=sim_time )
        layer_dic['GPE'].get_signal( 's', start=start_rec, stop=sim_time )
        layer_dic['STN'].get_signal( 's', start=start_rec, stop=sim_time )
        

        signal=layer_dic['SNR'].signals['spikes']
        m_hist=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
        r.append(numpy.mean(m_hist[1000:], axis=0))
        
        signal=layer_dic['GPE'].signals['spikes']
        m_hist=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
        r.append(numpy.mean(m_hist[1000:], axis=0))
        
        signal=layer_dic['STN'].signals['spikes']
        m_hist=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
        r.append(numpy.mean(m_hist[1000:], axis=0))
    
    r=numpy.array(r)
    r=numpy.array([r[0::3],r[1::3],r[2::3]])
    
    mr=numpy.mean(r,axis=1)
    mstd=numpy.std(r,axis=1)
    return mr, mstd

def simulate_sensitivity(load, save_result_at, n_exp, params_msn_d1, params_msn_d2, 
                         params_stn, synapse_models,model_params, sim_time, 
                         start_rec):
    
    p_weights=numpy.ones(17)
    p_weights_ch_names=[r'$g^{CTX-STN}_{0}$','$g^{GPE-STN}_{0}$','$g^{MSN_{D2}-GPe}_0$',
                         '$g^{STN-GPe}_{0}$','$g^{GPe-GPe}_{0}$','$g^{MSN_{D1}-SNr}_0$',
                         '$g^{STN-SNr}_{0}$','$g^{GPe-SNr}_{0}$']
    
    #p_weights_ch_names=['$g^{STN-SNr}_{0}$','$g^{GPe-SNr}_{0}$']
    p_weights_ch=[6, 8, 9, 10, 12, 13, 14, 16]
    #p_weights_ch=[14,16]
    p_conn=numpy.ones(7)
    p_conn_names=[r'$N_{MSN_{D1}-SNr}$', r'$N_{MSN_{D2}-GPe}$', r'$N_{STN-GPe}$',
                   r'$N_{GPe-GPe}$', r'$N_{GPe-STN}$', r'$N_{STN-SNr}$', 
                   r'$N_{GPe-SNr}$']
    #p_conn_names=[r'$N_{STN-SNr}$',  r'$N_{GPe-SNr}$']
    p_conn_ch=[0, 1, 2, 3, 4 ,5, 6]
    
    
    p_I_e=numpy.ones(3)
    #p_I_e_names=['$I_{In vivo}^{SNr}$', '$I_{In vivo}^{GPe}$', '$I_{In vivo}^{STN}$']
    #p_I_e_ch=[0,1,2]
    
    d=[]
    prop_ch=0.2 # Change each parameter up and down 0.2
    seed=0
    save_result_at_tmp=save_result_at+'weights'
    if not load[0]:
        dd=[]
        seed=2
        mr, mstd=simulate_sensitivity_fun(n_exp, params_msn_d1, params_msn_d2, params_stn,
                              synapse_models, sim_time, seed,
                              {}, threads, 
                              start_rec, model_params, p_weights, p_conn, p_I_e)
        dd=[[-1,0,mr[0],mr[1],mr[2], mstd[0],mstd[1],mstd[2], 0]]
        
        for p in p_weights_ch:
            for prop in [1+prop_ch,1-prop_ch]:
               seed+=1
                
               p_weights_tmp=copy.deepcopy(p_weights)
               p_weights_tmp[p]*=prop
               mr, mstd=simulate_sensitivity_fun(n_exp,params_msn_d1, params_msn_d2, params_stn,
                              synapse_models, sim_time, seed,
                              {}, threads, 
                              start_rec,model_params, p_weights_tmp, p_conn, p_I_e)
               
               dd.append([ 0,0, mr[0], mr[1], mr[2], mstd[0], mstd[1], mstd[2], 0])
        misc.pickle_save(dd,save_result_at_tmp)   
    else:
        dd=misc.pickle_load(save_result_at_tmp)   
    d.extend(dd)
    save_result_at_tmp=save_result_at+'conn'
    if not load[1]:    
        dd=[]
        for p in p_conn_ch:
            for prop in [1+prop_ch,1-prop_ch]:
               seed+=1
               p_conn_tmp=copy.deepcopy(p_conn)
               p_conn_tmp[p]*=prop
               mr, mstd=simulate_sensitivity_fun(n_exp,params_msn_d1, params_msn_d2, params_stn,
                              synapse_models, sim_time, seed,
                              {}, threads, 
                              start_rec, model_params, p_weights, p_conn_tmp, p_I_e)
               
               dd.append([1,0,mr[0],mr[1],mr[2], mstd[0],mstd[1],mstd[2], 0])
            
        misc.pickle_save(dd, save_result_at_tmp)   
    else:
        dd=misc.pickle_load(save_result_at_tmp)
    d.extend(dd)
#    save_result_at_tmp=save_result_at+'I_e'
#    if not load[2]:    
#        d=[]
#        for p in p_I_e_ch:
#            for prop in [1+prop_ch,1-prop_ch]:
#               seed+=1
#               p_I_e_tmp=copy.deepcopy(p_I_e)
#               p_I_e_tmp[p]*=prop
#               
#               mr, mstd=simulate_sensitivity_fun(n_exp,params_msn_d1, params_msn_d2, params_stn,
#                                                      synapse_models, sim_time, seed,
#                                                      {}, threads, start_rec, model_params, 
#                                                      p_weights, p_conn, p_I_e_tmp)
#               
#               d.append([2,0,mr[0],mr[1],mr[2], mstd[0],mstd[1],mstd[2], 0])
#                
#        misc.pickle_save(d,save_result_at_tmp)     
#    else:
#        d.extend(misc.pickle_load(save_result_at_tmp) ) 
        
    
       
    d=numpy.array(d)
    br=d[0,2:5] 
    bstd=d[0,5:8] 
    
    up=d[1::2,2:5]
    down=d[2::2,2:5]
    upstd=d[1::2,5:8]
    downstd=d[2::2,5:8]    
    dp=numpy.abs(up-down)/2.
    
    
    # get as precent change
    for i in range(up.shape[0]):
        up[i,:]-=br
        up[i,:]/=br
        up[i,:]*=100.
        down[i,:]-=br
        down[i,:]/=br    
        down[i,:]*=100.
        
        upstd[i,:]-=bstd
        upstd[i,:]/=bstd
        upstd[i,:]*=100.
        downstd[i,:]-=bstd
        downstd[i,:]/=bstd    
        downstd[i,:]*=100.
          
    data_rate_change=[[up, upstd],[down, downstd]]      
    return numpy.array(d), data_rate_change, dp, p_weights_ch_names + p_conn_names #+ p_I_e_names


threads=4
N_MSN=15000.0
p_mod_msn=0.01
p_mod_msn_d2=0.04
bg_rate=0.1*500.0
params_msn_d1={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [],
            'mod_times':[],  'n_mod':0, 
            'bg_rate':[bg_rate]}    
params_msn_d2={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [],
            'mod_times':[],  'n_mod':0, 'focus':False, 
            'skip':1, 'bg_rate':[bg_rate]} 
params_stn={'rate':219., 'mod':False,'mod_rate':0., 'mod_times':[]} 
start_rec=1000.
sim_time=4000.
n_exp=2
#synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p_stoc']
synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p']
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

save_result_at=OUTPUT_PATH+'/simulate_sensitivity.plk'
d, data_rate_change,dp, names=simulate_sensitivity([1, 1 ], save_result_at,n_exp,
                                           params_msn_d1, params_msn_d2, params_stn,
                                           synapse_models, model_params, sim_time,
                                           start_rec ) 
names=names[0:15]
model_params['conns'].update({'STN_SNR':{'syn':'STN_SNR_ampa_p3'}})
save_result_at=OUTPUT_PATH+'/simulate_sensitivity_STN_plast.plk'
d2, data_rate_change2, dp2, names2=simulate_sensitivity([1, 1], save_result_at,n_exp,
                                           params_msn_d1, params_msn_d2, params_stn,
                                           synapse_models, model_params, sim_time,
                                           start_rec )

 
 #Inspect results
plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0+275.0, fontsize=16)
font_size_text = 8
fig = pylab.figure( facecolor = 'w' )
ax_list = []
ax_list.append( MyAxes(fig, [ .075, .37, .135, .26 ] ) )    # text box
ax_list.append( MyAxes(fig, [ .26,  .7,  .4, .2 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .53,  .7,  .165, .2 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .8,   .7,  .165, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .26,  .3,  .4, .2 ] ) )    #     
#ax_list.append( MyAxes(fig, [ .53,  .4,  .165, .2 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .8,   .4,  .165, .2 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .26,  .1,  .165, .2 ] ) )    #     
#ax_list.append( MyAxes(fig, [ .53,  .1,  .165, .2 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .8,   .1,  .165, .2 ] ) )    # 


ax=ax_list[1]

plot_data(ax, data_rate_change[0], 'Static STN-SNr', names, show_label=False)
ax=ax_list[1]
plot_data(ax, data_rate_change[1], '', names, alpha=0.65)

ax=ax_list[2]
plot_data(ax, data_rate_change2[0], 'Depressing STN-SNr', names)
ax=ax_list[2]
plot_data(ax, data_rate_change2[1], '', names, alpha=0.65)


pylab.show()


# dpi does not matter since svg and pdf are both vectorbased
fig.savefig( picture_dir + '/' + SNAME  + '.svg', format = 'svg') 
fig.savefig( picture_dir + '/' + SNAME  + '.pdf', format = 'pdf')