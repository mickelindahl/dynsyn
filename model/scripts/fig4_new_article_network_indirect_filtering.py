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
from simulation_utils import simulate_network

SELECTION_THR=5.
model_name=os.getcwd().split('/')[-2]
picture_dir='/'.join(os.getcwd().split('/')[0:-3]) + '/pictures/'+model_name 
SNAME  = sys.argv[0].split('/')[-1].split('.')[0]

OUTPUT_PATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]

LIM_SYN_EVENTS=10000.

def plot_selection_vs_neurons(ax, GPEmeanRates, SNRmeanRates):
    colors=['b', 'b' ,'g', 'g' ] 

    linestyles=['-','--','-','--']
    
    
    for i, color in enumerate(colors):
        ax.plot(GPEmeanRates[i,:],SNRmeanRates[i,:],**{'color':color,
                                                       'linestyle':linestyles[i]})  
    
    
    vec=[0,100]

    ax.set_ylabel('Rate SNr (spikes/s)') 
    ax.set_xlabel(r'Rate $MSN_{D2}$ (spikes/s)')
    ax.my_set_no_ticks( yticks=8, xticks=6 )
    
    line1=ax.plot(1,1,'-k')
    line2=ax.plot(1,1,'--k')
    pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend
    leg=ax.legend([line1, line2],['$dep^{GPe}$', '$ref_{32 Hz}^{GPe}$'], loc='best')
    frame  = leg.get_frame() 
    frame.set_visible(False) 
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    pylab.setp(frame, edgecolor='w') 
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    pylab.setp(ltext, fontsize=pylab.rcParams['text.fontsize']-2, backgroundcolor='w')
    
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,100]) 
    
    labels=[ 'Static \nSTN-SNr', 'Depressing \nSTN-SNr']
    coords=[[0.55, 0.3], [ 0.3, 0.6]]
    colors=['b','g'] 
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=font_size_text, 
                 **{'color': color})
    ax.set_xlim(misc.adjust_limit([0,4]))
    ax.set_ylim(misc.adjust_limit([25,150]))          
    
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

def simulate_filtering(load, save_at, interval,syn_STN, params_msn_d1, params_msn_d2, params_stn,
                           synapse_model, sim_time=2000., seed=1,threads=4, 
                            start_rec=0):
    

    freq_filter=numpy.linspace(0.1,4.0,15)
    paused=numpy.zeros(len(freq_filter))
    seed=range(len(freq_filter))
    model_params={'misc':{'N_MSN':N_MSN},
                   'conns':{ 'MSN_D2_GPE':{ 'lines':False},
                               'STN_SNR' :{'syn':syn_STN}},
                   'neurons':{'MSN_D1':{'n':N_MSN},
                              'MSN_D2':{'n':N_MSN},
                              'GPE': {'paused':0}}}      
    
    if not load:
        mr=[]
        for f, pa, se in zip(freq_filter, paused,seed):
            
            params_msn_d2.update({ 'base_rates':[f]})
            
            params_msn_d2.update({'focus':False,
                                   'skip':1})
            
            model_params['neurons']['GPE']['paused']=pa
            
            layer_dic=simulate_network(params_msn_d1, params_msn_d2, params_stn,
                               synapse_models, sim_time=sim_time, seed=se,
                               I_e_add={}, threads=threads,
                               start_rec=start_rec,
                               model_params=model_params) 
        
            signal=layer_dic['SNR'].signals['spikes']
            spk_mean=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
            mr.append(numpy.mean(spk_mean[interval[0]:interval[1]]))
                
            signal=layer_dic['GPE'].signals['spikes']
            spk_mean=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
            mr.append(numpy.mean(spk_mean[interval[0]:interval[1]]))
        
        
            signal=layer_dic['STN'].signals['spikes']
            spk_mean=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
            mr.append(numpy.mean(spk_mean[interval[0]:interval[1]]))
            
        mr=numpy.array(mr)
        misc.pickle_save(mr, save_at)     
    else:
        mr=misc.pickle_load(save_at)  
    
    mr_SNR=mr[0::3]
    mr_GPE=mr[1::3]
    mr_STN=mr[2::3]
    
    d=numpy.array([freq_filter, mr_SNR, mr_GPE ,mr_STN])

    return d

params_msn_d1={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, 20., 0.1],
            'mod_times':[1,1000, 1000+500],  'n_mod':0}    
params_msn_d2={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, 0.1, 0.1],
            'mod_times':[1.,1000., 1000.+500.],  'n_mod':1} 
params_stn={'rate':219., 'mod':False,'mod_rate':0., 'mod_times':[1000., 1000.+500.]} 


N_MSN=1500.0

save_result_at=OUTPUT_PATH+'/simulate_network_indirect_filtering_mean_rate.plk'
synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p']
sim_time=4000.
start_rec=1000.0
interval=[0,3000] #  since recording starts at 3000 ms
d=[]
for syn_STN in ['STN_SNR_ampa_s', 'STN_SNR_ampa_p3']:
    
    for syn_GPE in ['GPE_SNR_gaba_p', 'GPE_SNR_gaba_s_ref']:
        synapse_models=['MSN_SNR_gaba_p1',syn_GPE]
        save_result_at=OUTPUT_PATH+'/simulate_filtering'+str(N_MSN)+syn_STN+syn_GPE+'.plk'
        d.append(simulate_filtering(1, save_result_at, interval, syn_STN, params_msn_d1, 
                                          params_msn_d2, params_stn,
                                          synapse_models, sim_time=sim_time, 
                                          threads=1, 
                                          start_rec=start_rec))
    


# DISPLAY
plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0, fontsize=16)
font_size_text = 14
fig = pylab.figure(facecolor = 'w' )

ax_list = []
ax_list.append( MyAxes(fig, [ .075, .37, .135, .26 ] ) )    # text box
ax_list.append( MyAxes(fig, [ .26,  .6,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .53,  .6,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .6,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .26,  .1,  .165, .34 ] ) )    #   
ax_list.append( MyAxes(fig,  [ .53,  .1,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .1,  .165, .34 ] ) )    # 


ax=ax_list[1]
i=0
MSNmeanRates=numpy.array([d[0][i],d[1][i], d[2][i], d[3][i]])
i=1
SNRmeanRates=numpy.array([d[0][i],d[1][i], d[2][i], d[3][i]])
plot_selection_vs_neurons(ax, MSNmeanRates, SNRmeanRates)

pylab.show()

# dpi does not matter since svg and pdf are both vectorbased
fig.savefig( picture_dir + '/' + SNAME  + '.svg', format = 'svg') 
fig.savefig( picture_dir + '/' + SNAME  + '.pdf', format = 'pdf')