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

model_name=os.getcwd().split('/')[-2]
picture_dir='/'.join(os.getcwd().split('/')[0:-3]) + '/pictures/'+model_name 
SNAME  = sys.argv[0].split('/')[-1].split('.')[0]

OUTPUT_PATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]
SELECTION_THR=5.

def plot_filtering(ax, STNmeanRates, SNRmeanRates, linestyle='-'):


  
    colors=['b', 'b' ,'g', 'g' ] 
    linestyles=['-','--','-','--']
    
    
    for i, color in enumerate(colors):
        ax.plot(STNmeanRates[i,:],SNRmeanRates[i,:],**{'color':color,'linestyle':linestyles[i]})  
     
    ax.set_ylabel('Rate SNr (spikes/s)') 
    ax.set_xlabel('Rate STN (spikes/s)')
    ax.my_set_no_ticks( yticks=8, xticks = 6 )
    

    
    line1=ax.plot(1,1,'-k')
    line2=ax.plot(1,1,'--k')
    pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend
    leg=ax.legend([line1, line2],['$dep^{GPe}$', '$ref_{30 Hz}^{GPe}$'], loc='best')
    frame  = leg.get_frame() 
    frame.set_visible(False) 
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    pylab.setp(frame, edgecolor='w') 
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    pylab.setp(ltext, fontsize=pylab.rcParams['text.fontsize']-2, backgroundcolor='w')
    
    ax.set_xlim(misc.adjust_limit([10, 50]))
    ax.set_ylim(misc.adjust_limit([0,140]))
    #lines = ax.lines
    #for line in lines:
    #    misc.slice_line(line, xlim=[10,25])   
                                  
    labels=[ 'Static \nSTN-SNr', 'Depressing \nSTN-SNr']
    coords=[[0.55, 0.3], [ 0.3, 0.6]]
    colors=['b','g'] 
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=font_size_text, 
                 **{'color': color})
             
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

def get_interval_rates(interval, mean_rates):
    interval_rates=[]
    for mr in mean_rates:
        mr=numpy.array(mr)
        interval_rates.append(numpy.mean(mr[:,interval[0]:interval[1]], axis=1))
    
    
    return numpy.array(interval_rates)

def simulate_filtering(load, save_at, interval, syn_STN, freq, params_msn_d1, 
                       params_msn_d2, params_stn, synapse_models, sim_time, 
                        threads, start_rec,N_MSN):
    
    
    
    mr=[]
    seed=range(len(freq))
    model_params={'misc':{'N_MSN':N_MSN},
                   'conns':{ 'MSN_D2_GPE':{ 'lines':False},
                               'STN_SNR' :{'syn':syn_STN}},
                   'neurons':{'MSN_D1':{'n':N_MSN},
                              'MSN_D2':{'n':N_MSN},
                              'GPE': {'paused':0}}}    
    if not load:
        for f, se in zip(freq, seed):
            
            params_stn.update({'rate':250., 'mod':True,'mod_rate':f, 'mod_times':[1., sim_time]})
            
            layer_dic=simulate_network(params_msn_d1, params_msn_d2, params_stn,
                               synapse_models, sim_time, seed,{}, threads,
                               start_rec, model_params) 
        
            signal=layer_dic['SNR'].signals['spikes']
            r=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
            mr.append(numpy.mean(r[interval[0]:interval[1]]))
    
            signal=layer_dic['GPE'].signals['spikes']
            r=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
            mr.append(numpy.mean(r[interval[0]:interval[1]]))
            
            signal=layer_dic['STN'].signals['spikes']
            r=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
            mr.append(numpy.mean(r[interval[0]:interval[1]]))
     
        mr=numpy.array(mr)
        misc.pickle_save(mr, save_at)     
    else:
        mr=misc.pickle_load(save_at)  
    
    mr_SNR=mr[0::3]
    mr_GPE=mr[1::3]
    mr_STN=mr[2::3]
    
    d=numpy.array([freq, mr_SNR, mr_GPE ,mr_STN])       
    return d
        
N_MSN=1500

params_msn_d1={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [],
            'mod_times':[],  'n_mod':0}    
params_msn_d2={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [],
            'mod_times':[],  'n_mod':0, 'focus':False,'skip':1, 'lines':False} 
params_stn={'rate':219., 'mod':False,'mod_rate':0., 'mod_times':[]} 


sim_time=5000.
freq_burst=numpy.linspace(1,1500,10)
seed=range(len(freq_burst))
mod=numpy.ones(len(freq_burst))*1000*N_MSN/15000.


save_result_at=OUTPUT_PATH+'/simulate_filtering.plk'

synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p']
sim_time=4000.
start_rec=1000.0
threads=4
interval=[0,3000]
d=[]
mr_SNR,mr_STN=[], []
for syn_STN in ['STN_SNR_ampa_s', 'STN_SNR_ampa_p3']:
    
    for syn_GPE in ['GPE_SNR_gaba_p', 'GPE_SNR_gaba_s_ref']:
        synapse_models=['MSN_SNR_gaba_p1',syn_GPE]
        save_result_at=OUTPUT_PATH+'/simulate_filtering'+str(N_MSN)+syn_STN+syn_GPE+'.plk'
        d=simulate_filtering(1, save_result_at, interval, syn_STN, freq_burst, 
                                    params_msn_d1, params_msn_d2,  params_stn,
                                    synapse_models, sim_time, 
                                    threads, start_rec, N_MSN)
        mr_STN.append(d[3])
        mr_SNR.append(d[1])
        
mr_SNR=numpy.array(mr_SNR)
mr_STN=numpy.array(mr_STN)
 
# DISPLAY
plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0, fontsize=16)
font_size_text = 14
fig = pylab.figure( facecolor = 'w' )

ax_list = []
ax_list.append( MyAxes(fig, [ .075, .37, .135, .26 ] ) )    # text box
ax_list.append( MyAxes(fig, [ .26,  .6,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .53,  .6,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .6,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .26,  .1,  .165, .34 ] ) )    #   
ax_list.append( MyAxes(fig,  [ .53,  .1,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .1,  .165, .34 ] ) )    # 


ax=ax_list[1]
plot_filtering(ax,mr_STN, mr_SNR)




ax=ax_list[0]
#plot_text(ax, info_string=s)
pylab.show()
# dpi does not matter since svg and pdf are both vectorbased
fig.savefig( picture_dir + '/' + SNAME  + '.svg', format = 'svg') 
fig.savefig( picture_dir + '/' + SNAME  + '.pdf', format = 'pdf')
