import numpy
import pylab
import os
import sys
import time as ttime
import pprint
import nest
import random
import copy
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
from simulation_utils import simulate_network_poisson


model_name=os.getcwd().split('/')[-2]
picture_dir='/'.join(os.getcwd().split('/')[0:-3]) + '/pictures/'+model_name 
SNAME  = sys.argv[0].split('/')[-1].split('.')[0]

OUTPUT_PATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]
SELECTION_THR=5.
INTERMEDIATE_THR=10.0

def plot_filtering(ax,data):
    
    colors=['b','g','r']
    linestyles=['--','-']
    for i in range(len(data)):
        c=colors[i]
        for j in range(len(data[i])):
            y=data[i][j][0]
            x=data[i][j][2] #numpy.arange(1,len(y)+1)
            ax.plot(x, y, **{'color':c,'linestyle':linestyles[j]} )
   
    line1=ax.plot(1,1,'-k')
    line2=ax.plot(1,1,'--k')
    
    pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend
    leg=ax.legend([line1, line2],['Static \nSTN-SNr', 'Depressing \nSTN-SNr'], loc='best')
    frame  = leg.get_frame() 
    frame.set_visible(False) 
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    pylab.setp(frame, edgecolor='w') 
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    pylab.setp(ltext, fontsize=pylab.rcParams['text.fontsize']-8, backgroundcolor='w')
    
    colors=['b','g','r']   
    labels=[r'First 100 ms' , r'250-350 ms',  
            r'Last 100 ms']
    coords=[[ 0.55, 0.15], [0.07, 0.47], [0.15, 0.78]]    
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=font_size_text, 
                 **{'color': color})
        
    ax.set_ylabel('Rate SNr (spikes/s)') 
    ax.set_xlabel('Rate STN (spikes/s)')
    ax.my_set_no_ticks( yticks=8, xticks = 6 )
    
    ax.set_xlim(misc.adjust_limit([15, 55]))
    #ax.set_ylim(misc.adjust_limit(ylim))
    
def plot_short(ax, data, title):
    mr_new, std_new=data
    start=90
    for i in range(3):
        y=mr_new[i][50:120]
        x=numpy.arange(1, len(y)+1)
        ax.plot(x, y)
        #ax.plot(mr_new[i]-std_new[0],'--k')
        #ax.plot(mr_new[i]+std_new[0],'--k')
    ax.set_title(title, fontsize=font_size_text)
    #ax.set_xlim([90,150])
    pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend
    ax.set_ylabel('Rate SNr (spikes/s)') 
    ax.set_xlabel('Time (ms)')
    
    colors=['b','g','r']
    labels=[r'SNr' , r'GPe',  
            r'STN']
    coords=[[ 0.55, 0.15], [0.07, 0.47], [0.15, 0.78]]    
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=font_size_text, 
                 **{'color': color})
    ax.my_set_no_ticks( yticks=8, xticks = 6 )
    
def plot_short2(ax, data, data2, title):
    mr_new, std_new=data
    start=90
    for i in range(3):
        y=mr_new[i][50:120]
        x=numpy.arange(1, len(y)+1)
        ax.plot(x, y)
        #ax.plot(mr_new[i]-std_new[0],'--k')
        #ax.plot(mr_new[i]+std_new[0],'--k')
        
    mr_new, std_new=data2
    start=90

    y=mr_new[0][50:120]
    x=numpy.arange(1, len(y)+1)
    ax.plot(x, y,'--b')
        
    ax.set_title(title, fontsize=font_size_text)
    #ax.set_xlim([90,150])
    pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend
    ax.set_ylabel('Rate SNr (spikes/s)') 
    ax.set_xlabel('Time (ms)')
    
    colors=['b','g','r']
    labels=[r'SNr' , r'GPe',  
            r'STN']
    coords=[[ 0.55, 0.15], [0.07, 0.47], [0.15, 0.78]]    
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=font_size_text, 
                 **{'color': color})
    ax.my_set_no_ticks( yticks=8, xticks = 6 )
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

def simulate_hyperdirect_fun(params_msn_d1, params_msn_d2, params_stn,
                             synapse_models, I_e_add,sim_time, seed,
                             threads, start_rec, model_params, flag_bg, stn_syn='STN_SNR_gaba_s'):
    

    if flag_bg:
        N_MSN=model_params['neurons']['MSN_D1']['n']
        
        params_msn_d1.update({'bg_rate':[0.1*(500)]})
        params_msn_d2.update({'bg_rate':[0.1*(500)]})
        
        # Change paramters to use poisson background
        params={'conns':{'MSN_D1_SNR':{'p':0},
                         'MSN_D2_GPE':{'p':0,
                                       'lines':False}},
                      'neurons':{'MSN_D1':{'n':0},
                                 'MSN_D2':{'n':0},
                                 'MSN_D1_bg':{'n':300,
                                              'lesion':False},
                                 'MSN_D2_bg':{'n':300,
                                              'lesion':False},  
                                 'GPE': {'paused':False}}}
        model_params=misc.dict_merge(model_params, params)


    
    
    layer_dic=simulate_network_poisson(params_msn_d1, params_msn_d2, params_stn,
                       synapse_models, sim_time=sim_time, seed=seed,
                       I_e_add=I_e_add, threads=4,start_rec=start_rec,
                       model_params=model_params) 

    r=[]
    signal=layer_dic['SNR'].signals['spikes']
    r.append(numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0))

    signal=layer_dic['GPE'].signals['spikes']
    r.append(numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0))
    
    signal=layer_dic['STN'].signals['spikes']
    r.append(numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0))
    
    numpy.array(r)    
    return r 

def simulate_hyperdirect(load, N_MSN, save_at, type, burst_time, burst_rate, n_exp, 
                         res, threads, stn_syn='STN_SNR_ampa_s', 
                         gpe_syn='GPE_SNR_gaba_p', flag_bg=False):
            
    ra=random.random()*200.
    start_rec=950+ra #Hundred ms before burst
    delay=10.0
    params_msn_d1={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [],
                'mod_times':[],  'n_mod':0., 'bg_rate':0}    
    params_msn_d2={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [],
                'mod_times':[],  'n_mod':0, 'focus':False,'skip':1,
                'bg_rate':0} 
    params_stn={'rate':219., 'mod':True,'mod_rate':burst_rate, 'mod_times':[1000.+ra, 1000.+ra+burst_time]} 
    synapse_models=['MSN_SNR_gaba_p1', gpe_syn]    
    
    model_params= {'misc':{'N_MSN':N_MSN},
                   'conns':{ 'MSN_D2_GPE':{'lines':False},
                            'STN_SNR':{'syn':stn_syn}},
                   'neurons':{'MSN_D1':{'n':0},
                              'MSN_D2':{'n':0}}} 
    if type =='normal':
       pass 
    if type =='no STN-GPE':
        model_params['conns']= {'STN_GPE':{'lesion':True}}
    if type =='no GPE-STN':
        model_params['conns']= {'GPE_STN':{'lesion':True}}    
    if type=='no STN-SNR':
        model_params['conns']= {'STN_SNR':{'lesion':True}}
    if type=='no GPE-SNR':
        model_params['conns']= {'GPE_SNR':{'lesion':True}}

    
    if flag_bg: save_at=save_at+'_bg'
    
    sim_time=1300.+burst_time
       
    raw_r=[]
    conv_r=[]
    mean_conv_r=[]
    std_conv_r=[]
    inputs=[]
    i=0
    if not load:      
        mr=[]   
        for i in range(n_exp):  
            seed=i
            mr_tmp=simulate_hyperdirect_fun(params_msn_d1,  params_msn_d2, params_stn,
                                        synapse_models, {},sim_time, seed,threads, 
                                        start_rec, model_params, flag_bg)
            
            mr.append(mr_tmp)
        
        misc.pickle_save(mr, save_at)
    else: 
       mr=misc.pickle_load(save_at)
    
    
    r_new=[[],[],[]]
    for r in mr:
        r=numpy.array(r)
        for i in range(3):
            r_new[i].append(r[i::4,:]) #
                     
    
    mr_new=[]
    std_new=[]
    for r in r_new:
        r=numpy.array(r)
        r=misc.convolve(r[:,0,:], res, 'triangle',single=False)
        mr_new.append(numpy.mean(r,axis=0))
        std_new.append(numpy.std(r,axis=0))
    #print mr_new 
    return mr_new, std_new

def conv_data(raw_r, conv_r, mean_conv_r, std_conv_r, bin=100, kernel='rectangle'):
        conv_r.append(misc.convolve(raw_r[-1], bin, kernel,single=False))       
        mean_conv_r.append(numpy.mean(conv_r[-1],axis=0))
        std_conv_r.append(numpy.std(conv_r[-1],axis=0))
        return conv_r, mean_conv_r, std_conv_r

def simulate_filtering(load,intervals,  N_MSN, save_at, type, 
                       n_exp, res, threads, flag_bg=False):

    burst_time=500.0
    burst_rates=numpy.linspace(100,2000,5)
    d=[]
    for syn_STN in ['STN_SNR_ampa_s', 'STN_SNR_ampa_p3']:
        for syn_GPE in ['GPE_SNR_gaba_p']:
            d_SNR, d_GPE,d_STN=[], [], []
            synapse_models=['MSN_SNR_gaba_p1',syn_GPE]
            
            d_tmp=[]
            for br in burst_rates:
                save_at=OUTPUT_PATH+'/simulate_filtering'+str(N_MSN)+syn_STN+syn_GPE+str(br)+'.plk'
                d_tmp=simulate_hyperdirect(load, N_MSN, save_at, type, burst_time, br, n_exp, 
                         res, threads, syn_STN, flag_bg=False)

                # Do not use std
                d_SNR.append(d_tmp[0][0])
                print d_tmp[0][0]
                d_GPE.append(d_tmp[0][1])
                d_STN.append(d_tmp[0][2])
            
            
            d.append([d_SNR, d_GPE,d_STN])
    
    data=[]
    for inter in intervals:
        data.append(mean_rates(d, inter))
    
    #ax=pylab.subplot(111)
    return data

def mean_rates(r, interval):    
    
    start_burst=50.0
    r=copy.deepcopy(r)
    for i_syn in range(len(r)):
        for i_neuron in range(len(r[i_syn])):
            for i_rate  in range(len(r[i_syn][i_neuron])):
#                    pylab.plot(r[i_syn][2][i_rate])
#                    
#                    if i_neuron==2:                        
#                        pylab.show()    
                    r_tmp=r[i_syn][i_neuron][i_rate][ start_burst+interval[0]:start_burst+interval[1]]
                    print r_tmp
                    r[i_syn][i_neuron][i_rate]=numpy.mean(r_tmp, axis=0)
                    print r[i_syn][i_neuron][i_rate], i_rate

    return r     
                                


N_MSN=15000
burst_time, burst_rate=3, 7500.0
n_exp, threads=20, 4
res=1
type= 'normal'
#stn_syn, gpe_syn='STN_SNR_ampa_s', 'GPE_SNR_gaba_p'
#stn_syn, gpe_syn='STN_SNR_ampa_s', 'GPE_SNR_gaba_p'
data=[]
stn_syn, gpe_syn='STN_SNR_ampa_p3', 'GPE_SNR_gaba_p'
save_at=OUTPUT_PATH+'/simulate_network_hyperdirect'+str(N_MSN)+str(n_exp)+type+stn_syn+gpe_syn+'.plk'
data.append(simulate_hyperdirect(0, N_MSN, save_at,type, burst_time, burst_rate, n_exp, 
                         res, threads, stn_syn, gpe_syn, flag_bg=True))

type='no STN-SNR'
stn_syn, gpe_syn='STN_SNR_ampa_p3', 'GPE_SNR_gaba_p'
save_at=OUTPUT_PATH+'/simulate_network_hyperdirect'+str(N_MSN)+str(n_exp)+type+stn_syn+gpe_syn+'.plk'
data.append(simulate_hyperdirect(0, N_MSN, save_at,type, burst_time, burst_rate, n_exp, 
                         res, threads, flag_bg=True))

type='no STN-GPE'
stn_syn, gpe_syn='STN_SNR_ampa_p3', 'GPE_SNR_gaba_p'
save_at=OUTPUT_PATH+'/simulate_network_hyperdirect'+str(N_MSN)+str(n_exp)+type+stn_syn+gpe_syn+'.plk'
data.append(simulate_hyperdirect(0, N_MSN, save_at,type, burst_time, burst_rate, n_exp, 
                         res, threads, flag_bg=True))

N_MSN=15000
n_exp, threads=5, 4
res=100
intervals=[[0,100],[250,350],[400,500]]
type= 'normal'
d=simulate_filtering(1, intervals, N_MSN, save_at, type,
                      n_exp, res, threads, flag_bg=True)


# DISPLAY
 #Inspect results
plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0+320.0, fontsize=16)
font_size_text = 16
fig = pylab.figure( facecolor = 'w' )
ax_list = []
ax_list.append( MyAxes(fig, [ .075, .37, .135, .26 ] ) )    # text box
ax_list.append( MyAxes(fig, [ .26,  .7,  .165, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .53,  .7,  .165, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .7,  .165, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .26,  .4,  .165, .2 ] ) )    #     
ax_list.append( MyAxes(fig, [ .53,  .4,  .165, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .4,  .165, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .26,  .1,  .165, .2 ] ) )    #     
#ax_list.append( MyAxes(fig, [ .53,  .1,  .165, .2 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .8,   .1,  .165, .2 ] ) )    # 
ax=ax_list[1]
plot_short(ax, data[0],'Control')


ax=ax_list[2]
plot_short(ax, data[1],  'STN to SNr lesion')

ax=ax_list[3]
plot_short(ax, data[2],'STN to GPe lesion')

ax=ax_list[4]
plot_filtering(ax,d)

pylab.show()


# dpi does not matter since svg and pdf are both vectorbased
fig.savefig( picture_dir + '/' + SNAME  + '.svg', format = 'svg') 
fig.savefig( picture_dir + '/' + SNAME  + '.pdf', format = 'pdf')



