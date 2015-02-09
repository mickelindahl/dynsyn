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

LIM_SYN_EVENTS=900.

def plot_filtering(ax, MSNmeanRates, SNRmeanRates):
    colors=['b','g','m']   
    linestyle=['--','--','-'] 
    labels=[r'$ref_{init}^{MSN_{D1}}$' , r'$ref_{max}^{MSN_{D1}}$',  
            r'$fac^{MSN_{D1}}$']
    coords=[[0.4, 0.42], [ 0.02, 0.001], [0.15, 0.78]]
    
    
    for i, color in enumerate(colors):
        ax.plot(MSNmeanRates,SNRmeanRates[i,:],**{'color':color, 'linestyle':linestyle[i]})  

   
    ax.plot([0,2.5],[SELECTION_THR]*len([0,2]),**{'color':'k',
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
    ax.text( 0.8, 0.20,'Thr' , fontsize=pylab.rcParams['font.size']-2, transform=ax.transAxes, **{ 'color' : 'k'})
    
    ax.set_ylabel('Rate SNr (spikes/s)') 
    ax.set_xlabel('Rate $MSN_{D1}$ (spikes/s)')
    ax.my_set_no_ticks( yticks=8, xticks = 6 )
    
    ax.set_xlim(misc.adjust_limit([0, 2.5]))
    ax.set_ylim(misc.adjust_limit([-3,34]))
    
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,2.6])   
                                  
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['font.size'], 
                 **{'color': color})
    pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend

def plot_filtering_vs_burst(ax, syn_events, mr):
    
    colors=['k','k']   
    linesyles=['-','-.']
    labels=[ r'$fac^{MSN_{D1}}$']
    coords=[[0.04, 0.17]]   
    
    pl_lines=[]
    for ls, color, x, y in zip(linesyles,colors,syn_events, mr):
        pl_lines.append(ax.plot(x, y,**{'color':color, 'linestyle':ls})) 

    
    
    vec=[50,LIM_SYN_EVENTS]
    ax.plot(vec,[SELECTION_THR]*len(vec),**{'color':'k', 'label':'', 
                                                        'linestyle':'--'})
    
    
    
    ax.text( 0.8, 0.18,'Thr' , fontsize=pylab.rcParams['font.size']-2,
             transform=ax.transAxes, **{ 'color' : 'k'})
    ax.set_ylabel('Rate SNr (spikes/s)') 
    ax.set_xlabel('Synaptic input SNr (events/s)')
    ax.my_set_no_ticks( yticks=8, xticks = 5 ) 
    
    x_arrow=400./LIM_SYN_EVENTS
    y_arrow=2./50.
    
    ax.arrow(x_arrow, y_arrow, 0, 0.02, transform=ax.transAxes,
             width=0.01, head_width=0.03, head_starts_at_zero=True,
            head_length=0.02,**{ 'color':'k'})
    
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['font.size'], 
                 **{'color': color})
    
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,LIM_SYN_EVENTS]) 
    
    pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend
    leg=ax.legend(lines,['Increased size \nof bursting \nsubpopulation', 
                            'Increased activity \nin all $MSN_{D1}$'], 
                  loc='best')
    frame  = leg.get_frame() 
    frame.set_visible(False) 
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    pylab.setp(frame, edgecolor='w') 
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    pylab.setp(ltext, fontsize=14, backgroundcolor='w')
    
    ax.set_xlim(misc.adjust_limit([50,LIM_SYN_EVENTS]))
    ax.set_ylim(misc.adjust_limit([0,40]))

def plot_const_syn_events(ax, nMSN_range, SNRmeanRates):

    colors=['b','g','m' ] 
    linestyle=['--','--','-'] 
    labels=[r'$ref_{init}^{MSN_{D1}}$' , r'$ref_{max}^{MSN_{D1}}$', 
            r'$fac^{MSN_{D1}}$']
    coords=[[0.1, 0.68], [ 0.05, 0.12], [ 0.35, 0.45]]
    
    for i, color in enumerate(colors):
        ax.plot(nMSN_range*100.0,SNRmeanRates[i,:],**{'color':color, 'linestyle':linestyle[i]})  
    
    ax.set_ylabel('Rate SNr (spikes/s)') 
    ax.set_xlabel('Percent bursting $MSN_{D1}$ (%)')
    ax.my_set_no_ticks( yticks=7, xticks=6 )
    
    
    #lines = ax.lines
    #for line in lines:
    #    misc.slice_line(line, xlim=[0,N_MAX_BURSTING])   
    
    ax.plot([0,max(nMSN_range*100.0)],[SELECTION_THR]*len([0,2]),**{'color':'k',
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
        
    ax.text( 0.1, 0.33,'Thr' , transform=ax.transAxes, 
             fontsize=pylab.rcParams['font.size']-2,
             **{ 'color' : 'k'})
    
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['font.size'], 
                 **{'color': color})
    pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend
    ax.set_xlim(misc.adjust_limit([0, 4.2]))      
    ax.set_ylim(misc.adjust_limit([-3,23]))
      
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

def simulate_filtering_fun(freq, params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, sim_time, seed, I_e_add, threads, 
                           start_rec, model_params):
    
    
    
    mean=[]
    for f, se in zip(freq, seed):
        
        params_msn_d1.update({ 'base_rates':[f]})
        
        params_msn_d2.update({'focus':False,
                               'skip':1})
        
        layer_dic=simulate_network(params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, sim_time=sim_time, seed=se,
                           I_e_add=I_e_add, threads=threads,
                           start_rec=start_rec,
                           model_params=model_params) 
    
        layer_dic['SNR'].get_signal( 's', start=start_rec, stop=sim_time )
        signal=layer_dic['SNR'].signals['spikes']
        spk_mean=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
        mean.append(numpy.mean(spk_mean))
        
    
    return mean   

def simulate_filtering(load, save_at, N_MSN, params_msn_d1, params_msn_d2, models_msn, sim_time, start_rec):
    
    
    freq_filter=numpy.linspace(0.1,2.6,10)
    seed=range(len(freq_filter))
    
    model_params={'misc':{'N_MSN':N_MSN},
                  'neurons':{'MSN_D1':{'n':N_MSN},
                             'MSN_D2':{'n':N_MSN}}}
    
    if not load:
        mr=[]
        
        
        for syn in models_msn:
           synapse_models=[syn, 'GPE_SNR_gaba_p']
           mr.append(simulate_filtering_fun(freq_filter, params_msn_d1, params_msn_d2, 
                                            params_stn,  synapse_models, sim_time, seed,
                                                    {}, threads, start_rec, model_params))
        
        mr=numpy.array(mr)
        misc.pickle_save(mr, save_at)
    
    else: 
        mr=misc.pickle_load(save_at) 
    
    syn_ev=freq_filter*N_MSN*500./N_MSN
    
    # Row one in mr is s min, then s max and finally the plastic synapse
    data=numpy.array([syn_ev, freq_filter,mr[0,:],mr[1,:],mr[2,:]]) 
    return data 

def simulate_filterting_burst_fun(mod, freq, base_rates, params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, sim_time=2000., seed=1,
                           I_e_add={}, threads=4, start_rec=0,
                           model_params={}):
    
 
    mr=[]   
    for f, m, se, br in zip(freq, mod, seed, base_rates):
        
        params_msn_d1.update({ 'mod_rates': [0.1, f, 0.1],
                               'mod_times':[1,1000, 1000+500],  
                               'n_mod':int(m)})
        params_msn_d1.update({'base_rates':[0.1, br, 0.1],'base_times':[1,1000, 1000+500]})
        
        params_msn_d2.update({'focus':False,
                               'skip':1})
        
        layer_dic=simulate_network(params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, sim_time=sim_time, seed=se,
                           I_e_add=I_e_add, threads=4,start_rec=start_rec,
                           model_params=model_params) 
    
        layer_dic['SNR'].get_signal( 's', start=start_rec, stop=sim_time )
        signal=layer_dic['SNR'].signals['spikes']

        mr.append(numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0))
         
    return numpy.array(mr)   

def simulate_filterting_burst(load, save_at,interval, N_MSN, params_msn_d1, 
                              params_msn_d2, models_msn, sim_time, start_re):
    
    N_MSN_syn_on_SNR=500.
    max_base_rate=0.9
    max_syn_events_burst=1000  # Max for burst, then need to add on the contribution of
                        # background, See at return statement.
    msn_burst_rate=20.
    msn_base_rate=0.1

    n_max_bursting=(max_syn_events_burst)/msn_burst_rate-1
    
    n_burst_per_SNR=numpy.arange(1,n_max_bursting,2)
    prop_mod=n_burst_per_SNR/N_MSN_syn_on_SNR
    mod_const_syn=prop_mod*N_MSN
    
    
    mod=prop_mod*float(N_MSN)
    mod=numpy.array([int(m) for m in mod])
    
    freq=numpy.ones(len(mod))*msn_burst_rate
    syn_events_burst=numpy.array(mod*freq+(N_MSN-mod)*0.1)*500./N_MSN
    
    model_params={'misc':{'N_MSN':N_MSN},
                  'conns':{ 'MSN_D2_GPE':{ 'lines':False}},
                  'neurons':{'MSN_D1':{'n':N_MSN},
                         'MSN_D2':{'n':N_MSN}}}
    seed=range(len(mod))
    base_rates=numpy.ones(len(mod))*0.1
    
    if not load:
        mr=simulate_filterting_burst_fun(mod, freq, base_rates, params_msn_d1, 
                                             params_msn_d2, params_stn, synapse_models, 
                                             sim_time, seed, {}, threads, 
                                             start_rec, model_params)
        misc.pickle_save(mr, save_at)
    else: 
        mr=misc.pickle_load(save_at)
    
    mrb=numpy.mean(mr[:,interval[0]:interval[1]], axis=1)
    syn_ev_b=numpy.array(mod*freq+(N_MSN-mod)*0.1)*500./N_MSN
    data=numpy.array([syn_ev_b, mrb])
    
    return data

def simulate_const_syn_event(load, save_at,interval, N_MSN, params_msn_d1, 
                                        params_msn_d2, models_msn, sim_time, 
                                        start_rec):
    # Solve (500-n)*x + 20*n=600, where 500 is total number of MSNs, 20 is burst
    # activation, x is MSN mean rate and n is number of bursters. 
    # Then x=(600-20*n)/(500-n)
    
    N_MSN_syn_on_SNR=500.
    ratio_neurons_syns=N_MSN/N_MSN_syn_on_SNR
    
    # Max based rate, that is with no bursters, determines max number of 
    # synaptic events.
    max_base_rate=0.9
    n_syn_events=max_base_rate*N_MSN_syn_on_SNR  
    msn_burst_rate=20.
    
    n_max_bursting=n_syn_events/msn_burst_rate-1
    
    n_burst_per_SNR=numpy.arange(1,n_max_bursting,2)
    prop_mod=n_burst_per_SNR/N_MSN_syn_on_SNR
    mod_const_syn=prop_mod*N_MSN
    
    # Base rates, max is "n_syn_events/N_MSN_syn_on_SNR", we then remove
    # the rate contribution buy mod "msn_burst_rate*prop_mod"
    base_rates=(n_syn_events- n_burst_per_SNR*msn_burst_rate)/N_MSN_syn_on_SNR
    
    
    synapse_models_msn=['MSN_SNR_gaba_s_min','MSN_SNR_gaba_s_max','MSN_SNR_gaba_p1']
    freq_const_syn=numpy.ones(len(mod_const_syn))*20.0
    model_params={'misc':{'N_MSN':N_MSN},
                  'conns':{ 'MSN_D2_GPE':{ 'lines':False}},
                  'neurons':{'MSN_D1':{'n':N_MSN},
                             'MSN_D2':{'n':N_MSN}}}
    
    seed=range(len(mod_const_syn))
    if not load:
        r=[]
        
        
        for syn in synapse_models_msn:
           synapse_models=[syn, 'GPE_SNR_gaba_p']

           
           r.append(simulate_filterting_burst_fun(mod_const_syn, freq_const_syn, 
                                                           base_rates, params_msn_d1, 
                                             params_msn_d2, params_stn, synapse_models, 
                                             sim_time, seed, {}, threads, 
                                             start_rec, model_params))
        
        r=numpy.array(r)
        m1=numpy.mean(r[0,:,0:500], axis=0)
        m2=numpy.mean(r[0,:,500:1000], axis=0)
        misc.pickle_save(r, save_at)
    
    else: 
        r=misc.pickle_load(save_at) 
    
    mr=[]    

    for i, syn in enumerate(synapse_models_msn):   
        mr.append(numpy.mean(numpy.array(r[i,:,interval[0]:interval[1]]), axis=1))
    
    
    mr=numpy.array(mr)
    
    data=numpy.array([prop_mod ,mr[0,:],mr[1,:],mr[2,:]]) 
    
    return data
        
params_msn_d1={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [],
            'mod_times':[],  'n_mod':0}    
params_msn_d2={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [],
            'mod_times':[],  'n_mod':0} 
params_stn={'rate':219., 'mod':False,'mod_rate':0., 'mod_times':[]} 
models_msn=['MSN_SNR_gaba_s_min','MSN_SNR_gaba_s_max','MSN_SNR_gaba_p1']


N_MSN, sim_time, threads, start_rec=15000, 4000., 4, 500.
save_at=OUTPUT_PATH+'/simulate_filtering'+str(N_MSN)+'.plk'
data=simulate_filtering(1, save_at, N_MSN, params_msn_d1, params_msn_d2, models_msn,
                        sim_time, start_rec)


N_MSN, sim_time, threads, start_rec=15000, 2000., 4, 1000.
interval=[100,500]
synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p']
save_at=OUTPUT_PATH+'/simulate_filtering_burst'+str(N_MSN)+'.plk'
data2=simulate_filterting_burst(1, save_at,interval, N_MSN, params_msn_d1, 
                                        params_msn_d2, models_msn, sim_time, 
                                        start_rec)


N_MSN, sim_time, threads, start_rec=15000, 2000., 4, 1000.
interval=[250., 500.] # To have threshold crossing as wiht direct example
save_at=OUTPUT_PATH+'/simulate_const_syn_event'+str(N_MSN)+'.plk'
data3=simulate_const_syn_event(0, save_at,interval, N_MSN, params_msn_d1, 
                                        params_msn_d2, models_msn, sim_time, 
                                        start_rec)
   

# DISPLAY
plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0, fontsize=16)
font_size_text = 12
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
plot_filtering(ax, data[1,:], data[2:5,:])

ax=ax_list[2]
plot_filtering_vs_burst(ax, [data[0,:],data2[0,:]] , [data[4,:],data2[1,:]])

ax=ax_list[3]
plot_const_syn_events(ax, data3[0,:] , data3[1:4,:])


ax=ax_list[0]
#plot_text(ax, info_string=s)
pylab.show()

# dpi does not matter since svg and pdf are both vectorbased
fig.savefig( picture_dir + '/' + SNAME  + '.svg', format = 'svg') 
fig.savefig( picture_dir + '/' + SNAME  + '.pdf', format = 'pdf')