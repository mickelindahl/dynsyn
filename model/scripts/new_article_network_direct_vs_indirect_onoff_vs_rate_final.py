import numpy
import pylab
import os
import sys
import time as ttime
import pprint
import nest
import random
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
from simulation_utils import simulate_network_poisson, simulate_network_direct_onoff_vs_rate, simulate_network_direct_indirect_onoff_vs_rate, simulate_network_hyper_direct_onoff_vs_rate


model_name=os.getcwd().split('/')[-2]
picture_dir='/'.join(os.getcwd().split('/')[0:-3]) + '/pictures/'+model_name 
SNAME  = sys.argv[0].split('/')[-1].split('.')[0]

OUTPUT_PATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]
SELECTION_THR=5.
INTERMEDIATE_THR=10.0


def plot_selection(ax, data, legend=True):
    #! Colors for spike voltage traces fo AMPA, NMDA and GANA (color 0, 1 and 2 in
    #! ``color_list``).
    color_map = 'Set1'
    no_colors = 6
    color_list = misc.make_N_colors(cmap_name=color_map, N=no_colors)
    
    
    sel_dic={0:[],1:[],2:[]}
    
    for d in zip(*data):
        sel_dic[d[4]].append((d[0],d[1]))
    
    pprint.pprint(sel_dic)
    
    plot_prop=[(sel_dic[2],'k.','No selection'),
               (sel_dic[0],'ob','Selection'),
               (sel_dic[1],'ys','Intermediate')]
    
    #! Create figure where figsize(width,height) and figure dimenstions window 
    #! width = figsize(width) x dpi and window hight = figsize(hight) x dpi
    #fig=pylab.figure(figsize=(11, 10), dpi=50)
    
    #! [left, bottom, width, hight]
    #ax=pylab.axes([0.05, .25, .9, .7])
    iter=0
    for cords, color_marker, sel  in plot_prop:
        if cords !=[]:
            cords_zip=zip(*cords)
            ax.plot(cords_zip[0], cords_zip[1], color_marker,  markersize=20, label=sel, zorder=2)
            iter+=1
      
    # beautify
    #ax.axis([0.0, 48.0, 0, 48])
    #ax.set_xlim(misc.adjust_limit([15, 80]))
    #ax.set_ylim(misc.adjust_limit([15,50]))
    #ax.set_aspect('equal', 'box')
    #ax.set_axes([0.1, 0.5, 0.7,0.1])
    #ax.set_xticks((0, 4, 8, 12 ,16, 20, 24, 28,32,36,40))
    #ax.set_yticks((0, 4, 8, 12 ,16, 20, 24, 28,32,36,40))
    ax.grid(True)
    ax.set_ylabel('Burst rate (Spikes/s)')
    ax.set_xlabel('Bursting $MSN_{D1}$ (%)')
    if legend:
        ax.legend(numpoints=1, ncol=1, loc=[-1.3,0.5], borderpad=.75)


    ax.my_set_no_ticks( yticks=8, xticks = 4 )
        
def plot_selection_vs_neurons_direct2(ax, data, title, ylim=[0,34]):
    colors=['k','k','k','k']   
    labels=[r'1 % bursting' , r'2 % bursting',  
            r'4 % bursting']
    linestyles=['-','--','-.','.']
    coords=[[0.4, 0.42], [ 0.02, 0.001], [0.15, 0.78]]


    for i, d in enumerate(data):

            c='k'
            lb=r'%i %s bursting' % (d[0][0],'%')
            ls='--'
            
            x,y_mean,y_std=d.take([1,2,3],axis=0)
            ax.plot(x,y_mean,**{'color':c, 'label': lb, 'linestyle':ls})  
            
            
            ax.fill_between(x,y_mean-y_std, y_mean+y_std, facecolor=c, alpha=0.5)
            
        #x,y,dy=fit_pol(MSNmeanRates,SNRmeanRates[i,:])    
        #ax.plot(x,y,**{'color':color})  
    leg=ax.legend(loc='best')
    frame  = leg.get_frame() 
    frame.set_visible(False)
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    pylab.setp(ltext, fontsize=10.) 
    '''
    ax.plot([0,2.5],[SELECTION_THR]*len([0,2]),**{'color':[0.6,0.6,0.6],
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
    ax.text( 0.8, 0.20,'Thr' , transform=ax.transAxes, **{ 'color' : [0.6,0.6,0.6]})
    
    '''
    ax.set_ylabel('Rate SNr (spikes/s)') 
    ax.set_xlabel('Burst rate $MSN_{D1}$ (spikes/s)')
    ax.my_set_no_ticks( yticks=8, xticks = 6 )
    ax.set_title(title, fontsize=12)
    ax.set_xlim(misc.adjust_limit([15, 50]))
    ax.set_ylim(misc.adjust_limit(ylim))
    
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,50])   
    
def plot_selection_vs_neurons_indirect2(ax, data, title, ylim=[0,34]):
    colors=['k','k','k','k']   
    labels=[r'1 % bursting' , r'2 % bursting',  
            r'4 % bursting']
    linestyles=['-','--','-.','.']
    coords=[[0.4, 0.42], [ 0.02, 0.001], [0.15, 0.78]]


    for i, d in enumerate(data):
        
            c='k'
            lb=r'%i %s bursting' % (d[0][0],'%')
            ls='--'
            
            x,y_mean,y_std=d.take([1,2,3],axis=0)
            ax.plot(x,y_mean,**{'color':c, 'label': lb, 'linestyle':ls})  
            
            
            ax.fill_between(x,y_mean-y_std, y_mean+y_std, facecolor=c, alpha=0.5)
            
        #x,y,dy=fit_pol(MSNmeanRates,SNRmeanRates[i,:])    
        #ax.plot(x,y,**{'color':color})  
    leg=ax.legend(loc='best')
    frame  = leg.get_frame() 
    frame.set_visible(False)
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    pylab.setp(ltext, fontsize=10.) 
    '''
    ax.plot([0,2.5],[SELECTION_THR]*len([0,2]),**{'color':[0.6,0.6,0.6],
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
    ax.text( 0.8, 0.20,'Thr' , transform=ax.transAxes, **{ 'color' : [0.6,0.6,0.6]})
    
    '''
    ax.set_ylabel('Rate SNr (spikes/s)') 
    ax.set_xlabel('Burst rate $MSN_{D1}$ (spikes/s)')
    ax.my_set_no_ticks( yticks=8, xticks = 6 )
    ax.set_title(title, fontsize=12)
    ax.set_xlim(misc.adjust_limit([15, 50]))
    ax.set_ylim(misc.adjust_limit(ylim))
    
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,50])   
def plot_selection_vs_neurons_direct_linear_fit_and_spread(ax, mods, slope, ylim=[-0.65,0]):
    colors=['b','g','r']   
    labels=[r'First 100 ms' , r'250-350 ms',  
            r'Last 100 ms']
    coords=[[ 0.55, 0.15], [0.07, 0.47], [0.15, 0.78]]
    
    ps=[]
    for i, color in enumerate(colors):
        ax.plot(mods*100.0,slope[i,:],**{'color':color,  'linestyle':'-'})  
        #x,y,dy=fit_pol(MSNmeanRates,SNRmeanRates[i,:])    
        #ax.plot(x,y,**{'color':color})    
    '''
    ax.plot([0,2.5],[SELECTION_THR]*len([0,2]),**{'color':[0.6,0.6,0.6],
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
    
    ax.text( 0.8, 0.20,'Thr' , transform=ax.transAxes, **{ 'color' : [0.6,0.6,0.6]})
    '''
    ax.set_ylabel('$\Delta_{SNr}/ \Delta_{MSN_{D1}}$') 
    ax.set_xlabel('Bursting $MSN_{D1}$ (%)')
    ax.my_set_no_ticks( yticks=8, xticks = 6 )
    
    #ax.set_xlim(misc.adjust_limit([17, 50]))
    ax.set_ylim(misc.adjust_limit(ylim))
    
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,50])   
                                  
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=12, 
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


def get_data_and_slopes(inputs,mean_conv_r,std_conv_r,resolution, N_MSN):
    
    data, data_split, slopes={}, {}, []
    for p in points:
        x, y, z_mean, z_std, sel=[], [], [], [], []

        for i, m, s in zip(inputs,mean_conv_r,std_conv_r):
            x.append(i[0]/float(N_MSN)*100) # As percent
            y.append(i[1])
            z_mean.append(m[p])
            z_std.append(s[p])
            
            if m[p]<SELECTION_THR: sel.append(0) 
            elif m[p]<INTERMEDIATE_THR: sel.append(1)
            else: sel.append(2)
        
        data[p]=numpy.array([x,y,z_mean, z_std, sel])        
        data_split[p]=numpy.split(data[p],resolution, axis=1)
        
        sl=[]
        for ds in data_split[p]:
            x,y, dy, p=fit_pol(ds[1],ds[2]) 
            sl.append(p[0])
        
        slopes.append(sl)
        
    slopes=numpy.array(slopes)  
    return data, data_split, slopes

def simulate_direct(load, N_MSN, save_at, threads, flag_bg=False):
            
    params_msn_d1={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, 20, 0.1],
                'mod_times':[1,1000, 1000+500],  'n_mod':0, 'bg_rate':0}    
    params_msn_d2={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, 20, 0.1],
                'mod_times':[1,1000, 1000+500],  'n_mod':0, 'focus':False,'skip':1,
                'bg_rate':0} 
    params_stn={'rate':250., 'mod':False,'mod_rate':0., 'mod_times':[1000., 1000.+500.]} 
    
    model_params={'neurons':{'MSN_D1':{'n':N_MSN},
                             'MSN_D2':{'n':N_MSN}}}    
    
    synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p']    
    
    if flag_bg: save_at=save_at+'_bg'
    
    sim_time=2000.    
    resolution=10
    n_exp=5
    
    burst_rate   =numpy.linspace(15 ,50 ,resolution)
    proportions  =numpy.linspace(0.01, 0.11,resolution) #arange(1,7,1)*150. 
    mods=proportions*N_MSN
    
    raw_r, conv_r, mean_conv_r, std_conv_r, inputs=[], [], [], [], []
    i=0
    if not load:
        for m in mods:
            for r in burst_rate:
                tmp_rates=[]                
                for e in range(n_exp):
                    seed=i
                    rates =simulate_network_direct_onoff_vs_rate(m, r, params_msn_d1,  params_msn_d2, params_stn, synapse_models, 
                                                                 {'SNR':280, 'STN':0,'GPE':20}, sim_time=sim_time, seed=seed, 
                                                                 threads=threads, start_rec=500.,model_params=model_params, flag_bg=flag_bg)
                    tmp_rates.append(list(rates))
                    i+=1
                
                raw_r.append(numpy.array(tmp_rates))
                conv_r.append(misc.convolve(raw_r[-1], 100, 'rectangle',single=False))      
                mean_conv_r.append(numpy.mean(conv_r[-1],axis=0))
                std_conv_r.append(numpy.std(conv_r[-1],axis=0))
                inputs.append((m,r))
                                
                misc.pickle_save([raw_r, conv_r, mean_conv_r, std_conv_r, inputs], save_at)
    else: 
        raw_r, conv_r, mean_conv_r, std_conv_r, inputs=misc.pickle_load(save_at)
    
    
    data, data_split, slopes=get_data_and_slopes(inputs, mean_conv_r,std_conv_r,resolution, N_MSN)
                 
    return data, data_split, slopes, burst_rate, proportions

def simulate_indirect_fun(mod, freq, params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, I_e_add,sim_time=2000., seed=1,
                            threads=4, start_rec=0,
                           model_params={}, flag_bg=False,max_mod=0):
    
    
    
    params_msn_d2.update({ 'mod_rates': [0.1, freq, 0.1],
                           'mod_times':[1,1000, 1000+500],  
                           'n_mod':int(mod)})
    if flag_bg:
        N_MSN=model_params['neurons']['MSN_D1']['n']
        
        params_msn_d1.update({'bg_rate':[0.1*(500-500*mod/float(N_MSN))]})
        params_msn_d2.update({'bg_rate':[0.1*(500)]})
        
        # Change paramters to use poisson background
        model_params={'misc':{'N_MSN':N_MSN},
                      'conns':{'MSN_D1_SNR':{'p':500./float(N_MSN)},
                               'MSN_D2_GPE':{'p':500./float(N_MSN),
                                'lines':False}},
                      'neurons':{'MSN_D1':{'n':0},
                                 'MSN_D2':{'n':max_mod},
                                 'MSN_D1_bg':{'n':300,
                                              'lesion':False},
                                 'MSN_D2_bg':{'n':300,
                                              'lesion':False},  
                                 'GPE': {'paused':False}}}

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

def simulate_indirect(load, N_MSN, save_at, threads,resolution=10, 
                      max_prop=1.1, flag_bg=False, skip=1, lines=False):
            
    params_msn_d1={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [],
                   'mod_times':[],  'n_mod':0, 'bg_rate':0}    
    params_msn_d2={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [],
                   'mod_times':[],  'n_mod':0, 'focus':False,'skip':skip, 'bg_rate':0} 
    params_stn={'rate':250., 'mod':False,'mod_rate':0., 'mod_times':[1000., 1000.+500.]} 

    synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p']
    
    
        
    model_params={'misc':{'N_MSN':N_MSN},
                  'conns':{ 'MSN_D2_GPE':{ 'lines':lines}},
                  'neurons':{'MSN_D1':{'n':N_MSN},
                             'MSN_D2':{'n':N_MSN}}}     
       
    if flag_bg: save_at=save_at+'_bg'
    
    sim_time=2000.
    n_exp=5
    
    burst_rate=numpy.linspace(15 ,50 ,resolution)
    proportions  =numpy.linspace(0.01,max_prop,resolution) #arange(1,7,1)*150. 
    mods=proportions*N_MSN
        
    raw_r, conv_r, mean_conv_r, std_conv_r, inputs=[], [], [], [], []
    i=0
    if not load:
        for m in mods:
            
            for r in burst_rate:
                tmp_rates=[]                
                for e in range(n_exp):
                    seed=i   
                    rates_SNR, rates_GPE =simulate_indirect_fun(m, r,  params_msn_d1, params_msn_d2, params_stn, synapse_models, 
                                                                {'SNR':280, 'STN':0,'GPE':20}, sim_time=sim_time, seed=seed, 
                                                                threads=threads, start_rec=500.,model_params=model_params,flag_bg=flag_bg, max_mod=max(mods))
                    
                    tmp_rates.append(list(rates_SNR))
                    
                    i+=1

                raw_r.append(numpy.array(tmp_rates))
                conv_r.append(misc.convolve(raw_r[-1], 100, 'rectangle',single=False))       
                mean_conv_r.append(numpy.mean(conv_r[-1],axis=0))
                std_conv_r.append(numpy.std(conv_r[-1],axis=0))
                inputs.append((m,r))
                
                misc.pickle_save([raw_r, conv_r, mean_conv_r, std_conv_r, inputs], save_at)
    else: 
        raw_r, conv_r, mean_conv_r, std_conv_r, inputs=misc.pickle_load(save_at)

    
    data, data_split, slopes=get_data_and_slopes(inputs, mean_conv_r,std_conv_r,resolution, N_MSN)

    return data, data_split, slopes, burst_rate, proportions
    
def simulate_direct_vs_indirect(load, N_MSN, burst_rate, save_at, threads,resolution=10,  flag_bg=False):
            


    params_msn_d1={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, burst_rate, 0.1],
                'mod_times':[1,1000, 1000+500],  'n_mod':0, 'bg_rate':0}    
    params_msn_d2={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, burst_rate, 0.1],
                'mod_times':[1,1000, 1000+500],  'n_mod':0, 'focus':False,'skip':1, 'bg_rate':0} 
    params_stn={'rate':250., 'mod':False,'mod_rate':0., 'mod_times':[1000., 1000.+500.]} 
    
    
    model_params= {'conns':{ 'MSN_D2_GPE':{'lines':False}},
                   'neurons':{'MSN_D1':{'n':N_MSN},
                              'MSN_D2':{'n':N_MSN}}} 
    
    synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p']    
    
    if flag_bg: save_at=save_at+'_bg'
    
    sim_time=2000.
    

    n_exp=2

    proportions_d1  =numpy.linspace(0.01,0.15,resolution) #arange(1,7,1)*150. 
    proportions_d2  =numpy.linspace(0.01,0.15,resolution) #arange(1,7,1)*150. 
    mods_d1=proportions_d1*N_MSN
    mods_d2=proportions_d2*N_MSN    
    
    raw_r=[]
    conv_r=[]
    mean_conv_r=[]
    std_conv_r=[]
    

    inputs=[]
    i=0
    if not load:
        for m_d1 in mods_d1:
            for m_d2 in mods_d2:
                tmp_rates=[]                
                for e in range(n_exp):
                    seed=i
                    
                    
                    rates_SNR, rates_GPE =simulate_network_direct_indirect_onoff_vs_rate(m_d1, m_d1, 
                                                           params_msn_d1, 
                                                           params_msn_d2, 
                                                           params_stn,
                                synapse_models, {'SNR':280, 'STN':0,'GPE':20},
                                sim_time=sim_time, seed=seed, 
                               threads=threads, start_rec=500.,model_params=model_params, flag_bg=flag_bg)
                    
                    tmp_rates.append(list(rates_SNR))
                    
                    i+=1

                raw_r.append(numpy.array(tmp_rates))
                conv_r.append(misc.convolve(raw_r[-1], 100, 'rectangle',single=False))       
                mean_conv_r.append(numpy.mean(conv_r[-1],axis=0))
                std_conv_r.append(numpy.std(conv_r[-1],axis=0))
                inputs.append((m_d1,m_d2))
                
                misc.pickle_save([raw_r, conv_r, mean_conv_r, std_conv_r, inputs], save_at)
    else: 
        raw_r, conv_r, mean_conv_r, std_conv_r, inputs=misc.pickle_load(save_at)

    
    data, data_split, slopes=get_data_and_slopes(inputs, mean_conv_r,std_conv_r,resolution, N_MSN)

    return data, data_split, slopes, burst_rate, proportions_d1, proportions_d2

def simulate_hyperdirect(load, N_MSN, save_at, threads, flag_bg=False):
            
    ra=random.random()*200.
    start_rec=900+ra
    delay=10.0
    params_msn_d1={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, 5000.0, 0.1],
                'mod_times':[1,1000.+ra+delay, 1000.+ra+2.+delay],  'n_mod':0., 'bg_rate':0}    
    params_msn_d2={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, 0.1, 0.1],
                'mod_times':[1,1000, 1000+500],  'n_mod':0, 'focus':False,'skip':1,
                'bg_rate':0} 
    params_stn={'rate':300., 'mod':True,'mod_rate':0., 'mod_times':[1000.+ra, 1000.+ra+2.]} 
    
    
    model_params= {'conns':{ 'MSN_D2_GPE':{'lines':False},
                              'STN_GPE':{'lesion':True},
                             'GPE_STN':{'lesion':True}},
                   'neurons':{'MSN_D1':{'n':N_MSN},
                              'MSN_D2':{'n':N_MSN}}} 
    
    synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p']    
    
    if flag_bg: save_at=save_at+'_bg'
    
    sim_time=1300.
    
    resolution=5
    n_exp=5

    proportions_d1  =numpy.linspace(0.01,0.1,resolution) #arange(1,7,1)*150. 
    proportions_d2  =numpy.linspace(0.01,0.1,resolution) #arange(1,7,1)*150. 
    mods_d1=proportions_d1*N_MSN
    mods_d2=proportions_d2*N_MSN    
    
    raw_r=[]
    conv_r=[]
    mean_conv_r=[]
    std_conv_r=[]
    m_d1=1000
    m_d2=0
    inputs=[]
    i=0
    if not load:      
        tmp_rates=[]   
        for e in range(n_exp):
            seed=i
            
            
            rates_SNR, rates_GPE =simulate_network_direct_indirect_onoff_vs_rate(m_d1, m_d1, 
                                                   params_msn_d1, 
                                                   params_msn_d2, 
                                                   params_stn,
                        synapse_models, {'SNR':280, 'STN':0,'GPE':20},
                        sim_time=sim_time, seed=seed, 
                       threads=threads, start_rec=start_rec,model_params=model_params, flag_bg=flag_bg)
            
            tmp_rates.append(list(rates_SNR))
            raw_r.append(numpy.array(tmp_rates))

            i+=1

        raw_r.append(numpy.array(tmp_rates))

        inputs.append((m_d1,m_d2))
        
        misc.pickle_save([raw_r, conv_r, mean_conv_r, std_conv_r, inputs], save_at)
    else: 
        raw_r, conv_r, mean_conv_r, std_conv_r, inputs=misc.pickle_load(save_at)
    conv_r=[]
    mean_conv_r=[]
    std_conv_r=[]

    conv_r, mean_conv_r, std_conv_r=conv_data(raw_r,conv_r, mean_conv_r, std_conv_r,
                                                       bin=1, kernel='triangle')
    
    
#    pylab.plot(mean_conv_r[-1])
#    pylab.plot(conv_r[-1][0])
#    pylab.plot(conv_r[-1][1])
#    pylab.plot(mean_conv_r[-1]-std_conv_r[-1],'--k')
#    pylab.plot(mean_conv_r[-1]+std_conv_r[-1],'--k')
#    pylab.show()

    return raw_r, conv_r

def conv_data(raw_r, conv_r, mean_conv_r, std_conv_r, bin=100, kernel='rectangle'):
        conv_r.append(misc.convolve(raw_r[-1], bin, kernel,single=False))       
        mean_conv_r.append(numpy.mean(conv_r[-1],axis=0))
        std_conv_r.append(numpy.std(conv_r[-1],axis=0))
        return conv_r, mean_conv_r, std_conv_r

def fit_pol(x,y):
    p=numpy.polyfit(x, y, 1, rcond=None, full=False)

    #f=lambda x: p[0]*x**3+p[1]*x**2+p[2]*x+p[3]
    #df=lambda x: 3*p[0]*x**2+2*p[1]*x+p[2]
    f=lambda x: p[0]*x+p[1]
    df=lambda x: p[0]    
    xr=numpy.linspace(min(x),max(x), 50)
    
    return xr, f(xr), df(xr), p
N_MSN=15000            
points=list(numpy.array([550, 650, 950])+10)


t=misc.timer()
save_result_at=OUTPUT_PATH+'/simulate_network_direct_onoff_vs_rate_direct'+str(N_MSN)+'.plk'   
d_d1, ds_d1, s_d1, br_d1, p_d1=simulate_direct(1, N_MSN, save_result_at, 4, flag_bg=False )
t=misc.timer(t)


N_MSN=15000
resolution=10
max_prop=0.15
save_result_at=OUTPUT_PATH+'/simulate_network_indirect_onoff_vs_rate'+str(N_MSN)+'.plk'   
data2=simulate_indirect(1, N_MSN, save_result_at, 4, resolution,max_prop, flag_bg=False)


N_MSN=15000
resolution=4
max_prop=0.2
save_result_at=OUTPUT_PATH+'/simulate_network_indirect_onoff_vs_rate_non_diff'+str(N_MSN)+str(resolution)+'.plk'   
data3=simulate_indirect(1, N_MSN, save_result_at, 4, resolution, max_prop, flag_bg=True,  skip=1, 
                      lines=True)

#N_MSN=15000
#burst_rate=30.
#save_result_at=OUTPUT_PATH+'/simulate_network_direct_indirect_onoff_vs_rate'+str(N_MSN)+'.plk'   
#d_d1_d2, d_split_d1_d2, s_d1_d2, br_d1_d2, p1_d1_d2, p2_d1_d2=simulate_direct_vs_indirect(1, N_MSN,burst_rate, save_result_at, 4,flag_bg=True )
#
#
#N_MSN=15000
#save_result_at=OUTPUT_PATH+'/simulate_network_hyperdirect'+str(N_MSN)+'.plk'   
#raw_r, conv_r = simulate_hyperdirect(1, N_MSN, save_result_at, 4, flag_bg=True)

# DISPLAY
 #Inspect results
plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0+320.0, fontsize=16)
font_size_text = 8
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
ax_list.append( MyAxes(fig, [ .53,  .1,  .165, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .1,  .165, .2 ] ) )    # 
ax=ax_list[1]
plot_selection_vs_neurons_direct2(ax, ds_d1[points[0]], title='First 100 ms of the burst')

ax=ax_list[2]
plot_selection_vs_neurons_direct2(ax, ds_d1[points[2]], title='Last 100 ms of the burst')

ax=ax_list[3]
plot_selection_vs_neurons_direct_linear_fit_and_spread(ax, p_d1, s_d1)

ax=ax_list[4]
plot_selection_vs_neurons_indirect2(ax, data2[1][points[0]], title='First 100 ms of the burst', ylim=[30,130])

ax=ax_list[5]
plot_selection_vs_neurons_indirect2(ax, data2[1][points[2]], title='Last 100 ms of the burst',ylim=[30,130])

ax=ax_list[6]
plot_selection_vs_neurons_direct_linear_fit_and_spread(ax,data2[4], data2[2], ylim=[0,1.5])


ax=ax_list[7]
#plot_selection(ax, d_d1[points[0]])
plot_selection_vs_neurons_indirect2(ax, data3[1][points[0]], title='First 100 ms of the burst', ylim=[30,130])

ax=ax_list[8]
plot_selection_vs_neurons_indirect2(ax, data3[1][points[2]], title='Last 100 ms of the burst',ylim=[30,130])
#plot_selection(ax, d_d1[points[2]], legend=False)

ax=ax_list[9]
plot_selection_vs_neurons_direct_linear_fit_and_spread(ax, data3[4], data3[2], ylim=[0,1.5])

# #Inspect results
#plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0+320.0, fontsize=16)
#font_size_text = 8
#fig = pylab.figure( facecolor = 'w' )
#ax_list = []
#ax_list.append( MyAxes(fig, [ .075, .37, .135, .26 ] ) )    # text box
#ax_list.append( MyAxes(fig, [ .26,  .7,  .165, .2 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .53,  .7,  .165, .2 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .8,   .7,  .165, .2 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .26,  .4,  .165, .2 ] ) )    #     
#ax_list.append( MyAxes(fig, [ .53,  .4,  .165, .2 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .8,   .4,  .165, .2 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .26,  .1,  .165, .2 ] ) )    #     
#ax_list.append( MyAxes(fig, [ .53,  .1,  .165, .2 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .8,   .1,  .165, .2 ] ) )    # 
#
#ax=ax_list[1]
#plot_selection(ax, d_d1_d2[points[0]])
#
#ax=ax_list[2]
#plot_selection(ax, d_d1_d2[points[2]], legend=False)
pylab.show()


# dpi does not matter since svg and pdf are both vectorbased
fig.savefig( picture_dir + '/' + SNAME  + '.svg', format = 'svg') 
fig.savefig( picture_dir + '/' + SNAME  + '.pdf', format = 'pdf')



print t