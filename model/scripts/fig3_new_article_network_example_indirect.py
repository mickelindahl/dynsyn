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
from matplotlib import rc


model_name=os.getcwd().split('/')[-2]
picture_dir='/'.join(os.getcwd().split('/')[0:-3]) + '/pictures/'+model_name 
SNAME  = sys.argv[0].split('/')[-1].split('.')[0]

OUTPUT_PATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]
SELECTION_THR=5.
ADJUST_XDATA_MS=500.

def plot_example_firing_rate_MSN(ax, layer, name,  color='b', ylim=[]):
    time_bin=20


    signal=layer.signals['spikes']
    #signal.my_firing_rate(bin=time_bin, display=ax,
    #                      kwargs={'color':color})

    print name, 'CV:', numpy.mean(signal.cv_isi())
    print name, 'mean:', signal.mean_rate()
    print name, 'std:', signal.mean_rate_std()
    hist=signal.spike_histogram(time_bin=1, normalized=True)
    spk_mean=numpy.mean(hist, axis=0)
    spk_mean=misc.convolve(spk_mean, 100, 'triangle',single=True)[0]
    time=numpy.arange(1,len(spk_mean)+1)
    ax.plot(time,spk_mean)
    
    spk_mean=numpy.mean(hist[numpy.array(layer.id_mod)-layer.ids[0],:], axis=0)
    spk_mean=misc.convolve(spk_mean, 100, 'triangle',single=True)[0]
    time=numpy.arange(1,len(spk_mean)+1)
    ax.plot(time,spk_mean,'r')

    pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend
    ax.set_ylabel('Rate '+name+' (spikes/s)') 
    ax.set_xlabel('Time (ms)')    
    ax.my_set_no_ticks( yticks=6, xticks = 5 ) 
    ax.set_xlim(misc.adjust_limit([0,1500]))
    ax.set_ylim(misc.adjust_limit(ylim)) 

    ax.text( 0.05, 0.85, 'Bursting \n$MSN_{D2}$' , backgroundcolor='w',
             transform=ax.transAxes, fontsize=font_size_text+2, **{ 'color' : 'r' }) 
    ax.text( 0.05, 0.75, 'Average \nover all' , backgroundcolor='w',
             transform=ax.transAxes, fontsize=font_size_text+2, **{ 'color' : 'b' })  

def plot_SNR_rate(ax, d_SNR, d_SNR2, name='SNr'):
    
    d_SNR.extend(d_SNR2)
    colors=['m', 'm','b','b']
    linstyles=['--','-','--','-']
    labels=[ 'Diffuse $MSN_{D2}$  \nwith $ref_{30 Hz}^{GPe}$' ,
             'Diffuse $MSN_{D2}$ \nwith $dep^{GPe}$',
             'Non-diffuse $MSN_{D2}$ \nwith $ref_{32 Hz}^{GPe}$']
    for d, c, ls in zip(d_SNR, colors, linstyles):

        y_mean=d[0][ADJUST_XDATA_MS:]
        y_std=d[1][ADJUST_XDATA_MS:]
        x=numpy.arange(1,len(y_mean)+1)
        ax.plot(x,y_mean,**{'color':c, 'linestyle':ls})  
        ax.fill_between(x,y_mean-y_std, y_mean+y_std, facecolor=c, alpha=0.25)
    
    ax.text( 0.05, 0.75, 'Diffuse $MSN_{D2}$ to \nGPe projection' , 
             fontsize=font_size_text,
             backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : colors[0] }) 
    ax.text( 0.05, 0.85,'Non-Diffuse $MSN_{D2}$ to \nGPe projection',
             fontsize=font_size_text,
              backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : colors[2]})    
    
    line1=ax.plot(1,1,'-k')
    line2=ax.plot(1,1,'--k')
    pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend
    leg=ax.legend([line1, line2],['$dep^{GPe}$', '$ref_{30 Hz}^{GPe}$'], 
                  loc='best')
    frame  = leg.get_frame() 
    frame.set_visible(False) 
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    pylab.setp(frame, edgecolor='w') 
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    pylab.setp(ltext, fontsize=pylab.rcParams['font.size'], backgroundcolor='w')
    
    
    ax.set_ylim(misc.adjust_limit([10,80]))
    ax.set_xlim(misc.adjust_limit([0,1500]))
    
    ax.set_ylabel('Rate '+name+' (spikes/s)') 
    ax.set_xlabel('Time (ms)')    
    ax.my_set_no_ticks( yticks=6, xticks = 8 ) 

def plot_GPE(ax, data, data2, name='GPe'):
       
    data.extend(data2)
    colors=['b','b','b','m']
    linstyles=['--','-','-','-']
    for d, c, ls in zip(data, colors, linstyles):

        y_mean=d[0][ADJUST_XDATA_MS:]
        y_std=d[1][ADJUST_XDATA_MS:]
        x=numpy.arange(1,len(y_mean)+1)
        ax.plot(x,y_mean,**{'color':c, 'linestyle':ls})  
        ax.fill_between(x,y_mean-y_std, y_mean+y_std, facecolor=c, alpha=0.25)
    
    
    ax.text( 0.05, 0.9, 'Non-Diffuse $MSN_{D2}$ to \nGPe projection', fontsize=font_size_text,
             backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : 'b' } )      
    ax.text( 0.55, 0.75, 'Average response, \nnon-diffuse (dotted)', fontsize=font_size_text,
             backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : 'k' }) 
    ax.text( 0.05, 0.6, 'Disinhibited \nsubpopulation', fontsize=font_size_text,
             backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : 'k' })   
    ax.text( 0.55, 0.45, 'Directly inhibited \nsubpopulation', fontsize=font_size_text,
             backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : 'k' }) 
    ax.text( 0.05, 0.2, 'Diffuse $MSN_{D2}$ to \nGPe projection', fontsize=font_size_text,
            backgroundcolor='w',
             transform=ax.transAxes, **{ 'color' : 'm' }) 
    #ax.text( 0.55, 0.2, 'Average response \nin all ', fontsize=font_size_text,
    #         backgroundcolor='w',
    #         transform=ax.transAxes, **{ 'color' : 'k' }) 
  
    pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend
    ax.set_ylim(misc.adjust_limit([0,60]))
    ax.set_xlim(misc.adjust_limit([0,1500]))      
    
    ax.set_ylabel('Rate '+name+' (spikes/s)') 
    ax.set_xlabel('Time (ms)')    
    ax.my_set_no_ticks( yticks=6, xticks = 8 ) 
    
def plot_example_raster(ax, layer, name):
    global ADJUST_XDATA_MS
    
    ids_base=numpy.array(list(set(layer.ids).difference(layer.id_mod)))
    mod_at=10000.
    ids_base1=ids_base[0:mod_at]
    ids_base2=ids_base[mod_at:]
    n_mod=len(layer.id_mod)
    
    layer.signals['spikes'].my_raster_plot(display=ax, id_list=ids_base1, reduce=6,
                                      kwargs={'color':'b', 'zorder':1})  
    layer.signals['spikes'].my_raster_plot(display=ax, id_list=ids_base2, reduce=6,
                                      kwargs={'color':'b', 'zorder':1}, id_start=mod_at+n_mod,)  
    if layer.id_mod:
        layer.signals['spikes'].my_raster_plot(display=ax, id_list=layer.id_mod, reduce=6,
                                               id_start=mod_at,
                                      kwargs={'color':'r', 'zorder':1})  
        
    pylab.rcParams.update({'path.simplify':False} )
    lines = ax.lines
    ax.set_ylabel(name+' id')
    ax.my_set_no_ticks( yticks=6, xticks = 5 )
    ax.set_ylim([0,15000])
    #ax.set_ylim([layer.ids[0],layer.ids[-1]])
    
    ax.text( 0.05, 0.05, 'Non-bursting $MSN_{D2}$' , backgroundcolor='w',
             transform=ax.transAxes,fontsize=16, **{ 'color' : 'b' })  
    ax.text( 0.05, 0.15, 'Bursting $MSN_{D2}$' , backgroundcolor='w',
             transform=ax.transAxes, fontsize=16, **{ 'color' : 'r' }) 
    
    for line in lines:
        line.set_xdata(line.get_xdata()-1000)
    
    ax.set_xlim(misc.adjust_limit([0,1500]))

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

def simulate_diffuse(load, save_at, n_exp, mod, res=100.):
    #synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p_stoc']
    synapse_models_mod=[ 'GPE_SNR_gaba_s_ref', 'GPE_SNR_gaba_p']
    
    model_params= {'misc':{'N_MSN':N_MSN},
                   'conns':{ 'MSN_D2_GPE':{ 'lines':False}},
                   'neurons':{'MSN_D1':{'n':N_MSN},
                              'MSN_D2':{'n':N_MSN}}} 
    
    params_msn_d2={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, 20, 0.1],
            'mod_times':[1,1500, 1500+500],  'n_mod':mod, 'focus':False, 
            'skip':1} 
    if not load:
        d={}
        for model,seed in zip(synapse_models_mod, [1,1]):
            synapse_models=['MSN_SNR_gaba_p1', model]
            spk_mean=[]
            for i in range(n_exp):
                seed=i
                layer_dic=simulate_network(params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, sim_time=sim_time, seed=seed,
                           I_e_add={'SNR':300, 'STN':0,'GPE':30}, threads=4, 
                           start_rec=start_rec, model_params=model_params) 
                
                layer_dic['SNR'].get_signal( 's', start=start_rec, stop=sim_time )
                signal=layer_dic['SNR'].signals['spikes']
                spk_mean.append(numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0))
                
                layer_dic['GPE'].get_signal( 's', start=start_rec, stop=sim_time )
                signal=layer_dic['GPE'].signals['spikes']
                spk_mean.append(numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0))
   
            d[model]=spk_mean
            
            misc.pickle_save(d,save_at)     
    else:
        d=misc.pickle_load(save_at)  
        
    
    d_SNR,d_GPE=[],[]
    for r in d.values():
        r=numpy.array(r)
        r=misc.convolve(r, res, 'triangle',single=False)
        d_SNR.append([numpy.mean(r[0::2,:], axis=0),
                      numpy.std(r[0::2,:], axis=0)]) # SNR
        d_GPE.append([numpy.mean(r[1::2,:], axis=0),
                      numpy.std(r[1::2,:], axis=0)]) # SNR
        
    return d, d_SNR, d_GPE

def simulate_non_diffuse(load, save_at, n_exp, mod, skip, res=100):
    #synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p_stoc']
    synapse_models_mod=[ 'GPE_SNR_gaba_s_ref', 'GPE_SNR_gaba_p']
    

    model_params= {'misc':{'N_MSN':N_MSN},
                   'conns':{ 'MSN_D2_GPE':{ 'lines':True}},
                   'neurons':{'MSN_D1':{'n':N_MSN},
                              'MSN_D2':{'n':N_MSN}}}
    params_msn_d2={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, 20.0, 0.1],
            'mod_times':[1,1500, 1500+500],  'n_mod':mod, 'focus':False, 
            'skip':skip}
    
    if not load:
        d={}
        for model,seed in zip(synapse_models_mod, [1,1]):
            synapse_models=['MSN_SNR_gaba_p1', model]
            spk_mean=[]
            for i in range(n_exp):
                seed=i
                layer_dic=simulate_network(params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, sim_time=sim_time, seed=seed,
                           I_e_add={}, threads=4, 
                           start_rec=start_rec, model_params=model_params) 
                
                signal=layer_dic['SNR'].signals['spikes']
                spk_mean.append(numpy.mean( signal.spike_histogram(time_bin=1, normalized=True),axis=0) )
            
                signal=layer_dic['GPE'].signals['spikes']
                hist=signal.spike_histogram(time_bin=1, normalized=True)
                spk_mean.append(numpy.mean( hist, axis=0) )
               
                down, up=[],[]
                for i in range(hist.shape[0]):
                    if numpy.mean(hist[i,1000:1500])<numpy.mean(hist[i,0:1000]):
                        down.append(hist[i,:])
                    else:
                        up.append(hist[i,:])
                
                spk_mean.append(numpy.mean(numpy.array(down), axis=0))
                spk_mean.append(numpy.mean(numpy.array(up), axis=0)) 
            
            d[model]=spk_mean
            
            misc.pickle_save(d,save_at)     
    else:
        d=misc.pickle_load(save_at)  
    
    d_SNR, d_GPE=[],[]
    for r in d.values():
        r=numpy.array(r)
        r=misc.convolve(r, res, 'triangle',single=False)
        d_SNR.append([numpy.mean(r[0::4,:], axis=0),
                      numpy.std(r[0::4,:], axis=0)]) # SNR

              
    d_GPE.append([numpy.mean(r[1::4,:], axis=0),
                  numpy.std(r[1::4,:], axis=0)]) # GPE
    d_GPE.append([numpy.mean(r[2::4,:], axis=0),
                  numpy.std(r[2::4,:], axis=0)]) # GPE down
    d_GPE.append([numpy.mean(r[3::4,:], axis=0),
                  numpy.std(r[3::4,:], axis=0)]) # GPE up 
        
        
    return d, d_SNR, d_GPE
font_size_text = 12
N_MSN=15000.
params_msn_d1={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [],
               'mod_times':[],  'n_mod':0}    
params_stn={'rate':219., 'mod':False,'mod_rate':0., 'mod_times':[]} 

sim_time=3000.
n_exp=5.
start_rec=500.

mod1=int(0.05*N_MSN)
save_at=OUTPUT_PATH+'/simulate_network_diffuse'+str(N_MSN)+'_'+str(n_exp)+'.plk'
d, d_SNR, d_GPE=simulate_diffuse(1, save_at, n_exp, mod1)

N_MSN=15000.
# =0.05 vs .65 skip=9
skip=9
mod2=int(0.065*N_MSN)
save_at=OUTPUT_PATH+'/simulate_network_non_diffuse'+str(N_MSN)+'.plk'
d, d_SNR2, d_GPE2=simulate_non_diffuse(1, save_at, n_exp, mod2, skip)
#plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0, fontsize=12)
#ax=pylab.subplot(211)
#plot_GPE(ax ,d_GPE2, d_GPE)
##
#ax=pylab.subplot(212)
#plot_SNR_rate(ax, d_SNR, d_SNR2)
#pylab.show()


N_MSN=15000.
#synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p_stoc']
synapse_models_mod=[ 'GPE_SNR_gaba_s_ref', 'GPE_SNR_gaba_p']
model_params= {'misc':{'N_MSN':N_MSN},
               'conns':{ 'MSN_D2_GPE':{ 'lines':False}},
               'neurons':{'MSN_D1':{'n':N_MSN},
                          'MSN_D2':{'n':N_MSN}}} 
params_msn_d1={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [],
            'mod_times':[],  'n_mod':0}    
params_stn={'rate':250., 'mod':False,'mod_rate':0., 'mod_times':[]} 
params_msn_d2={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, 20.0, 0.1],
            'mod_times':[1,2000, 2000+500],  'n_mod':mod1, 'focus':False, 
            'skip':1}

save_result_at=OUTPUT_PATH+'/simulate_network_example_indirect'+str(N_MSN)+'.plk'
if 0:

    synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p']
    layer=simulate_network(params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, sim_time=sim_time, seed=1,
                           I_e_add={}, threads=4, 
                           start_rec=1500., model_params=model_params)   
    misc.pickle_save(layer, save_result_at)  
else:
    layer=misc.pickle_load(save_result_at)  

 
#Inspect results
plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0+275.0, fontsize=16)
font_size_text = 14
fig = pylab.figure( facecolor = 'w' )
ax_list = []
ax_list.append( MyAxes(fig, [ .075, .37, .135, .26 ] ) )    # text box
ax_list.append( MyAxes(fig, [ .26,  .7,  .165, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .53,  .7,  .165, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .7,  .165, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .26,  .4,  .165*2.312, .2 ] ) )    #     
#ax_list.append( MyAxes(fig, [ .53,  .4,  .165, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .4,  .165, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .26,  .1,  .165*2.312, .2 ] ) )    #     
#ax_list.append( MyAxes(fig, [ .53,  .1,  .165, .2 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .1,  .165, .2 ] ) )    # 



ax=ax_list[1]
plot_example_raster(ax, layer['MSN_D2'], r'$MSN_{D2}$',)

ax=ax_list[2]
plot_example_firing_rate_MSN(ax, layer['MSN_D2'],r'$MSN_{D2}$',ylim=[0,25])

ax=ax_list[4]
#plot_GPE_diffuse(ax, d_GPE)
plot_GPE(ax, d_GPE2, d_GPE)

ax=ax_list[6]
plot_SNR_rate(ax, d_SNR, d_SNR2)


pylab.show()

# dpi does not matter since svg and pdf are both vectorbased
fig.savefig( picture_dir + '/' + SNAME  + '.svg', format = 'svg') 
fig.savefig( picture_dir + '/' + SNAME  + '.pdf', format = 'pdf')
