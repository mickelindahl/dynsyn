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
    print spk_mean.shape
    time=numpy.arange(1,len(spk_mean)+1)
    ax.plot(time,spk_mean)
    
    shift=0
    if layer.id_mod: 
        print numpy.array(layer.id_mod)-layer.ids[0]-1
        spk_mean=numpy.mean(hist[numpy.array(layer.id_mod)-layer.ids[0],:]-shift, axis=0)
        spk_mean=misc.convolve(spk_mean, 50, 'triangle',single=True)[0]
        print spk_mean.shape
        time=numpy.arange(1,len(spk_mean)+1)
        ax.plot(time,spk_mean,'r')

    pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend
    
    ax.set_ylabel('Rate  '+name+' (spikes/s)') 
    ax.set_xlabel('Time (ms)')    
    ax.my_set_no_ticks( yticks=6, xticks = 5 ) 
    ax.set_xlim(misc.adjust_limit([0,1500]))
    ax.set_ylim(misc.adjust_limit(ylim))    

    ax.text( 0.05, 0.65, 'Bursting \n$MSN_{D1}$' , backgroundcolor='w',
             transform=ax.transAxes, fontsize=font_size_text, **{ 'color' : 'r' }) 
    ax.text( 0.05, 0.85, 'Average \nover all' , backgroundcolor='w',
             transform=ax.transAxes, fontsize=font_size_text, **{ 'color' : 'b' })  
    
    for line in ax.lines:
        line.set_xdata(line.get_xdata()-500)
        
def plot_SNR_rate(ax, d_SNR, name='SNr',ylim=[0,50]):
    
    coords=[[0.1, 0.22], [0.2, 0.74], [0.65, 0.45]]
    labels=['$ref_{max}^{MSN_{D1}}$', '$ref_{init}^{MSN_{D1}}$',  '$fac^{MSN_{D1}}$']
    colors=['g','b','m']
    linestyles=['--','--','-']
    for d, c, ls in zip(d_SNR, colors, linestyles):

        y_mean=d[0][500:2000]
        y_std=d[1][500:2000]
        x=numpy.arange(1,len(y_mean)+1)
        ax.plot(x,y_mean,**{'color':c, 'linestyle':ls})  
        ax.fill_between(x,y_mean-y_std, y_mean+y_std, facecolor=c, alpha=0.5)
    
    
    for coord, label, color in zip(coords, labels, colors):
        ax.text( coord[0], coord[1], label , transform=ax.transAxes, 
                 fontsize=pylab.rcParams['font.size']+2,
                 **{ 'color' : color})  
    
    
    c='k'#[0.6,0.6,0.6]
    ax.plot([0, 1490],[SELECTION_THR]*len([0, 1490]),**{'color':c,
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
    ax.text( 0.9, 0.15,'Thr' , transform=ax.transAxes, **{ 'color' :c}) 
    
    ax.set_ylabel('Rate '+name+' (spikes/s)') 
    ax.set_xlabel('Time (ms)')    
    ax.my_set_no_ticks( yticks=6, xticks = 8 ) 
    
    pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend
    ax.set_xlim(misc.adjust_limit([0,1500]))
    ax.set_ylim(misc.adjust_limit(ylim)) 
        
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

 
    lines = ax.lines
    ax.set_ylabel(name+' id')
    ax.my_set_no_ticks( yticks=6, xticks = 5 )
    ax.set_ylim([0,15000])
    
    ax.text( 0.05, 0.05, 'Non-bursting $MSN_{D1}$' , backgroundcolor='w',
             transform=ax.transAxes, fontsize=font_size_text, **{ 'color' : 'b' })  
    ax.text( 0.05, 0.15, 'Bursting $MSN_{D1}$' , backgroundcolor='w',
             transform=ax.transAxes, fontsize=font_size_text, **{ 'color' : 'r' }) 
    
    for line in lines:
        line.set_xdata(line.get_xdata()-ADJUST_XDATA_MS)
    
    ax.set_xlim(misc.adjust_limit([0,1500]))
    for line in ax.lines:
        line.set_xdata(line.get_xdata()-500)
        
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

def simluate(load, save_at, n_exp, synapse_models_mod, mod, res=100. ):
    model_params= {'misc':{'N_MSN':N_MSN},
               'neurons':{'MSN_D1':{'n':N_MSN},
                          'MSN_D2':{'n':N_MSN}}}
    params_msn_d1={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, 20, 0.1],
            'mod_times':[1,1500, 1500+500],  'n_mod':mod}    
    if not load:
        d={}
        for model in synapse_models_mod:
            synapse_models=[model, 'GPE_SNR_gaba_p']
            spk_mean=[]
            for i in range(n_exp):
                seed=i
                layer_dic=simulate_network(params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, sim_time=sim_time, seed=seed,
                           I_e_add={}, threads=4, 
                           start_rec=start_rec, model_params=model_params) 
                
                layer_dic['SNR'].get_signal( 's', start=start_rec, stop=sim_time )
                signal=layer_dic['SNR'].signals['spikes']
                spk_mean.append(numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0))
            
            d[model]=spk_mean
            
            misc.pickle_save(d,save_at)     
    else:
        d=misc.pickle_load(save_at)  
        
    
    d_SNR,d_GPE=[],[]
    for r in d.values():
        r=numpy.array(r)
        r=misc.convolve(r, res, 'triangle',single=False)
        d_SNR.append([numpy.mean(r[0::1,:], axis=0),
                      numpy.std(r[0::1,:], axis=0)]) # SNR

        
    return d, d_SNR

prop_mod=0.04
N_MSN=15000.0

params_msn_d2={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [],
            'mod_times':[],  'n_mod':0, 'focus':True, 
            'skip':1} 
params_stn={'rate':219., 'mod':False,'mod_rate':0., 'mod_times':[]} 
sim_time=2700.
start_rec=500.0


#synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p_stoc']
synapse_models_mod=['MSN_SNR_gaba_s_min','MSN_SNR_gaba_s_max','MSN_SNR_gaba_p1']

n_exp=5
save_at=OUTPUT_PATH+'/simulate_network_exp'+str(N_MSN)+'_'+str(n_exp)+'.plk1'
d, d_SNR=simluate(1, save_at, n_exp, synapse_models_mod, int(N_MSN*prop_mod), res=100.)
#ax=pylab.subplot(111)

#pylab.show()


N_MSN=15000.0
params_msn_d1={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, 20, 0.1],
            'mod_times':[1,1500, 1500+500],  'n_mod':int(N_MSN*prop_mod)}    
model_params= {'misc':{'N_MSN':N_MSN},
               'neurons':{'MSN_D1':{'n':N_MSN},
                          'MSN_D2':{'n':N_MSN}}}
save_result_at=OUTPUT_PATH+'/simulate_network_example_direct.plk'
if 0:
        synapse_models=['MSN_SNR_gaba_s_min', 'GPE_SNR_gaba_p']
        layer=simulate_network(params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, sim_time=sim_time, seed=1,
                           I_e_add={}, threads=4, 
                           start_rec=500., model_params=model_params)    
        misc.pickle_save(layer, save_result_at)  
else:
    layer=misc.pickle_load(save_result_at)  
 
 #Inspect results
plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0, fontsize=16)
font_size_text = 16
fig = pylab.figure( facecolor = 'w' )
ax_list = []
ax_list.append( MyAxes(fig, [ .075, .37, .135, .26 ] ) )    # text box
ax_list.append( MyAxes(fig, [ .26,  .6,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .53,  .6,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .26,  .1,  .165*2.312, .34 ] ) )    #     
ax_list.append( MyAxes(fig, [ .8,   .6,  .165, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .53,  .1,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .1,  .165, .34 ] ) )    # 



ax=ax_list[1]
plot_example_raster(ax, layer['MSN_D1'], r'$MSN_{D1}$',)
ax=ax_list[2]
plot_example_firing_rate_MSN(ax, layer['MSN_D1'], r'$MSN_{D1}$',ylim=[0,25])
ax=ax_list[3]
plot_SNR_rate(ax, d_SNR)

pylab.show()

# dpi does not matter since svg and pdf are both vectorbased
fig.savefig( picture_dir + '/' + SNAME  + '.svg', format = 'svg') 
fig.savefig( picture_dir + '/' + SNAME  + '.pdf', format = 'pdf')