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
from simulation_utils import simulate_network_direct_sep_freq, simulate_network_direct_only_snr, simulate_network


model_name=os.getcwd().split('/')[-2]
picture_dir='/'.join(os.getcwd().split('/')[0:-3]) + '/pictures/'+model_name 
SNAME  = sys.argv[0].split('/')[-1].split('.')[0]
SELECTION_THR=5.
OUTPUT_PATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]
def plot_example_snr(ax, SNR):
    time_bin=20
    signal=SNR.signals['spikes']
    colors = ['k']
    
    sim_time=2500
    
    hist=signal.spike_histogram(time_bin=1, normalized=True)
    spk_mean=numpy.mean(hist, axis=0)
    spk_mean=misc.convolve(spk_mean, 100, 'triangle',single=True)[0]
    time=numpy.arange(1,len(spk_mean)+1)
    ax.plot(time,spk_mean,'k')


    ax.my_set_no_ticks( yticks=6, xticks=7 ) 
    lines = ax.lines

    lines = ax.lines
    for line in lines:

        misc.slice_line(line, xlim=[50,2450])
    
    ax.plot([0, sim_time],[SELECTION_THR]*len([0, sim_time]),**{'color':[0.6,0.6,0.6],
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
    ax.text( 0.9, 0.14,'Thr' , transform=ax.transAxes, **{ 'color' : [0.6,0.6,0.6]})   
    misc.slice_line(lines[0], xlim=[0,sim_time])
    
    ax.plot([1000,1500],[39,39],color='k', marker='|')
    ax.text( 0.47, 0.84,'stop' , transform=ax.transAxes, **{ 'color' : 'k' }) 

    ax.plot([500,700],[30,30],color='k', marker='|')
    ax.text( 0.23, 0.67,'1st' , transform=ax.transAxes, **{ 'color' : 'k' })
    
    ax.plot([1500,1700],[30,30],color='k', marker='|') 
    ax.text( 0.60, 0.67, '2nd' , transform=ax.transAxes, **{ 'color' : 'k' }) 

    ax.set_xlim(misc.adjust_limit([0,sim_time]))
    ax.set_ylim(misc.adjust_limit([0,45]))
    
    ax.set_ylabel('Firing rate SNr (spikes/s)')

def plot_rate_first_and_second_bursts_full(ax, x, data):
    colors = ['k']
    max_delay=3100
    
    ax.text( 0.1, 0.4, r'$\delta_{fac}^{MSN}$' , 
             fontsize=pylab.rcParams['text.fontsize'], 
             transform=ax.transAxes, 
             **{ 'color' : colors[0]})
     
    line1=ax.plot(1,1,'-k')
    line2=ax.plot(1,1,'-.k')
    leg=ax.legend([line1, line2],['1st', '2nd'], loc='best')
    frame  = leg.get_frame() 
    frame.set_visible(False)
    
    ax.plot([0, max_delay],[SELECTION_THR]*len([0, max_delay]),**{'color':[0.6,0.6,0.6],
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
    ax.text( 0.8, 0.3,'Thr' , transform=ax.transAxes, **{ 'color' : [0.6,0.6,0.6]})   

    
    ax.plot(x,data[0],**{'label':'1st', 'color': colors[0]})
    ax.plot(x,data[1],**{'label':'2nd', 'color': colors[0], 
                                    'linestyle':'-.'})
    ax.set_ylabel('Firing rate SNr (spikes/s)') 
    ax.set_xlabel('Burst stop (ms)')
    ax.my_set_no_ticks( yticks=6, xticks = 5 )

    lines = ax.lines
    misc.slice_line(lines[2], xlim=[0,max_delay])
    misc.slice_line(lines[3], xlim=[0,max_delay])  

    ax.set_xlim(misc.adjust_limit([0,max_delay]))
    ax.set_ylim(misc.adjust_limit([3,10]))
    
        
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


params_msn_d1={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, 20, 0.1, 20, 0.1],
            'mod_times':[1,1000, 1000+500, 2000, 2500],  'n_mod':60*10}    
params_msn_d2={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, 20, 0.1],
            'mod_times':[1,1000, 1000+500],  'n_mod':0, 'focus':True, 
            'skip':1} 
params_stn={'rate':350., 'mod':False,'mod_rate':0., 'mod_times':[1000., 1000.+500.]} 
sim_time=3000.

#synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p_stoc']
synapse_models_mod=['MSN_SNR_gaba_p1']

save_result_at=OUTPUT_PATH+'/simulate_network_example_robust_pause.plk'
if 0:
    layer_dics=[]
    for model,seed in zip(synapse_models_mod, [1]):
        synapse_models=[model, 'GPE_SNR_gaba_p']
        layer_dics.append(simulate_network(params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, sim_time=sim_time, seed=seed,
                           I_e_add={'SNR':300, 'STN':0,'GPE':30}, threads=4, 
                           start_rec=500.))    
    misc.pickle_save(layer_dics, save_result_at)  
else:
    layer_dics=misc.pickle_load(save_result_at)  
 

params_msn_d1={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, 20, 0.1, 20, 0.1],
            'mod_times':[1,1000, 1000+500],  'n_mod':60}    
params_msn_d2={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, 20, 0.1],
            'mod_times':[1,1000, 1000+500],  'n_mod':0,'focus':False, 'skip':1} 
params_stn={'rate':350., 'mod':False,'mod_rate':0., 'mod_times':[1000., 1000.+500.]} 
synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p']

N_MSN=15000.
seed=[1]
sim_time=5600.

#synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p_stoc']
delays=numpy.linspace(100, 3000, 11)
save_result_at=OUTPUT_PATH+'/simulate_network_direct_robust_pause.plk'
model_params={'neurons':{'MSN_D1':{'n':N_MSN},
                         'MSN_D2':{'n':N_MSN}}}
if 0:
    spk=[]
    for delay in delays:
        params_msn_d1.update({'mod_times':[1,1000, 1500,1500+delay, delay+2000]})
        spk.append(simulate_network_direct_only_snr(params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, sim_time=sim_time, seed=[int(delay)],
                           I_e_add={'SNR':200, 'STN':0,'GPE':30}, threads=4, start_rec=500.,
                           model_params=model_params))
    misc.pickle_save(spk, save_result_at)
else: 
    spk=misc.pickle_load(save_result_at)
 
#pylab.plot(numpy.array(spk)[0,0,:])
#pylab.show() 
mean_rates_b1=numpy.mean(numpy.array(spk)[:,0,510:710], axis=1)
mean_rates_b2=[]
for i in range(len(delays)):
    delay=delays[i]
    mean_rates_b2.append(numpy.mean(numpy.array(spk)[i,0,1010+delay:delay+1210], axis=0))
mean_rates_b2=numpy.array(mean_rates_b2)

 #Inspect results
 #Inspect results
plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0, fontsize=12)
font_size_text = 8
fig = pylab.figure( facecolor = 'w' )
ax_list = []
ax_list.append( MyAxes(fig, [ .075, .37, .135, .26 ] ) )    # text box
ax_list.append( MyAxes(fig, [ .26,  .6,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .53,  .6,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .26,  .1,  .165, .34 ] ) )    #     
ax_list.append( MyAxes(fig, [ .8,   .6,  .165, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .53,  .1,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .1,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .1,  .165, .2 ] ) )    # 


ax=ax_list[1]
plot_rate_first_and_second_bursts_full(ax, delays, [mean_rates_b1, mean_rates_b2])

ax=ax_list[3]
plot_example_snr(ax, layer_dics[0]['SNR'])


ax=ax_list[0]
#plot_text(ax, info_string=s)
pylab.show()

# dpi does not matter since svg and pdf are both vectorbased
fig.savefig( picture_dir + '/' + SNAME  + '.svg', format = 'svg') 
fig.savefig( picture_dir + '/' + SNAME  + '.pdf', format = 'pdf')