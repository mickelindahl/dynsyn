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
from simulation_utils import simulate_network

model_name=os.getcwd().split('/')[-2]
picture_dir='/'.join(os.getcwd().split('/')[0:-3]) + '/pictures/'+model_name 
SNAME  = sys.argv[0].split('/')[-1].split('.')[0]

OUTPUT_PATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]
SELECTION_THR=5.
ADJUST_XDATA_MS=0.

def plot_example_SNR(ax, data, flag=0, title=''):
    time_bin=20
    
    if flag==1:
        colors=['b','g','r','m']#,misc.make_N_colors('jet',5)
    if flag==0:
    
        colors=['b','r','m']#,misc.make_N_colors('jet',5)
    
    #colors=['m','c',colors[1]]   
    #labels=[r'$\delta_{fac}^{MSN}$' , r'$\delta_{dep}^{GPe}$',  
    #        r'$\delta_{fac}^{MSN}$+$\delta_{dep}^{GPe}$']
    coords=[[0.05, 0.2], [ 0.05, 0.65], [0.05, 0.8]]


    for color, y_mean, y_std in zip(colors, data[0], data[1]):
        y_mean=y_mean[500:]
        y_std=y_std[500:]
        #spk_mean=numpy.mean(sp, axis=0)
        #spk_mean=misc.convolve(spk_mean, 100, 'triangle',single=True)[0]
        x=numpy.arange(1,len(y_mean)+1)
        ax.plot(x,y_mean,**{'color':color})  
        ax.fill_between(x,y_mean-y_std, y_mean+y_std, facecolor=color, alpha=0.5)


       
    
    lines = ax.lines
    for line in lines:
        line.set_xdata(line.get_xdata()-ADJUST_XDATA_MS)
        misc.slice_line(line, xlim=[0,1490])
    pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend   
    
    if flag==0:
        leg=ax.legend(lines,[ 'Burst 30 Hz in 3 % $MSN_{D1}$' , 
            'Added burst STN static syn','Added burst STN depressed syn'], loc=2)
    if flag==1:
       leg=ax.legend(lines,['Burst 30 Hz in 3 % $MSN_{D1}$' , 'Added burst 30 Hz in 4 % $MSN_{D2}$',  
            'Added burst STN static syn','Added burst STN depressed syn'], loc=2)
    
    frame  = leg.get_frame() 
    frame.set_visible(False) 
    pylab.setp(frame, edgecolor='w') 
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    pylab.setp(ltext, fontsize=12., backgroundcolor='w') 
   
    ax.plot([0, 1490],[SELECTION_THR]*len([0, 1490]),**{'color':'k',
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--',
                                                        'linewidth':2})
    ax.text( 0.9, 0.11,'Thr' , fontsize=14, transform=ax.transAxes, **{ 'color' : 'k'})   
       
    #ax.set_title('bin=%i'%(time_bin),**{'fontsize':12})
    ax.set_ylabel('Rate SNr (spikes/s)') 
    ax.set_xlabel('Time (ms)')
    ax.my_set_no_ticks( yticks=7, xticks=8 ) 

    ax.set_xlim(misc.adjust_limit([0,1500]))
    ax.set_ylim(misc.adjust_limit([0,40]))
    ax.set_title(title, fontsize=16)  
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

def simulate(load, save_at, n_exp, res, params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, sim_time, seed, start_rec,
                           I_e_add, threads, model_params, dis_conn_GPE_STN=False ):
    
    

    if not load:
        r=[]
        
        for i in range(n_exp):
            seed=i
            p=numpy.random.randint(1)
            start_rec+=p
            sim_time+=p
            layer_dic=simulate_network(params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, sim_time, seed, I_e_add, threads, 
                           start_rec, model_params, dis_conn_GPE_STN=dis_conn_GPE_STN) 
    
            layer_dic['SNR'].get_signal( 's', start=start_rec, stop=sim_time )
            signal=layer_dic['SNR'].signals['spikes']
            r.append(numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0))

        
        numpy.array(r)
        misc.pickle_save(r, save_at)
    else: 
        r=misc.pickle_load(save_at)    
    r=numpy.array(r)
    r=misc.convolve(r, res, 'triangle',single=False)    
    mr=numpy.mean(r, axis=0)
    mstd=numpy.std(r, axis=0)
    d=[mr, mstd]
    return d

params_msn_d1={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, 30, 0.1],
            'mod_times':[1,1500, 1500+500],  'n_mod':0}    
params_msn_d2={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, 30, 0.1],
            'mod_times':[1,1500, 1500+500],  'n_mod':0, 'focus':False, 
            'skip':1} 

synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p']

N_MSN=15000
sim_time=3000.
n_exp=5
threads=4
start_rec=500.0
res=100 #convolve kernel size 100
d1_precent=0.03
d2_precent=0.04
stn_percent=1.

types=['Only_D1','D1_and_D2','Add_STN_static','Add_STN_plastic']
mod_d1=[N_MSN*d1_precent, N_MSN*d1_precent,N_MSN*d1_precent, N_MSN*d1_precent]
mod_d2=[ 0,N_MSN*d2_precent,0, 0 ]
rate_stn=[0., 0., 189.*stn_percent, 189.*stn_percent]
stn_syn=['STN_SNR_ampa_s','STN_SNR_ampa_s','STN_SNR_ampa_s','STN_SNR_ampa_p3']

mr, std=[],[]
for mod1, mod2, ssyn, rstn, type in zip(mod_d1, mod_d2, stn_syn, rate_stn, types):
        
        save_at=OUTPUT_PATH+'/simulate'+type+str(n_exp)+str(N_MSN)+'.plk'
        
        model_params= {'misc':{'N_MSN':N_MSN},
                       'conns':{ 'MSN_D2_GPE':{ 'lines':False},
                        'STN_SNR':{'syn':ssyn}}, 
                       'neurons':{'MSN_D1':{'n':N_MSN},
                                  'MSN_D2':{'n':N_MSN}}}
        params_msn_d1.update({'n_mod':int(mod1)})
        params_msn_d2.update({'n_mod':int(mod2)})
        params_stn={'rate':219., 'mod':True, 'mod_rate':rstn, 'mod_times':[1500., 1500.+500.]}
                
        data=simulate(1, save_at, n_exp, res, params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, sim_time, 1, start_rec,
                            {}, threads, model_params)

        mr.append(data[0])
        std.append(data[1])


N_MSN=15000
types=['Only_D1', 'Add_STN_static','Add_STN_plastic']
mod_d1=[N_MSN*d1_precent,N_MSN*d1_precent, N_MSN*d1_precent]
mod_d2=[0.,0., 0. ]
rate_stn=[0., 189.*stn_percent, 189.*stn_percent]
stn_syn=['STN_SNR_ampa_s', 'STN_SNR_ampa_s','STN_SNR_ampa_p3']
dis_conn=['GPE','GPE','GPE']
n_exp=5
mr2, std2=[],[]
for mod1, mod2, ssyn, rstn, type , ds in zip(mod_d1, mod_d2, stn_syn, rate_stn, types, dis_conn):
        
        save_at=OUTPUT_PATH+'/simulate_dis_conn_GPE_STN'+type+str(n_exp)+str(N_MSN)+'.plk'
        
        model_params= {'misc':{'N_MSN':N_MSN},
                       'conns':{ 'MSN_D2_GPE':{ 'lines':False},
                        'STN_SNR':{'syn':ssyn}}, 
                       'neurons':{'MSN_D1':{'n':N_MSN},
                                  'MSN_D2':{'n':N_MSN}}}
        params_msn_d1.update({'n_mod':int(mod1)})
        params_msn_d2.update({'n_mod':int(mod2)})
        params_stn={'rate':219., 'mod':True, 'mod_rate':rstn, 'mod_times':[1500., 1500.+500.]}
                
        data2=simulate(1, save_at, n_exp, res, params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, sim_time, 1, start_rec,
                            {}, threads, model_params, dis_conn_GPE_STN=ds)

        mr2.append(data2[0])
        std2.append(data2[1])

 #Inspect results
# DISPLAY
plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0, fontsize=16)
font_size_text = 10
fig = pylab.figure( facecolor = 'w' )

ax_list = []
ax_list.append( MyAxes(fig, [ .075, .37, .135, .26 ] ) )    # text box
ax_list.append( MyAxes(fig, [ .26,  .6,  .165*2.312, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .26,  .1,  .165*2.312, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .8,   .6,  .165, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .26,  .1,  .165, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .53,  .1,  .165, .34 ] ) )    # 
ax=ax_list[1]
plot_example_SNR(ax, [mr, std], flag=1,title='STN and GPe converge in SNr')
ax=ax_list[2]
plot_example_SNR(ax, [mr2, std2], flag=0, title='STN and GPe do not converge in SNr')

ax=ax_list[0]
#plot_text(ax, info_string=s)
pylab.show()

# dpi does not matter since svg and pdf are both vectorbased
fig.savefig( picture_dir + '/' + SNAME  + '.svg', format = 'svg') 
fig.savefig( picture_dir + '/' + SNAME  + '.pdf', format = 'pdf')