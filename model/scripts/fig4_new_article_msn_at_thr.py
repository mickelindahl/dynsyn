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

model_name=os.getcwd().split('/')[-2]
picture_dir='/'.join(os.getcwd().split('/')[0:-3]) + '/pictures/'+model_name 
SNAME  = sys.argv[0].split('/')[-1].split('.')[0]

# Add model, code and current directories to python path
sys.path.append(os.getcwd())  
sys.path.append(model_dir)
sys.path.append(code_dir+'/nest_toolbox') 
from src import my_nest, misc, my_topology, plot_settings
from src.my_axes import MyAxes 
import nest.topology as tp
from simulation_utils import simulate_network, simulate_network_poisson
import scipy.optimize as opt

OUTPUT_PATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]

def plot_thr_prop(ax, thr_prop, rate_MSN):    
    colors=['m','m', 'b','g']
    linstyles=['-','-','--','--']   
    
    labels=[r'$ref_{init}^{MSN_{D1}}$' , r'$ref_{max}^{MSN_{D1}}$',  
            r'$fac_{MSN_{D1}}$']  
    
    coords=[[0.03, 0.8], [ 0.05, 0.07], [0.03, 0.50]]   
    
    
    for y, ls, c in zip(thr_prop, linstyles, colors):
        
        ax.plot(rate_MSN, y*100.0, **{'color':c,'linestyle':ls}) 
    
    ax.set_xlabel(r'Rate $MSN_{D1}$ (spikes/s)') 
    ax.set_ylabel(r'Bursting at thr (%)')
    
    
    colors=['b','g','m']   
    labels=[r'$ref_{init}^{MSN_{D1}}$' , r'$ref_{max}^{MSN_{D1}}$',  
            r'$fac^{MSN_{D1}}$', ]
    coords=[[0.4, 0.52], [ 0.05, 0.22], [0.1, 0.74]]
    
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['text.fontsize'], 
                 **{'color': color})
        
    colors=['k','k']
    labels=['First 100 ms of the \nburst','Last 100 ms of the \nburst']
    coords=[[0.3, 0.54], [0.3, 0.34]]
    
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=font_size_text, 
                 **{'color': color})
        
    #line1=ax.plot(1,1,'-k')
    #line2=ax.plot(1,1,'--k')
    pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend
    #leg=ax.legend([line1, line2],[, ], loc='best')
    #frame  = leg.get_frame() 
    #frame.set_visible(False) 
    #pylab.setp(frame, edgecolor='w') 
    #ltext  = leg.get_texts()  # all the text.Text instance in the legend
    #pylab.setp(ltext, fontsize=font_size_text, backgroundcolor='w') 
    
  
    ax.set_xlim(misc.adjust_limit([6,48]))
    ax.set_ylim(misc.adjust_limit([-3,30]))   

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

def error_fun(x, sim_time, interval, syn, burst_rate, N_MSN):
        
        p_mod_msn=x[0]
        target_thr=5.0 
        n_mod=int(p_mod_msn*N_MSN)
        n_exp=1
        params_msn_d1={'base_rates':[0.1], 'base_times':[1.], 'mod_rates': [0.1, float(burst_rate), 0.1],
                       'mod_times':[1,1000, 1000+500],  'n_mod':n_mod, 
                       'bg_rate':[0.1*(500-500*p_mod_msn)]}    
        
        model_params={'misc':{'N_MSN':N_MSN},
                      'conns':{'MSN_D1_SNR':{'p':500./float(N_MSN)},
                               'MSN_D2_GPE':{'p':500./float(N_MSN),
                                     'lines':False}},
                     'neurons':{'MSN_D1':{'n':n_mod},
                         'MSN_D2':{'n':0.},
                         'MSN_D1_bg':{'n':300,
                                  'lesion':False},
                         'MSN_D2_bg':{'n':300,
                                      'lesion':False},  
                         'GPE': {'paused':False}}}
        seed=0
        
        synapse_models=[syn, 'GPE_SNR_gaba_s_ref']
        
        if n_mod>=1:
            e_tmp=[]
            for n in range(n_exp):
                layer_dic=simulate_network_poisson(params_msn_d1, params_msn_d2, params_stn,
                           synapse_models, sim_time, seed, {}, threads, 
                           start_rec, model_params)  
            
        
                layer_dic['SNR'].get_signal( 's', start=start_rec, stop=sim_time )
                signal=layer_dic['SNR'].signals['spikes']
                m_hist=numpy.mean(signal.spike_histogram(time_bin=1, normalized=True), axis=0)
                e_tmp.append(numpy.mean(m_hist[interval[0]:interval[1]], axis=0)-target_thr)
            e=sum(e_tmp)/len(e_tmp)
            
        else:
            
            e=(-n_mod*10)**2

             
        print e**2
        return e

def my_opt(x, sim_time, interval, syn, burst_rate, N_MSN, maxiter=20):
        upper=1
        lower=0
        go=True
        i=0
    
        e_rel=[]
        e_vec=[]
        x_vec=[]
        while go:
            e=error_fun([x], sim_time, interval, syn, burst_rate, N_MSN)
            
            if e>0:
                lower=x
                x1=x+(upper-lower)/2.
        
            if e<0:
                upper=x
                x1=x-(upper-lower)/2.
                
            if (maxiter<i):
                go=False
                if i>1:
                    if numpy.abs(e_rel[-1])<0.01:
                        go=False
                
            else:
                x=x1
                i+=1
            
            e_vec.append(e)
            x_vec.append(x)
            
            if i>1:
                e_rel.append(sum(e_vec[0:-1])/len(e_vec[0:-1])-sum(e_vec)/len(e_vec))
                
        
        allvecs=numpy.array([e_vec, x_vec])   
        return x, e, i, allvecs
    
def my_fmin(load, save_at, x0, interval, syn, N_MSN, burst_rate):
    
    x=x0
    if not load:

        x, e, i, allvecs=  my_opt(x, sim_time, interval, syn, burst_rate, N_MSN, maxiter=10)  
        misc.pickle_save([x,e,i, allvecs], save_at)
    else:
        [x,e,i, allvecs]=misc.pickle_load(save_at)        
    return x, e 

def fmin(load, save_at, x0, interval, syn, N_MSN, burst_rate):
    

     #[current, w_GPE_STN]
    args=(sim_time, interval,  syn, burst_rate, N_MSN )
    if not load:
        [xopt,fopt, iter, funcalls , warnflag, allvecs] = opt.fmin(error_fun, x0, args=args, maxiter=20, maxfun=10,full_output=1, retall=1)

        misc.pickle_save([xopt,fopt, iter, funcalls , warnflag, allvecs], save_at)
    else:
        [xopt,fopt, iter, funcalls , warnflag, allvecs]=misc.pickle_load(save_at)        
    return xopt,  fopt 

params_msn_d2={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, 0.1, 0.1],
            'mod_times':[1.,1500., 1500.+500.],  'n_mod':0, 'focus':False, 
            'skip':1, 'bg_rate':[0.1*(500-500*0)]} 
params_stn={'rate':219., 'mod':False,'mod_rate':400., 'mod_times':[1000., 1000.+500.]} 


threads=1

intervals=numpy.array([ [0,100.],[400,500.],[0,500.],[0,500.]])+10.
N_MSN=15000.
burst_rates=numpy.arange(6,50,1)
synapse_models=['MSN_SNR_gaba_p1', 'MSN_SNR_gaba_p1','MSN_SNR_gaba_s_min','MSN_SNR_gaba_s_max',]
sim_time=2000.
start_rec=1000.

load_to=50
d=[]
x0=0.25*0.5 # Inital guess

for br in burst_rates:
    for inter, syn in zip(intervals,synapse_models):
        if load_to>br:load=1
        else: load=0
        print br
        save_at=OUTPUT_PATH+'/simulate_network_fmin'+str(br)+syn+str(inter[0])+'-'+str(inter[0])+'.plk' 
        #x, e=fmin(0,save_at, x0, inter, syn, N_MSN, br)
        x, e=my_fmin(load,save_at, x0, inter, syn, N_MSN, br)
        d.append([x,e, inter[0],inter[1]])
        
        
        print x,e
        x0=x*0.5

d=numpy.array(d)
thr_prop=numpy.array([d[0::4,0], d[1::4,0], d[2::4,0], d[3::4,0]])


f_MSN=10.
N_MSN=int(1500*f_MSN)
p_mod_msn=0.01
p_mod_msn_d2=0.04


 #Inspect results
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
plot_thr_prop(ax, thr_prop, burst_rates)

pylab.show()

# dpi does not matter since svg and pdf are both vectorbased
fig.savefig( picture_dir + '/' + SNAME  + '.svg', format = 'svg') 
fig.savefig( picture_dir + '/' + SNAME  + '.pdf', format = 'pdf')
