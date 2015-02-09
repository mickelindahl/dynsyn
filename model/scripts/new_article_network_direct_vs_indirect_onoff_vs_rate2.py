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
from simulation_utils import simulate_network_direct_onoff_vs_rate, simulate_network_indirect_onoff_vs_rate, simulate_network_hyper_direct_onoff_vs_rate


model_name=os.getcwd().split('/')[-2]
picture_dir='/'.join(os.getcwd().split('/')[0:-3]) + '/pictures/'+model_name 
SNAME  = sys.argv[0].split('/')[-1].split('.')[0]

OUTPUT_PATH  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]
SELECTION_THR=5.
def plot_selection_vs_neurons_direct(ax, MSNmeanRates, SNRmeanRates):
    colors=['b','g','m']   
    labels=[r'$\delta_{weak}^{MSN}$' , r'$\delta_{strong}^{MSN}$',  
            r'$\delta_{fac}^{MSN}$']
    coords=[[0.4, 0.42], [ 0.02, 0.001], [0.15, 0.78]]
    
    
    for i, color in enumerate(colors):
        ax.plot(MSNmeanRates,SNRmeanRates[i,:],**{'color':color, 'marker':'.', 'linestyle':''})  
        #x,y,dy=fit_pol(MSNmeanRates,SNRmeanRates[i,:])
        #x,y,dy=fit_exp(MSNmeanRates,SNRmeanRates[i,:])
        #ax.plot(x,y,**{'color':color})  
   
    ax.plot([0,2.5],[SELECTION_THR]*len([0,2]),**{'color':[0.6,0.6,0.6],
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
    ax.text( 0.8, 0.20,'Thr' , transform=ax.transAxes, **{ 'color' : [0.6,0.6,0.6]})
    
    ax.set_ylabel('Firing rate SNr (spikes/s)') 
    ax.set_xlabel('Firing rate MSN (spikes/s)')
    ax.my_set_no_ticks( yticks=8, xticks = 6 )
    
    ax.set_xlim(misc.adjust_limit([0, 50]))
    ax.set_ylim(misc.adjust_limit([0,34]))
    
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,50])   
                                  
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['text.fontsize'], 
                 **{'color': color})
        
        
def plot_selection_vs_neurons_direct2(ax, MSNmeanRates, SNRmeanRates, title):
    colors=['k','k','k']   
    labels=[r'1 % bursting' , r'2 % bursting',  
            r'4 % bursting']
    linestyle=['-','--','-.']
    coords=[[0.4, 0.42], [ 0.02, 0.001], [0.15, 0.78]]
    

    for i, color, ls, lb in zip([0,1,3], colors, linestyle, labels):
        ax.plot(MSNmeanRates,SNRmeanRates[i,:],**{'color':color, 'label': lb, 'linestyle':ls})  
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
    ax.set_ylabel('Firing rate SNr (spikes/s)') 
    ax.set_xlabel('Burst firing rate $MSN_{D1}$ (spikes/s)')
    ax.my_set_no_ticks( yticks=8, xticks = 6 )
    ax.set_title(title, fontsize=12)
    ax.set_xlim(misc.adjust_limit([15, 50]))
    ax.set_ylim(misc.adjust_limit([0,34]))
    
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,50])   
    

def plot_selection_vs_neurons_direct_linear_fit_and_spread(ax, mods, slope):
    colors=['b','g','r']   
    labels=[r'First 100 ms' , r'250-350 ms',  
            r'Last 100 ms']
    coords=[[ 0.55, 0.15], [0.07, 0.47], [0.15, 0.78]]
    
    ps=[]
    for i, color in enumerate(colors):
        ax.plot(mods,slope[i,:],**{'color':color,  'linestyle':'-'})  
        #x,y,dy=fit_pol(MSNmeanRates,SNRmeanRates[i,:])    
        #ax.plot(x,y,**{'color':color})    
    '''
    ax.plot([0,2.5],[SELECTION_THR]*len([0,2]),**{'color':[0.6,0.6,0.6],
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
    
    ax.text( 0.8, 0.20,'Thr' , transform=ax.transAxes, **{ 'color' : [0.6,0.6,0.6]})
    '''
    ax.set_ylabel('$\Delta_{SNr}/ \Delta_{MSN_{D1}}$') 
    ax.set_xlabel('Percent bursting $MSN_{D1}$ (%)')
    ax.my_set_no_ticks( yticks=8, xticks = 6 )
    
    #ax.set_xlim(misc.adjust_limit([17, 50]))
    ax.set_ylim(misc.adjust_limit([-0.65,0]))
    
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,50])   
                                  
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=12, 
                 **{'color': color})

def plot_selection_vs_neurons_direct2_dy(ax, MSNmeanRates, SNRmeanRates):
    colors=['b','g','r']   
    labels=[r'$First 100 ms$' , r'250-350 ms$',  
            r'$Last 100 ms$']
    coords=[[0.4, 0.42], [ 0.02, 0.001], [0.15, 0.78]]
    
    
    for i, color in enumerate(colors):
   #     ax.plot(MSNmeanRates,SNRmeanRates[i,:],**{'color':color, 'marker':'.', 'linestyle':''})  
        x,y,dy=fit_pol(MSNmeanRates,SNRmeanRates[i,:])
        ax.plot(x,dy,**{'color':color})  
   
    ax.plot([0,2.5],[SELECTION_THR]*len([0,2]),**{'color':[0.6,0.6,0.6],
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
    ax.text( 0.8, 0.20,'Thr' , transform=ax.transAxes, **{ 'color' : [0.6,0.6,0.6]})
    
    ax.set_ylabel('Firing rate SNr (spikes/s)') 
    ax.set_xlabel('Firing rate MSN (spikes/s)')
    ax.my_set_no_ticks( yticks=8, xticks = 6 )
    
    ax.set_xlim(misc.adjust_limit([0, 50]))
    #ax.set_ylim(misc.adjust_limit([0,34]))
    
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,50])   
                                  
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['text.fontsize'], 
                 **{'color': color})
        
        
def plot_selection_vs_neurons_indirect(ax, MSNmeanRates, SNRmeanRates):
    colors=['r','c' ] 
    labels=[r'$\delta_{ref}^{GPe}$' , r'$\delta_{dep}^{GPe}$']
    coords=[[0.4, 0.42], [ 0.02, 0.001]]
    
    
    for i, color in enumerate(colors):
        ax.plot(MSNmeanRates,SNRmeanRates[i,:],**{'color':color})  

   
    ax.plot([0,2.5],[SELECTION_THR]*len([0,2]),**{'color':[0.6,0.6,0.6],
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
    ax.text( 0.8, 0.20,'Thr' , transform=ax.transAxes, **{ 'color' : [0.6,0.6,0.6]})
    
    ax.set_ylabel('Firing rate SNr (spikes/s)') 
    ax.set_xlabel('Firing rate MSN (spikes/s)')
    ax.my_set_no_ticks( yticks=8, xticks = 6 )
    
    #ax.set_xlim(misc.adjust_limit([0, 50]))
    ax.set_ylim(misc.adjust_limit([0,100]))
    
    #lines = ax.lines
    #for line in lines:
    #    misc.slice_line(line, xlim=[0,50])   
                                  
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=pylab.rcParams['text.fontsize'], 
                 **{'color': color})
      
def plot_selection_vs_neurons_indirect2(ax, MSNmeanRates, SNRmeanRates):
    colors=['b','g','r']   
    labels=[r'First 100 ms' , r'250-350 ms',  
            r'Last 100 ms']
    coords=[[0.4, 0.42], [ 0.02, 0.001], [0.15, 0.78]]
    
    
    for i, color in enumerate(colors):
        ax.plot(MSNmeanRates,SNRmeanRates[i,:],**{'color':color, 'marker':'.', 'linestyle':'--'})  
        #x,y,dy=fit_pol(MSNmeanRates,SNRmeanRates[i,:])    
        #ax.plot(x,y,**{'color':color})  
   
    ax.plot([0,2.5],[SELECTION_THR]*len([0,2]),**{'color':[0.6,0.6,0.6],
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
    ax.text( 0.8, 0.20,'Thr' , transform=ax.transAxes, **{ 'color' : [0.6,0.6,0.6]})
    
    ax.set_ylabel('Firing rate SNr (spikes/s)') 
    ax.set_xlabel('Firing rate $MSN_{D2}$ (spikes/s)')
    ax.my_set_no_ticks( yticks=8, xticks = 6 )
    
    ax.set_xlim(misc.adjust_limit([17, 50]))
    ax.set_ylim(misc.adjust_limit([34, 110]))
    
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,50])   
                                  
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=12, 
                 **{'color': color})

def plot_selection_vs_neurons_indirect2(ax, MSNmeanRates, SNRmeanRates, title):
    colors=['k','k','k']   
    labels=[r'1 % bursting' , r'3 % bursting',  
            r'6 % bursting']
    linestyle=['-','--','-.']
    coords=[[0.4, 0.42], [ 0.02, 0.001], [0.15, 0.78]]
    
    
    for i, color, ls, lb in zip([0,2,5], colors, linestyle, labels):
        ax.plot(MSNmeanRates,SNRmeanRates[i,:],**{'color':color, 'label': lb, 'linestyle':ls})    
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
    ax.set_ylabel('Firing rate SNr (spikes/s)') 
    ax.set_xlabel('Burst firing rate $MSN_{D2}$ (spikes/s)')
    ax.my_set_no_ticks( yticks=8, xticks = 6 )
    ax.set_title(title, fontsize=12)
    ax.set_xlim(misc.adjust_limit([15, 50]))
    ax.set_ylim(misc.adjust_limit([34,105]))
    
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,50])   

def plot_selection_vs_neurons_indirect_linear_fit_and_spread(ax, mods, slope):
    colors=['b','g','r']   
    labels=[r'First 100 ms' , r'250-350 ms',  
            r'Last 100 ms']
    coords=[[ 0.55, 0.15], [0.07, 0.47], [0.15, 0.78]]
    
    ps=[]
    for i, color in enumerate(colors):
        ax.plot(mods,slope[i,:],**{'color':color,  'linestyle':'-'})  
        #x,y,dy=fit_pol(MSNmeanRates,SNRmeanRates[i,:])    
        #ax.plot(x,y,**{'color':color})    
    '''
    ax.plot([0,2.5],[SELECTION_THR]*len([0,2]),**{'color':[0.6,0.6,0.6],
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
    
    ax.text( 0.8, 0.20,'Thr' , transform=ax.transAxes, **{ 'color' : [0.6,0.6,0.6]})
    '''
    ax.set_ylabel('$\Delta_{SNr}/ \Delta_{MSN_{D2}}$') 
    ax.set_xlabel('Percent bursting $MSN_{D2}$ (%)')
    ax.my_set_no_ticks( yticks=8, xticks = 6 )
    
    #ax.set_xlim(misc.adjust_limit([15, 50]))
    ax.set_ylim(misc.adjust_limit([0, 1.3]))
    
    lines = ax.lines
    for line in lines:
        misc.slice_line(line, xlim=[0,50])   
                                  
    for label, coord, color in zip(labels,coords,colors):
        ax.text( coord[0], coord[1], label, transform=ax.transAxes, 
                 fontsize=12, 
                 **{'color': color})
        
def plot_selection_vs_neurons_hyperdirect2(ax, MSNmeanRates, SNRmeanRates,linestyle='-'):
    colors=['b','g','r']   
    labels=[r'First 100 ms' , r'250-350 ms',  
            r'Last 100 ms']
    coords=[[0.55, 0.05],  [0.35, 0.7], [ 0.1, 0.45]]
    
    
    for i, color in enumerate(colors):
        ax.plot(MSNmeanRates, SNRmeanRates[i,:],**{'color':color, 'marker':'.', 'linestyle':linestyle})  
        #x,y,dy=fit_pol(MSNmeanRates,SNRmeanRates[i,:])    
        #ax.plot(x,y,**{'color':color})  
    '''
    ax.plot([0,2.5],[SELECTION_THR]*len([0,2]),**{'color':[0.6,0.6,0.6],
                                                        'label':'', 'linewidth':1,
                                                        'linestyle':'--'})
    ax.text( 0.8, 0.20,'Thr' , transform=ax.transAxes, **{ 'color' : [0.6,0.6,0.6]})
    '''
        
    line1=ax.plot(1,1,'-k')
    line2=ax.plot(1,1,'--k')
    leg=ax.legend([line1, line2],[r'Static STN $\rightarrow$ SNr', r'Depressing STN $\rightarrow$ SNr'], loc='best')
    frame  = leg.get_frame() 
    #frame.set_visible(False) 
    pylab.setp(frame, edgecolor='w') 
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    pylab.setp(ltext, fontsize=10., backgroundcolor='w') 
        
    ax.set_ylabel('Firing rate SNr (spikes/s)') 
    ax.set_xlabel('Firing rate STN (spikes/s)')
    ax.my_set_no_ticks( yticks=8, xticks = 6 )
    
    ax.set_xlim(misc.adjust_limit([10, 36]))
    ax.set_ylim(misc.adjust_limit([25, 80]))
    
    lines = ax.lines
    #for line in lines:
    #    misc.slice_line(line, xlim=[0,50])   
                                  
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


def get_interval_rates(interval, mean_rates):
    interval_rates=[]
    for mr in mean_rates:
        mr=numpy.array(mr)
        interval_rates.append(numpy.mean(mr[:,interval[0]:interval[1]], axis=1))
    
    
    return numpy.array(interval_rates)

def fit_pol(x,y):
    p=numpy.polyfit(x, y, 1, rcond=None, full=False)

    #f=lambda x: p[0]*x**3+p[1]*x**2+p[2]*x+p[3]
    #df=lambda x: 3*p[0]*x**2+2*p[1]*x+p[2]
    f=lambda x: p[0]*x+p[1]
    df=lambda x: p[0]    
    xr=numpy.linspace(min(x),max(x), 50)
    
    return xr, f(xr), df(xr), p

def fit_exp(x,y):
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit
    
    f=lambda x:(sum(numpy.exp(-dist[source]/x))-n)**2
    fmin(f[source + '_' + target],numpy.array([0.1]))    
    
    def func(x, a, b, c):
        return a * np.exp(-b * x) + c
    def dfunc(x,a,b,c):
        return -b*a * np.exp(-b * x) 
    popt, pcov = curve_fit(func, x, y)
    xr=numpy.linspace(min(x),max(x), 50)
    return xr, func(xr, *popt), dfunc(xr, *popt)



def simulate_direct(load, N_MSN, save_at, threads, points=[550, 650, 950]):
            


    params_msn_d1={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, 20, 0.1],
                'mod_times':[1,1000, 1000+500],  'n_mod':0}    
    params_msn_d2={'base_rates':[0.1], 'base_times':[1], 'mod_rates': [0.1, 20, 0.1],
                'mod_times':[1,1000, 1000+500],  'n_mod':0, 'focus':False,'skip':1} 
    params_stn={'rate':300., 'mod':False,'mod_rate':0., 'mod_times':[1000., 1000.+500.]} 
    
    
    model_params={'neurons':{'MSN_D1':{'n':N_MSN},
                             'MSN_D2':{'n':N_MSN}}}
    
    synapse_models=['MSN_SNR_gaba_p1', 'GPE_SNR_gaba_p']    
    
    sim_time=2000.
    
    resolution=2
    n_exp=2
    
    burst_rate=numpy.linspace(15 ,50 ,resolution)
    proportions  =numpy.linspace(0.01,0.07,resolution) #arange(1,7,1)*150. 
    mods=proportions*N_MSN
        
    
    pop_raw_rates=[]
    pop_conv_rates=[]
    mean_pop_conv_rates=[]
    std_pop_conv_rates=[]
    

    inputs=[]
    i=0
    if not load:
        for m in mods:
            for r in burst_rate:
                tmp_rates=[]                
                for e in range(n_exp):
                    seed=i
                    
                    
                    rates =simulate_network_direct_onoff_vs_rate(m, r, 
                                                           params_msn_d1, 
                                                           params_msn_d2, 
                                                           params_stn,
                                synapse_models, {'SNR':280, 'STN':0,'GPE':20},
                                sim_time=sim_time, seed=seed, 
                               threads=threads, start_rec=500.,model_params=model_params)
                    
                    tmp_rates.append(list(rates))
                    
                    i+=1
                
                tmp_rates=numpy.array(tmp_rates)
                pop_raw_rates.append(tmp_rates)
                pop_conv_rates.append(misc.convolve(tmp_rates, 100, 'rectangle',single=False))
                
                inputs.append((m,r))
                
                
                
                mean_pop_conv_rates.append(numpy.mean(pop_conv_rates[-1],axis=0))
                std_pop_conv_rates.append(numpy.std(pop_conv_rates[-1],axis=0))
                
#                mean=mean_pop_conv_rates[-1]
#                std=std_pop_conv_rates[-1]
#            
#                for i in range(n_exp):
#                    pylab.plot(pop_conv_rates[-1][i,:])
#                pylab.plot(mean,'k')
#                pylab.plot(mean-std,'--k')
#                pylab.plot(mean+std,'--k')
#                pylab.show()
                
                
                misc.pickle_save([pop_raw_rates, 
                                  pop_conv_rates, 
                                  mean_pop_conv_rates, 
                                  std_pop_conv_rates,
                                  inputs], save_result_at)
    else: 
        pop_raw_rates, pop_conv_rates, mean_pop_conv_rates, std_pop_conv_rates, inputs=misc.pickle_load(save_result_at)


    data={}
    for p in points:
        
        x=[]
        y=[]
        z_mean=[]
        z_std=[]
        for i, m, s in zip(inputs,mean_pop_conv_rates,std_pop_conv_rates):
            x.append(i[0])
            y.append(i[1])
            z_mean.append(m[p])
            z_std.append(s[p])
        
        data[p]=numpy.array([x,y,z_mean, z_std]) 
            
        spl=numpy.split(data[p],2, axis=1)
#        
#        for a in spl:
#            pylab.plot(a[1,:],a[2,:])
#            pylab.plot(a[1,:],a[2,:]+a[3,:],'--k')
#            pylab.plot(a[1,:],a[2,:]-a[3,:],'--k')
#    
#    pylab.show()
    

    return data

N_MSN=1500            
save_result_at=OUTPUT_PATH+'/simulate_network_direct_onoff_vs_rate_direct'+str(N_MSN)+'.plk'   
data=simulate_direct(1, N_MSN, save_result_at, 4 )


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

interval_rates_direct=[]
slopes_direct=[]
for start, stop  in [[510,610], [660,760], [910,1010]]:
    ird=get_interval_rates([start,stop],mean_rates)
    interval_rates_direct.append(ird)
    sl=[]
    for d in ird:
        x,y, dy, p=fit_pol(freq_burst1,d) 
        sl.append(p[0])
    slopes_direct.append(sl)

ax=ax_list[1]
plot_selection_vs_neurons_direct2(ax, data[550], title='First 100 ms of the burst')

ax=ax_list[2]
plot_selection_vs_neurons_direct2(ax, data[550], title='First 100 ms of the burst')

ax=ax_list[3]
plot_selection_vs_neurons_direct_linear_fit_and_spread(ax, mods_direct_percent, slopes_direct)
#plot_selection_vs_neurons_direct2(ax, freq_burst1, interval_rates_direct[2])
#plot_selection_vs_neurons_direct(ax, freq_burst1, interval_rates_direct_3)

pylab.show()
slopes_direct=numpy.array(slopes_direct)   
interval_rates_direct=numpy.array(interval_rates_direct)



sim_time=2000.
freq_burst2=numpy.linspace(15,50,6)
seed=range(len(freq_burst2))
mod=numpy.ones(len(freq_burst2))*1000*N_MSN/15000.

model_params= {'conns':{ 'MSN_D2_GPE':{ 'lines':False}},
               'neurons':{'MSN_D1':{'n':N_MSN},
                          'MSN_D2':{'n':N_MSN}}} 

save_result_at=OUTPUT_PATH+'/simulate_network_indirect_onoff_vs_rate'+str(N_MSN)+'.plk'
synapse_models=['MSN_SNR_gaba_p1','GPE_SNR_gaba_p']

#mods=[150, 300, 600, 900]
mods=numpy.arange(1,11,1)*150.
mods_indirect_percent=N_MSN/15000.*numpy.array(mods)/N_MSN*100.
if 0:
    snr_mr, gpe_mr=[],[]
    for m in mods:
       mod_msn=numpy.ones(len(freq_burst2))*N_MSN/15000.*m 
       snr_m, gpe_m=simulate_network_indirect_onoff_vs_rate(mod_msn, freq_burst2, 
                                                       params_msn_d1, 
                                                       params_msn_d2, 
                                                       params_stn,
                           synapse_models, sim_time=sim_time, seed=seed, 
                           threads=4, start_rec=500.,model_params=model_params)
       snr_mr.append(snr_m)
       gpe_mr.append(gpe_m)
    misc.pickle_save([snr_mr, gpe_mr], save_result_at)
else: 
    snr_mr, gpe_mr=misc.pickle_load(save_result_at)

interval_rates_indirect=[]
slopes_indirect=[]
for start, stop  in [[510,610], [660,760], [910,1010]]:
    ird=get_interval_rates([start,stop],snr_mr)
    interval_rates_indirect.append(ird)
    sl=[]
    for d in ird:
        x,y, dy, p=fit_pol(freq_burst1,d) 
        sl.append(p[0])
    slopes_indirect.append(sl)
    
slopes_indirect=numpy.array(slopes_indirect)   
interval_rates_indirect=numpy.array(interval_rates_indirect)


interval_rates_indirect_1=get_interval_rates([510,610],snr_mr)
interval_rates_indirect_2=get_interval_rates([660,760],snr_mr)
interval_rates_indirect_3=get_interval_rates([910,1010],snr_mr)
interval_rates_indirect_12=get_interval_rates([510,600],gpe_mr)
interval_rates_indirect_22=get_interval_rates([600,700],gpe_mr)
interval_rates_indirect_32=get_interval_rates([900,1000],gpe_mr)  
        
#pylab.plot(numpy.mean(snr_mr[0],axis=0))

sim_time=2000.
freq_burst3=numpy.linspace(1,1000,5)
seed=range(len(freq_burst3))
mod=numpy.ones(len(freq_burst3))*1000*N_MSN/15000.

model_params= {'conns':{ 'MSN_D2_GPE':{ 'lines':False}},
               'neurons':{'MSN_D1':{'n':N_MSN},
                          'MSN_D2':{'n':N_MSN}}} 

synapse_models_gpe=['GPE_SNR_gaba_s_ref', 'GPE_SNR_gaba_p']
save_result_at=OUTPUT_PATH+'/simulate_network_hyper_direct_onoff_vs_rate'+str(N_MSN)+'.plk'
if 0:
    snr_mr3, gpe_mr3, stn_mr3=[],[],[]
    for syn in synapse_models_gpe:
       synapse_models=['MSN_SNR_gaba_p1',syn]
       snr_m3, gpe_m3, stn_m3=simulate_network_hyper_direct_onoff_vs_rate(freq_burst3, 
                                                       params_msn_d1, 
                                                       params_msn_d2, 
                                                       params_stn,
                           synapse_models, sim_time=sim_time, seed=seed, 
                           threads=4, start_rec=500.,model_params=model_params)
       snr_mr3.append(snr_m3)
       gpe_mr3.append(gpe_m3)
       stn_mr3.append(stn_m3)
    misc.pickle_save([snr_mr3, gpe_mr3, stn_mr3], save_result_at)
else: 
    snr_mr3, gpe_mr3, stn_mr3=misc.pickle_load(save_result_at)

interval_rates_hyperdirect_1=get_interval_rates([510,610],snr_mr3)
interval_rates_hyperdirect_2=get_interval_rates([660,760],snr_mr3)
interval_rates_hyperdirect_3=get_interval_rates([910,1010],snr_mr3)
interval_rates_hyperdirect_12=get_interval_rates([500,600],gpe_mr3)
interval_rates_hyperdirect_22=get_interval_rates([600,700],gpe_mr3)
interval_rates_hyperdirect_32=get_interval_rates([900,1000],gpe_mr3)  
interval_rates_hyperdirect_13=get_interval_rates([500,600],stn_mr3)
interval_rates_hyperdirect_23=get_interval_rates([600,700],stn_mr3)
interval_rates_hyperdirect_33=get_interval_rates([900,1000],stn_mr3)  
x_hyperdirect=get_interval_rates([510,1010],stn_mr3) 


model_params= {'conns':{ 'MSN_D2_GPE':{ 'lines':False},
                        'STN_SNR':{'syn':'STN_SNR_ampa_p3'}},  
               'neurons':{'MSN_D1':{'n':N_MSN},
                          'MSN_D2':{'n':N_MSN}}}

synapse_models_gpe=['GPE_SNR_gaba_s_ref', 'GPE_SNR_gaba_p']
save_result_at=OUTPUT_PATH+'/simulate_network_hyper_direct2_onoff_vs_rate.plk'
if 0:
    snr_mr, gpe_mr, stn_mr=[],[],[]
    for syn in synapse_models_gpe:
       synapse_models=['MSN_SNR_gaba_p1',syn]
       snr_m, gpe_m, stn_m=simulate_network_hyper_direct_onoff_vs_rate(freq_burst3, 
                                                       params_msn_d1, 
                                                       params_msn_d2, 
                                                       params_stn,
                           synapse_models, sim_time=sim_time, seed=seed, 
                           threads=4, start_rec=500.,model_params=model_params)
       snr_mr.append(snr_m)
       gpe_mr.append(gpe_m)
       stn_mr.append(stn_m)
    misc.pickle_save([snr_mr, gpe_mr, stn_mr], save_result_at)
else: 
    snr_mr, gpe_mr, stn_mr=misc.pickle_load(save_result_at)

interval2_rates_hyperdirect_1=get_interval_rates([510,610],snr_mr)
interval2_rates_hyperdirect_2=get_interval_rates([660,760],snr_mr)
interval2_rates_hyperdirect_3=get_interval_rates([910,1010],snr_mr)
interval2_rates_hyperdirect_12=get_interval_rates([500,600],gpe_mr)
interval2_rates_hyperdirect_22=get_interval_rates([600,700],gpe_mr)
interval2_rates_hyperdirect_32=get_interval_rates([900,1000],gpe_mr)  
interval2_rates_hyperdirect_13=get_interval_rates([500,600],stn_mr)
interval2_rates_hyperdirect_23=get_interval_rates([600,700],stn_mr)
interval2_rates_hyperdirect_33=get_interval_rates([900,1000],stn_mr)  
x_hyperdirect2=get_interval_rates([510,1010],stn_mr) 







interval_plast_indirect=numpy.array([interval_rates_indirect_1[1,:], 
                                     interval_rates_indirect_2[1,:], 
                                     interval_rates_indirect_3[1,:]])

interval_plast_hyperdirect=numpy.array([interval_rates_hyperdirect_1[1,:], 
                                     interval_rates_hyperdirect_2[1,:], 
                                     interval_rates_hyperdirect_3[1,:]])

interval2_plast_hyperdirect=numpy.array([interval2_rates_hyperdirect_1[1,:], 
                                     interval2_rates_hyperdirect_2[1,:], 
                                     interval2_rates_hyperdirect_3[1,:]])
ax=ax_list[1]
plot_selection_vs_neurons_direct2(ax, freq_burst1, interval_rates_direct[0], title='First 100 ms of the burst')
#plot_selection_vs_neurons_direct(ax, freq_burst1, interval_rates_direct_1)
ax=ax_list[2]
plot_selection_vs_neurons_direct2(ax, freq_burst1, interval_rates_direct[2],title='Last 100 ms of the burst')
#plot_selection_vs_neurons_direct2_dy(ax, freq_burst1, interval_plast)
#plot_selection_vs_neurons_direct(ax, freq_burst1, interval_rates_direct_2)
ax=ax_list[3]
plot_selection_vs_neurons_direct_linear_fit_and_spread(ax, mods_direct_percent, slopes_direct)
#plot_selection_vs_neurons_direct2(ax, freq_burst1, interval_rates_direct[2])
#plot_selection_vs_neurons_direct(ax, freq_burst1, interval_rates_direct_3)
ax=ax_list[4]
plot_selection_vs_neurons_indirect2(ax, freq_burst1, interval_rates_indirect[0], title='First 100 ms of the burst')
#plot_selection_vs_neurons_indirect2(ax, freq_burst2, interval_plast_indirect)
#plot_selection_vs_neurons_indirect(ax, freq_burst2, interval_rates_indirect_1)
#plot_selection_vs_neurons_indirect(ax, freq_burst2, interval_rates_indirect_12)
ax=ax_list[5]
plot_selection_vs_neurons_indirect2(ax, freq_burst1, interval_rates_indirect[2], title='Last 100 ms of the burst')

#plot_selection_vs_neurons_indirect(ax, freq_burst2, interval_rates_indirect_2)
#plot_selection_vs_neurons_indirect(ax, freq_burst2, interval_rates_indirect_22)
ax=ax_list[6]
plot_selection_vs_neurons_indirect_linear_fit_and_spread(ax, mods_indirect_percent, slopes_indirect)
#plot_selection_vs_neurons_indirect(ax, freq_burst2, interval_rates_indirect_3)
#plot_selection_vs_neurons_indirect(ax, freq_burst2, interval_rates_indirect_32)


ax=ax_list[7]
plot_selection_vs_neurons_hyperdirect2(ax, x_hyperdirect[0], interval_plast_hyperdirect)
#plot_selection_vs_neurons_indirect(ax, freq_burst3, interval_rates_hyperdirect_1)
#plot_selection_vs_neurons_indirect(ax, freq_burst3, interval_rates_hyperdirect_12)
#plot_selection_vs_neurons_indirect(ax, freq_burst3, interval_rates_hyperdirect_13)

#ax=ax_list[8]
plot_selection_vs_neurons_hyperdirect2(ax, x_hyperdirect2[0], interval2_plast_hyperdirect, linestyle='--')
#plot_selection_vs_neurons_indirect(ax, freq_burst3, interval_rates_hyperdirect_2)
#plot_selection_vs_neurons_indirect(ax, freq_burst3, interval_rates_hyperdirect_22)
#plot_selection_vs_neurons_indirect(ax, freq_burst3, interval_rates_hyperdirect_23)
ax=ax_list[9]
plot_selection_vs_neurons_indirect(ax, freq_burst3, interval_rates_hyperdirect_3)
plot_selection_vs_neurons_indirect(ax, freq_burst3, interval_rates_hyperdirect_32)
plot_selection_vs_neurons_indirect(ax, freq_burst3, interval_rates_hyperdirect_33)

ax=ax_list[0]
#plot_text(ax, info_string=s)

pylab.show()

# dpi does not matter since svg and pdf are both vectorbased
fig.savefig( picture_dir + '/' + SNAME  + '.svg', format = 'svg') 
fig.savefig( picture_dir + '/' + SNAME  + '.pdf', format = 'pdf')