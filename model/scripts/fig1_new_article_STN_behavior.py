#! Imports
import numpy
import pylab
import os
import sys


# Add directories to python path
sys.path.append(os.getcwd())                            
parent_dir='/'.join(os.getcwd().split('/')[0:-1])       
                   
model_dir=   '/'.join(os.getcwd().split('/')[0:-1])        
code_dir=  '/'.join(os.getcwd().split('/')[0:-2]) 

model_name=os.getcwd().split('/')[-2]
picture_dir='/'.join(os.getcwd().split('/')[0:-3]) + '/pictures/'+model_name 
                
sys.path.append(model_dir) 
sys.path.append(code_dir+'/nest_toolbox') 
spath  = os.getcwd()+'/output/' + sys.argv[0].split('/')[-1].split('.')[0]

from model_params import models                               # Then import models     
from src import misc, my_nest, my_signals, plot_settings
from src.my_population import MyGroup 
from src.my_axes import MyAxes 

NEURON_MODELS=['STN_75_aeif']

# 270 in curr gives 62.5 Hz in fring frequency
I_E=0.0
I_E2=21

def plot_example(ax, STN):
    STN.signals['V_m'].plot(display=ax, kwargs={'color':'k'})
    
    ax.set_xlim(misc.adjust_limit([0,2000]))
    ax.set_ylim(misc.adjust_limit([-80,30]))
    
    ax.my_set_no_ticks( yticks=8, xticks = 6 )  

def plot_example_inh_current(ax, STN):
    STN.signals['V_m'].plot(display=ax, kwargs={'color':'k'})
    ax.set_xlim(misc.adjust_limit([0,1000]))
    ax.set_ylim(misc.adjust_limit([-80,30]))
    ax.my_set_no_ticks( yticks=8, xticks = 6 )  
    ax.set_ylabel('STN potential (mV)')

def plot_IV(ax, current, voltage):
    color='r'
    #lower_limit=-65
    #upper_limit=-80
    #current=current[voltage<lower_limit]
    #voltage=voltage[voltage<lower_limit]
    #current=current[voltage>upper_limit]
    #voltage=voltage[voltage>upper_limit]
    current=current-current[-1]
    ax.plot(current, voltage, **{'color':color})
    
    line=ax.lines

    ax.my_set_no_ticks( yticks=8, xticks = 6 )  
    ax.set_xlabel('Current (pA)') 
    ax.set_ylabel('Potential (mV)')

    df=my_nest.GetDefaults(NEURON_MODELS[0]) 
    
    speed=numpy.diff(voltage)/numpy.diff(current/1000.) 
    mean_resistance=numpy.mean(speed)
    mean_membrane_time_constant=mean_resistance*df['C_m']*0.001
    
    


    point=-15
    ax.text( 0.50, 0.75,str(round(speed[point]))+r' M$\Omega$', backgroundcolor='w',
             transform=ax.transAxes,fontsize=16,  **{ 'color' : 'k'}) 
    ax.plot(current[point], voltage[point],color=color, marker='.',markersize=10)
    x=[current[point], current[point]-100]
    y=[voltage[point], voltage[point]-(100)*speed[point]/1000.]
    ax.plot(x,y,color=color, linestyle='--')
    
    ax.set_xlim(misc.adjust_limit([-200,0]))
    ax.set_ylim(misc.adjust_limit([-85,-55]))
    pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend
    
def plot_IF(ax, I_vec, fIsi, mIsi, lIsi):
    color='r'
    ax.plot(I_vec, 1000./fIsi, **{ 'color':color,
                                  'linestyle':'--'})
    ax.plot(I_vec, 1000./lIsi, **{'color':color})

    print 1000./lIsi[I_vec<100]
    print I_vec[I_vec<100]
    ax.my_set_no_ticks( yticks=8, xticks = 6 )  
    ax.set_xlabel('Current (pA)') 
    ax.set_ylabel('Rate (spikes/s)')
    ax.set_xlim(misc.adjust_limit([-100,300]))
    ax.set_ylim(misc.adjust_limit([0,200]))
    
    
    
    line1=ax.plot(1,1,'--k')
    line2=ax.plot(1,1,'-k')
    pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend
    leg=ax.legend([line1, line2],['First ISI','Steady-state ISI'], loc='upper left')
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    llines = leg.get_lines()  # all the lines.Line2D instance in the legend
    frame  = leg.get_frame() 
    pylab.setp(ltext, fontsize=14) 
    frame.set_visible(False)
    
    
    speed=numpy.diff(1000./lIsi)/numpy.diff(I_vec) 
    #speed=speed[speed>0]
    #ax.text( 0.3, 0.02, 'Mean slope steady-state\nISI '+str(round(numpy.mean(speed),2))+' Hz/pA', backgroundcolor='w',
    #         transform=ax.transAxes,fontsize=10,  **{ 'color' : color }) 
    point=numpy.nonzero(I_vec==81)[0][0]
    ax.text( 0.6, 0.2, 'STN', backgroundcolor='w',
             transform=ax.transAxes,fontsize=16,  **{ 'color' : color }) 
    ax.text( 0.5, 0.02, str(round(speed[point],2))+' Hz/pA', backgroundcolor='w',
             transform=ax.transAxes,fontsize=16,  **{ 'color' : 'k' }) 
    ax.plot(I_vec[point], 1000./lIsi[point],color=color, marker='.',markersize=10)
 
    ax.set_xlim(misc.adjust_limit([-100,300]))
    ax.set_ylim(misc.adjust_limit([0,260]))


def plot_example_rebound_spike(ax, STN):
    STN.signals['V_m'].plot(display=ax, kwargs={'color':'k'})   
    ax.set_xlim(misc.adjust_limit([0,2500]))
    ax.set_ylim(misc.adjust_limit([-100,30]))
    line=ax.lines
    misc.slice_line(line[0], xlim=[0,2500])

    ax.my_set_no_ticks( yticks=3, xticks = 6 )  
    ax.set_ylabel('STN potential (mV)')

def plot_ahp(ax,  I_vec, delays ):

    ax.plot(I_vec, numpy.array(delays)/1000.0, **{'label':'Steady-state ISI', 'color':'k','marker':'s'})
    ax.my_set_no_ticks( yticks=3, xticks = 6 )  
    ax.set_xlabel('Current (pA)') 
    ax.set_ylabel('Duration of AHP (s)')
    ax.set_xlim(misc.adjust_limit([0,350]))
    #ax.set_ylim(misc.adjust_limit([0.2,0.8]))

def plot_text(ax, infoString=''):
    
    my_nest.ResetKernel()
    model_list, model_dict=models()
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    
    STN = MyGroup( NEURON_MODELS[0], 1, mm_dt = 0.1)
    statusSTN = my_nest.GetStatus( STN[:] )[0]
    
    tb = ''     
    tb = tb + infoString
    
    tb = tb + '\n'
    for key, val in statusSTN.iteritems():
        if not key in ['vp','state','t_spike','local','parent','Delta_T',
                       'tau_minus_triplet','address', 't_ref','thread',
                       'frozen','archiver_length', 'global_id', 'local_id',
                       'recordables','receptor_types']:
            tb = tb + ' %s %5s %3s \n' % ( key+':', str(val), 
                                   '--' )

    ax.text( 0.85, 0.5, tb , fontsize= font_size_text,
             horizontalalignment='right',
             verticalalignment='center',
             transform=ax.transAxes,     # to define coordinates in right scale
             **{ 'fontname' : 'monospace' })                           
    
    ax.my_remove_axis( xaxis=True, yaxis=True )
    ax.my_remove_spine(left=True,  bottom=True, right=True, top=True)

def simulate_example(I_e):
    
    simTime  = 3000.  # ms
    
    my_nest.ResetKernel()
    model_list, model_dict=models()
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    
    I_e0=my_nest.GetDefaults(NEURON_MODELS[0])['I_e']
    STN = MyGroup( NEURON_MODELS[0], 1, sd=True, mm=True, mm_dt = 1., 
                   params={'I_e':I_e+I_e0} )
    '''
    scg = my_nest.Create( 'step_current_generator',n=1 )  
    rec=my_nest.GetStatus(STN[:])[0]['receptor_types']
    

    my_nest.SetStatus(scg, {'amplitude_times':[280.,1700.],
                                'amplitude_values':[100.0,0.]})
    my_nest.Connect( scg, STN.ids, 
                         params = { 'receptor_type' : rec['CURR'] } )
    
    '''
    my_nest.MySimulate(simTime)
    STN.get_signal( 'v','V_m', stop=simTime ) # retrieve signal
    
    STN.get_signal( 's') # retrieve signal
    meanRate=round(STN.signals['spikes'].mean_rate(0,1000),1)
    print STN.signals['spikes'].isi()
    STN.signals['V_m'].my_set_spike_peak( 21, spkSignal= STN.signals['spikes'] )
    
    s='\n'
    s = s + 'Example:\n'
    s = s + ' %s %5s %3s %s %5s %3s \n' % ( 'Mean rate:', meanRate,  'Hz', 
                                            'I_e', I_e,'pA' )

    infoString=s
    
    return STN, infoString

def simulate_example_inh_current(I_vec):
    
    simTime  = 1000.  # ms
    my_nest.ResetKernel()
    model_list, model_dict=models()
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    df=my_nest.GetDefaults(NEURON_MODELS[0])
    n=len(I_vec)
    
    STN = MyGroup( NEURON_MODELS[0], n, sd=True,  mm=True, mm_dt = 1.0,
                   record_from=['V_m','u'] )
    I_e0=my_nest.GetStatus(STN[:])[0]['I_e']
    my_nest.SetStatus(STN[:], params={'I_e':I_e0+I_E+50}) # Set I_e
    
    I_e = my_nest.GetStatus(STN.ids,'I_e')[0]
    scg = my_nest.Create( 'step_current_generator',n=n )  
    rec=my_nest.GetStatus(STN[:])[0]['receptor_types']
    
    for source, target, I in zip(scg, STN[:], I_vec):
        my_nest.SetStatus([source], {'amplitude_times':[280.,700.],
                                'amplitude_values':[float(I),0.]})
        my_nest.Connect( [source], [target], 
                         params = { 'receptor_type' : rec['CURR'] } )
    
    
    my_nest.MySimulate(simTime)
    STN.get_signal( 'v','V_m', stop=simTime ) # retrieve signal   
    STN.get_signal( 's') # retrieve signal
    STN.signals['V_m'].my_set_spike_peak( 21, spkSignal= STN.signals['spikes'] )
    
    e = my_nest.GetStatus(STN.mm)[0]['events']    # get events 
    #pylab.plot(e['u'])
    #pylab.show()
    meanRate=round(STN.signals['spikes'].mean_rate(0,500),1)

    s='\n'
    s =s + 'Example inhibitory current:\n'
    s = s + ' %s %5s %3s %s %5s %3s \n' % ( 'Mean rate:', meanRate,  'Hz', 
                                            'I_e', I_e,'pA' )
    s = s + 'Steps:\n'
    s = s + ' %5s %3s \n' % ( I_vec,  'pA' )
    infoString=s
    
    return STN, infoString

def simulate_IV(I_vec):
    
    I_e = 0.+I_E
    
    my_nest.ResetKernel()
    model_list, model_dict=models()
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    
    STN = MyGroup( NEURON_MODELS[0], 1, sd=True, mm=True, mm_dt = 0.1)
    
    I_e0=my_nest.GetStatus(STN[:])[0]['I_e']
    my_nest.SetStatus(STN[:], params={'I_e':I_e0+I_E}) # Set I_e 
    I_e = my_nest.GetStatus(STN.ids,'I_e')[0]
    
    I_vec, voltage = STN.IV_I_clamp(I_vec)   
    
    #current=current
    speed=numpy.diff(voltage)/numpy.diff(I_vec)*1000.
    speed=speed[speed>0]
    s='\n'
    s =s + 'IV:\n'
    s = s + ' %s %5s %3s \n' % ( 'I_e:', I_e,  'pA' )
    '''
    s = s + ' %s %4s %s %4s %s %4s\n' % ( 'Speed (mV/pA=MOhm), min:', 
                                          str(min(speed))[0:4],  
                                          'max',str(max(speed))[0:4],
                                          'mean', 
                                          str(sum(speed)/len(speed))[0:4])
    '''#
    infoString=s

    I_vec=numpy.array(I_vec)
    voltage=numpy.array(voltage)
    
    return I_vec, voltage, infoString

def simulate_IF(I_vec):
    
    tStim = 700+1300
    
    my_nest.ResetKernel()
    model_list, model_dict=models()
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    
    STN = MyGroup( NEURON_MODELS[0], 1, mm=True, sd=True, mm_dt = 0.1 )
    
    I_e0=my_nest.GetStatus(STN[:])[0]['I_e']
    my_nest.SetStatus(STN[:], params={'I_e':I_e0+I_E}) # Set I_e 
    I_e = my_nest.GetStatus(STN.ids,'I_e')[0]

    I_vec, fIsi, mIsi, lIsi = STN.IF(I_vec, tStim=tStim)   
    
    #ax=pylab.subplot(111)
    #for isi in isi_list:
    #    if not isi[0]==1000000.:
    #        ax.plot(1000/isi[:])
    #pylab.show()
    
    speed_f=numpy.diff(1000.0/fIsi)/numpy.diff(I_vec)
    speed_l=numpy.diff(1000.0/lIsi)/numpy.diff(I_vec)
    speed_f=speed_f[speed_f>0]
    speed_l=speed_l[speed_l>0]
    s='\n'
    s =s + 'IF:\n'
    s = s + ' %s %5s %3s \n' % ( 'First to Last ISI:', tStim,  'ms' )
    s = s + ' %s %5s %3s \n' % ( 'Added I_e:', I_e,  'pA' )
    s = s + ' %s %4s %s %4s %s %4s\n' % ( 'Speed first ((Hz/pA), min:', 
                                          str(min(speed_f))[0:4],  
                                          'max',str(max(speed_f))[0:4],
                                          'mean', 
                                          str(sum(speed_f)/len(speed_f))[0:4])
    s = s + ' %s %4s %s %4s %s %4s\n' % ( 'Speed last (Hz/pA), min:', 
                                          str(min(speed_l))[0:4],  
                                          'max',str(max(speed_l))[0:4],
                                          'mean', 
                                          str(sum(speed_l)/len(speed_l))[0:4])
    infoString=s

    return I_vec, fIsi, mIsi, lIsi, infoString

def simulate_example_rebound_spike(I_vec, time):
    
    simTime  = 3000.  # ms
    my_nest.ResetKernel()
    model_list, model_dict=models()
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    
    n=len(I_vec)
    
    STN = MyGroup( NEURON_MODELS[0], n, sd=True,  mm=True, mm_dt = 1.0 )
    I_e0=my_nest.GetStatus(STN[:])[0]['I_e']
    #my_nest.SetStatus(STN[:], params={'I_e':-10.}) # Set I_e
    my_nest.SetStatus(STN[:], params={'I_e':.5}) # Set I_e
    I_e = my_nest.GetStatus(STN.ids,'I_e')[0]
    scg = my_nest.Create( 'step_current_generator',n=n )  
    rec=my_nest.GetStatus(STN[:])[0]['receptor_types']
    
    for source, target, I in zip(scg, STN[:], I_vec):
        my_nest.SetStatus([source], {'amplitude_times':[500.,time+500.],
                                     'amplitude_values':[float(I),0.]})
        my_nest.Connect( [source], [target], 
                         params = { 'receptor_type' : rec['CURR'] } )
    
    
    my_nest.MySimulate(simTime)
    STN.get_signal( 'v','V_m', stop=simTime ) # retrieve signal
    STN.get_signal( 's') # retrieve signal
    STN.signals['V_m'].my_set_spike_peak( 21, spkSignal= STN.signals['spikes'] )

    
    
    meanRate=round(STN.signals['spikes'].mean_rate(0,500),1)

    s='\n'
    s =s + 'Example inhibitory current:\n'
    s = s + ' %s %5s %3s %s %5s %3s \n' % ( 'Mean rate:', meanRate,  'Hz', 
                                            'I_e', I_e,'pA' )
    s = s + 'Steps:\n'
    s = s + ' %5s %3s \n' % ( I_vec,  'pA' )
    infoString=s
    
    return STN, infoString

def simulate_ahp(I_vec):
    
    simTime  = 3000.  # ms
    my_nest.ResetKernel()
    model_list, model_dict=models()
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    
    n=len(I_vec)
    
    STN = MyGroup( NEURON_MODELS[0], n, sd=True,  mm=True, mm_dt = 1.0 )
    I_e0=my_nest.GetStatus(STN[:])[0]['I_e']
    #my_nest.SetStatus(STN[:], params={'I_e':-10.}) # Set I_e
    my_nest.SetStatus(STN[:], params={'I_e':1.0}) # Set I_e
    I_e = my_nest.GetStatus(STN.ids,'I_e')[0]
    scg = my_nest.Create( 'step_current_generator',n=n )  
    rec=my_nest.GetStatus(STN[:])[0]['receptor_types']
    
    for source, target, I in zip(scg, STN[:], I_vec):
        my_nest.SetStatus([source], {'amplitude_times':[500.,1000.],
                                'amplitude_values':[float(I),0.]})
        my_nest.Connect( [source], [target], 
                         params = { 'receptor_type' : rec['CURR'] } )
    
    
    my_nest.MySimulate(simTime)
    
    STN.get_signal( 's') # retrieve signal
    STN.signals['spikes']= STN.signals['spikes'].time_slice(700,2000)
    
    
    delays=[]
    for i, curr in enumerate(I_vec):
      delays.append(max(numpy.diff(STN.signals['spikes'].spiketrains[i+1.0].spike_times)));
    

    meanRate=round(STN.signals['spikes'].mean_rate(0,500),1)

    s='\n'
    s =s + 'Example inhibitory current:\n'
    s = s + ' %s %5s %3s %s %5s %3s \n' % ( 'Mean rate:', meanRate,  'Hz', 
                                            'I_e', I_e,'pA' )
    s = s + 'Steps:\n'
    s = s + ' %5s %3s \n' % ( I_vec,  'pA' )
    infoString=s
    
    return I_vec, delays

print 'Simulation'   
    
# SIMULATION
infoString=''

I_vec_ahp, delays=simulate_ahp(numpy.arange(0, 350,20) )


# example
STN_example, s=simulate_example(0)
infoString=infoString+s

# IV relation 
#current, voltage, s  = simulate_IV(numpy.arange(-250,-80,10) )
current, voltage, s  = simulate_IV(numpy.arange(-200,0,2) )
infoString=infoString+s

# example inh curr
STN_example_inh_curr, s = simulate_example_inh_current([3,3])
#STN_example_inh_curr, s = simulate_example_inh_current([-74-I_E,-71-I_E])
infoString=infoString+s
#numpy.arange(-85.-I_E,-75.-I_E,5)
# IF relation 
#I_vec=numpy.arange(-30-I_E, -I_E,1)
#I_vec=numpy.append(I_vec,numpy.arange(1, 600,10) )
I_vec=numpy.arange(-99, 301,10)
I_vec, fIsi, mIsi, lIsi, s = simulate_IF(I_vec )
infoString=infoString+s

# example rebound spik
#STN_example_inh_curr, s = simulate_example_inh_current(numpy.arange(-140.,-130.,5))
STN_example_rebound_spike, s   = simulate_example_rebound_spike([-70.], time=300.)
STN_example_rebound_spike_2, s = simulate_example_rebound_spike([-70.], time=450.)
STN_example_rebound_spike_3, s = simulate_example_rebound_spike([-70.], time=600.)
STN_example_rebound_spike_4, s = simulate_example_rebound_spike([-40.], time=300.)
STN_example_rebound_spike_5, s = simulate_example_rebound_spike([-70.], time=300.)
STN_example_rebound_spike_6, s = simulate_example_rebound_spike([-100.],time=300.)


# DISPLAY
plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0, fontsize=16)
font_size_text = 8
fig = pylab.figure( facecolor = 'w' )


ax_list = []
ax_list.append( MyAxes(fig, [ .075, .37, .135, .26 ] ) )    # text box
ax_list.append( MyAxes(fig, [ .26,  .6,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .53, .6,  .165, .34 ] ) )    #  
ax_list.append( MyAxes(fig, [ .8,   .6,  .165, .075 ] ) )    # 
ax_list.append( MyAxes(fig, [ .26,  .3,  .165, .075 ] ) )    # 
ax_list.append( MyAxes(fig, [ .53, .3,  .165, .075 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .3,  .165, .075 ] ) )    # 
ax_list.append( MyAxes(fig, [ .26,  .1,  .165, .075 ] ) )    # 
ax_list.append( MyAxes(fig, [ .53, .1,  .165, .075 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .1,  .165, .075 ] ) )    # 
# Text
ax=ax_list[0]
plot_text(ax, infoString)

# Example
#ax=ax_list[1]
#plot_example(ax, STN_example)

# Example inh current
ax=ax_list[3]
plot_ahp(ax,  I_vec_ahp, delays )

# IV
ax=ax_list[1]
plot_IV(ax, current, voltage)

# IF
ax=ax_list[2]
plot_IF(ax, I_vec, fIsi, mIsi, lIsi )

# Example inh current
ax=ax_list[4]
plot_example_rebound_spike(ax, STN_example_rebound_spike)

ax=ax_list[5]
plot_example_rebound_spike(ax, STN_example_rebound_spike_2)


ax=ax_list[6]
plot_example_rebound_spike(ax, STN_example_rebound_spike_3)


ax=ax_list[7]
plot_example_rebound_spike(ax, STN_example_rebound_spike_4)

ax=ax_list[8]
plot_example_rebound_spike(ax, STN_example_rebound_spike_5)


ax=ax_list[9]
plot_example_rebound_spike(ax, STN_example_rebound_spike_6)


plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0, fontsize=16)
font_size_text = 8
fig2 = pylab.figure( facecolor = 'w' )


ax=MyAxes(fig2, [ .1, .1, .8, .8 ] )
plot_example_rebound_spike(ax, STN_example)
pylab.show()

name = sys.argv[0].split('/')[-1].split('.')[0]
fig.savefig( picture_dir + '/' + name  + '.svg', format = 'svg',transparent=True)
fig.savefig( picture_dir + '/' + name  + '.pdf', format = 'pdf',transparent=True)