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

NEURON_MODELS=['GPE_aeif']

# 270 in curr gives 62.5 Hz in fring frequency
I_E=0.0
I_E2=63
def plot_example(ax, GPE):
    GPE.signals['V_m'].plot(display=ax, kwargs={'color':'k'})
    
    ax.set_xlim(misc.adjust_limit([0,2000]))
    ax.set_ylim(misc.adjust_limit([-80,30]))
    
    ax.my_set_no_ticks( yticks=8, xticks = 6 )  

def plot_example_inh_current(ax, GPE, x_lim=[0,1000]):
    GPE.signals['V_m'].plot(display=ax, kwargs={'color':'k'})
    ax.set_xlim(misc.adjust_limit(x_lim))
    ax.set_ylim(misc.adjust_limit([-80,30]))
    ax.my_set_no_ticks( yticks=3, xticks = 6 )  
    ax.set_ylabel('GPe potential (mV)')
    #pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend
    
def plot_IV(ax, current, voltage):
    color='g'
    '''
    lower_limit=-55
    upper_limit=-100
    current=current[voltage<lower_limit]
    voltage=voltage[voltage<lower_limit]
    current=current[voltage>upper_limit]
    voltage=voltage[voltage>upper_limit]
    '''
    current=current-current[-1]    
    ax.plot(current, voltage, **{'color':color})
    
    line=ax.lines



    ax.my_set_no_ticks( yticks=8, xticks = 6 )  
    ax.set_xlabel('Current (pA)') 
    ax.set_ylabel('Potential (mV)')
    #ax.set_xlim(misc.adjust_limit([-220-I_E,-110-I_E]))
    #ax.set_ylim(misc.adjust_limit([-80,limit]))
    
    df=my_nest.GetDefaults(NEURON_MODELS[0]) 
    
    speed=numpy.diff(voltage)/numpy.diff(current/1000.) 
    mean_resistance=numpy.mean(speed)
    mean_membrane_time_constant=mean_resistance*df['C_m']*0.001

    point=-15
    ax.text( 0.50, 0.75, str(round(speed[point]))+r' M$\Omega$', backgroundcolor='w',
             transform=ax.transAxes,fontsize=16,  **{ 'color' : 'k' }) 
    ax.plot(current[point], voltage[point],color=color, marker='.',markersize=10)
    x=[current[point], current[point]-100]
    y=[voltage[point], voltage[point]-(100)*speed[point]/1000.]
    ax.plot(x,y,color=color, linestyle='--')
    
    ax.set_xlim(misc.adjust_limit([-200,0]))
    ax.set_ylim(misc.adjust_limit([-85,-55]))
    
    pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend
def plot_IF(ax, I_vec, fIsi, mIsi, lIsi):
    color='g'
    ax.plot(I_vec, 1000./fIsi, **{'label':'First ISI', 'color':color, 'linestyle':'--'})
    ax.plot(I_vec, 1000./lIsi, **{'label':'Steady-state ISI', 'color':color})

    
    print 1000./lIsi[I_vec<100]
    print I_vec[I_vec<100]
    ax.my_set_no_ticks( yticks=8, xticks = 6 )  
    ax.set_xlabel('Current (pA)') 
    ax.set_ylabel('Rate (spikes/s)')

    leg=ax.legend(numpoints=1, loc='upper left')
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    llines = leg.get_lines()  # all the lines.Line2D instance in the legend
    frame  = leg.get_frame() 
    pylab.setp(ltext, fontsize=14) 
    frame.set_visible(False)
    
    
    speed=numpy.diff(1000./lIsi)/numpy.diff(I_vec) 

    #ax.text( 0.3, 0.02, 'Mean slope steady-state\nISI '+str(round(numpy.median(speed),2))+' Hz/pA', backgroundcolor='w',
    #        transform=ax.transAxes,fontsize=10,  **{ 'color' : color }) 
    point=numpy.nonzero(I_vec==101)[0][0]
    ax.text( 0.6, 0.2, 'GPe', backgroundcolor='w',
             transform=ax.transAxes,fontsize=16,  **{ 'color' : color }) 
    ax.text( 0.6, 0.02, str(round(speed[point],2))+' Hz/pA', backgroundcolor='w',
             transform=ax.transAxes,fontsize=16,  **{ 'color' : 'k' }) 
    ax.plot(I_vec[point], 1000./lIsi[point],color=color, marker='.',markersize=10)

    ax.set_xlim(misc.adjust_limit([-100,300]))
    ax.set_ylim(misc.adjust_limit([0,260]))
    
    pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend
    
def plot_example_rebound_spike(ax, GPE):
    GPE.signals['V_m'].plot(display=ax, kwargs={'color':'k'})   
    ax.set_xlim(misc.adjust_limit([300,1000]))
    ax.set_ylim(misc.adjust_limit([-85,30]))
    line=ax.lines
    misc.slice_line(line[0], xlim=[300,1000])

    ax.my_set_no_ticks( yticks=3, xticks = 6 )  
    ax.set_ylabel('GPe potential (mV)')
    
    #pylab.setp(ax.lines, linewidth=2.0) # Need to pu ti before generating legend
def plot_text(ax, infoString=''):
    
    my_nest.ResetKernel()
    model_list, model_dict=models()
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    
    GPE = MyGroup( NEURON_MODELS[0], 1, mm_dt = 0.1)
    statusGPE = my_nest.GetStatus( GPE[:] )[0]
    
    tb = ''     
    tb = tb + infoString
    
    tb = tb + '\n'
    for key, val in statusGPE.iteritems():
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
    
    simTime  = 2000.  # ms
    
    my_nest.ResetKernel()
    model_list, model_dict=models()
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    
    I_e0=my_nest.GetDefaults(NEURON_MODELS[0])['I_e']
    GPE = MyGroup( NEURON_MODELS[0], 1, sd=True, mm=True, mm_dt = 1., params={'I_e':I_e+I_e0}  )
    
    my_nest.MySimulate(simTime)
    GPE.get_signal( 'v','V_m', stop=simTime ) # retrieve signal
    
    GPE.get_signal( 's') # retrieve signal
    meanRate=round(GPE.signals['spikes'].mean_rate(500,2000),1)
    print GPE.signals['spikes'].isi()
    GPE.signals['V_m'].my_set_spike_peak( 21, spkSignal= GPE.signals['spikes'] )
    
    s='\n'
    s = s + 'Example:\n'
    s = s + ' %s %5s %3s %s %5s %3s \n' % ( 'Mean rate:', meanRate,  'Hz', 
                                            'I_e', I_e,'pA' )
    
    infoString=s
    
    return GPE, infoString

def simulate_example_inh_current(I_vec):
    
    simTime  = 1000.  # ms
    my_nest.ResetKernel()
    model_list, model_dict=models()
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    
    n=len(I_vec)
    
    GPE = MyGroup( NEURON_MODELS[0], n, sd=True,  mm=True, mm_dt = 1.0 )
    I_e0=my_nest.GetStatus(GPE[:])[0]['I_e']
    my_nest.SetStatus(GPE[:], params={'I_e':I_e0+I_E}) # Set I_e
    
    I_e = my_nest.GetStatus(GPE.ids,'I_e')[0]
    scg = my_nest.Create( 'step_current_generator',n=n )  
    rec=my_nest.GetStatus(GPE[:])[0]['receptor_types']
    
    for source, target, I in zip(scg, GPE[:], I_vec):
        my_nest.SetStatus([source], {'amplitude_times':[280.,700.],
                                'amplitude_values':[float(I),0.]})
        my_nest.Connect( [source], [target], 
                         params = { 'receptor_type' : rec['CURR'] } )
    
    
    my_nest.MySimulate(simTime)
    GPE.get_signal( 'v','V_m', stop=simTime ) # retrieve signal
    GPE.get_signal( 's') # retrieve signal
    GPE.signals['V_m'].my_set_spike_peak( 15, spkSignal= GPE.signals['spikes'] )

    
    
    meanRate=round(GPE.signals['spikes'].mean_rate(0,500),1)

    s='\n'
    s =s + 'Example inhibitory current:\n'
    s = s + ' %s %5s %3s %s %5s %3s \n' % ( 'Mean rate:', meanRate,  'Hz', 
                                            'I_e', I_e,'pA' )
    s = s + 'Steps:\n'
    s = s + ' %5s %3s \n' % ( I_vec,  'pA' )
    infoString=s
    
    return GPE, infoString

def simulate_example_irregular_firing(I_vec=[0]):
    
    simTime  = 2000.  # ms
    my_nest.ResetKernel()
    model_list, model_dict=models()
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    
    n=len(I_vec)
    
    GPE = MyGroup( NEURON_MODELS[0], n, sd=True,  mm=True, mm_dt = 1.0 )
    I_e=I_vec[0]
    my_nest.SetStatus(GPE[:], params={'I_e':I_e}) # Set I_e

    
    
    scg = my_nest.Create( 'step_current_generator',n=n )  
    noise=my_nest.Create('noise_generator', params={'mean':0.,'std':10.})
    rec=my_nest.GetStatus(GPE[:])[0]['receptor_types']
    
    for source, target, I in zip(scg, GPE[:], I_vec):
        #I=5.
        my_nest.SetStatus([source], {'amplitude_times':[1., simTime],
                                'amplitude_values':[-5.,float(I)]})
        my_nest.Connect( [source], [target], 
                         params = { 'receptor_type' : rec['CURR'] } )
        my_nest.Connect( noise, [target], 
                         params = { 'receptor_type' : rec['CURR'] } )
    
    my_nest.MySimulate(simTime)
    GPE.get_signal( 'v','V_m', stop=simTime ) # retrieve signal
    GPE.get_signal( 's') # retrieve signal
    GPE.signals['V_m'].my_set_spike_peak( 15, spkSignal= GPE.signals['spikes'] )
    
    
    #a=GPE.signals['V_m'].analog_signals[1].signal
    #pylab.plot(a)
#   #a=a[500:]
#   pylab.subplot(211).plot(a, 'r')
    #pylab.show()
#    
#    a=a    
#    a=a-numpy.mean(a)
#    
#    numpy.savetxt("foo.csv", a, delimiter=",")
#
#    ff=numpy.abs(numpy.fft.fft(a))
#    #pylab.plot(ff)
#    c=numpy.correlate(a, a, mode='full')
#    pylab.subplot(212).plot(c)
#    pylab.show()

    meanRate=round(GPE.signals['spikes'].mean_rate(0,500),1)

    

    s='\n'
    s =s + 'Example inhibitory current:\n'
    s = s + ' %s %5s %3s %s %5s %3s \n' % ( 'Mean rate:', meanRate,  'Hz', 
                                            'I_e', I_e,'pA' )
    s = s + 'Steps:\n'
    s = s + ' %5s %3s \n' % ( I_vec,  'pA' )
    infoString=s
    
    return GPE, infoString

def simulate_IV(I_vec):
    
    I_e = 0.+I_E
    
    my_nest.ResetKernel()
    model_list, model_dict=models()
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    
    GPE = MyGroup( NEURON_MODELS[0], 1, sd=True, mm=True, mm_dt = 0.1)
    
    I_e0=my_nest.GetStatus(GPE[:])[0]['I_e']
    my_nest.SetStatus(GPE[:], params={'I_e':I_e0+I_E}) # Set I_e 
    I_e = my_nest.GetStatus(GPE.ids,'I_e')[0]
    
    I_vec, voltage = GPE.IV_I_clamp(I_vec)   
    
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
    
    tStim = 700
    
    my_nest.ResetKernel()
    model_list, model_dict=models()
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    
    GPE = MyGroup( NEURON_MODELS[0], 1, mm=True, sd=True, mm_dt = 0.1 )
    
    I_e0=my_nest.GetStatus(GPE[:])[0]['I_e']
    my_nest.SetStatus(GPE[:], params={'I_e':I_e0+I_E}) # Set I_e 
    I_e = my_nest.GetStatus(GPE.ids,'I_e')[0]

    I_vec, fIsi, mIsi, lIsi = GPE.IF(I_vec, tStim=tStim)   
    
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

def simulate_example_rebound_spike(I_vec):
    
    simTime  = 5000.  # ms
    my_nest.ResetKernel()
    model_list, model_dict=models()
    my_nest.MyLoadModels( model_list, NEURON_MODELS )
    
    n=len(I_vec)
    
    GPE = MyGroup( NEURON_MODELS[0], n, sd=True,  mm=True, mm_dt = 1.0 )
    I_e0=my_nest.GetStatus(GPE[:])[0]['I_e']
    #my_nest.SetStatus(GPE[:], params={'I_e':-10.}) # Set I_e
    my_nest.SetStatus(GPE[:], params={'I_e':-10.}) # Set I_e
    I_e = my_nest.GetStatus(GPE.ids,'I_e')[0]
    scg = my_nest.Create( 'step_current_generator',n=n )  
    rec=my_nest.GetStatus(GPE[:])[0]['receptor_types']
    
    for source, target, I in zip(scg, GPE[:], I_vec):
        my_nest.SetStatus([source], {'amplitude_times':[500.,700.],
                                'amplitude_values':[float(I),0.]})
        my_nest.Connect( [source], [target], 
                         params = { 'receptor_type' : rec['CURR'] } )
    
    
    my_nest.MySimulate(simTime)
    GPE.get_signal( 'v','V_m', stop=simTime ) # retrieve signal
    GPE.get_signal( 's') # retrieve signal
    GPE.signals['V_m'].my_set_spike_peak( 21, spkSignal= GPE.signals['spikes'] )

    
    
    meanRate=round(GPE.signals['spikes'].mean_rate(0,500),1)

    s='\n'
    s =s + 'Example inhibitory current:\n'
    s = s + ' %s %5s %3s %s %5s %3s \n' % ( 'Mean rate:', meanRate,  'Hz', 
                                            'I_e', I_e,'pA' )
    s = s + 'Steps:\n'
    s = s + ' %5s %3s \n' % ( I_vec,  'pA' )
    infoString=s
    
    return GPE, infoString

print 'Simulation'   
    
# SIMULATION
infoString=''
# example
GPE_example, s=simulate_example(0.)
infoString=infoString+s

# IV relation 
#current, voltage, s  = simulate_IV(numpy.arange(-250,-80,10) )
current, voltage, s  = simulate_IV(numpy.arange(-210,-10,2) )
infoString=infoString+s

# example inh curr
#GPE_example_inh_curr, s = simulate_example_inh_current(numpy.arange(-140.,-130.,5))
GPE_example_inh_curr, s = simulate_example_inh_current([-80-I_E])

infoString=infoString+s
#numpy.arange(-85.-I_E,-75.-I_E,5)
# IF relation 
I_vec, fIsi, mIsi, lIsi, s = simulate_IF(numpy.arange(-99, 299,10) )
infoString=infoString+s

# example rebound spike
#GPE_example_inh_curr, s = simulate_example_inh_current(numpy.arange(-140.,-130.,5))
GPE_example_rebound_spike, s = simulate_example_rebound_spike([-80.])


GPE_example_irregular_firing, s=simulate_example_irregular_firing([5.])
GPE_example_irregular_firing2, s=simulate_example_irregular_firing([10.])
GPE_example_irregular_firing3, s=simulate_example_irregular_firing([40.])


# DISPLAY
plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0, fontsize=16)
font_size_text = 8
fig = pylab.figure( facecolor = 'w' )

ax_list = []
ax_list.append( MyAxes(fig, [ .075, .37, .135, .26 ] ) )    # text box
ax_list.append( MyAxes(fig, [ .26,  .6,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .53, .6,  .165, .34 ] ) )    #  
ax_list.append( MyAxes(fig, [ .8,   .6,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .26,  .1+.2,  .165, .075 ] ) )    # 
ax_list.append( MyAxes(fig, [ .53, .1+.2,  .165, .075 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .1+.2,  .165, .075 ] ) )    # 
ax_list.append( MyAxes(fig, [ .26,  .1,  .165, .075 ] ) )    # 
ax_list.append( MyAxes(fig, [ .53, .1,  .165, .075 ] ) )    # 
ax_list.append( MyAxes(fig, [ .8,   .1,  .165, .075 ] ) )    # 
# Text
ax=ax_list[0]
plot_text(ax, infoString)

# Example
#ax=ax_list[1]
#plot_example(ax, GPE_example)



# IV
ax=ax_list[1]
plot_IV(ax, current, voltage)

# IF
ax=ax_list[2]
plot_IF(ax, I_vec, fIsi, mIsi, lIsi )

# Example inh current
ax=ax_list[3]
plot_example(ax, GPE_example)

# Example inh current
ax=ax_list[4]
plot_example_rebound_spike(ax, GPE_example_rebound_spike)

ax=ax_list[5]
plot_example_inh_current(ax, GPE_example_irregular_firing,x_lim=[0,2000])

ax=ax_list[6]
plot_example_inh_current(ax, GPE_example_irregular_firing2,x_lim=[0,2000])
ax=ax_list[7]
plot_example_inh_current(ax, GPE_example_irregular_firing3,x_lim=[0,2000])

pylab.show()

name = sys.argv[0].split('/')[-1].split('.')[0]
fig.savefig( picture_dir + '/' + name  + '.svg', format = 'svg',transparent=True)
fig.savefig( picture_dir + '/' + name  + '.pdf', format = 'pdf',transparent=True)