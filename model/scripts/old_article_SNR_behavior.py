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

neuronModels=['SNR_izh']
# 270 in curr gives 62.5 Hz in fring frequency
I_E=0.0-30

def plot_example(ax, SNR):
    SNR.signals['V_m'].plot(display=ax, kwargs={'color':'k'})
    
    ax.set_xlim(misc.adjust_limit([0,2000]))
    ax.set_ylim(misc.adjust_limit([-80,30]))
    
    ax.my_set_no_ticks( yticks=5, xticks = 5 )  

def plot_example_inh_current(ax, SNR):
    SNR.signals['V_m'].plot(display=ax, kwargs={'color':'k'})
    ax.set_xlim(misc.adjust_limit([0,2000]))
    ax.set_ylim(misc.adjust_limit([-80,30]))
    ax.my_set_no_ticks( yticks=6, xticks=4 )
    ax.set_ylabel('Potential (mV)')

def plot_IV(ax, current, voltage):
    ax.plot(current, voltage, **{'color':'k'})
    ax.my_set_no_ticks( yticks=5, xticks=5 ) 
    ax.set_xlabel('Current (pA)') 
    ax.set_ylabel('Potential (mV)')
    ax.set_xlim(misc.adjust_limit([-150-I_E,-90-I_E]))
    ax.set_ylim(misc.adjust_limit([-80,-64]))

def plot_IF(ax, I_vec, fIsi, mIsi, lIsi):

    ax.plot(I_vec, 1000./fIsi, **{'label':'$First_{ISI}$', 'color':'k'})
    ax.plot(I_vec, 1000./lIsi, **{'label':'$Last_{ISI}$', 'color':'k',
                                  'linestyle':'--'})

    print 1000./lIsi[I_vec<100]
    print I_vec[I_vec<100]
    ax.my_set_no_ticks( yticks=5, xticks = 5 ) 
    ax.set_xlabel('Current (pA)') 
    ax.set_ylabel('Frequency (Hz)')
    ax.set_xlim(misc.adjust_limit([-100,300]))
    ax.set_ylim(misc.adjust_limit([0,120]))
    leg=ax.legend(numpoints=1, loc='upper left')
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    llines = leg.get_lines()  # all the lines.Line2D instance in the legend
    frame  = leg.get_frame() 
    pylab.setp(ltext, fontsize=14) 
    frame.set_visible(False)
    
    #for line in llines:        
    #    line.set_xdata(line.get_xdata()-400)    
    #    line.set_ydata(line.get_ydata()+68.) 
    
def plot_text(ax, infoString=''):
    
    my_nest.ResetKernel()
    model_list=models()
    my_nest.MyLoadModels( model_list, neuronModels )
    
    SNR = MyGroup( neuronModels[0], 1, mm_dt = 0.1)
    statusSNR = my_nest.GetStatus( SNR[:] )[0]
    
    tb = ''     
    tb = tb + infoString
    
    tb = tb + '\n'
    for key, val in statusSNR.iteritems():
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
    
    simTime  = 1000.  # ms
    
    my_nest.ResetKernel()
    model_list=models()
    my_nest.MyLoadModels( model_list, neuronModels )
    
    SNR = MyGroup( neuronModels[0], 1, sd=True, mm=True, mm_dt = 1., 
                   params={'I_e':I_e} )
    
    my_nest.MySimulate(simTime)
    SNR.get_signal( 'v','V_m', stop=simTime ) # retrieve signal
    
    SNR.get_signal( 's') # retrieve signal
    meanRate=round(SNR.signals['spikes'].mean_rate(0,1000),1)
    print SNR.signals['spikes'].isi()
    SNR.signals['V_m'].my_set_spike_peak( 21, spkSignal= SNR.signals['spikes'] )
    
    s='\n'
    s = s + 'Example:\n'
    s = s + ' %s %5s %3s %s %5s %3s \n' % ( 'Mean rate:', meanRate,  'Hz', 
                                            'I_e', I_e,'pA' )
    
    infoString=s
    
    return SNR, infoString

def simulate_example_inh_current(I_vec):
    
    simTime  = 2000.  # ms
    my_nest.ResetKernel()
    model_list=models()
    my_nest.MyLoadModels( model_list, neuronModels )
    
    n=len(I_vec)
    
    SNR = MyGroup( neuronModels[0], n, sd=True,  mm=True, mm_dt = 1.0 )
    I_e0=my_nest.GetStatus(SNR[:])[0]['I_e']
    my_nest.SetStatus(SNR[:], params={'I_e':I_e0+I_E}) # Set I_e
    
    I_e = my_nest.GetStatus(SNR.ids,'I_e')[0]
    scg = my_nest.Create( 'step_current_generator',n=n )  
    rec=my_nest.GetStatus(SNR[:])[0]['receptor_types']
    
    for source, target, I in zip(scg, SNR[:], I_vec):
        my_nest.SetStatus([source], {'amplitude_times':[500.,1500.],
                                'amplitude_values':[float(I),0.]})
        my_nest.Connect( [source], [target], 
                         params = { 'receptor_type' : rec['CURR'] } )
    
    
    my_nest.MySimulate(simTime)
    SNR.get_signal( 'v','V_m', stop=simTime ) # retrieve signal
    SNR.get_signal( 's') # retrieve signal
    SNR.signals['V_m'].my_set_spike_peak( 21, spkSignal= SNR.signals['spikes'] )

    
    
    meanRate=round(SNR.signals['spikes'].mean_rate(0,500),1)

    s='\n'
    s =s + 'Example inhibitory current:\n'
    s = s + ' %s %5s %3s %s %5s %3s \n' % ( 'Mean rate:', meanRate,  'Hz', 
                                            'I_e', I_e,'pA' )
    s = s + 'Steps:\n'
    s = s + ' %5s %3s \n' % ( I_vec,  'pA' )
    infoString=s
    
    return SNR, infoString

def simulate_IV(I_vec):
    
    I_e = 0.+I_E
    
    my_nest.ResetKernel()
    model_list=models()
    my_nest.MyLoadModels( model_list, neuronModels )
    
    SNR = MyGroup( neuronModels[0], 1, sd=True, mm=True, mm_dt = 0.1)
    
    I_e0=my_nest.GetStatus(SNR[:])[0]['I_e']
    my_nest.SetStatus(SNR[:], params={'I_e':I_e0+I_E}) # Set I_e 
    I_e = my_nest.GetStatus(SNR.ids,'I_e')[0]
    
    I_vec, voltage = SNR.IV_I_clamp(I_vec)   
    
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

    return I_vec, voltage, infoString

def simulate_IF(I_vec):
    
    tStim = 500
    
    my_nest.ResetKernel()
    model_list=models()
    my_nest.MyLoadModels( model_list, neuronModels )
    
    SNR = MyGroup( neuronModels[0], 1, mm=True, sd=True, mm_dt = 0.1 )
    
    I_e0=my_nest.GetStatus(SNR[:])[0]['I_e']
    my_nest.SetStatus(SNR[:], params={'I_e':I_e0+I_E}) # Set I_e 
    I_e = my_nest.GetStatus(SNR.ids,'I_e')[0]

    I_vec, fIsi, mIsi, lIsi = SNR.IF(I_vec, tStim=tStim)   
    
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
                                          str(min(speed_f))[0:4],  
                                          'max',str(max(speed_f))[0:4],
                                          'mean', 
                                          str(sum(speed_f)/len(speed_f))[0:4])
    infoString=s

    return I_vec, fIsi, mIsi, lIsi, infoString

print 'Simulation'   
    
# SIMULATION
infoString=''
# example
SNR_example, s=simulate_example(81.+I_E)
infoString=infoString+s

# IV relation 
#current, voltage, s  = simulate_IV(numpy.arange(-250,-80,10) )
current, voltage, s  = simulate_IV(numpy.arange(-150,-30,10) )
infoString=infoString+s

# example inh curr
#SNR_example_inh_curr, s = simulate_example_inh_current(numpy.arange(-140.,-130.,5))
SNR_example_inh_curr, s = simulate_example_inh_current(numpy.arange(-86.-I_E,-76.-I_E,5))
infoString=infoString+s

# IF relation 
I_vec, fIsi, mIsi, lIsi, s = simulate_IF(numpy.arange(-200-I_E, 300-I_E,10) )
infoString=infoString+s

# DISPLAY
plot_settings.set_mode(pylab, mode='by_fontsize', w = 1100.0, h = 450.0, fontsize=12)
font_size_text = 8
fig = pylab.figure( facecolor = 'w' )

ax_list = []
ax_list.append( MyAxes(fig, [ .075, .37, .135, .26 ] ) )    # text box
ax_list.append( MyAxes(fig, [ .26,  .6,  .165, .34 ] ) )    # 
ax_list.append( MyAxes(fig, [ .53, .6,  .165, .34 ] ) )    #  
ax_list.append( MyAxes(fig, [ .8,   .6,  .165, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .26,  .1,  .165, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .53, .1,  .165, .34 ] ) )    # 
#ax_list.append( MyAxes(fig, [ .8,   .1,  .165, .34 ] ) )    # 

# Text
ax=ax_list[0]
plot_text(ax, infoString)

# Example
#ax=ax_list[1]
#plot_example(ax, SNR_example)

# Example inh current
ax=ax_list[3]
plot_example_inh_current(ax, SNR_example_inh_curr)

# IV
ax=ax_list[1]
plot_IV(ax, current, voltage)

# IF
ax=ax_list[2]
plot_IF(ax, I_vec, fIsi, mIsi, lIsi )

pylab.show()

name = sys.argv[0].split('/')[-1].split('.')[0]
fig.savefig( picture_dir + '/' + name  + '.svg', format = 'svg',transparent=True)
fig.savefig( picture_dir + '/' + name  + '.pdf', format = 'pdf',transparent=True)