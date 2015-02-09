from  nest_toolbox.my_population import MyGroup
import nest_toolbox.my_signals as my_signals
import nest_toolbox.my_nest as my_nest
import nest_toolbox.plot_settings as ps
import numpy
import os
import pylab

my_nest.ResetKernel()
       
my_nest.CopyModel('static_synapse', 'syn_ex', )
my_nest.CopyModel('static_synapse', 'syn_in', )


# display recordables for illustration
print 'MSN recordables: ', my_nest.GetDefaults('aeif_cond_exp')['recordables']

n = MyGroup( 'aeif_cond_exp', n=1, params = {'I_e': 700.}, 
           mm_dt=.1, sname='test')

# Create spike generators and connect
gex = my_nest.Create('spike_generator',
                  params = {'spike_times': numpy.arange(10.0, 1000.0, 50.0)})
gin = my_nest.Create('spike_generator',
                  params = {'spike_times': numpy.arange(35.0, 1000.0, 50.0)})

my_nest.Connect(gex, n, params={'weight': 10.0}, model='syn_ex' ) # excitatory
my_nest.Connect(gin, n, params={'weight': -10.0}, model='syn_in' ) # inhibitory

# simulate
my_nest.Simulate(1000)

n.get_signal('v', 'V_m')
n.get_signal('g', 'g_ex')
n.get_signal('g', 'g_in')
n.get_signal('s', stop=1000)

ps.set_mode(mode='my_medium', w = 900.0, h = 500.0)

fig = pylab.figure( facecolor = 'w' )
fig.suptitle('Example simulation aeif_cond_exp stimulated with excitatory and inhibitory spikes', fontsize=15)

ax_list=[]
ax_list.append(pylab.axes([0.1, 0.6, 0.35, 0.3]))
ax_list.append(pylab.axes([0.6, 0.6, 0.35, 0.3]))
ax_list.append(pylab.axes([0.1, 0.1, 0.35, 0.3]))
ax_list.append(pylab.axes([0.6, 0.1, 0.35, 0.3]))

ax=ax_list[0]
n.signals['V_m'].plot(display=ax)

ax=ax_list[1]
n.signals['g_ex'].plot(display=ax,kwargs={'label':'Exc'})
n.signals['g_in'].plot(display=ax, kwargs={'label':'Inh'})
ax.legend()

ax=ax_list[2]
n.signals['spikes'].firing_rate( bin=200,display=ax)
ax.set_title('Firing rate from histogram (bin=%i)'%(200), fontsize=15)

ax=ax_list[3]
n.signals['spikes'].firing_rate_sliding_window(bin=200, display=ax, step=10)
ax.set_title('Firing rate from sliding window (bin=%i)'%(200), fontsize=15)

pylab.show()







