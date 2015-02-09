from my_group import MyGroup
import my_signals as msi
import my_nest
import numpy
import os
import pylab
import plot_settings as ps
my_nest.Install('/home/lindahlm/activity-phd/programs/NEST/ml_module/build-pre-ver100725/ml_module')

my_nest.ResetKernel()

# MSN params
params = {'a':0.01, 'b_1':-20., 'b_2':-20., 'c':-55., 'C_m':15.2, 'd':91., 'k':1., 
          'V_b':-80., 'V_peak':40., 'V_r':-80., 'V_t':-29.7, 'I_e':0., 
          'V_m' :-80.,'AMPA_Tau_decay'  : 6., 'GABAA_1_Tau_decay' : 11.4 ,
          'GABAA_1_E_rev'     : -64.,}                                      
my_nest.CopyModel('izhik_cond_exp', 'MSN', params)

rec = my_nest.GetDefaults('izhik_cond_exp')['receptor_types']                  
my_nest.CopyModel('static_synapse', 'syn_ex', 
                  { 'delay':10., 'receptor_type' : rec[ 'AMPA' ] } )
my_nest.CopyModel('static_synapse', 'syn_in', 
                  { 'delay':10., 'receptor_type' : rec[ 'GABAA_1' ] } )


# display recordables for illustration
print 'MSN recordables: ', my_nest.GetDefaults('MSN')['recordables']

n = MyGroup( 'MSN', n=1, params = {'I_e': 240.}, 
           mm_dt=.1, sname='test')

# Create spike generators and connect
gex = my_nest.Create('spike_generator',
                  params = {'spike_times': numpy.arange(10.0, 3000.0, 50.0)})
gin = my_nest.Create('spike_generator',
                  params = {'spike_times': numpy.arange(35.0, 3000.0, 50.0)})

my_nest.Connect(gex, n, params={'weight': 1.0}, model='syn_ex' ) # excitatory
my_nest.Connect(gin, n, params={'weight': 1.0},  model='syn_in') # inhibitory

# simulate
my_nest.Simulate(3000)

n.get_signal('v', 'V_m')
n.get_signal('g', 'g_AMPA')
n.get_signal('g', 'g_GABAA_1')
n.get_signal('s', stop=3000)

ps.set_mode(mode='my_medium', w = 900.0, h = 500.0)

fig = pylab.figure( facecolor = 'w' )
fig.suptitle('Example simulation MSN stimulated with excitatory and inhibitory spikes', fontsize=15)

ax_list=[]
ax_list.append(pylab.axes([0.1, 0.6, 0.35, 0.3]))
ax_list.append(pylab.axes([0.6, 0.6, 0.35, 0.3]))
ax_list.append(pylab.axes([0.1, 0.1, 0.35, 0.3]))
ax_list.append(pylab.axes([0.6, 0.1, 0.35, 0.3]))

ax=ax_list[0]
n.signals['V_m'].plot(display=ax)

ax=ax_list[1]
n.signals['g_AMPA'].plot(display=ax,kwargs={'label':'Exc'})
n.signals['g_GABAA_1'].plot(display=ax, kwargs={'label':'Inh'})
ax.legend()

ax=ax_list[2]
n.signals['spikes'].firing_rate( bin=500,display=ax)
ax.set_title('Firing rate from histogram (bin=%i)'%(500), fontsize=15)

ax=ax_list[3]
n.signals['spikes'].firing_rate_sliding_window(bin=500, display=ax, step=10)
ax.set_title('Firing rate from sliding window (bin=%i)'%(500), fontsize=15)

pylab.show()







