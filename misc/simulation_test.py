import nest
import nest.topology as tp
import numpy
import pylab
import time 

print nest.GetDefaults('iaf_neuron')
nest.SetDefaults('iaf_neuron',{'I_e':16.*250.0/10.0})

n_neurons=100
sim_time=5000.0
start_rec=1000.0
res=100.0
weight=-10.0
delay=1.0
dev=0.5
eps=0.0001
l = tp.CreateLayer({'rows' : n_neurons , 'columns' : 1 , 'elements' : 'iaf_neuron',
                    'edge_wrap': True, 'extent':[1.,1.], 'center' : [0.,0.] } )

projections = {'allow_multapses': False,
               'allow_autapses' : False,
               'delays':{'uniform':{'min':(1-dev)*delay,'max':(1.+dev)*delay}}, 
                'weights':{'uniform':{'min':(1+dev)*weight,'max':(1-dev)*weight}},
               'connection_type' : 'convergent',
               'mask':{'rectangular' : { 'lower_left' : [ -0.5 , -0.5 ], 
                                        'upper_right' : [ 0.5, 0.5 ]}},
               'kernel': 0.5,
               'sources': {'model' : 'iaf_neuron' },
               'targets': {'model' : 'iaf_neuron' }}

tp.ConnectLayers(l, l, projections)

if nest.version()=='NEST 2.0.0':
    sources=nest.GetChildren(l)
if nest.version()=='NEST 2.2.1':
    sources=nest.GetChildren(l)[0]
    
sd=nest.Create('spike_detector')
nest.SetStatus(sd, {"withgid": True,'start':start_rec })
nest.ConvergentConnect(sources, sd)

mm=nest.Create('multimeter')
recodables=['V_m']
nest.SetStatus(mm, {'interval': 0.1, 'record_from': recodables, 'start':start_rec})
nest.Connect(mm, [sources[0]])

nest.Simulate(sim_time)

x=nest.GetStatus(sd)[0]['events']['times']
y=nest.GetStatus(sd)[0]['events']['senders']

hist, bin_edges=numpy.histogram(x, sim_time-start_rec ,[start_rec, sim_time])
hist=numpy.array(hist)
hist=hist/float(n_neurons)*1000.0

kernel=numpy.arange(1,numpy.ceil(res/2.)+1)
kernel=numpy.append(kernel, numpy.arange(numpy.floor(res/2.),0,-1))
kernel=kernel/float(numpy.sum(kernel))
hist_conv=numpy.convolve(kernel, hist, mode='same')

ax=pylab.subplot(211)
ax.plot(hist_conv)
ax.set_title(nest.version()+', m_rate:'+str(numpy.mean(hist)))
ax=pylab.subplot(212)
st_mm=nest.GetStatus(mm)[0]
ax.plot(st_mm['events']['V_m'])
pylab.show()

