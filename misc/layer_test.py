import nest
import nest.topology as tp
import time 
l_dummy1 = tp.CreateLayer({'rows' : 15000 , 'columns' : 1 , 'elements' : 'iaf_neuron',
                    'edge_wrap': True, 'extent':[1.,1.], 'center' : [0.,0.] } )
l_dummy2 = tp.CreateLayer({'rows' : 15000 , 'columns' : 1 , 'elements' : 'iaf_neuron',
                    'edge_wrap': True, 'extent':[1.,1.], 'center' : [0.,0.] } )
l = tp.CreateLayer({'rows' : 100 , 'columns' : 1 , 'elements' : 'iaf_neuron',
                    'edge_wrap': True, 'extent':[1.,1.], 'center' : [0.,0.] } )
projections = {'allow_multapses': False,
               'allow_autapses' : False,
               'connection_type' : 'convergent',
               'mask':{'rectangular' : { 'lower_left' : [ -0.5 , -0.5 ], 
                                        'upper_right' : [ 0.5 , 0.5 ]}},
               'kernel': 0.5,
               'sources': {'model' : 'iaf_neuron' },
               'targets': {'model' : 'iaf_neuron' }}

tp.ConnectLayers(l, l, projections)
tp.ConnectLayers(l_dummy1, l, projections)
tp.ConnectLayers(l_dummy2, l, projections)

        

if nest.version()=='NEST 2.0.0':
    sources=nest.GetChildren(l)
if nest.version()=='NEST 2.2.1':
    sources=nest.GetChildren(l)[0]

t=time.time()
targets=tp.GetTargetNodes(sources, l)
print time.time()-t


