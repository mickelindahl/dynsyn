'''
Mikael Lindahl 2010


Module:
misc


'''

# Imports
import matplotlib
import nest
import numpy
import os
import pylab
from NeuroTools.stgen import StGen
from numpy import array, concatenate
from numpy.random import poisson 

from scipy.cluster.vq import vq, kmeans, kmeans2, whiten

'''
Method for method to evenly split up a colormap into an RGB colors array
'''

     
def kmean_cluster( data, k, iter = 10):
    '''
    The k-means algorithm takes as input the number of clusters to generate, 
    k, and a set of observation vectors to cluster. It returns a set of 
    centroids, one for each of the k clusters. An observation vector is
    classified with the cluster number or centroid index of the centroid  
    closest to it.
    
    Inputs:
        data          - A M by N array of M observations in N dimensions or a 
                         length M array of M one-dimensional observations.
        k             - The number of clusters to form as well as the number 
                         of centroids to generate. If minit initialization 
                         string is 'matrix', or if a ndarray is given instead, 
                         it is interpreted as initial cluster to use instead.
        iter          - Number of iterations of the k-means algrithm to run. 
                        Note that this differs in meaning from the iters 
                        parameter to the kmeans function.
    
    Returns: 
        codebook      - A k by N array of k centroids. The i'th centroid 
                        codebook[i] is represented with the code i. The 
                        centroids and codes     generated represent the lowest 
                        distortion seen, not necessarily the globally minimal
                         distortion.
        distortion    - The distortion between the observations passed and the 
                        centroids generated.
        code          - A length N array holding the code book index for each 
                        observation.
        dist          - The distortion (distance) between the observation and 
                        its nearest code.
                       
    Examples:
        >>> from misc import kmean_cluster
        >>> features  = array([[ 1.9,2.3],
        ...                    [ 1.5,2.5],
        ...                    [ 0.8,0.6],
        ...                    [ 0.4,1.8],
        ...                    [ 0.1,0.1],
        ...                    [ 0.2,1.8],
        ...                    [ 2.0,0.5],
        ...                    [ 0.3,1.5],
        ...                    [ 1.0,1.0]])
        >>> kmean_cluster( features, 2, iter = 10):
        (array([[ 2.3110306 ,  2.86287398],
           [ 0.93218041,  1.24398691]]), 0.85684700941625547)
    '''

    whitened = whiten(data)
    codebook, distortion = kmeans( whitened, k, iter = iter)
    code, dist    =  vq( whitened, codebook )

    return codebook, distortion, code, dist   
    
def kmean_image( data, code, times = False, ids = False, display=False, kwargs={} ):
        '''
        Plots the kmean data accoringly to clusters
        Inputs:
            data      - A M by N array of M observations in N dimensions or a 
                        length M array of M one-dimensional observations.
            code       - A length N array holding the code book index for each 
                        observation.
            times      - vector with the times of the sliding window
        '''
        
        if not display: ax = pylab.axes()
        else:           ax = display
        
        
        sorted_index = numpy.argsort(code)
        sorted_data  = numpy.array([data[i] for i in sorted_index])
               
        kwargs.update( { 'origin' : 'lower', } )
        if any(times) and any(ids): image = ax.imshow(sorted_data, extent=[times[0],times[-1], ids[0], ids[-1]], **kwargs)
        else: image = ax.imshow(sorted_data, **kwargs)
        ax.set_xlabel('Time (ms)')
        ax.set_ylabel('Neuron #')
        ax.set_aspect(aspect='auto')

        return image
      
def kmean_raster_plot( data_spk, code, times = False, ids = False, display=False, kwargs={} ):
        '''
        Plots the kmean data accoringly to clusters
        Inputs:
            data_spk   - List with spike data for each id
            code       - A length N array holding the code book index for each 
                        observation.
            times      - vector with the times of the sliding window
        '''
        
        if not display: ax = pylab.axes()
        else:           ax = display
        
        
        sorted_index = numpy.argsort(code)
        sorted_data  = [data_spk[i] for i in sorted_index]
           
        ids = []
        spike_times = [] 
   
        for i, d in enumerate(sorted_data):
            ids.extend( [ i for j in range( len( d ) ) ] )
            spike_times.extend(d)
                
               
        if len(spike_times) > 0:
           ax.plot(spike_times, ids, ',', **kwargs)
        ax.set_xlabel('Time (ms)')
        ax.set_ylabel('Neuron #')
        ax.set_ylim(min(ids), max(ids))
    
def make_N_colors(cmap_name, N):
     #! `cmap_name` is a string chosen from the values in cm.datad.keys().You 
     #! will want N to be one minus number of colors you want.
     cmap = matplotlib.cm.get_cmap(cmap_name, N)
     
     
     #! Return list with rgb colors. With return cmap(np.arange(N) a list of 
     #! RGBA values is returned, actually an Nx4 ndarray.
     return cmap(numpy.arange(N))[:, :-1]
  
def IF_curve( id = None, sim_time = None, step_currents = None ):    
    '''
    Function that creates I-F curve
    Inputs:
            id             - id of neuron to use for calculating I-F relation
            sim_time       - lenght of each step current injection in miliseconds
            step_currents  - step currents to inject
    
    Returns: 
            f_isi          - first interspike interval 
            m_isi          - mean interspike interval
            l_isi          - last interspike interval
            trace_isi      - interspike interval traces used to calculate I-F relation  
            
    Examples:
            >> n  = nest.Create('izhik_cond_exp')
            >> sc = [ float( x ) for x in range( 10, 270, 50 ) ]
            >> f_isi, m_isi, l_isi = IF_curve( id = n, sim_time = 500, step_currents = sc ):
    '''
    
    if isinstance( id, int ): id =[ id ] 
    
    traces_isi           = []                                                       # sacve isi traces
    f_isi, m_isi, l_isi  = [], [], []                                               # first, mean and last isi
    dt                   = 2000.0 + sim_time                                        # each I-F trace in ms
    dt_sub               = 2000.0                                                   # time no current is injected ( current is injected 
                                                                                    # between dt_sub and dt )
    T                    = 0                                                        # accumulated simulation time
    
    scg = nest.Create( 'step_current_generator' )                                   # step current generator
    sd  = nest.Create( 'spike_detector' )                                             # spike detector
    nest.SetStatus( sd, { "withgid" : True } )

    rec = nest.GetStatus( id )[0]['receptor_types']

    nest.Connect( scg, id, params = { 'receptor_type' : rec['CURR'] } )             # connect step current generator
    nest.Connect( id,  sd )                                                         # connect spike detector

    
    for curr in step_currents:
        
        nest.SetStatus( scg, params = { 'amplitude_times'  : [ T + dt_sub, T + dt ],       # amplitude times 
                                        'amplitude_values' : [ curr      , 0.     ] } )    # amplitude values
        
        nest.Simulate( dt )                                                         # simulate
        spikes      = nest.GetStatus( sd )[ 0 ][ 'events' ][ 'times' ]              # retrieve spikes
        spike_index = numpy.nonzero( spikes > T + dt_sub, )[ 0 ]                 # retrieve index in spikes of spikes after dt*i +dt_sub  
        
        isi = numpy.diff( spikes[spike_index[ 0 ] : -1 ] , n =  1)                  # get isi 
        f_isi.append( isi[ 0 ] )                                                    # retrieve first isi
        m_isi.append( numpy.mean( isi ) )                                           # retrieve mean isi
        l_isi.append( isi[ -1 ] )                                                   # retrieve last isi
    
        traces_isi.append(isi)                                                      # save trace
    
        T += dt                                                                     # update accumulated time
    return traces_isi, f_isi, m_isi, l_isi

def IV_I_clamp( id = None, sim_time = None, step_currents = None, dt_pre = None, dt_post = None ):    
    '''
    Function that creates I-V by injecting hyperpolarizing currents and then 
    measuring steady-state membrane (current clamp). Each trace is preceded and
    followed by 1/5 of the simulation time of no stimulation
    Inputs:
            id             - id of neuron to use for calculating I-F relation
            sim_time       - lenght of each step current stimulation in ms
            step_currents  - step currents to inject
    
    Returns: 
            v_ss           - steady state voltage 
            traces_v       - voltage traces   
            trace_t        - time points from 0 to trace simulation time as many
                             as recorded data points for each voltage trace.
            
    Examples:
            >> n  = nest.Create('izhik_cond_exp')
            >> sc = [ float( x ) for x in range( -300, 100, 50 ) ] 
            >> tr_t, tr_v, v_ss = IV_I_clamp( id = n, sim_time = 500, step_currents = sc ):
    '''
    
    if isinstance( id, int ): id =[ id ] 
    
    if not dt_pre: dt_pre   = 1/5.*sim_time                                            # preceding time without stimulation
    dt_stim  = sim_time                                                 # stimulation time in ms
    if not dt_post: dt_post  = 1/5.*sim_time                                            # post time without stimulation
    dt_dead  = 1000.0                                                   # dead time between runs such that membrane potential settles
    
    traces_v = []                                                       # save voltage traces
    v_ss     = []                                                       # steady state voltage
                                                                                    # between dt_sub and dt )
    dt       = dt_pre + dt_stim + dt_post + dt_dead
    T        = 0                                                        # accumulated simulation time
    
    scg = nest.Create( 'step_current_generator' )                                   # step current generator
    mm  = nest.Create( 'multimeter' )                                               # voltage detector
    dt_mm = 0.1
    nest.SetStatus( mm,  {'interval': .1, 'record_from': ['V_m']})
    
    rec = nest.GetStatus( id )[0]['receptor_types']

    nest.Connect( scg, id, params = { 'receptor_type' : rec['CURR'] } )             # connect step current generator
    nest.Connect( mm,  id)                                                          # connect miltimeter
    
    for curr in step_currents:
        
        nest.SetStatus( scg, params = { 'amplitude_times'  : [ T + dt_dead + dt_pre, T + dt_dead + dt_pre + dt_stim ],       # amplitude times 
                                        'amplitude_values' : [ curr,       0.,                  ] } )    # amplitude values
        
        nest.Simulate( dt )                                                         # simulate
        print nest.GetStatus( mm )[ 0 ]
        V_m          = nest.GetStatus( mm )[ 0 ][ 'events' ][ 'V_m' ]               # retrieve membrane potential trace
        times        = nest.GetStatus( mm )[ 0 ][ 'events' ][ 'times' ]             # retrieve recording times
        v_ss_index   = numpy.nonzero( times == T +  dt_dead + dt_pre + dt_stim, )[ 0 ]         # retrieve index in spikes of spikes after dt*i +dt_sub  
        
        v_ss.append( V_m[ v_ss_index[ 0 ] ] )                                       # steady state voltage
        
        traces_v.append( V_m[ ( T + dt_dead)/dt_mm : (T + dt)/dt_mm ] )                                        # voltage trace from T : T + dt)
        print V_m[ T : T + dt ]
        T += dt                                                                     # update accumulated time
    
    print len( traces_v[0] )
    trace_t              = numpy.linspace(0, dt, len( traces_v[0] ) ) 
    
    return trace_t, traces_v, v_ss

def PRC():
    '''
    Function that calculates excitatory and inhibitory phase response curves
    
    '''
    
    if isinstance( id, int ): id =[ id ] 
    
    dt_pre   = 1/5.*sim_time                                            # preceding time without stimulation
    dt_stim  = sim_time                                                 # stimulation time in ms
    dt_post  = 1/5.*sim_time                                            # post time without stimulation
    dt_dead  = 1000.0                                                   # dead time between runs such that membrane potential settles
    
    traces_v = []                                                       # save voltage traces
    v_ss     = []                                                       # steady state voltage
                                                                                    # between dt_sub and dt )
    dt       = dt_pre + dt_stim + dt_post + dt_dead
    T        = 0                                                        # accumulated simulation time
   
def inh_poisson_spikes(rate, t, t_stop, n_rep=1 ,seed=None):   
    '''
    Returns a SpikeTrain whose spikes are a realization of an inhomogeneous 
    poisson process (dynamic rate). The implementation uses the thinning 
    method, as presented in the references (see neurotools).
    
    From Professor David Heeger 2000 - "Poisson Model of Spike Generation" 
    
    Generating Poisson Spike Trains
    There are two commonly used procedures for numerically generating Poisson 
    spike trains. The first approach is based on the approximation in Eq. 2 
   
    P{1 spike during the interval(t-dt/2,t+dt/2)}=r(t)dt                     (2) 
    
    for the probability of a spike occurring during a short time interval. For 
    the homogeneous Poisson process, this expression can be rewritten
    (removing the time dependence) as
    
    P{1 spike during dt} ~ rdt
    
    This equation can be used to generate a Poisson spike train by first 
    subdividing time into a bunch of short intervals, each of duration dt. Then 
    generate a sequence of random numbers x[i] , uniformly distributed between 
    0 and 1. For each interval, if x[i] <=rdt, generate a spike. Otherwise, 
    no spike is generated. This procedure is appropriate only when dt is very 
    small, i.e, only when rdt<<1. Typically, dt = 1 msec should suffice. The 
    problem with this approach is that each spike is assigned a discrete 
    time bin, not a continuous time value.
    
    The second approach for generating a homogeneous Poisson spike train, 
    that circumvents this problem, is simply to choose interspike intervals 
    randomly from the exponential distribution. Each successive spike time is 
    given by the previous spike time plus the randomly drawn interspike 
    interval . Now each spike is assigned a continuous time value instead of a 
    discrete time bin. However, to do anything with the simulated spike train 
    (e.g., use it to provide synaptic input to another simulated neuron), it is
     usually much more convenient to discretely sample the spike train
    (e.g., in 1 msec bins), which makes this approach for generating the spike 
    times equivalent to the first approach described above.
    
    I use neruotools 
    
    
    Inputs:
             rate   - an array of the rates (Hz) where rate[i] is active on interval 
                     [t[i],t[i+1]] (if t=t.append(t_stop)
            t      - an array specifying the time bins (in milliseconds) at which to 
                     specify the rate
            t_stop - length of time to simulate process (in ms)
            n_rep  - number times to repeat spike pattern  
            seed   - seed random generator
    
    Returns: 
            spike  - array of spikes
            
    Examples:
          spike_train_dyn = misc.inh_poisson_generator(rate = np.array([50., 80., 30.]),
                                                  t = [0., 1000., 2000.],
                                                  t_stop = 2.5,
                                                  array = False)

      
    '''
    t=numpy.array(t)
    rate=numpy.array(rate)

    times=[]
    rates=[]
    for i in range(n_rep):
        times=concatenate((times,t + t_stop*i/n_rep) ,1)
        rates=concatenate((rates,rate),1)

    stgen=StGen(seed=seed) 
    spikes=stgen.inh_poisson_generator(rate = rates, t=times, t_stop=t_stop, array=True)
    
    return spikes
    
    
    
def read_f_name(data_path, contain_string = None):
    '''
    read in file names in a directory. Files can be filtered by extension
    with extension
    '''
    
    # Get file names in directory    
    file_list=os.listdir(data_path)
    
    # Get files with certain extension 
    if contain_string:
        new_file_list=[]
        for file_path in file_list:
            
            if contain_string in file_path.split('/')[-1]:
                
                new_file_list.append(file_path)
        
        file_list = new_file_list        
            
    return file_list        
    
def wrap(ids, start, stop):
    if stop > len(ids): return ids[ start : ] + ids[ : stop - len(ids) ] 
    else:               return ids[ start : stop ] 
