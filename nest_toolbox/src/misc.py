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
import pickle
from NeuroTools.stgen import StGen
from numpy import array, concatenate
from numpy.random import poisson 
import scipy.signal as signal
from copy import deepcopy
import time
from scipy.cluster.vq import vq, kmeans, kmeans2, whiten


def adjust_limit(limit, percent_x=0.03, percent_y=0.04):
    ''' Adjust limit by withdrawing percent of diff from start and
        adding percent of diff to start.  
    
        Input:
            limit      - the limit to manipulate
            percent    - percent margin to add
    
    '''
    
    x=limit[0]
    y=limit[1]
    diff=y-x
    coords=[x-diff*percent_x, y+diff*percent_y]
    
    return coords
    

'''
Method for method to evenly split up a colormap into an RGB colors array
'''
def autocorrelation(spikes, bin=1, max_time=1000):
    N=len(spikes)
    
    N_order_isi=[]
    N_order_isi_histogram=[] # save isi histogram order 1 to n, see Perkel 1967
    top_bins=numpy.arange(1, max_time, bin)
    bottom_bins=numpy.arange(-max_time+1, 0, bin)
    middle_bins=array([bottom_bins[-1],0.001,0.001,top_bins[0]])
    bins=array([])
    bins=numpy.append(bins, bottom_bins)
    bins=numpy.append(bins, middle_bins)
    bins=numpy.append(bins, top_bins)
    isi_matrix=spikes-numpy.transpose(array([spikes]))*numpy.ones([1,N])
    
    # Zero order isi is on the diagonal, first on the first sub diagonal 
    # above an below the diagonal third on the next diagonals, and n 
    # order is the right top value and left bottom value.
    hist_autocorr,xaxis =numpy.histogram(isi_matrix, bins=bins, normed=False)
    
    # Set zero interval bin count to zero
    hist_autocorr[hist_autocorr==max(hist_autocorr)]=0
    
    hist_autocorr=hist_autocorr/float((N*bin))*1000
    
    xaxis=xaxis[0:-1]+bin/2
    
    return hist_autocorr, xaxis

   
    
def convolve(binned_data, bin_extent, kernel_type, axis=0, single=False):
    ''' Convolve data with a specific kernel. (Low pass filtering)
        Inputs:
            binned_data - binned spike count
            bin_extent  - extent of kernel
            kernel      - kernel to filter by    
    '''
    
    if kernel_type=='triangle':
        kernel=numpy.arange(1,numpy.ceil(bin_extent/2.)+1)
        kernel=numpy.append(kernel, numpy.arange(numpy.floor(bin_extent/2.),0,-1))
        kernel=kernel/float(numpy.sum(kernel))
    
    if kernel_type=='rectangle':
        kernel=numpy.ones(bin_extent+1)
        kernel=kernel/float(numpy.sum(kernel))
    
    conv_data=[]

    if not single:
        for i in range(binned_data.shape[0]):
            conv_data.append(numpy.convolve(kernel, binned_data[i,:], mode='same'))     

    else:
        
        conv_data.append(numpy.convolve(kernel, binned_data, mode='same'))     

    conv_data=numpy.array(conv_data)
    
    return conv_data
    

def dict_merge(a, b):
    '''recursively merges dict's. not just simple a['key'] = b['key'], if
    both a and bhave a key who's value is a dict then dict_merge is called
    on both values and the result stored in the returned dictionary.'''
    if not isinstance(b, dict):
        return b
    result = deepcopy(a)
    for k, v in b.iteritems():
        if k in result and isinstance(result[k], dict):
                result[k] = dict_merge(result[k], v)
        else:
            result[k] = deepcopy(v)
    return result           

def sigmoid(p,x):
        x0,y0,c,k=p
        y = c / (1 + np.exp(-k*(x-x0))) + y0
        return y
    
def fit_sigmoid():
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.optimize


    def residuals(p,x,y):
        return y - sigmoid(p,x) 

    def resize(arr,lower=0.0,upper=1.0):
        arr=arr.copy()
        if lower>upper: lower,upper=upper,lower
        arr -= arr.min()
        arr *= (upper-lower)/arr.max()
        arr += lower
        return arr
    
    
    
    
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
     '''       

    times=[]
    rates=[]
    t=numpy.array(t)
    for i in range(n_rep):
        times=concatenate((times,t + t_stop*i/n_rep) ,1)
        rates=concatenate((rates,rate),1)

    stgen=StGen(seed=seed) 
    spikes=stgen.inh_poisson_generator(rate = rates, t=times, t_stop=t_stop, array=True)
    
    return spikes   
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


def fano_factor(d):
    
    var=numpy.var(d, axis=0)
    mean=numpy.mean(d, axis=0)
    ff=var/mean
    return ff
    

def mutual_information(count_sr_list):
    '''
    
    Arguments:
        p_sr        - list with the joint probability
                      count of stimulus and response.
                      Each position in count_sr_list[i] contains a matrix 
                      n x m matrix where n is number of stimulus 
                      and m is the number of responses.

    '''    
    mi=[]
        
    for p_sr in count_sr_list:
        p_sr=numpy.array(p_sr)/float(numpy.sum(numpy.sum((p_sr)))) # Probabilities needs to sum to one
        p_s=numpy.sum(p_sr, axis=1) #Marginal distribution p_s
        p_r=numpy.sum(p_sr, axis=0) #Marginal distribution p_r
        
        
        
        tmp=0
        for i in range(p_sr.shape[0]):
            for j in range(p_sr.shape[1]):
                if p_sr[i,j]>0:
                    tmp+=p_sr[i,j]*numpy.log2(p_sr[i,j]/(p_s[i]*p_r[j]))
                if numpy.isnan(tmp):
                    print tmp
        mi.append(tmp)
    
    return numpy.array(mi)

def timer(t=0):
    if not t:
        return time.time() 
    
    else:
        stop = time.time() 
    
        s = stop - t
        m = s // 60
        s = s - m*60
        return 'Simulation time: %i minutes, %i seconds' % ( m, s )    
            

def text_save(data, fileName):
    '''
    
    Arguments:
        data        - data to be pickled
        fileName    - full path or just file name
    '''
    f=open(fileName, 'w')
    f.write(data)
    f.close()
    
def text_load(fileName):
    '''
    
    Arguments:
        data        - data to be pickled
        fileName    - full path or just file name
    '''
    f=open(fileName, 'r')
    data=f.read()
    f.close()
    return data

def pickle_save(data, fileName):
    '''
    
    Arguments:
        data        - data to be pickled
        fileName    - full path or just file name
    '''
    try:
        f=open(fileName, 'wb') #open in binary mode
    except:
        parts=fileName.split('/')
        os.mkdir('/'.join(parts[0:-1]))    
        f=open(fileName, 'wb')    
    
    pickle.dump(data, f)
    f.close()

def pickle_save_groups(group_list, fileName):
    fileName=fileName+'-'+str(nest.Rank())
    pickle_save(group_list, fileName)
    
def pickle_load(fileName):
    '''
    
    Arguments:
        fileName    - full path or just file name
    '''
    f=open(fileName, 'rb') # make sure file are read in binary mode
    data=pickle.load(f)
    f.close()
    return data


def pickle_load_groups(file_name):
    path='/'.join(file_name.split('/')[0:-1])
    name=file_name.split('/')[-1]
    fileNames = read_f_name( path, contain_string=name )
    
    check_parts=pickle_load(path+'/'+fileNames[0])
    parts=[]   
    for name in fileNames:
        if isinstance(check_parts[0],list):
            parts.append(pickle_load(path+'/'+name))
        else:
            # Put it in a list
            parts.append([pickle_load(path+'/'+name)])
    
    # To put merge parts groups in 
    groups=parts[0]   

    # Find what data is recorded from    
    recorded_from=parts[0][0][0].signals.keys() 
  
    # Iterate over groups
    for j in range(len(parts[0])):
        for k in range(len(parts[0][0])) :   
            group=parts[0][j][k]
                       
            # For each recordable
            for recordable in recorded_from:                                           
                
                if 'spikes' == recordable:
                    
                    spk_signal=group.signals[recordable]
                    
                    # Iterate over parts and merge spike data
                    # Only happens for mpi data
                    
                    for i_part in range(1, len(parts[1:])+1):
                        group_k = parts[i_part][j][k]
                        add_spk_signal = group_k.signals[recordable]
                        spk_signal.merge(add_spk_signal)
                        
                    group.signals[recordable] = spk_signal
                
                # Must be analog signal    
                else:
                    
                    ag_signal=group.signals[recordable]
                    
                    # Iterate over parts and merge analog data
                    # Only happens for mpi data
                    for i_part in range(1, len(parts[1:])+1):
                        group_k = parts[i_part][j][k]
                        ag_signal_k = group_k.signals[recordable]
                        
                        for id, signal in ag_signal_k.analog_signals.iteritems():
                            ag_signal.append(id, signal)
                    
                    group.signals[recordable] = ag_signal
                    
                groups[j][k]=group
    
    return groups         
            
        

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
    
    
def slice_line(line, xlim=False, ylim=False):
    
    if xlim:
        x=line.get_xdata()
        y=line.get_ydata()
        
        y=y[x>=xlim[0]]
        x=x[x>=xlim[0]]
        
        y=y[x<=xlim[1]]
        x=x[x<=xlim[1]]
        line.set_xdata(x)
        line.set_ydata(y)
        
    if ylim:
        x=line.get_xdata()
        y=line.get_ydata()
        
        x=x[y>=ylim[0]]
        y=y[x>=ylim[0]]
        
        x=x[y<=ylim[1]]
        y=y[y<=ylim[1]]    
        
        line.set_ydata(y)
        line.set_ydata(x)



def time_resolved_rate(binned_data, bin_extent, kernel_type, res):
    '''
        The rate of at each bin is calculated by a kernel which uses 
        nearby data points to estimate the bin rate. 
                
        Inputs:
            binned_data - binned spike data, n x m where
                          n is number of traces and m is 
                          number of samples
            bin_extent  - extent of kernel
            kernel      - kernel to filter by
            res         - binned data resolution
    
    '''   
    
    if kernel_type=='triangle':
        kernel=numpy.arange(1,numpy.ceil(bin_extent/2.)+1)
        kernel=numpy.append(kernel, numpy.arange(numpy.floor(bin_extent/2.),0,-1))
        kernel=kernel/float(numpy.sum(kernel))
    
    rate_traces=[]
    for i in range(binned_data.shape[0]):
        rate_traces.append(signal.lfilter(kernel, 1/res, binned_data[i,:])*1000.0)     
    rate_traces=numpy.array(rate_traces)
    return rate_traces
    
def wrap(ids, start, stop):
    if stop > len(ids): return ids[ start : ] + ids[ : stop - len(ids) ] 
    else:               return ids[ start : stop ] 
