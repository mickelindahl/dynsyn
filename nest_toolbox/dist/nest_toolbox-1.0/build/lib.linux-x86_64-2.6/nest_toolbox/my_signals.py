# coding:latin
'''
Mikael Lindahl 2010


File:
my_signals

Wrapper for signals module in NeuroTools.  Additional functions can be added 
here. In self.list neurotool object is created through which it can be 
directly accessed 


A collection of functions to create, manipulate and play with analog signals and
spikes 

Classes
-------
VmList           - AnalogSignalList object used for Vm traces
ConductanceList  - AnalogSignalList object used for conductance traces
CurrentList      - AnalogSignalList object used for current traces
SpikeList        - SpikeList object used for spike trains
'''


import numpy
import pylab


# Import StandardPickleFile for saving of spike object
from NeuroTools.io import StandardPickleFile

# For setting number of ticks
from matplotlib.ticker import MaxNLocator


import plot_settings as ps

from NeuroTools import signals
from NeuroTools.signals import ConductanceList
from NeuroTools.signals import CurrentList
from NeuroTools.signals import VmList
from NeuroTools.signals import SpikeList
from NeuroTools.plotting import get_display

class MyConductanceList(ConductanceList):
    '''
    MyConductanceList(signals, id_list, dt, t_start=0, t_stop=None, dims=None )
    inherit from ConductanceList which has base class AnalogSignalList.  
    
    Arguments:
        signal  - the vector with the data of the AnalogSignal
        dt      - the time step between two data points of the sampled analog signal
        t_start - begining of the signal, in ms.
        t_stop  - end of the SpikeList, in ms. If None, will be inferred from the data
    '''
    
    def __init__(self, signals, id_list, dt, t_start=0, t_stop=None, dims=None ):
        ''' 
        Constructor 
        
        Inherited attributes:
        self.t_start        = float(t_start)
        self.t_stop         = t_stop
        self.dt             = float(dt)
        self.dimensions     = dims
        self.analog_signals = {}
        
        New attributes:
        self.ids = sorted( id_list )     # sorted id list
        
        '''
        
        # About the super() function 
        # (from python programming - Michael Dawson, page 277) 
        # Incorporate the superclass ConductanceList method's functionality. 
        # To add a new attribute ids i need to override the inherited 
        # constructor method from ConductanceList. I also want my new 
        # constructor to create all the attributes from ConductanceList. 
        # This can be done with the function super(). It lets you invoke the 
        # method of a base class(also called a superclass). The first argument 
        # in the function call, 'MyConductanceList', says I want to invoke a 
        # method of the superclass (or base class) of MyConductanceList which 
        # is ConductanceList. The next argument. seöf, passes a reference to 
        # the object so that ConductanceList can get to the object and add 
        # its attributes to it. The next part of the statement __init__(
        # signals, id_list, dt, t_start, t_stop, dims) tells python I want to
        # invoke the constructor method of ConductanceList and a want to pass 
        # it the values of signals, id_list, dt, t_start, t_stop and dims.
        
        # Invoke __init__ of base class ConductanceList which has 
        # AnalogSignalList as base class where the method exist.
        super( MyConductanceList, self ).__init__( signals, id_list, dt, 
                                                   t_start, t_stop, dims)
        # Add new attribute
        self.ids = sorted( id_list )     # sorted id list
        
    def plot(self,  id_list=None, display=True, kwargs={} ):
        """
        Plot all cells in the AnalogSignalList defined by id_list. Right now
        here to exemplify use of super() on super class method.
        
        Inputs:
            id_list - can be a integer (and then N cells are randomly selected) or a 
                      list of ids. If None, we use all the ids of the SpikeList
            display - axes handle
            kwargs  - dictionary contening extra parameters that will be sent to the plot 
                      function
        
        Examples:
            >> z = subplot(221)
            >> aslist.plot(5, display=z, kwargs={'color':'r'})
        """
        
        # Invoke plot function in ConductanceList
        super(MyConductanceList, self).plot(id_list, None, display, kwargs )
                 
    def save(self, userFileName):
        '''
        Save analog list

        Inputs:
            userFileName    - name of file to save
            
        Examples:
            >> userFileName = /home/savename.dat
            >> aslist.save(userFileName)           
        '''
        userFile = StandardPickleFile( userFileName )      # create user file 
       
        # Invoke save function of base class ConductanceList which has 
        # AnalogSignalList as base class where the method exist.
        super(MyConductanceList, self).save( userFile )
                  
class MyCurrentList(CurrentList):
    ''' 
    MyCurrentList(signals, id_list, dt, t_start=0, t_stop=None, dims=None )
    inherit from CurrentList which has base class AnalogSignalList.
    
    Arguments:
        signal  - the vector with the data of the AnalogSignal
        dt      - the time step between two data points of the sampled analog signal
        t_start - begining of the signal, in ms.
        t_stop  - end of the SpikeList, in ms. If None, will be inferred from the data
    '''  
    
    def __init__( self, signals, id_list, dt, t_start=0, t_stop=None, dims=None ):
        ''' 
        Constructor 
        
        Inherited attributes:
        self.t_start        = float(t_start)
        self.t_stop         = t_stop
        self.dt             = float(dt)
        self.dimensions     = dims
        self.analog_signals = {}
        
        New attributes:
        self.ids = sorted( id_list )     # sorted id list
        
        '''
        
        # Invoke __init__ of base class ConductanceList which has 
        # AnalogSignalList as base class where the method exist.
        super( MyCurrentList, self ).__init__( signals, id_list, dt, 
                                                   t_start, t_stop, dims)
        # Add new attribute
        self.ids = sorted( id_list )     # sorted id list
    
    def plot(self, id_list = None, display = True, kwargs = {} ):
        """
        Plot all cells in the AnalogSignalList defined by id_list
        
        Inputs:
            id_list - can be a integer (and then N cells are randomly selected) or a 
                      list of ids. If None, we use all the ids of the SpikeList
            display - if True, a new figure is created. Could also be a subplot
            kwargs  - dictionary contening extra parameters that will be sent to the plot 
                      function
        
        Examples:
            >> z = subplot(221)
            >> aslist.plot(5, display=z, kwargs={'color':'r'})
        """
        
        # Invoke plot function in ConductanceList
        super(MyCurrentList, self).plot(id_list, None, display, kwargs )


    def save(self, userFileName):
        '''
        Save analog list

        Inputs:
            userFileName    - name of file to save
            
        Examples:
            >> userFileName = /home/savename.dat
            >> aslist.save(userFileName)           
        '''
        userFile = StandardPickleFile( userFileName )      # create user file 
       
        # Invoke save function of base class ConductanceList which has 
        # AnalogSignalList as base class where the method exist.
        super(MyCurrentList, self).save( userFile )
         
class MyVmList(VmList):
    ''' 
    MyVmList(signals, id_list, dt, t_start=0, t_stop=None, dims=None )
    inherit from VmList which has base class AnalogSignalList.
    
    Arguments:
        signal  - the vector with the data of the AnalogSignal
        dt      - the time step between two data points of the sampled analog signal
        t_start - begining of the signal, in ms.
        t_stop  - end of the SpikeList, in ms. If None, will be inferred from the data
    '''  
    
    def __init__(self, signals, id_list, dt, t_start=0, t_stop=None, dims=None ):
        ''' 
        Constructor 
        
        Inherited attributes:
        self.t_start        = float(t_start)
        self.t_stop         = t_stop
        self.dt             = float(dt)
        self.dimensions     = dims
        self.analog_signals = {}
        
        New attributes:
        self.ids = sorted( id_list )     # sorted id list
        
        '''
        
        # Invoke __init__ of base class ConductanceList which has 
        # AnalogSignalList as base class where the method exist.
        super( MyVmList, self ).__init__( signals, id_list, dt, 
                                                   t_start, t_stop, dims)
        # Add new attribute
        self.ids = sorted( id_list )     # sorted id list
    
    def plot(self, id_list = None, display = True, kwargs = {} ):
        """
        Plot all cells in the AnalogSignalList defined by id_list
        
        Inputs:
            id_list - can be a integer (and then N cells are randomly selected) or a 
                      list of ids. If None, we use all the ids of the SpikeList
            display - if True, a new figure is created. Could also be a subplot
            kwargs  - dictionary contening extra parameters that will be sent to the plot 
                      function
        
        Examples:
            >> z = subplot(221)
            >> aslist.plot(5, display=z, kwargs={'color':'r'})
        """
        
        # Invoke plot function in ConductanceList
        super(MyVmList, self).plot(id_list, None, display, kwargs )

    
    def raw_data_id_order(self): 
        '''
        Return data matrix with raw data for each id on row ordered in descending 
        id order. 
        '''    

        values = numpy.concatenate( [self.analog_signals[id].signal for id in self.ids]) 
        #ids    = numpy.concatenate( [id*numpy.ones( len( obj.analog_signals[id].signal ), int ) for id in self.ids] ) 

        data = values.reshape( len( self.ids ),-1 )
        
        return data
    
    def save(self, userFileName):
        '''
        Save analog list

        Inputs:
            userFileName    - name of file to save
            
        Examples:
            >> userFileName = /home/savename.dat
            >> aslist.save(userFileName)           
        '''
        userFile = StandardPickleFile( userFileName )      # create user file 
       
        # Invoke save function of base class ConductanceList which has 
        # AnalogSignalList as base class where the method exist.
        super(MyVmList, self).save( userFile )


    def image_plot(self, display = None, kwargs = {}):
        
        if not display: ax = pylab.axes()
        else:           ax = display

        
        data_Vm = numpy.array( self.raw_data() )
        data_Vm = numpy.array( zip( *data_Vm ) )
        data_Vm = data_Vm[ 0 ].reshape( len( self.ids ), -1 )  # Transpose
        image = ax.imshow( data_Vm, origin = 'lower', **kwargs)
        ax.set_xlabel('Membrane potential (mV)')
        ax.set_ylabel('Neuron #')  
        ax.set_aspect(aspect='auto')
        return image
             
class MySpikeList(SpikeList):
    """
    MySpikeList(spikes, id_list, t_start=None, t_stop=None, dims=None)
    
    Inherits from SpikeList
    
    Return a SpikeList object which will be a list of SpikeTrain objects.

    Inputs:
        spikes  - a list of (id,time) tuples (id being in id_list)
        id_list - the list of the ids of all recorded cells (needed for silent cells)
        t_start - begining of the SpikeList, in ms. If None, will be infered from the data
        t_stop  - end of the SpikeList, in ms. If None, will be infered from the data
        dims    - dimensions of the recorded population, if not 1D population
    
    t_start and t_stop are shared for all SpikeTrains object within the SpikeList
    
    Examples:
        >> sl = SpikeList([(0, 0.1), (1, 0.1), (0, 0.2)], range(2))
        >> type( sl[0] )
            <type SpikeTrain>
    
    See also
        load_spikelist
    """
    
    def __init__(self, spikes, id_list, t_start=None, t_stop=None, dims=None):

        # Invoke __init__ of base class ConductanceList which has 
        # AnalogSignalList as base class where the method exist.
        super( MySpikeList, self ).__init__( spikes, id_list, t_start, t_stop, dims)
        
        # Add new attribute
        self.ids = sorted( id_list )     # sorted id list
    
    
    def firing_rate_sliding_window(self, bin=100, display=True, id_list=[], step=1, stop=None, kwargs={}):
        ''' 
        Calculate spike rates at ``sample_step`` using s sliding rectangular 
        window.
        
        Arguments:
        bin           Bin size of sliding window
        display       If True, a new figure is created. Could also be a subplot
        id_list       List of ids to calculate firing rate for
        step          Step size for moving sliding window (ms)
        stop          End of spike train
        kwargs        Additional plot arguments
          
        Here the window is centered over over each time point at sampe_step with 
        window size equalling bin_size. Takes spike times, number of 
        neurons no_neurons, bin_size and sample_step as input argument.
        '''  
        
        ax = get_display(display)
        
        if stop is None: stop = self.t_stop
        if not any(id_list): id_list = self.ids    
        
        spikes = {}     # dictionary with spike times   
        for id in id_list: 
            spikes[ id ] = self.spiketrains[ id ].spike_times.copy()                   
        
        n  = int( (stop - bin )/step )          # number of windows
        f = bin/2                               # first window at bin/2
        l = step*n + bin/2                      # last window at bin*n - bin/2 
 
        # Time axis for sliding windows, n_win + 1 due to end points 
        # [ 0, 1, 2] -> two windows and tre timings
        timeAxis = numpy.linspace(f, l, n + 1)           
        
        
        firingRates = []                          # sliding time window data  
        dataSpk   = []
        
        #! Calculate number of spikes in ms bins and then sliding window rates 
        for id in id_list:
                    
            spk = spikes[id] 
            dataSpk.append(spk)

            i = 0
            j = bin/2
            j_max  = stop
            rates = []
                 
            #! For each ms i in 0 to stop at step intervals 
            for tPoint in timeAxis:                                    
                
                sum = numpy.sum( ( tPoint - bin/2 <= spk ) * ( spk < tPoint + bin/2 ) )
                rates.append( 1000.0 * sum/float( bin ) )
                  

            firingRates.append(rates)
        
        firingRates = numpy.array(firingRates)        # Convert to numpy array
   
        meanFiringRates = numpy.mean( firingRates, axis = 0 )
       
        ax.plot(timeAxis, meanFiringRates,**kwargs)
        ax.set_xlabel( 'Time (ms)' )
        ax.set_ylabel( 'Frequency (spike/s)' )
        
        return timeAxis, firingRates, dataSpk   
        
 
    def image_firing_rate_slide_window(self, bin=100, id_list=[], step=1, stop=None, display=True, kwargs={}):
        '''
        Function that create an image with bin size sliding window firing rate as 
        calculated at step intervals. kwargs - dictionary contening extra 
        parameters that will be sent to the plot function    
        
        Arguments:
        bin           Bin size of sliding window
        display       If True, a new figure is created. Could also be a subplot
        id_list       List of ids to calculate firing rate for
        step          Step size for moving sliding window (ms)
        stop          End of spike train
        kwargs        Additional plot arguments
        
        '''
        ax = get_display(display)
        
        if not any(id_list): id_list = self.ids
        
        t, r, spk = self.firing_rate_slide_window(bin, display, id_list, 
                                                    step, stop, kwargs)
               
        kwargs.update( { 'origin' : 'lower', } )
        image = ax.imshow(r, extent=[t[0],t[-1],self.ids[0],self.ids[-1]], 
                          **kwargs)
        ax.set_xlabel('Time (ms)')
        ax.set_ylabel('Neuron #')
        ax.set_aspect(aspect='auto')

        return image

    def spike_histogram_n_rep(self, bin=10, id_list=[], n_rep=1, normalized=True):
        '''
        Generate an array with all the spike_histograms of all the SpikeTrains. 
        Each time series should be n_rep repeated simulations. Histogram is then
        generated for one simulation period repeated n_rep times. It is thus
        possible to supply several similar simulations and get the 
        mean firing statistics for all of them. 
        
        Inputs:
            bin        - the time bin used to gather the data
            id_list    - list with ids to use for creation of histogram
            n_rep      - number of experimental repetitions with same 
                         stimulation. E.g. if n_runs=3
            normalized - if True, the histogram are in Hz (spikes/second), 
                         otherwise they are in spikes/bin

        
        In neurophysiology, this is similar to a peristimulus time histogram 
        and poststimulus time histogram, both abbreviated PSTH or PST histogram, 
        are histograms of the times at which neurons fire. These histograms are 
        used to visualize the rate and timing of neuronal spike discharges in 
        relation to an external stimulus or event. The peristimulus time 
        histogram is sometimes called perievent time histogram, and 
        post-stimulus and peri-stimulus are often hyphenated.
        
        To make a PSTH, a spike train recorded from a single neuron is aligned 
        with the onset, or a fixed phase point, of an identical stimulus 
        repeatedly presented to an animal. The aligned sequences are 
        superimposed in time, and then used to construct a histogram.[2] 
        
        Construction procedure
        For each spike train i
        
        1. Align spike sequences with the onset of stimuli that repeated n 
           times. For periodic stimuli, wrap the response sequence back to 
           time zero after each time period T, and count n (n_rep) as the 
           total number of periods of data.
        2. Divide the stimulus period or observation period T into N bins 
           of size Δ seconds ( bin ).
        3. Count the number of spikes k_i_j from all n sequences that 
           fall in the bin j.
        4. Calculate  i_j given by k_i_j/( n*Δ ) in units of estimated 
           spikes per second at time j * Δ.

        The optimal bin size Δ is a minimizer of the formula, (2k-v)/Δ2, where
        k and v are mean and variance of k_i. [3]  
                                                                       
        '''
        
        if not any(id_list): id_list = self.ids
        
        spkList = self.id_slice( id_list )      # Retrieve spike trains 
        
        # Short cuts
        start    = spkList.t_start 
        period   = spkList.t_stop/n_rep 
        
        spike_train_hist = {}       # Create data table for histogram
           
        # Add spikes to data
        for id, spk in spkList.spiketrains.iteritems():
            
            
            # First take modulus t_period of spike times ending up with a vector
            # with spike times between 0 and period    
            spikes_mod_period = spk.spike_times % period                                        
            
            # Create sequence of bins edges, histogram between 0 and period
            bin_sequence = numpy.arange( start, period + bin, bin )    
            
            # Create histogram over bin sequence and convert to spikes per bin 
            hist_n_rep, edges = numpy.histogram( spikes_mod_period, bin_sequence)
            hist_n_rep        = hist_n_rep/float( n_rep )     
        
            spike_train_hist[ id ] = hist_n_rep     # Add to dictionary                                              
 
        # Create numpy array with all spike trains histograms   
        histograms = numpy.array( spike_train_hist.values() )                                      
        
        # Convert data to spike rates
        if normalized: histograms *= 1000.0/bin                                            

        # Center time axis over bin mid points 
        timeAxis = numpy.arange(start+bin/2., period, bin)
        return timeAxis, histograms
                                           
    def image_spike_histogram(sself, bin=10, display=True, id_list=[], n_rep=1, normalized=True, kwargs = {} ):
        '''
        Plot an image of all the spike_histograms generated by 
        spike_histogram_n_rep
    
        
        Arguments:
        bin        - the time bin used to gather the data
        display    - If True, a new figure is created. Could also be a subplot
        id_list    - list with ids to use for creation of histogram
        n_rep      - Number of experimental repetitions with same stimulation. 
                     E.g. if n_runs=3 then it is assumed that three 
                     experimental runs is recorded and the data will be 
                     sliced up into three time intervals.
        normalized - if True, the histogram are in Hz (spikes/second), 
                     otherwise they are in spikes/bin
        kwargs     . Additional plot arguments
        '''        
        
        ax = get_display(display)
        
        timeAxis, histograms = self.spike_histogram_n_rep( bin, id_list, 
                                                            n_rep, normalized )
             
        kwargs.update( { 'origin' : 'lower', } )
        image = ax.imshow( histograms, **kwargs )
        ax.set_xlabel( 'Time (ms)'     )
        ax.set_ylabel( 'Neuron #'      )
        ax.set_aspect( aspect = 'auto' )
        
        
        n_points=len(timeAxis)
        xticks = numpy.arange( 0, n_points, n_points*0.2)
        xticklabels = numpy.arange( 0, timeAxis[-1], timeAxis[-1]*0.2)
        ax.set_xticks( xticks )
        ax.set_xticklabels( [ str( l ) for l in xticklabels ] )
        
        n_ids=len(id_list)
        yticks =  numpy.arange( 0, n_ids, n_ids*0.2)
        ax.set_yticks( yticks)
        ax.set_yticklabels( id_list[ yticks ] )
    
        return image
           
    def plot_spike_histogram(self, bin=10, display=True, id_list=[], n_rep=1, normalized=True, kwargs = {} ):    
        '''
        Plot an mean spike histogram for for all ids generated by 
        spike_histogram_n_rep
        
        Arguments:
        bin        - the time bin used to gather the data
        display    - If True, a new figure is created. Could also be a subplot
        id_list    - list with ids to use for creation of histogram
        n_rep      - Number of experimental repetitions with same stimulation. 
                     E.g. if n_runs=3 then it is assumed that three 
                     experimental runs is recorded and the data will be 
                     sliced up into three time intervals.
        normalized - if True, the histogram are in Hz (spikes/second), 
                     otherwise they are in spikes/bin
    
        '''        
        ax = get_display(display)
        
        timeAxis, histograms = self.spike_histogram_n_rep( bin, id_list, 
                                                            n_rep, normalized )
              
        ax.bar(left=t_vec, height=data, width=0.8, 
               bottom=None, hold=None, **kwargs)
              
    def firing_rate(self, bin=200, display=True, id_list=[], n_rep=1, kwargs={}):          
        '''
        Plot all the instantaneous firing rates along time (in Hz) 
        over all ids in id_list from histogram generated by spike_histogram_n_rept.
        
        Inputs:
            bin        - the time bin used to gather the data
            display    - If True, a new figure is created. Could also be a subplot
            id_list    - list with ids to use for creation of histogram
            n_rep      - Number of experimental repetitions with same stimulation. 
                     E.g. if n_runs=3 then it is assumed that three 
                     experimental runs is recorded and the data will be 
                     sliced up into three time intervals.
            kwargs     - dictionary contening extra parameters that will be sent to the plot 
                         function
        '''
        
        ax = get_display(display)
        
        timeAxis, histograms = self.spike_histogram_n_rep( bin, id_list, 
                                                            n_rep, True )
        
        firingRates = numpy.mean( histograms, axis = 0 )
       
        ax.plot(timeAxis, firingRates, **kwargs)
        ax.set_xlabel( 'Time (ms)' )
        ax.set_ylabel( 'Frequency (spike/s)' )
        
          
    def raster_plot_cluster(self, id_list=[], t_start=None, t_stop=None, display=True, clusters = [], kwargs={}):
        """
        (functional?) Generate a raster plot for the SpikeList in a subwindow of interest,
        defined by id_list, t_start and t_stop. 
        
        Inputs:
            id_list - can be a integer (and then N cells are randomly selected) or a list of ids. If None, 
                      we use all the ids of the SpikeList
            t_start - in ms. If not defined, the one of the SpikeList object is used
            t_stop  - in ms. If not defined, the one of the SpikeList object is used
            display - if True, a new figure is created. Could also be a subplot
            clusters - vector containing code for cluster belonging of each spike 
                       train opatain from clustering analysis. Plotted accordingly.
            kwargs  - dictionary contening extra parameters that will be sent to the plot 
                      function
        
        Examples:
            >> z = subplot(221)
            >> spikelist.raster_plot(display=z, kwargs={'color':'r'})
        
        See also
            SpikeTrain.raster_plot
        """
        
        ax = get_display(display)
        
        spk = self.list
        
        if t_start is None: t_start = spk.t_start
        if t_stop is None:  t_stop  = spk.t_stop
        
        ids, spike_times = spk.convert(format="[ids, times]")
        idx = numpy.where((spike_times >= t_start) & (spike_times <= t_stop))[0]
        
        sorted_index = numpy.argsort( clusters )                            # Sort spike trains accoringly to clusters
        for i, id in enumerate(self.ids):
            ids[ids==id] = -self.ids[ sorted_index[ i ] ]
        ids = abs(ids)
                 
        
        if len(spike_times) > 0:
            ax.plot(spike_times, ids, ',', **kwargs)
        xlabel = "Time (ms)"
        ylabel = "Neuron #"
        set_labels(ax, xlabel, ylabel)
        
        min_id = numpy.min(spk.id_list())
        max_id = numpy.max(spk.id_list())
        length = t_stop - t_start
        set_axis_limits(ax, t_start-0.05*length, t_stop+0.05*length, min_id-2, max_id+2)
        pylab.draw()
                   

    def save(self, userFileName):
        '''
        Save spike list

        Inputs:
            userFileName    - name of file to save
            
        Examples:
            >> userFileName = /home/savename.dat
            >> aslist.save(userFileName)           
        '''
        userFile = StandardPickleFile( userFileName )      # create user file 
       
        # Invoke save function of base class ConductanceList which has 
        # AnalogSignalList as base class where the method exist.
        super(MyVmList, self).save( userFile )
       
def load(userFileName, dataType):  
    
    """
    load(userFileName, dataType)
    Convenient data loader for results saved as NeuroTools StandardPickleFile. 
    Return the corresponding NeuroTools object. Datatype argument may become 
    optionnal in the future, but for now it is necessary to specify the type 
    of the recorded data. To have a better control on the parameters of the 
    NeuroTools objects, see the load_*** functions.
    
    Inputs:
        userFileName - the user file name
        datatype - A string to specify the type od the data in
                    's' : spikes
                    'g' : conductances
                    'v' : membrane traces
                    'c' : currents
    """
    userFile=StandardPickleFile(userFileName)
       
    if dataType in ('s', 'spikes'):
        return signals.load_spikelist(userFile)
    elif dataType == 'v':
        
        #  Need t_start to be None, othervice NeuroTools overwrite loaded 
        # t_start with default value 0
        return signals.load_vmlist(userFile, t_start=None)  
    elif dataType == 'c':
        return signals.load_currentlist(userFile, t_start=None)
    elif dataType == 'g':
        return signals.load_conductancelist(userFile,t_start=None)
    else:
        raise Exception("The datatype %s is not handled ! Should be 's','g','c' or 'v'" %datatype)
    

    
       