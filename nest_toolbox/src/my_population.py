''''
Mikael Lindahl August 2011


 Creates a group of neurons and attaches recording devices spike detector and
 multimeter. 

 Usage:

 Declare a mygroup: 
   mygroup = MyGroup(models = 'iaf_neuron', n=1, params = {}, mm_dt=1.0, sname='', spath='', sd=True, mm=True)


 When simulation is finished conductance, current and voltage lists can be
 obtained by calling the function get_signal() with parameter dataype, 
 ('g' for conductance, 'c' for current, 'v' for voltage and 's' or 'spike' for
 spike date) and my_nest recordable type (except for spike data). With 
 save_signal() and load_signal() data can instead be written and read from disk.
 
 Function can take several modes and number of neurons per model. Then the models
 has to have the same recordables. Also if parameters is provided in params these
 are set for all the models. 
'''

from NeuroTools import signals 

import numpy
import os
import random # Random generator
import my_signals # Own wrapper class for signal module in Neurotools
from my_signals import MyConductanceList, MyCurrentList, MyVmList, MySpikeList
import my_nest
import my_topology
import misc
import copy
from numpy.random import random_integers

# Really scary to start using randomstate. Funny stuff can happend.
# Just now by restarting eclips I got rid of oscillations in the network... 
#from numpy.random import RandomState
#random_integers  = RandomState(3).random_integers 
        
class MyGroup(object):
    '''
    MyGroup(self, model = 'iaf_neuron', n=1, params = {}, mm_dt=1.0, sname='', spath='', sd=True, mm=True)
    
    Arguments:
        model      my_nest model type, can be a list
        n           number of model to create, can be a list
        params      common parameters for model to be set
        mm_dt       multimeter recording precision
        sname       file basename
        spath       Path to save file at 
        sd          boolean, True if spikes should me recorded
        mm          boolean, True if mulitmeter should record  
        
    ''' 
    def __init__(self, model = 'iaf_neuron', n=1, params = {}, mm_dt=1.0, 
                 sname='', spath='', sname_nb=0, sd=False, sd_params={},
                 mm=False, record_from=[], ids=[] ):
        '''
        Constructor
        
        Arguments:
            ids         provide ids if nodes already created
            model       my_nest model type, can be a list
            n           number of model to create, can be a list
            params      common parameters for model to be set
            mm_dt       multimeter recording precision
            sname       file basename
            spath       Path to save file at 
            sd          boolean, True if spikes should me recorded
            mm          boolean, True if mulitmeter should record  
        '''         
        
        self.connections    = {}        # Set after network has been built with FindConnections
        self.ids            = []
        self.local_ids      = []
        self.mm             = []        # Id of multimeter
        self.mm_dt          = 0         # Recording interval multimeter
        self.model          = model
        self.params         = []
        self.record_from    = []
        self.receptor_types = {}
        self.recordables    = {}
        self.sd             = []        # Id of spike detector
        self.sd_params      = sd_params
        self.sname_nb       = sname_nb        # number for sname string
        self.sname          = ''        # Specific file basename
        self.spath          = ''        # Path to save file at 
        self.signals        = {}        # dictionary with signals for current, conductance, voltage or spikes    
        
        
        if not self.sname: 
            self.sname=model+'-'+str(sname_nb)+'-'
        else:
            self.sname=self.snam+'-'+str(sname_nb)+'-'
            
        # If no spath is provided current path plus data_tmp is set to 
        # spath.
        if spath is '': 
            self.spath = os.getcwd()+'/output_tmp'    
        else: 
            self.spath=spath
                    
        # Create save dir if it do not exist
        try: 
            msg=os.system('mkdir ' + self.spath + ' 2>/dev/null' )                                   
        except: 
            pass    
        
        if ids:
            self.ids=ids             
        else:
            self.ids=my_nest.Create(model, n, params) # Create models

        # Get local ids on this processor. Necessary to have for mpi run.
        for id in self.ids:           
            nodetype=my_nest.GetStatus([id])[0]['model']        
            if nodetype != 'proxynode':
                 self.local_ids.append(id)
        
        self.params = my_nest.GetStatus(self.ids)
        
        
        # Pick out recordables and receptor types using first model.
        try: 
            self.recordables = my_nest.GetDefaults(model)['recordables']
        except: 
            pass 
        try:    
            self.receptor_types = my_nest.GetDefaults(model)['receptor_types']
        except:    
            pass     
        
        if self.recordables:
            if any(record_from): 
                self.record_from=record_from
            else:                
                self.record_from=self.recordables
                
        
        # Add spike detector
        if sd: 
            self.sd = my_nest.Create("spike_detector")       
            self.sd_params.update({"withgid": True })
            my_nest.SetStatus(self.sd, self.sd_params)
            my_nest.ConvergentConnect(self.ids, self.sd)
        
        # Record with multimeter from first neuron 
        
        if mm: 
            self.mm = my_nest.Create("multimeter")    
            self.mm_dt = mm_dt # Recording interval
            my_nest.SetStatus(self.mm, {'interval': self.mm_dt, 'record_from': self.record_from})       
            my_nest.DivergentConnect(self.mm, self.ids)
    
    def __getitem__(self, key):
        ''' 
        Calling self with one index
        '''
        ids=numpy.array(self.ids) 
        ids=ids[key]

        if type(ids).__name__ in ['int', 'int64']:
            return ids
        else:
            return list(ids)
      
    def __getslice__(self, i,j):
        ''' 
        Calling self by its own 'self' or as slice object 'self[i,j] 
        '''
         
        return self.ids[i:j]

    def __len__(self):
        ''' 
        Return lenght of ids list. Neccesary to have to be able to call 
        self[1:-1] where -1 inforce length lookup
        '''
        return len(self.ids)   
             
    #def __repr__(self):
    #    return self.ids
    
    def __str__(self):
        '''
        Function called when printing object.
        '''
        return str(self.ids)
    
    def _create_signal_object(self, dataType, recordable='spikes', start=None, stop=None ):
        '''
        -_create_signal_object(self, self, dataType, recordable='times', stop=None )
        Creates NeuroTool signal object for the recordable simulation data.  

        
        Arguments:
        dataType        type of data. 's' or 'spikes' for spike data, 
                        'g' for conductance data, 'c' for current data and 
                        'v' for voltage data
        recordable      Need to be supplied for conductance, current and 
                        voltage data. It is the name of my_nest recorded data with
                        multimeter, e.g. V_m, I_GABAA_1, g_NMDA.
        stop            end of signal in ms
        '''
        
        # Short cuts
        ids=self.local_ids      # Eacj processor has its set of ids
        mm_dt=self.mm_dt
        spath=self.spath
        sname=self.sname
        
        
        # File to save to
        extension= '-' + recordable + '-' +  str(my_nest.Rank()) + '.dat'
        fileName = spath + '/' + sname + extension
        
        # Spike data
        if dataType in ['s', 'spikes']:     
            e  = my_nest.GetStatus(self.sd)[0]['events']   # get events 
            s = e['senders']                            # get senders
            t = e['times']                              # get spike times
        
            if stop: s, t = s[t<stop], t[t<stop]    # Cut out data
            signal  = zip( s, t )                   # create signal 
            
            
        # Mulitmeter data, conductance, current or voltage data    
        elif dataType in ['g', 'c', 'v']:     
            e = my_nest.GetStatus(self.mm)[0]['events']    # get events 
            v = e[recordable]                           # get analog value
            s = e['senders']                            # get senders
            t = e['times']                              # get spike times 
            
            
                   
            if stop: 
                s, v = s[t<stop], v[t<stop]    # Cut out data
                start = stop - len(s)/len(ids)*float(mm_dt)
            else:
                start = t[0]        # start time for NeuroTools  
                
            signal  = zip( s, v )                   # create signal  
            
              
             
            
        if dataType in ['s', 'spikes']: list = MySpikeList( signal, ids, start, 
                                                            stop)             
        if dataType in ['g']: list = MyConductanceList(signal, ids, mm_dt, 
                                                       start,stop)
        if dataType in ['c']: list = MyCurrentList(signal, ids, mm_dt, start,
                                                   stop)
        if dataType in ['v']: list = MyVmList(signal, ids, mm_dt, start, 
                                              stop)    
        
        return list        
    
    def count_afferents( self, connecting_group ):
        ''' 
        Calculated number off afferents from connecting_group onto each 
        neuron in self.
        '''
        print 'Counting affarents from', self.models, 'onto', connecting_group.models
        
        connecting_group.FindConnections()
        
        d = {}
        for id in self.ids: d[ id ] = 0
        for source_id, targets in connecting_group.connections.iteritems():
            for target_id in targets:
                if target_id in self.ids: d[ target_id ] += 1                      # Check that it is self that is the target
                    
                
        return d        
                  
    def find_connections(self):
        '''
        FindConnections(self)
        Find connections for each node in layer
        '''
        
        # Clear 
        self.connections={}
        
        for node in self.ids:
            
            self.connections[str(node)] = [target for target in 
                                           my_nest.GetStatus(my_nest.FindConnections([node]), 'target') 
                                           if target not in self.sd + self.mm]
    
    def gaussian_conn_par(self, type='delay', std_rel=0.1, sead=None):
        '''
        GaussianConnPar(self, type='delay', std_rel=0.1, sead=None)
        Make connection parameters gaussian distributed
        '''
               
        #! Used identical sead. 
        if sead:
            random.seed(sead)
        
        if not self.connections:
            self.find_connections() 
        
        for source, targets in self.connections.iteritems():
            
            for target in targets: 
                
                conn  = my_nest.FindConnections([int(source)], [target])        
                val = my_nest.GetStatus(conn, type)[0] 
                my_nest.SetStatus(conn, params={type :  numpy.abs( random.gauss( val, std_rel * val ) ) })
            
    def gaussian_model_par(self, types=['C_m'], std_rel=0.1, sead=None):
        '''
        GaussianModelPar(self, types=['C_m'], std_rel=0.1, sead=None)
        Make model parameters gaussian distributed
        '''
        
        #! Used identical sead. 
        if sead:
            random.seed(sead)
            
        for par in types:
            
      
            for node in self.ids:
                
                # OBS Does not work when GetStatus is taken on all nodes. Simulation
                # hang. Need to reproduce this on the side and post on my_nest forum.
                # sta=my_nest.GetStatus([nodes[0]])[0]
                
                # If node type is proxynode then do nothing. The node is then a
                # a shadow node for mpi process as I understand it.
                nodetype=my_nest.GetStatus([node])[0]['model']
               
                if nodetype != 'proxynode':
  
                    val = my_nest.GetStatus([node],par)[0]
                    rand = numpy.abs( random.gauss( val, std_rel * val ) )
                    if par in ['V_t', 'V_r']: rand = -rand
                    my_nest.SetStatus([node], { par : rand })                    
    
    def get_conn_par(self, type='delay'): 
        '''
        Get all connections parameter values for type
        '''
        
        conn_par={}
        
        if not self.connections: self.Find_connections() 
        
        for source, targets in self.connections.iteritems():
        
            conn_par[source]=[my_nest.GetStatus(my_nest.FindConnections([int(source)], 
                                                                  [target]), 'delay')[0] 
                              for target in targets] 
        
        return conn_par        
    
    def get_model_par(self, type='C_m'): 
         '''
         Retrieve one or several parameter values from all nodes
         '''
         model_par=[my_nest.GetStatus([node],type)[0] for node in self.ids]

         return model_par
                   
    def get_signal(self, dataType, recordable='spikes', start=None, stop=None ):
        '''
        get_signal(self, self, dataType, recordable='spikes', stop=None )
        Sets group NeuroTool signal object for the recordable simulation data.  

        
        Arguments:
        dataType        type of data. 's' or 'spikes' for spike data, 
                        'g' for conductance data, 'c' for current data and 
                        'v' for voltage data
        recordable      Need to be supplied for conductance, current and 
                        voltage data. It is the name of my_nest recorded data with
                        multimeter, e.g. V_m, I_GABAA_1, g_NMDA.
        stop            end of signal in ms
        '''
        list=self._create_signal_object(dataType, recordable, start, stop)
        self.signals[recordable]=list
           
    def load_signal(self, dataType, recordable='spikes'):
        '''
        load_signals(self, dataType, recordable)
        Loads simulation data. Needs data type and name of recordable to load.
        
        Arguments
        dataType        type of data. 's' or 'spikes' for spike data, 
                        'g' for conductance data, 'c' for current data and 
                        'v' for voltage data
        recordable      Need to be supplied for conductance, current and 
                        voltage data. It is the name of my_nest recorded data with
                        multimeter, e.g. V_m, I_GABAA_1, g_NMDA. 
        ''' 
        # Short cuts
        spath=self.spath
        sname=self.sname
            
        fileNames = misc.read_f_name( spath, 
                                      contain_string=sname + recordable )
                            
        # First load first spike list
        list  = my_signals.my_load(spath+'/'+fileNames[0], dataType)
        
        # Then load the rest of spike list and merge. ( happens 
        # only with data from mpi run)
        for name in fileNames[1:]:
            if dataType in ['s', 'spikes']: 
                list.merge(my_signals.my_load(spath+'/'+name, dataType))
            else:
                # Since each processor has its own ids (local ids) append can
                # be used here.
                list2=my_signals.my_load(spath+'/'+name, dataType)
                for id, signal in list2.analog_signals.iteritems():
                    list.append(id, signal)
                
        # Add to signals    
        self.signals[recordable] = list
    
    def IF( self, I_vec, id = None, tStim = None ):    
        '''
        Function that creates I-F curve
        Inputs:
                id      - id of neuron to use for calculating 
                                 I-F relation
                tStim   - lenght of each step current injection 
                                in miliseconds
                I_vec   - step currents to inject
        
        Returns: 
                fIsi          - first interspike interval 
                mIsi          - mean interspike interval

                
        Examples:
                >> n  = nest.Create('izhik_cond_exp')
                >> sc = [ float( x ) for x in range( 10, 270, 50 ) ]
                >> f_isi, m_isi, l_isi = IF_curve( id = n, sim_time = 500, I_vec = sc ):
        '''
        if not id: id=self.ids[0]
        if isinstance( id, int ): id =[ id ] 

        
        fIsi, mIsi, lIsi  = [], [], []  # first, mean and last isi
        if not tStim: tStim = 500.0 
        tAcum = 0
    
        I_e0=my_nest.GetStatus(id)[0]['I_e'] # Retrieve neuron base current
        
        isi_list=[]
        for I_e in I_vec:
            
            my_nest.SetStatus( id, params = { 'I_e': float(I_e+I_e0) } )   
            my_nest.SetStatus( id, params = { 'V_m': float(-61) } )            
            #my_nest.SetStatus( id, params = { 'w': float(0) } )
            my_nest.SetStatus( id, params = { 'u': float(0) } )
            simulate=True
            tStart=tAcum
            while simulate:
                my_nest.Simulate( tStim )
                tAcum+=tStim    
                
                self.get_signal('s', start=tStart, stop=tAcum)
                                
                signal=self.signals['spikes'].time_slice(tStart, tAcum)
                
                
                if signal.mean_rate()>0.1 or tAcum>20000:
                    simulate=False
                    
            isi=signal.isi()[0]      
            if not any(isi): isi=[1000000.] 
              
            fIsi.append( isi[ 0 ] )            # retrieve first isi
            mIsi.append( numpy.mean( isi ) )   # retrieve mean isi
            if len(isi)>100:lIsi.append( isi[ -1 ] )
            else:lIsi.append( isi[ -1 ] )           # retrieve last isi
            isi_list.append(numpy.array(isi))
            
        fIsi=numpy.array(fIsi)
        mIsi=numpy.array(mIsi)
        lIsi=numpy.array(lIsi)
        
        I_vec=numpy.array(I_vec)
            
        return I_vec, fIsi, mIsi, lIsi


    def I_PSE(self, I_vec, synapse_model, id=0,  receptor='I_GABAA_1'):
        '''
        Assure no simulations has been run before this (reset kernel).
        Function creates relation between maz size of postsynaptic 
        event (current, conductance, etc). The type is set by receptor. 
        
        Inputs:
            I_vec          - step currents to clamp at
            id             - id of neuron to use for calculating I-F relation
                             If not providet id=ids[0]
     
        Returns: 
            v_vec          - voltage clamped at
            size_vec       - PSE size at each voltage 
                
        Examples:
                >> n  = my_nest.Create('izhik_cond_exp')
                >> sc = [ float( x ) for x in range( -300, 100, 50 ) ] 
                >> tr_t, tr_v, v_ss = IV_I_clamp( id = n, tSim = 500, 
                I_vec = sc ):
        '''
        
        vSteadyState=[]
        
        if not id: id=self.ids[0]
        if isinstance( id, int ): id =[ id ]
        
        simTime  = 700. # ms
        spikes_at = numpy.arange(500., len(I_vec)*simTime,simTime) # ms
        
    
        voltage=[] # mV
        pse=[] # post synaptic event
        
        sg = my_nest.Create('spike_generator', params={'spike_times':spikes_at} )
    
        my_nest.Connect(sg, id, model=synapse_model)
        
        simTimeTot=0
        for I_e in I_vec:
            
            my_nest.SetStatus(self[:], params={'I_e':float(I_e)})
            my_nest.MySimulate(simTime)
            simTimeTot+=simTime
        
        
        self.get_signal( receptor[0].lower(),receptor, 
                         stop=simTimeTot ) # retrieve signal
        
        simTimeAcum=0
        
        
        for I_e in I_vec:
            
            
            size=[]
   
            signal=self.signals[receptor].my_time_slice(400+simTimeAcum, 
                                                        700+simTimeAcum)
            simTimeAcum+=simTime     
            # First signal object at position 1    
            clamped_at  = signal[1].signal[999]
            minV=min(signal[1].signal)
            maxV=max(signal[1].signal)
            if abs(minV-clamped_at)<abs(maxV-clamped_at):        
                size.append(max(signal[1].signal)-clamped_at)
                
            else:
                size.append(min(signal[1].signal)-clamped_at)
        
            voltage.append(clamped_at)
            pse.append(size[0])
                                   
    
        return voltage, pse
         
    def IV_I_clamp(self, I_vec, id = None, tStim = 2000):    
        '''
        Assure no simulations has been run before this (reset kernel).
        Function that creates I-V by injecting hyperpolarizing currents 
        and then measuring steady-state membrane (current clamp). Each 
        trace is preceded and followed by 1/5 of the simulation time 
        of no stimulation
        Inputs:
            I_vec          - step currents to inject
            id             - id of neuron to use for calculating I-F relation
            tStim          - lenght of each step current stimulation in ms

        Returns: 
            I_vec          - current for each voltage
            vSteadyState   - steady state voltage 
                
        Examples:
                >> n  = my_nest.Create('izhik_cond_exp')
                >> sc = [ float( x ) for x in range( -300, 100, 50 ) ] 
                >> tr_t, tr_v, v_ss = IV_I_clamp( id = n, tSim = 500, 
                I_vec = sc ):
        '''
        
        vSteadyState=[]
        
        if not id: id=self.ids[0]
        if isinstance( id, int ): id =[ id ] 
        
                                                                           
        tAcum  = 1    # accumulated simulation time, step_current_generator
                      # recuires it to start at t>0    
        
        scg = my_nest.Create( 'step_current_generator' )  
        rec=my_nest.GetStatus(id)[0]['receptor_types']
        my_nest.Connect( scg, id, params = { 'receptor_type' : rec['CURR'] } )             
        
        ampTimes=[]
        ampValues=[]
        for I_e in I_vec:

            ampTimes.extend([ float(tAcum) ])
            ampValues.extend([ float(I_e) ])
            tAcum+=tStim
            
        my_nest.SetStatus( scg, params = { 'amplitude_times':ampTimes,       
                                           'amplitude_values' :ampValues } )   
        my_nest.Simulate(tAcum)    
        
        self.get_signal( 'v','V_m', stop=tAcum ) # retrieve signal
        self.get_signal('s')
        if 0 < self.signals['spikes'].mean_rate():
            print 'hej'
        tAcum  = 1 
        for I_e in I_vec:
            
            if 0>=self.signals['spikes'].mean_rate(tAcum+10, tAcum+tStim):
                signal=self.signals['V_m'].my_time_slice(tAcum+10, tAcum+tStim)
                vSteadyState.append( signal[1].signal[-1] )
            tAcum+=tStim
    
            
        I_vec=I_vec[0:len(vSteadyState)]         
        return I_vec, vSteadyState
            
    def mean_weights(self):
        '''
        Return a dictionary with mean_weights for each synapse type with mean weight
        and receptor type
        '''
        
        print 'Calculating mean weights', self.models
        
        syn_dict = {}                                                                 # container weights per synapse type
        rev_rt   = {}                                                                 # receptor type number dictionary
        rt_nb    = []                                                                 # receptor type numbers

            

        for source in self.ids:                                                     # retrieve all weights per synapse type
            for conn in my_nest.GetStatus( my_nest.FindConnections( [ source ] ) ):
                st = conn[ 'synapse_type' ] 

                if syn_dict.has_key( st ):
                    syn_dict[ st ]['weights'].append( conn['weight'] )
                else:
                    syn_dict[ st ] = {}
                    syn_dict[ st ]['weights']       = [ conn['weight'] ]
                    syn_dict[ st ]['receptor_type'] = { st : my_nest.GetDefaults( st )[ 'receptor_type' ]  } 
                    
        
        
        mw = {}                                                                     # container mean weights per synapse type
        for key, val in syn_dict.iteritems():
            if len( val[ 'weights' ] ) > 0: 
                syn_dict[ key ]['mean_weight'] = sum( val[ 'weights' ] )/len( val[ 'weights' ] )   # calculate mean weight             
                syn_dict[ key ]['weights'] = numpy.array( syn_dict[ key ]['weights'] )    
                syn_dict[ key ]['nb_conn'] = len( val[ 'weights' ] )
            else:
                syn_dict[ key ]['mean_weight'] = 0
                syn_dict[ key ]['nb_conn']     = 0
            
        
        return syn_dict   
    
    def merge(self, group2):
        self.ids.extend(group2.ids) 
        self.ids.sort()
        
        # Merge signals
        if self.signals:
            for key in self.signals.keys():
                if key in ['spikes']:
                    self.signals['spikes'].merge(group2.signals['spikes'])
 

        
        
    def slice(self, ids ):
        group=copy.deepcopy(self)    
        
        group.ids=ids
        # Slice signals
        if group.signals:
            for key in group.signals.keys():
                if key in ['spikes']:
                    group.signals['spikes']=self.signals['spikes'].id_slice(group.ids)
         
        return group
     
    def print_neuron(self):
        '''
        PrintNeuron(self)
        Print layer info 
        '''
        print ' '
        print 'Model: ' + self.models
        print 'Ids: ' + str(self.ids)
        print 'recordables: ' +  str(self.recordables) 
        print 'sd: ' + str(self.sd) 
        print 'mm: ' + str(self.mm) 
        
        print 'Params:'
        for key, val in self.params.iteritems():
            print key + ': ' + str(val)
    
    def print_connections(self):
        '''
        PrintConnections(self)
        Print connections for each node in the layer
        '''
        print 'Node targets:'

        if not self.connections:
            self.FindConnections() 
        
        for key, value in self.connections.iteritems():
            print key + ' ' +str(value) 
     
    def save_group(self):
        '''
        Not complete
        '''
                  
        group = {}
    
        group["connections"] = self.connections
        group["models"] = self.models
        group["ids"] = self.ids
        group["recordables"] = self.recordables
        group["receptor_types"] = self.receptor_types
        group["sname"] = self.sname
        group["GID"] = str(my_nest.GetStatus(self.sd)[0]['global_id'])
        group["mm_dt"] = self.mm_dt
        
        output = open(self.spath + '/' +self.sname + '-' +  str(my_nest.Rank()) + '.' + 'pickle', 'w')

        # Pickle dictionary 
        pickle.dump(group, output)

        output.close() 
       
    def save_signal(self, dataType, recordable='spikes', start=None, stop=None):  
        '''
        save_signal(self, self, dataType, recordable='spikes', stop=None )
        Save NeuroTool signal object for the recordable simulation data.  

        
        Arguments:
        dataType        type of data. 's' or 'spikes' for spike data, 
                        'g' for conductance data, 'c' for current data and 
                        'v' for voltage data
        recordable      Need to be supplied for conductance, current and 
                        voltage data. It is the name of my_nest recorded data with
                        multimeter, e.g. V_m, I_GABAA_1, g_NMDA.
        stop            end of signal in ms
        '''
        
        fileName = (self.spath + '/' + self.sname + recordable + '-'+ 
                    str(my_nest.Rank()) + '.dat')
        
        # Create save dir if it do not exist
        try: 
            os.system('mkdir ' + self.spath + ' 2>/dev/null')
        except:
            pass

        try:  
            os.system('rm ' + fileName  + ' 2>/dev/null' )  
        except:
            pass
            
                
        if recordable in self.signals.keys(): 
            list=self.signals[recordable]
        else:
            list = self._create_signal_object(dataType, recordable, start, stop)
        list.my_save(fileName)     # Save signal
    
    
class MyPoissonInput(MyGroup):
        
    def __init__(self,model = 'spike_generator', n=1, params = {}, mm_dt=1.0, 
                 sname='', spath='', sname_nb='', sd=False, sd_params={},
                 mm=False, record_from=[], ids=[]):
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
        # is ConductanceList. The next argument. self, passes a reference to 
        # the object so that ConductanceList can get to the object and add 
        # its attributes to it. The next part of the statement __init__(
        # signals, id_list, dt, t_start, t_stop, dims) tells python I want to
        # invoke the constructor method of ConductanceList and a want to pass 
        # it the values of signals, id_list, dt, t_start, t_stop and dims.
        
        # Invoke __init__ of base class ConductanceList which has 
        # AnalogSignalList as base class where the method exist.
        super( MyPoissonInput, self ).__init__(model, n, params, mm_dt, 
                                                   sname, spath, sname_nb, sd, 
                                                   sd_params, mm, record_from, ids)
        self._init_extra_attributes()
    

        
        
    def _init_extra_attributes(self):
        # Add new attribute
        self.rates={}
        self.times={}
        self.t_stop={}
        
        for id in self.ids:
            self.rates[id]=[]    # rates for id, id as key
            self.times[id]=[]    # times for id, id as key
            self.t_stop[id]=0   # t_stop for id, id as key
    
#    def set_spike_times(self, rates=[], times=[], t_stop=None, ids=None, seed=None):
#        
#        
#        if ids is None:
#            for id in self.ids:
#                seed=random_integers(0,10**5)
#                      
#                spikeTimes=misc.inh_poisson_spikes( rates, times, t_stop=t_stop, 
#                                            n_rep=1, seed=seed)
#                if any(spikeTimes):
#                    my_nest.SetStatus([id], params={'spike_times':spikeTimes})
#        
#        else:
#            for id in ids:
#                seed=random_integers(0,10**5)
#                      
#                spikeTimes=misc.inh_poisson_spikes( rates, times, t_stop=t_stop, n_rep=1, seed=seed)
#                if any(spikeTimes):
#                    my_nest.SetStatus([id], params={'spike_times':spikeTimes})

    def set_spike_times(self, rates=[], times=[], t_stop=None, ids=None, seed=None, idx=None):
        
        df=my_nest.GetDefaults(self.model)

        if ids is None and (not idx is None):            
            tmp_ids=numpy.array(self.ids)
            ids=list(tmp_ids[idx])
        if ids is None:            
            ids=self.ids 
            
        #print df.keys()
        #print df['model']      
        # Spike generator
        if 'spike_generator' == df['model']:
            
            
            
            for id in ids:
                seed=random_integers(0,10**5)
                          
                spikeTimes=misc.inh_poisson_spikes( rates, times, t_stop=t_stop, n_rep=1, seed=seed)
                if any(spikeTimes):
                    rs=my_nest.GetKernelStatus('resolution')
                    n_dec=int(-numpy.log10(rs))
                    spikeTimes=numpy.round(spikeTimes,n_dec)
                    my_nest.SetStatus([id], params={'spike_times':spikeTimes})
                    
        
        # MIP
        elif 'mip_generator' == df['model']:
            c=df['p_copy']
            
            seed=random_integers(0,10**6)
            new_ids=[]
            t_starts=times
            t_stops=times[1:]+[t_stop]
            for id in ids:
                i=0
                for r, start, stop in rates, t_starts, t_stops:               
                   r_mother=r/c
                   params={'rate':r_mother, 'start':start, 'stop':stop,
                           'p_copy':c, 'mother_seed':seed}
                   if i==0:
                       nest.SetStatus(id, params)
                   else:                   
                        new_id=my_nest.Create('mip_generator', 1, params)
                        
                   new_ids.append(new_id)
            self.ids.append(new_ids)  
        
        # Poisson generator
        elif 'poisson_generator' == df['model']:

            t_starts=times
            t_stops=list(times[1:])+list([t_stop])
            for i in range(len(ids)):
                
                id=ids[i]
                model=my_nest.GetStatus([id],'model')[0]
                
                if model=='parrot_neuron':                  
                    j=1
                else:
                    ids[i]=my_nest.Create('parrot_neuron')[0]
                    my_nest.Connect([id], [ids[i]])
                    #self.ids_mul_visited[id]=True
                    j=0
                for r, start, stop in zip(rates, t_starts, t_stops):
                   params={'rate':r, 'start':start, 'stop':stop}
                   if j==0:
                       #print my_nest.GetStatus([id])
                       my_nest.SetStatus([id], params)
                   else:    
                        if params['rate']!=0:       
                            #print params       
                            new_id=my_nest.Create('poisson_generator', 1, params)
                            #self.ids_mul[ids[i]].append(new_id[0])
                            my_nest.Connect(new_id, [ids[i]])
                   j+=1
            
            if not idx is None:
                 tmp_ids=numpy.array(self.ids)
                 tmp_ids[idx]=list(ids)
                 self.ids=list(tmp_ids)
            else:
                 self.ids=list(ids)
            
            self.local_ids=list(self.ids) # Nedd to put on locals also
        
class MyLayerGroup(MyGroup):
    
    def __init__(self, layer_props={}, model = 'iaf_neuron', n=1, params = {}, mm_dt=1.0, 
                 sname='', spath='', sname_nb=0, sd=False, sd_params={},
                 mm=False, record_from=[]): 
                 
                layer=my_topology.CreateLayer(layer_props)
                ids=my_nest.GetLeaves(layer)[0]
                
                super( MyLayerGroup, self ).__init__(  model, n, params, mm_dt, 
                                                   sname, spath, sname_nb, sd, 
                                                   sd_params, mm, record_from,ids)
                self._init_extra_attributes_layer_group(layer, layer_props)
    


        
        
    def _init_extra_attributes_layer_group(self, layer, props):
        
        
        # Add new attribute
        self.layer_id=layer
        self.conn={}
        self.id_mod=[]
    def add_connection(self, source, type, props):
        ''' Store connection properties of an incomming synapse' ''' 
        name=str(source) + '_'+type
        self.conn[name]={}
        self.conn[name]['mask']=props['mask']
        self.conn[name]['kernel']=props['kernel']   
        
    def get_kernel(self, source, type):
        name=str(source) + '_'+type
        
        return self.conn[name]['kernel']
    
    def get_mask(self, source, type):
        name=str(source) + '_'+type
        
        return self.conn[name]['mask']
    
    def stat_connections(self, source_layer, type):
        import time
        import nest.topology as tp
        
        t=time.time()        
        targets=tp.GetTargetNodes(source_layer.ids, self.layer_id)
        
        m=[]
        for trg in targets:
            m.append(len(trg))
        print time.time()-t
        
        t=time.time()  
        m=[]
        for id in source_layer.ids:
            m.append(len(my_topology.GetTargetNodes([id], self.layer_id)[0]))
        print time.time()-t
        
        print len(m)
        mean=numpy.mean(numpy.array(m))
        std=numpy.std(numpy.array(m))
        return mean, std
            
            
            
    
    def plot(self, ax=None, nodecolor='b', nodesize=20):        
       my_topology.MyPlotLayer(self.layer_id, ax,nodecolor, nodesize)
    
    def plot_targets(self, type, source_layer, src_nrn=None, 
                     tgt_model=None, syn_type=None, ax=None,
                     src_color='red', src_size=50, tgt_color='red', tgt_size=20,
                     mask_color='red', kernel_color='red'):
        
        if not src_nrn:
            src_nrn = my_topology.FindCenterElement ( source_layer.layer_id)
        
        kernel=self.get_kernel(source=source_layer, type=type) 
        mask=self.get_mask(source=source_layer, type=type)
        
        ax, n_targets=my_topology.MyPlotTargets(src_nrn, self.layer_id, ax = ax, 
                   mask=mask , kernel=kernel, src_size=src_size , 
                   tgt_color=tgt_color , tgt_size=tgt_size , 
                   kernel_color =kernel_color)
        return n_targets
class MyLayerPoissonInput(MyPoissonInput):
    
    def __init__(self, layer_props={},  ids=[], model = 'spike_generator', n=1, params = {}, mm_dt=1.0, 
                 sname='', spath='', sname_nb='', sd=False, sd_params={},
                 mm=False, record_from=[]):
                 
                layer=my_topology.CreateLayer(layer_props)
                ids=my_nest.GetLeaves(layer)[0]
                
                super( MyLayerPoissonInput, self ).__init__(model, n, params, mm_dt, 
                                                   sname, spath, sname_nb, sd, 
                                                   sd_params, mm, record_from,ids)
                self._init_extra_attributes_poisson_input(layer, layer_props)
          
    def _init_extra_attributes_poisson_input(self, layer, props):
        
        
        # Add new attribute
        self.layer_id=layer
        self.conn={}
        self.id_mod=[]
    
    def sort_ids(self, pos=[[0,0]]):     
        print self.layer_id
        node=my_topology.FindCenterElement(self.layer_id)
        print node
        ids=self.ids
        d=numpy.array(my_topology.Distance(node*len(ids),ids))
        idx=sorted(range(len(d)), key=d.__getitem__, reverse=True)
        print d[idx]
        return idx 
    def plot(self, ax=None, nodecolor='b', nodesize=20):        
       my_topology.MyPlotLayer(self.layer_id, ax,nodecolor, nodesize)
    


              
         