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
 spike date) and nest recordable type (except for spike data). With 
 save_signal() and load_signal() data can instead be written and read from disk.
 
 Function can take several modes and number of neurons per model. Then the models
 has to have the same recordables. Also if parameters is provided in params these
 are set for all the models. 
'''

import nest
from NeuroTools import signals 

import numpy
import os
import random # Random generator
import my_signals # Own wrapper class for signal module in Neurotools
from my_signals import MyConductanceList, MyCurrentList, MyVmList, MySpikeList
import misc


class MyGroup:
    '''
    MyGroup(self, models = 'iaf_neuron', n=1, params = {}, mm_dt=1.0, sname='', spath='', sd=True, mm=True)
    
    Arguments:
        models      nest model type, can be a list
        n           number of models to create, can be a list
        params      common parameters for models to be set
        mm_dt       multimeter recording precision
        sname       file basename
        spath       Path to save file at 
        sd          boolean, True if spikes should me recorded
        mm          boolean, True if mulitmeter should record  
        
    ''' 
    def __init__(self, models = 'iaf_neuron', n=1, params = {}, mm_dt=1.0, sname='', spath='', sd=True, mm=True):
        '''
        Constructor
        
        Arguments:
            models      nest model type, can be a list
            n           number of models to create, can be a list
            params      common parameters for models to be set
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
        self.models         = []
        self.params         = []
        self.receptor_types = {}
        self.recordables    = {}
        self.sd             = []        # Id of spike detector
        self.sname          = ''        # Specific file basename
        self.spath          = ''        # Path to save file at 
        self.signals        = {}        # dictionary with signals for current, conductance, voltage or spikes    
        
        # If model is a string put it in a list
        if isinstance(models, str):
            models =[models]
        
        # If n is a integer put it in a list    
        if isinstance(n, int):
            n =[n]  
                
        self.models = models
        
        # If no sname is provided take first model neuron as save name.
        if sname is '': self.sname = self.models[0]
        else: self.sname = sname
        
        
        # If no spath is provided current path plus data_tmp is set to 
        # spath.
        if spath is '': self.spath = os.getcwd()+'/output_tmp'    
        else: self.spath = spath
                    
        # Create save dir if it do not exist
        try: os.system('mkdir ' + self.spath )                                   
        except: print self.spath  + ' already exists'
        
        print 'Group save path: '+ self.spath 
                       
        # For each model create it with i neurons.
        for model, i in zip(models,n):
        
            self.ids.extend(nest.Create(model, i, params))
  
        
        
        # Get local ids on this processor. Necessary to have for mpi run.
        for id in self.ids:
            
            nodetype=nest.GetStatus([id])[0]['model']
        
            if nodetype != 'proxynode':
            
                 self.local_ids.append(id)
        
        self.params = nest.GetStatus(self.ids)
        
        
        # Pick out recordables and receptor types using first model.
        try: self.recordables = nest.GetDefaults(models[0])['recordables']
        except: print 'No recordables'
        try:    self.receptor_types = nest.GetDefaults(models[0])['receptor_types']
        except:    print 'No receptor types'      
        
        
        # Add spike detector
        if sd: self.sd = nest.Create("spike_detector")
        
        nest.SetStatus(self.sd, {"withgid": True })
        nest.ConvergentConnect(self.ids, self.sd)
        
        # Record with multimeter from first neuron 
        
        if mm: 
            self.mm = nest.Create("multimeter")    
            self.mm_dt = mm_dt # Recording interval
            nest.SetStatus(self.mm, {'interval': self.mm_dt, 'record_from': self.recordables})
        
            nest.DivergentConnect(self.mm, self.ids)
      
    def __getitem__(self, key):
        ''' 
        Now self can be called as an list object returning the ids. 
        '''
         
        return self.ids[key]

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
    
    def _create_signal_object(self, dataType, recordable='spikes', stop=None ):
        '''
        -_create_signal_object(self, self, dataType, recordable='times', stop=None )
        Creates NeuroTool signal object for the recordable simulation data.  

        
        Arguments:
        dataType        type of data. 's' or 'spikes' for spike data, 
                        'g' for conductance data, 'c' for current data and 
                        'v' for voltage data
        recordable      Need to be supplied for conductance, current and 
                        voltage data. It is the name of nest recorded data with
                        multimeter, e.g. V_m, I_GABAA_1, g_NMDA.
        stop            end of signal in ms
        '''
        
        # Short cuts
        ids=self.local_ids
        mm_dt=self.mm_dt
        spath=self.spath
        sname=self.sname
        
        
        # File to save to
        extension= '-' + recordable + '-' +  str(nest.Rank()) + '.dat'
        fileName = spath + '/' + sname + extension
        
        # Spike data
        if dataType in ['s', 'spikes']:     
            e  = nest.GetStatus(self.sd)[0]['events']   # get events 
            s = e['senders']                            # get senders
            t = e['times']                              # get spike times
        
            if stop: s, t = s[t<stop], t[t<stop]    # Cut out data
            signal  = zip( s, t )                   # create signal 
            
            start=0
        # Mulitmeter data, conductance, current or voltage data    
        elif dataType in ['g', 'c', 'v']:     
            e = nest.GetStatus(self.mm)[0]['events']    # get events 
            v = e[recordable]                           # get analog value
            s = e['senders']                            # get senders
            t = e['times']                              # get spike times
                   
            if stop: s, v = s[t<stop], v[t<stop]    # Cut out data
            signal  = zip( s, v )                   # create signal  
            
            start = t[0]        # start time for NeuroTools    
             
            
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
                                           nest.GetStatus(nest.FindConnections([node]), 'target') 
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
            self.FindConnections() 
        
        for source, targets in self.connections.iteritems():
            
            for target in targets: 
                
                conn  = nest.FindConnections([int(source)], [target])        
                val = nest.GetStatus(conn, type)[0] 
                nest.SetStatus(conn, params={type :  numpy.abs( random.gauss( val, std_rel * val ) ) })
            
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
                # hang. Need to reproduce this on the side and post on nest forum.
                # sta=nest.GetStatus([nodes[0]])[0]
                
                # If node type is proxynode then do nothing. The node is then a
                # a shadow node for mpi process as I understand it.
                nodetype=nest.GetStatus([node])[0]['model']
               
                if nodetype != 'proxynode':
  
                    val = nest.GetStatus([node],par)[0]
                    rand = numpy.abs( random.gauss( val, std_rel * val ) )
                    if par in ['V_t', 'V_r']: rand = -rand
                    nest.SetStatus([node], { par : rand })                    
                    print par,rand
    
    def get_conn_par(self, type='delay'): 
        '''
        Get all connections parameter values for type
        '''
        
        conn_par={}
        
        if not self.connections: self.Find_connections() 
        
        for source, targets in self.connections.iteritems():
        
            conn_par[source]=[nest.GetStatus(nest.FindConnections([int(source)], 
                                                                  [target]), 'delay')[0] 
                              for target in targets] 
        
        return conn_par        
    
    def get_model_par(self, type='C_m'): 
         '''
         Retrieve one or several parameter values from all nodes
         '''
         model_par=[nest.GetStatus([node],type)[0] for node in self.ids]

         return model_par
                   
    def get_signal(self, dataType, recordable='spikes', stop=None ):
        '''
        get_signal(self, self, dataType, recordable='spikes', stop=None )
        Sets group NeuroTool signal object for the recordable simulation data.  

        
        Arguments:
        dataType        type of data. 's' or 'spikes' for spike data, 
                        'g' for conductance data, 'c' for current data and 
                        'v' for voltage data
        recordable      Need to be supplied for conductance, current and 
                        voltage data. It is the name of nest recorded data with
                        multimeter, e.g. V_m, I_GABAA_1, g_NMDA.
        stop            end of signal in ms
        '''
        list=self._create_signal_object(dataType, recordable, stop)
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
                        voltage data. It is the name of nest recorded data with
                        multimeter, e.g. V_m, I_GABAA_1, g_NMDA. 
        ''' 
        # Short cuts
        spath=self.spath
        sname=self.sname
            
        fileNames = misc.read_f_name( spath, 
                                      contain_string=sname + '-' + recordable )
                            
        # First load first spike list
        list  = my_signals.load(spath+'/'+fileNames[0], dataType)
        
        # Then load the rest of spike list and merge. ( happens 
        # only with data from mpi run)
        for name in fileNames[1:]:
            list.merge(my_signals.load(spath+'/'+name, dataType))
            
        # Add to signals    
        self.signals[recordable] = list
            
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
            for conn in nest.GetStatus( nest.FindConnections( [ source ] ) ):
                st = conn[ 'synapse_type' ] 

                if syn_dict.has_key( st ):
                    syn_dict[ st ]['weights'].append( conn['weight'] )
                else:
                    syn_dict[ st ] = {}
                    syn_dict[ st ]['weights']       = [ conn['weight'] ]
                    syn_dict[ st ]['receptor_type'] = { st : nest.GetDefaults( st )[ 'receptor_type' ]  } 
                    
        
        
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
        group["GID"] = str(nest.GetStatus(self.sd)[0]['global_id'])
        group["mm_dt"] = self.mm_dt
        
        output = open(self.spath + '/' +self.sname + '-' +  str(nest.Rank()) + '.' + 'pickle', 'w')

        # Pickle dictionary 
        pickle.dump(group, output)

        output.close() 
       
    def save_signal(self, dataType, recordable='spikes', stop=None):  
        '''
        save_signal(self, self, dataType, recordable='spikes', stop=None )
        Save NeuroTool signal object for the recordable simulation data.  

        
        Arguments:
        dataType        type of data. 's' or 'spikes' for spike data, 
                        'g' for conductance data, 'c' for current data and 
                        'v' for voltage data
        recordable      Need to be supplied for conductance, current and 
                        voltage data. It is the name of nest recorded data with
                        multimeter, e.g. V_m, I_GABAA_1, g_NMDA.
        stop            end of signal in ms
        '''
        list = self._create_signal_object(dataType, recordable, stop)
        list.save(fileName)     # Save signal
    
        
        
        
    
         