'''
Mikael Lindahl 2010


Module:
mynest

Here my own nest functions can be defined. For example connection functions
setting random weight or delay on connections.


'''

# Imports
import nest
import numpy
import numpy.random as rand
from nest import *



def C(pre, post, params = None, delay = None, model = "static_synapse" ):    
    
    '''
    As NEST Connect, Make one-to-one connections of type model between 
    the nodes in pre and the nodes in post. pre and post have to be lists 
    of the same length. If params is given (as dictionary or list of
    dictionaries), they are used as parameters for the connections. If
    params is given as a single float or as list of floats, it is used
    as weight(s), in which case delay also has to be given as float or
    as list of floats.
    '''
    
    if isinstance(model, str): model = [ model ]  
    if isinstance(pre, int):   pre   = [ pre   ]                                # since NEST needs list objects
    if isinstance(post, int):  post  = [ post  ]                                # since NEST needs list objects
     
    nest.Connect(pre, post, params = params, delay = delay, model = model[ 0 ] ) 
    
    if len( model ) > 1:                                                        # pre -> post with model[1], model[2] etc  
        for m in model[1:]:
            nest.Connect(pre, post, params = params, delay = delay, model = m ) 
                



def CC( pre,  post, 
        params = { 'd_mu'  : [], 'w_mu'  : [], 'd_rel_std' : 0.0, 'w_rel_std' : 0.0 }, 
        model = [ 'static_synapse' ] ):
    '''
    As NEST ConvergentConnect except that it can take severals models and make same
    connections for both models. This can be used to connect same source with 
    target with both AMPA and NMDA. Mean and standard deviation for delays and
    weights can be given. Mean as a list of values one for each model and
    standard deviation as single values applied to all synapse models.
     
    Inputs:
        pre                      - list of source ids       
        post                     - list of target ids  
    	params[ 'd_mu' ]         - mean delay for each model
    	params[ 'w_mu' ]       	 - mean weight for each model
        params[ 'd_rel_std' ]    - relative mean standard deviation delays
        params[ 'w_rel_std' ]    - relative mean standard deviation weights
        model                    - synapse model, can be a list of models
    '''
    
    for key, val in params.iteritems():
        exec '%s = %s' % (key, str(val))                                        # params dictionary entries as variables
    
    
    if isinstance(model, str): model = [model]                                  # if model is a string put into a list   
    if isinstance(pre, int):   pre   = [ pre   ]                                # since NEST needs list objects
    else:                      pre   = pre[ : ]                                 # else retreive list
    if isinstance(post, int):  post  = [ post  ]                                # since NEST needs list objects
    else:                      post   = post[ : ]                               # else retreive list
    
        
    n = len( pre )
    
    if not d_mu: d_mu = [ nest.GetDefaults( m )[ 'delay' ] for m in model ]     # get mean delays d_mu
    d_sigma  = [ d_rel_std*mu for mu in d_mu ]                                  # get delay standard deviation
    
    if not w_mu: w_mu = [ nest.GetDefaults( m )[ 'weight' ] for m in model ]    # get mean weights w_mu   
    w_sigma = [ w_rel_std*mu for mu in w_mu ]                                   # get weight standard deviation            
    
    delays  = [ rand.normal( mu, sigma, [ len( pre ), n ] )                     # calculate delays, randomized if sigma != 0
               if sigma else numpy.ones( ( len( pre ), n ) )*mu           
    		   for mu, sigma in zip( d_mu, d_sigma )   ]                                                                
      
    weights = [ rand.normal( mu, sigma, [ len( pre ), n ] )                     # calculate weights, randomized if sigma != 0
               if sigma else numpy.ones( ( len( pre ), n ) )*mu           
    		   for mu, sigma in zip( w_mu, w_sigma ) ]                                                                  
    
    
    for i, post_id in enumerate( post ):                                        # connect with model[0]                                                                            
        d = [ dj for dj in delays[  0 ][ i ] ] 
        w = [ wj for wj in weights[ 0 ][ i ] ]
        nest.ConvergentConnect(  pre,  [ post_id ], weight = w , delay = d, 
                                 model = model[ 0 ] )                        
          
        if len( model ) > 1:                                                    # pre -> post with model[1], model[2] etc  
            for m, dlays, wghts, in zip( model[1:], delays[1:], weights[1:] ):
            	d = [ dj for dj in dlays[ i ] ]
                w = [ wj for wj in wghts[ i ] ]
            	nest.ConvergentConnect( pre, [ post_id ], weight = w, 
                                        delay = d, model = m )                   
   
        
        
def DC( pre, post,
        params = { 'd_mu'  : [], 'w_mu'  : [], 'd_rel_std' : 0.0, 'w_rel_std' : 0.0 }, 
        model = [ 'static_synapse' ] ):        
    '''
    As NEST DivergentConnect except that it can take severals models and make same
    connections for both models. This can be used to connect same source with 
    target with both AMPA and NMDA.Mean and standard deviation for delays and
    weights can be given. Mean as a list of values one for each model and
    standard deviation as single values applied to all synapse models.
     
    Inputs:
        pre                      - list of source ids       
        post                     - list of target ids  
        params[ 'd_mu' ]         - mean delay for each model
        params[ 'w_mu' ]            - mean weight for each model
        params[ 'd_rel_std' ]    - relative mean standard deviation delays
        params[ 'w_rel_std' ]    - relative mean standard deviation weights
        model                    - synapse model, can be a list of models
    '''
  
    for key, val in params.iteritems():
        exec '%s = %s' % (key, str(val))                                        # params dictionary entries as variables  
  
    if isinstance(model, str): model = [model]                                  # if model is a string put into a list   
    if isinstance(pre, int):   pre   = [ pre   ]                                # since NEST needs list objects
    else:                      pre   = pre[ : ]                                 # else retreive list
    if isinstance(post, int):  post  = [ post  ]                                # since NEST needs list objects
    else:                      post   = post[ : ]                               # else retreive list

    n = len( post )    

    if not d_mu: d_mu =  [nest.GetDefaults( m )[ 'delay' ] for m in model ]     # get mean delays d_mu
    d_sigma  = [ d_rel_std*mu for mu in d_mu ]                                  # get delay standard deviation
    
    if not w_mu: w_mu = [ nest.GetDefaults( m )[ 'weight' ] for m in model ]    # get mean weights w_mu   
    w_sigma = [ w_rel_std*mu for mu in w_mu ]                                   # get weight standard deviation            
    
    delays  = [ rand.normal( mu, sigma, [ len( pre ), n ] )                     # calculate delays, randomized if sigma != 0
               if sigma else numpy.ones( ( len( pre ), n ) )*mu          
			   for mu, sigma in zip( d_mu, d_sigma )   ]                                                                
      
    weights = [ rand.normal( mu, sigma, [ len( pre ), n ] )                     # calculate weights, randomized if sigma != 0
               if sigma else numpy.ones( ( len( pre ), n ) )*mu           
			   for mu, sigma in zip( w_mu, w_sigma ) ]                                                                  

    
    for i, pre_id in enumerate( pre ):                                          # connect with model[0                                    
        d = [ dj for dj in delays[  0 ][ i ] ] 
        w = [ wj for wj in weights[ 0 ][ i ] ]
        nest.DivergentConnect( [ pre_id  ], post, weight = w , delay = d, 
                               model = model[ 0 ] )                          
               
        if len( model ) > 1:                                                    # pre -> post with model[1], model[2] etc   
            for m, dlays, wghts, in zip( model[1:], delays[1:], weights[1:] ):
            	d = [ dj for dj in dlays[ i ] ]
                w = [ wj for wj in wghts[ i ] ]
            	nest.DivergentConnect( [ pre_id  ], post, weight = w, 
                                       delay = d, model = m )                           
    
            
def RDC(pre, post, n, 
        params = { 'd_mu'  : [], 'w_mu'  : [], 'd_rel_std' : 0.0, 'w_rel_std' : 0.0 },  
        model = 'static_synapse' ):
    '''
    As NEST RandomDivergentConnect, except it can take severals models and 
    make same connections for both models. This can be used to connect same
    source with  target with both AMPA and NMDA.Mean and standard deviation 
    for delays and weights can be given. Mean as a list of values one for 
    each model and standard deviation as single values applied to all 
    synapse models.
     
    Inputs:
        pre                      - list of source ids       
        post                     - list of target ids  
        n                        - number of neurons each source neuron connects to
        params[ 'd_mu' ]         - mean delay for each model
        params[ 'w_mu' ]            - mean weight for each model
        params[ 'd_rel_std' ]    - relative mean standard deviation delays
        params[ 'w_rel_std' ]    - relative mean standard deviation weights
        model                    - synapse model, can be a list of models
    '''

    rn = rand.normal                                                            # copy function for generation of normal distributed random numbers

    for key, val in params.iteritems():
        exec '%s = %s' % (key, str(val))                                        # params dictionary entries as variables
    
    if isinstance(model, str): model = [ model ]                                # if model is a string put into a list   

    if not d_mu: d_mu = [nest.GetDefaults( m )[ 'delay' ] for m in model ]      # get mean delays d_mu
    d_sigma  = [ d_rel_std*mu for mu in d_mu ]                                  # get delay standard deviation
    
    if not w_mu: w_mu = [nest.GetDefaults( m )[ 'weight' ] for m in model ]     # get mean weights w_mu   
    w_sigma = [ w_rel_std*mu for mu in w_mu ]                                   # get weight standard deviation            
    
    delays  = [ rand.normal( mu, sigma, [ len( pre ), n ] )                     # calculate delays, randomized if sigma != 0
               if sigma else numpy.ones( ( len( pre ), n ) )*mu           
			   for mu, sigma in zip( d_mu, d_sigma )   ]                                                                
      
    weights = [ rand.normal( mu, sigma, [ len( pre ), n ] )                     # calculate weights, randomized if sigma != 0
               if sigma else numpy.ones( ( len( pre ), n ) )*mu            
			   for mu, sigma in zip( w_mu, w_sigma ) ]                                                                  
     
    for i, pre_id in enumerate( pre ):                                                                              
        j = 0
        
         
        if d_sigma[ j ]: d = rn( d_mu[ j ], d_sigma[ j ], [ 1, n ] )[ 0 ]       # if sigma len( targets ) randomized delays                               
        else:  d = numpy.ones( ( 1, n ) )[ 0 ]*d_mu[ j ]                        # else no randomized delays  
                                  
        if w_sigma[ j ]: w = rn( w_mu[ j ], w_sigma[ j ], [ 1, n ] )[ 0 ]       # if signa len( targets ) randomized weights                              
        else: w = numpy.ones( ( 1, n ) )[ 0 ]*w_mu[ j ]                         # else no randomized weights
             
        d, w = list( d ), list( w )                                             # as lists             
        nest.RandomDivergentConnect( [ pre_id ], post, n, weight = w ,          # connect with model[0]
                                     delay = d, model = model[ 0 ] )                  
        
        if len( model ) > 1:
            targets = [ conn['target'] for conn in                              # find connections made by RandomDivergentConnect 
                       nest.GetStatus(nest.FindConnections( [ pre_id ] ) ) 
                       if conn[ 'synapse_type' ] == model[ 0 ]]                                                           
            nt      = len( targets )
            
            for j, m in enumerate( model[ 1 : ], start = 1 ):                   # pre -> post with model[1], model[2] etc  
                if d_sigma[ j ]: d = rn( d_mu[ j ], d_sigma[ j ], [ 1, nt ] )[ 0 ] # if sigma len( targets ) randomized delays 
                else:            d = numpy.ones( ( 1, nt ) )[ 0 ]*d_mu[ j ]     # else no randomized delays   
                                         
                if w_sigma[ j ]: w =rn( w_mu[ j ], w_sigma[ j ], [ 1, nt ] )[ 0 ]  # if signa len( targets ) randomized delays 
                else:            w = numpy.ones( ( 1, nt ) )[ 0 ]*w_mu[ j ]      # else not   
                
                d, w = list( d ), list( w )         
                nest.RandomDivergentConnect( [ pre_id ], targets, n, weight = w, 
                                       delay = d, model = m )                           
                
                
def RCC(pre, post, n, 
        params = { 'd_mu'  : [], 'w_mu'  : [], 'd_rel_std' : 0.0, 'w_rel_std' : 0.0 },  
        model = 'static_synapse' ):
    '''
    As NEST RandomDivergentConnect, except it can take severals models and 
    make same connections for both models. This can be used to connect same
    source with  target with both AMPA and NMDA.Mean and standard deviation 
    for delays and weights can be given. Mean as a list of values one for 
    each model and standard deviation as single values applied to all 
    synapse models.
     
    Inputs:
        pre                      - list of source ids       
        post                     - list of target ids  
        n                        - number of neurons each source neuron connects to
        params[ 'd_mu' ]         - mean delay for each model
        params[ 'w_mu' ]            - mean weight for each model
        params[ 'd_rel_std' ]    - relative mean standard deviation delays
        params[ 'w_rel_std' ]    - relative mean standard deviation weights
        model                    - synapse model, can be a list of models
    '''

    rn = rand.normal                                                            # copy function for generation of normal distributed random numbers

    for key, val in params.iteritems():
        exec '%s = %s' % (key, str(val))                                        # params dictionary entries as variables
    
    if isinstance(model, str): model = [ model ]                                # if model is a string put into a list   

    if not d_mu: d_mu = [nest.GetDefaults( m )[ 'delay' ] for m in model ]      # get mean delays d_mu
    d_sigma  = [ d_rel_std*mu for mu in d_mu ]                                  # get delay standard deviation
    
    if not w_mu: w_mu = [nest.GetDefaults( m )[ 'weight' ] for m in model ]     # get mean weights w_mu   
    w_sigma = [ w_rel_std*mu for mu in w_mu ]                                   # get weight standard deviation            
    
    delays  = [ rand.normal( mu, sigma, [ len( pre ), n ] )                     # calculate delays, randomized if sigma != 0
               if sigma else numpy.ones( ( len( pre ), n ) )*mu           
               for mu, sigma in zip( d_mu, d_sigma )   ]                                                                
      
    weights = [ rand.normal( mu, sigma, [ len( pre ), n ] )                     # calculate weights, randomized if sigma != 0
               if sigma else numpy.ones( ( len( pre ), n ) )*mu            
               for mu, sigma in zip( w_mu, w_sigma ) ]                                                                  
     
    for i, pre_id in enumerate( pre ):                                                                              
        j = 0
        
         
        if d_sigma[ j ]: d = rn( d_mu[ j ], d_sigma[ j ], [ 1, n ] )[ 0 ]       # if sigma len( targets ) randomized delays                               
        else:  d = numpy.ones( ( 1, n ) )[ 0 ]*d_mu[ j ]                        # else no randomized delays  
                                  
        if w_sigma[ j ]: w = rn( w_mu[ j ], w_sigma[ j ], [ 1, n ] )[ 0 ]       # if signa len( targets ) randomized weights                              
        else: w = numpy.ones( ( 1, n ) )[ 0 ]*w_mu[ j ]                         # else no randomized weights
             
        d, w = list( d ), list( w )                                             # as lists             
        nest.RandomDivergentConnect( [ pre_id ], post, n, weight = w ,          # connect with model[0]
                                     delay = d, model = model[ 0 ] )                  
        
        if len( model ) > 1:
            targets = [ conn['target'] for conn in                              # find connections made by RandomDivergentConnect 
                       nest.GetStatus(nest.FindConnections( [ pre_id ] ) ) 
                       if conn[ 'synapse_type' ] == model[ 0 ]]                                                           
            nt      = len( targets )
            
            for j, m in enumerate( model[ 1 : ], start = 1 ):                   # pre -> post with model[1], model[2] etc  
                if d_sigma[ j ]: d = rn( d_mu[ j ], d_sigma[ j ], [ 1, nt ] )[ 0 ] # if sigma len( targets ) randomized delays 
                else:            d = numpy.ones( ( 1, nt ) )[ 0 ]*d_mu[ j ]     # else no randomized delays   
                                         
                if w_sigma[ j ]: w =rn( w_mu[ j ], w_sigma[ j ], [ 1, nt ] )[ 0 ]  # if signa len( targets ) randomized delays 
                else:            w = numpy.ones( ( 1, nt ) )[ 0 ]*w_mu[ j ]      # else not   
                
                d, w = list( d ), list( w )         
                nest.ConvergentConnect( [ pre_id ], targets, weight = w, 
                                       delay = d, model = m )                  
        
           
            
            
            
            
