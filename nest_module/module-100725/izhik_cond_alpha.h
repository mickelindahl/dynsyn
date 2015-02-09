/*
 *  izhik_cond_alpha.h
 *
 *  This file is part of ML_MODULE
 *
 *  Copyright (C) 2005-2009 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 *   Izhikivich simple model described in Dynamical systems in neuro science
 *   2007
 *  (C) 2010 Mikael Lindahl
 *
 */

#ifndef IZHIK_COND_ALPHA_H
#define IZHIK_COND_ALPHA_H

#include "config.h"

#ifdef HAVE_GSL

#include "nest.h"
#include "event.h"
#include "archiving_node.h"
#include "ring_buffer.h"
#include "connection.h"
#include "universal_data_logger.h"
#include "recordables_map.h"

// next line will disappear when all models support multimeter
//#include "analog_data_logger.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

/* BeginDocumentation
Name: izhik_cond_alpha - Izhikevich simple neuron

Description:
izhik_cond_alpha is an implementation of based on  Izhikevich simple model
presented in (Izhikevich 2007). The model has two equation one is the voltage
and the other is a slow recovery variable and alpha shaped synaptic conductancies.


Model equations:
C_mv'=k( V - V_r )( v - V_t ) - u - I
I = g_AMPA( V - E_AMPA) + g_NMDA( V - E_NMDA) + g_GABAA( V - E_GABAA) + I_e
if v_b < V
		b <- b_1 and p <- p_1
else
	  b <- b_2 and p <- p_2

u'=a( b( V - V_r )^p - u + u_const )
 If v > V_peak then V = c and u = u + d

Parameters: 
The following parameters can be set in the status dictionary.

Neuron parameters:
Voltage equation parameters
V_m        double - Membrane potential in mV 
V_r        double - Membrane resting potential
V_t			   double - Instantaneous activation threshold
C_m        double - Capacity of the membrane in pF
k					 double - factor determining instantaneous I-V relation
c					 double - Reset value membrane potential
V_peak     double - spike cut of value
I_e        double - Constant input current in pA.

Recovery variable equations
u          double - Recovery variable
V_b     	 double - Recovery variable voltage threshold
a			     double - Time constant slow dynamics
b_1				 double - Slope factor slow dynamics voltage interval [-inf v_b]
b_2				 double - Slope factor slow dynamics interval [v_b +inf]
p_1 			 double - Recovery variable polynomial voltage dependency factor  interval  [-inf v_b]
p_2 			 double - Recovery variable polynomial voltage dependency factor  interval  [v_b +inf]
d					 double - slow variable change at spike
u_const    double - Constant input current in pA to recovery variable u.

Synapse parameters
E_AMPA     double - AMPA reversal potential in mV.
E_NMDA     double - NMDA reversal potential in mV.
E_GABAA    double - GABAA reversal potential in mV.
tau_AMPA   double - Rise time of the AMPA synaptic alpha function in ms.
tau_NMDA   double - Rise time of the NMDA synaptic alpha function in ms.
tau_GABAA  double - Rise time of the GABAAsynaptic alpha function in ms.
NMDA_Sact  double - For voltage dependence of NMDA-synapse mV, see eq. above
NMDA_Vact  double - For voltage dependence of NMDA-synapse mV, see eq. above

Sends: SpikeEvent

Receives: SpikeEvent, CurrentEvent, DataLoggingRequest

Author: Mikael, Lindahl

SeeAlso: iaf_cond_exp, iaf_cond_alpha, iaf_cond_alpha_mc, ht_neuron
*/

namespace mynest
{
  /**
   * Function computing right-hand side of ODE for GSL solver.
   * @note Must be declared here so we can befriend it in class.
   * @note Must have C-linkage for passing to GSL. Internally, it is
   *       a first-class C++ function, but cannot be a member function
   *       because of the C-linkage.
   * @note No point in declaring it inline, since it is called
   *       through a function pointer.
   * @param void* Pointer to model neuron instance.
   */
  extern "C"
  int izhik_cond_alpha_dynamics (double, const double*, double*, void*);

  /**
   * Integrate-and-fire neuron model with two conductance-based synapses.
   *
   * @note Per 2009-04-17, this class has been revised to our newest 
   *       insights into class design. Please use THIS CLASS as a reference
   *       when designing your own models with nonlinear dynamics.
   *       One weakness of this class is that it distinguishes between
   *       inputs to the two synapses by the sign of the synaptic weight.
   *       It would be better to use receptor_types, cf izhik_cond_alpha_mc.
   */
  class izhik_cond_alpha : public nest::Archiving_Node
  {
    
    // Boilerplate function declarations --------------------------------

  public:
    
    izhik_cond_alpha();
    izhik_cond_alpha(const izhik_cond_alpha&);
    ~izhik_cond_alpha();

    /*
     * Import all overloaded virtual functions that we
     * override in this class.  For background information, 
     * see http://www.gotw.ca/gotw/005.htm.
     */
#ifndef IS_BLUEGENE
    using nest::Node::check_connection;
#endif
    using nest::Node::connect_sender;
    using nest::Node::handle;

    nest::port check_connection(nest::Connection&, nest::port);

    nest::port connect_sender(nest::SpikeEvent &, nest::port);
    nest::port connect_sender(nest::CurrentEvent &, nest::port);
    nest::port connect_sender(nest::DataLoggingRequest &, nest::port);

    void handle(nest::SpikeEvent &);
    void handle(nest::CurrentEvent &);
    void handle(nest::DataLoggingRequest &);
        
    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &);

  private:
    void init_node_(const Node& proto); // No nest:: here
    void init_state_(const Node& proto); // No nest:: here
    void init_buffers_();
    void calibrate();
    void update(nest::Time const &, const nest::long_t, const nest::long_t);

    // END Boilerplate function declarations ----------------------------

    // Enumerations and constants specifying structure and properties ----

    /**
     * Minimal spike receptor type.
     * @note Start with 1 so we can forbid port 0 to avoid accidental
     *       creation of connections with no receptor type set.
     */
    static const nest::port MIN_SPIKE_RECEPTOR = 1;

    /**
     * Spike receptors (SUP_SPIKE_RECEPTOR=3).
     */
    // Three spike receptors AMPA, NMDA and GABAA
    enum SpikeSynapseTypes { AMPA=MIN_SPIKE_RECEPTOR, NMDA, GABAA,
    	SUP_SPIKE_RECEPTOR };

    static const nest::size_t NUM_SPIKE_RECEPTORS = SUP_SPIKE_RECEPTOR - MIN_SPIKE_RECEPTOR;

    /**
     * Minimal current receptor type.
     *  @note Start with SUP_SPIKE_RECEPTOR to avoid any overlap and
     *        accidental mix-ups.
     */
    static const nest::port MIN_CURR_RECEPTOR = SUP_SPIKE_RECEPTOR;

    /**
     * Current receptors (SUP_CURR_RECEPTOR = 5).
     */
    enum CurrentSynapseTypes { I_CURR = MIN_CURR_RECEPTOR,
    	SUP_CURR_RECEPTOR };

    static const nest::size_t NUM_CURR_RECEPTORS = SUP_CURR_RECEPTOR - MIN_CURR_RECEPTOR;

    // Enumerations and constants specifying structure and properties ----


    // Friends --------------------------------------------------------

    // make dynamics function quasi-member. Add mynest since this is a function
    //in namespace mynest.
    friend int mynest::izhik_cond_alpha_dynamics(double, const double*, double*, void*);


    // The next two classes need to be friends to access the State_ class/member
    // Add nest:: since these are nest classes in namespace nest
    friend class nest::RecordablesMap<izhik_cond_alpha>;
    friend class nest::UniversalDataLogger<izhik_cond_alpha>;

  private:

    // Parameters class ------------------------------------------------- 

    //! Model parameters
    struct Parameters_ {

    	// Voltage equation parameters
    	nest::double_t V_m; 					//!< Membrane potential in mV
    	nest::double_t V_r;        //!< Membrane resting potential
    	nest::double_t V_t;			   //!< Instantaneous activation threshold
    	nest::double_t C_m;        //!< Capacity of the membrane in pF
    	nest::double_t k;					 //!< Factor determining instantaneous I-V relation
    	nest::double_t c;					 //!< Reset value membrane potential
    	nest::double_t V_peak;     //!< spike cut of value
    	nest::double_t I_e;        //!< Constant input current in pA.

    	// Recovery variable equations
    	nest::double_t u; 					//!< Recovery variable
     	nest::double_t V_b;			   //!< Recovery variable voltage threshold
    	nest::double_t a;			     //!< Time constant slow dynamics
    	nest::double_t b_1;				 //!< Slope factor slow dynamics voltage interval [-inf v_b]
     	nest::double_t b_2;				 //!< Slope factor slow dynamics interval [v_b +inf]
     	nest::double_t p_1;				 //!< Recovery variable polynomial voltage dependency factor  interval  [-inf v_b]
     	nest::double_t p_2;				 //!< Recovery variable voltage dependency factor interval 1 [v_b +inf]
     	nest::double_t d;					 //!< slow variable change at spike
    	nest::double_t u_const;    //!< Constant input current in pA to recovery variable u.

    	// Synaptic parameters
      nest::double_t E_AMPA;      //!< AMPA reversal Potential in mV
      nest::double_t E_NMDA;      //!< NMDA reversal Potential in mV
      nest::double_t E_GABAA;     //!< GABAA reversal Potential in mV
      nest::double_t tau_synAMPA;    //!< Synaptic Time Constant AMPA Synapse in ms
      nest::double_t tau_synNMDA;    //!< Synaptic Time Constant NMDA Synapse in ms
      nest::double_t tau_synGABAA;    //!< Synaptic Time Constant GABAA Synapse in ms
      nest::double_t NMDA_Vact;  //!< mV, inactive for V << Vact, inflection of sigmoid
      nest::double_t NMDA_Sact;  //!< mV, scale of inactivation
  
      Parameters_();        //!< Set default parameter values

      void get(DictionaryDatum&) const;  //!< Store current values in dictionary
      void set(const DictionaryDatum&);  //!< Set values from dicitonary
    };
    
    // State variables class -------------------------------------------- 

    /**
     * State variables of the model.
     * 
     * State variables consist of the state vector for the subthreshold
     * dynamics and the refractory count. The state vector must be a
     * C-style array to be compatible with GSL ODE solvers.
     *
     * @note Copy constructor and assignment operator are required because
     *       of the C-style array.
     */
    struct State_ {
      
      //! Symbolic indices to the elements of the state vector y
      enum StateVecElems_ { V_M = 0, u,
														DG_AMPA, G_AMPA,
														DG_NMDA, G_NMDA,
														DG_GABAA, G_GABAA,
														STATE_VEC_SIZE };

      //! state vector, must be C-array for GSL solver
      nest::double_t y[STATE_VEC_SIZE];

      //!< number of refractory steps remaining
      //nest::int_t    r;

      // Contructors, copy-constructors and destructors has to be defined in
      // .cpp in order to retrieve them.
      nest::double_t I_;  		 //!< Total input current from synapses and current injection
      nest::double_t I_AMPA_;  //!< AMPA current; member only to allow recording
      nest::double_t I_NMDA_;  //!<NMDA current; member only to allow recording
      nest::double_t I_GABAA_; //!< GABAA current; member only to allow recording

      // Taken from ht_neuron use to be only State_(const Parameters_&) and
      // State_(const State_&). Changed since compiler complained that I_AMPA,
      // I_NMDA and I_GABAA read-only structure
      State_();
      State_(const Parameters_& p);
      State_(const State_& s);
      ~State_();

      State_& operator=(const State_& s);

      void get(DictionaryDatum&) const;  //!< Store current values in dictionary

      /**
       * Set state from values in dictionary.
       * Requires Parameters_ as argument to, eg, check bounds.'
       */
      void set(const DictionaryDatum&, const Parameters_&);
    };    

    // Buffers class -------------------------------------------------------- 

    /**
     * Buffers of the model.
     * Buffers are on par with state variables in terms of persistence,
     * i.e., initalized only upon first Simulate call after ResetKernel
     * or ResetNetwork, but are implementation details hidden from the user.
     */
    struct Buffers_ {
      Buffers_(izhik_cond_alpha&); //!<Sets buffer pointers to 0
      Buffers_(const Buffers_&, izhik_cond_alpha&); //!<Sets buffer pointers to 0

      //! Logger for all analog data
      nest::UniversalDataLogger<izhik_cond_alpha> logger_;

      /** buffers and sums up incoming spikes/currents */
      // One ring buffer for each synapse. Need to register this in receptor
      // dictionary.
      nest::RingBuffer spikes_AMPA_;
      nest::RingBuffer spikes_NMDA_;
      nest::RingBuffer spikes_GABAA_;
      nest::RingBuffer currents_;

      /* GSL ODE stuff */
      gsl_odeiv_step*    s_;    //!< stepping function
      gsl_odeiv_control* c_;    //!< adaptive stepsize control function
      gsl_odeiv_evolve*  e_;    //!< evolution function
      gsl_odeiv_system   sys_;  //!< struct describing system
      
      // IntergrationStep_ should be reset with the neuron on ResetNetwork,
      // but remain unchanged during calibration. Since it is initialized with
      // step_, and the resolution cannot change after nodes have been created,
      // it is safe to place both here.
      nest::double_t step_;           //!< step size in ms
      double   IntegrationStep_;//!< current integration time step, updated by GSL
    };
    
    // Variables class ------------------------------------------------------- 
    
    /**
     * Internal variables of the model.
     * Variables are re-initialized upon each call to Simulate.
     */
    struct Variables_ { 
      /**
       * Impulse to add to DG_AMPA on spike arrival to evoke unit-amplitude
       * conductance excursion.
       */
    	nest::double_t PSConInit_AMPA;

      /**
    	* Impulse to add to DG_NMDA on spike arrival to evoke unit-amplitude
    	* conductance excursion.
    	*/
    	 nest::double_t PSConInit_NMDA;
      
      /**
       * Impulse to add to DG_GABAA on spike arrival to evoke unit-amplitude
       * conductance excursion.
       */
    	nest::double_t PSConInit_GABAA;
      
      //! refractory time in steps
    	//nest::int_t    RefractoryCounts;
      
      //! make external input current available to dynamics function
    	nest::double_t I_CURR;
    };
    
    // Access functions for UniversalDataLogger -------------------------------
    
    //! Read out state vector elements and currents, used by UniversalDataLogger
    template <State_::StateVecElems_ elem>
    nest::double_t get_y_elem_() const { return S_.y[elem]; }
    nest::double_t get_I_() const { return S_.I_; }
    nest::double_t get_I_AMPA_() const { return S_.I_AMPA_; }
    nest::double_t get_I_NMDA_() const { return S_.I_NMDA_; }
    nest::double_t get_I_GABAA_() const { return S_.I_GABAA_; }
    
    //! Read out remaining refractory time, used by UniversalDataLogger
    //nest::double_t get_r_() const { return nest::Time::get_resolution().get_ms() * S_.r; }
    
    // Data members ----------------------------------------------------------- 

    // keep the order of these lines, seems to give best performance
    Parameters_ P_;
    State_      S_;
    Variables_  V_;
    Buffers_    B_;

    //! Mapping of recordables names to access functions
    static nest::RecordablesMap<izhik_cond_alpha> recordablesMap_;
  };
  

  // Boilerplate inline function definitions ----------------------------------

  inline
  nest::port mynest::izhik_cond_alpha::check_connection(nest::Connection& c, nest::port receptor_type)
  {
  	nest::SpikeEvent e;
    e.set_sender(*this);
    c.check_event(e);
    return c.get_target()->connect_sender(e, receptor_type);
  }

  inline
  nest::port mynest::izhik_cond_alpha::connect_sender(nest::SpikeEvent&, nest::port receptor_type)
  {
  	// If receptor type is less than 1 =(MIN_SPIKE_RECEPTOR) or greater or equal to 4
  	// (=SUP_SPIKE_RECEPTOR) then provided receptor type is not a spike receptor.
	  if ( receptor_type < MIN_SPIKE_RECEPTOR || receptor_type >= SUP_SPIKE_RECEPTOR )
	  	// Unknown receptor type is less than 0 or greater than 6
	  	// (SUP_CURR_RECEPTOR).
		  if ( receptor_type < 0 || receptor_type >= SUP_CURR_RECEPTOR )
			  throw nest::UnknownReceptorType(receptor_type, get_name());
	    // Otherwise it is a current receptor or receptor 0 (data logging request
	    // not used here and therefore incompatible.
		  else
			  throw nest::IncompatibleReceptorType(receptor_type, get_name(), "SpikeEvent");
	  // If we arrive here the receptor type is a spike receptor and either 1, 2 or 3 e.i.
	  // greater or equal to MIN_SPIKE_RECEPTOR = 1, and less than SUP_SPIKE_RECEPTOR
	  // = 4. Then 0, 1, or 2 is returned.
	  return receptor_type - MIN_SPIKE_RECEPTOR;
  }
 
  inline
  nest::port mynest::izhik_cond_alpha::connect_sender(nest::CurrentEvent&, nest::port receptor_type)
  {
  	// If receptor type is less than 4 (MIN_CURR_RECEPTOR) or greater or equal
  	// to 5 (SUP_CURR_RECEPTOR) the provided receptor type is not current
  	// receptor.
	  if ( receptor_type < MIN_CURR_RECEPTOR || receptor_type >= SUP_CURR_RECEPTOR )
	  	// If receptor is not a current receptor but still a receptor type that is
	  	// the receptor type is greater or equal to 0 or less than 3
	  	// (MIN_CURR_RECEPTOR).
	  	if ( receptor_type >= 0 && receptor_type < MIN_CURR_RECEPTOR )
			  throw nest::IncompatibleReceptorType(receptor_type, get_name(), "CurrentEvent");
		  // Otherwise unknown receptor type.
	  	else
			  throw nest::UnknownReceptorType(receptor_type, get_name());
	  //MIN_CURR_RECEPTOR =4, If here receptor type equals 4  and 0 is returned.
	  return receptor_type - MIN_CURR_RECEPTOR;
  }
 
  inline
  nest::port mynest::izhik_cond_alpha::connect_sender(nest::DataLoggingRequest& dlr,
		  nest::port receptor_type)
  {
  	// If receptor type does not equal 0 then it is not a data logging request
  	// receptor.
    if ( receptor_type != 0 )
    	// If not a spike or current receptor that is less than 0 or greater or
    	//  equal to 4 (SUP_CURR_RECEPTOR).
      if ( receptor_type < 0 || receptor_type >= SUP_CURR_RECEPTOR )
        throw nest::UnknownReceptorType(receptor_type, get_name());
      // Otherwise it is a spike or current receptor type.
      else
        throw nest::IncompatibleReceptorType(receptor_type, get_name(), "DataLoggingRequest");
    B_.logger_.connect_logging_device(dlr, recordablesMap_);
    return 0;
  }

  inline
  void izhik_cond_alpha::get_status(DictionaryDatum &d) const
  {
    P_.get(d);
    S_.get(d);
    nest::Archiving_Node::get_status(d);

    (*d)[nest::names::recordables] = recordablesMap_.get_list();

    /**
     * @TODO dictionary construction should be done only once for
     * static member in default c'tor, but this leads to
     * a seg fault on exit, see #328
     */
    DictionaryDatum receptor_dict_ = new Dictionary();
    (*receptor_dict_)[Name("AMPA")]  = AMPA;
    (*receptor_dict_)[Name("NMDA")]  = NMDA;
    (*receptor_dict_)[Name("GABAA")]  = GABAA;
    (*receptor_dict_)[Name("I_CURR")] = I_CURR;

    (*d)[nest::names::receptor_types] = receptor_dict_;

  }

  inline
  void izhik_cond_alpha::set_status(const DictionaryDatum &d)
  {
    Parameters_ ptmp = P_;  // temporary copy in case of errors
    ptmp.set(d);                       // throws if BadProperty
    State_      stmp = S_;  // temporary copy in case of errors
    stmp.set(d, ptmp);                 // throws if BadProperty

    // We now know that (ptmp, stmp) are consistent. We do not 
    // write them back to (P_, S_) before we are also sure that 
    // the properties to be set in the parent class are internally 
    // consistent.
    nest::Archiving_Node::set_status(d);

    // if we get here, temporaries contain consistent set of properties
    P_ = ptmp;
    S_ = stmp;
  }

} // namespace


#endif //HAVE_GSL
#endif //IZHIK_COND_ALPHA_H
