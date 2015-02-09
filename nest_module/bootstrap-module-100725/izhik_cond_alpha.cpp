/*
 *  izhik_cond_alpha.cpp
 *
 *  This file is part of NEST
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
 */


#include "izhik_cond_alpha.h"

#ifdef HAVE_GSL

#include "exceptions.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
//#include "analog_data_logger_impl.h"
#include "universal_data_logger_impl.h"
#include <limits>

#include <iomanip>
#include <iostream>
#include <cstdio>

using namespace nest; //added

/* ---------------------------------------------------------------- 
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap<mynest::izhik_cond_alpha> mynest::izhik_cond_alpha::recordablesMap_;

namespace nest   // template specialization must be placed in namespace
{
  /*
   * Override the create() method with one call to RecordablesMap::insert_() 
   * for each quantity to be recorded.
   */
  template <>
  void RecordablesMap<mynest::izhik_cond_alpha>::create()
  {
  	// use standard names whereever you can for consistency!
  	// Recording current seems
  	insert_(names::V_m,
  			&mynest::izhik_cond_alpha::get_y_elem_<mynest::izhik_cond_alpha::State_::V_M>);
  	insert_(Name("u"),
  			&mynest::izhik_cond_alpha::get_y_elem_<mynest::izhik_cond_alpha::State_::u>);
  	insert_(Name("g_AMPA"),
  			&mynest::izhik_cond_alpha::get_y_elem_<mynest::izhik_cond_alpha::State_::G_AMPA>);
  	insert_(Name("g_NMDA"),
  			&mynest::izhik_cond_alpha::get_y_elem_<mynest::izhik_cond_alpha::State_::G_NMDA>);
  	insert_(Name("g_GABAA"),
  			&mynest::izhik_cond_alpha::get_y_elem_<mynest::izhik_cond_alpha::State_::G_GABAA>);
    insert_(Name("I"), &mynest::izhik_cond_alpha::get_I_);
  	insert_(Name("I_AMPA"), &mynest::izhik_cond_alpha::get_I_AMPA_);
    insert_(Name("I_NMDA"), &mynest::izhik_cond_alpha::get_I_NMDA_);
    insert_(Name("I_GABAA"), &mynest::izhik_cond_alpha::get_I_GABAA_);
    //insert_(names::t_ref_remaining,
	  //  &mynest::izhik_cond_alpha::get_r_);
  }
}

/* ---------------------------------------------------------------- 
 * Iteration function
 * ---------------------------------------------------------------- */

extern "C"
inline int mynest::izhik_cond_alpha_dynamics(double, const double y[], double f[], void* pnode)
{ 
  // some shorthands
  //typedef mynest::izhik_cond_alpha         N;
  typedef mynest::izhik_cond_alpha::State_ S;

  // get access to node so we can almost work as in a member class
  assert(pnode);
  mynest::izhik_cond_alpha& node =  *(reinterpret_cast<mynest::izhik_cond_alpha*>(pnode));

  // easier access to membrane potential and recovery variable
  const nest::double_t& V = y[S::V_M];
  const nest::double_t& u = y[S::u];

  // y[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.y[]. 
  
  // The following code is verbose for the sake of clarity. We assume that a
  // good compiler will optimize the verbosity away.
  const nest::double_t I_AMPA = - y[S::G_AMPA] * ( V - node.P_.E_AMPA );
  const nest::double_t I_NMDA = - y[S::G_NMDA] * ( V - node.P_.E_NMDA )
      / ( 1 + std::exp( (node.P_.NMDA_Vact - V)/node.P_.NMDA_Sact ) );
  const nest::double_t I_GABAA = - y[S::G_GABAA] * ( V - node.P_.E_GABAA );

  // Set state variable used for recording AMPA, NMDA and GABAA current
  // contributions
  node.S_.I_AMPA_ = I_AMPA;
  node.S_.I_NMDA_ = I_NMDA;
  node.S_.I_GABAA_ = I_GABAA;

  // Total input current from synapses and external input
  node.S_.I_ = I_AMPA + I_NMDA + I_GABAA + node.V_.I_CURR + node.P_.I_e;

  // Neuron dynamics
  // dV_m/dt
  f[0]= ( node.P_.k*( V - node.P_.V_r )*(  V - node.P_.V_t ) - u + node.S_.I_ )/node.P_.C_m;

  // du/dt
  if (node.P_.V_b < V)
  	f[1]= node.P_.a*( node.P_.b_1*std::pow( V - node.P_.V_b, node.P_.p_1 ) - u + node.P_.u_const);
  else
  	f[1]= node.P_.a*( node.P_.b_2*std::pow( V - node.P_.V_b, node.P_.p_2 ) - u + node.P_.u_const);


  // Synapse dynamics
  // d dg_AMPA/dt, dg_AMPA/dt
  f[2] = -y[S::DG_AMPA] / node.P_.tau_synAMPA;
  f[3] =  y[S::DG_AMPA] - (y[S::G_AMPA]/node.P_.tau_synAMPA);

  // d dg_NMDA/dt, dg_NMDA/dt
  f[4] = -y[S::DG_NMDA] / node.P_.tau_synNMDA;
  f[5] =  y[S::DG_NMDA] - (y[S::G_NMDA]/node.P_.tau_synNMDA);

  // d dg_GABAA/dt, dg_GABAA/dt
  f[6] = -y[S::DG_GABAA] / node.P_.tau_synGABAA;
  f[7] =  y[S::DG_GABAA] - (y[S::G_GABAA]/node.P_.tau_synGABAA);

  return GSL_SUCCESS;
 }

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */
    
mynest::izhik_cond_alpha::Parameters_::Parameters_()
  : V_m			( -50.0  ), 					//!< Membrane potential in mV
		V_r			( -50.0  ),        //!< Membrane resting potential
		V_t			( -40.0  ),			   //!< Instantaneous activation threshold
		C_m			(  100.0 ),        //!< Capacity of the membrane in pF
		k				(  0.7   ),					 //!< Factor determining instantaneous I-V relation
		c				(  -50.0 ),					 //!< Reset value membrane potential
		V_peak	(  30.0  ),     //!< spike cut of value
    I_e     (  0.0    ),   // pA
		u			  (   0.0  ), 					//!< Recovery variable pA
   	V_b     ( -50.0  ),			   //!< Recovery variable voltage threshold mV
		a				(  0.03  ),			     //!< Time constant slow dynamics
		b_1			( -2.0   ),					 //!< Slope factor 1 slow dynamics
		b_2			( -2.0   ),					 //!< Slope factor 2 slow dynamics
		p_1			( 1.0   ),					 //!< Polynomial voltage dependency factor 1
		p_2			( 1.0   ),					 //!< Polynomial voltage dependency factor 2
		u_const ( 0.0   ),					 //!< Constant current recovery variable
		d				(  100.0 ),					 //!< slow variable change at spike
    E_AMPA    (  0.0    ),  // mV
    E_NMDA    (  0.0    ),  // mV
    E_GABAA    (-70.0    ),  // mV
    tau_synAMPA(  3.0    ),  // ms
    tau_synNMDA(  100.0    ),  // ms
    tau_synGABAA(  3.0    ),  // ms
		NMDA_Vact      (-58.0 ), // mV
		NMDA_Sact      (  2.5 ) // mV
{
  recordablesMap_.create();
}

mynest::izhik_cond_alpha::State_::State_(const Parameters_& p)
  : I_(0.0),
		I_AMPA_(0.0),
    I_NMDA_(0.0),
    I_GABAA_(0.0)
{
  y[V_M] = p.V_r;  // initialize to resting potential
  y[u] = p.u;   // initialize to u
  for ( size_t i = 2 ; i < STATE_VEC_SIZE ; ++i )
    y[i] = 0;
}

mynest::izhik_cond_alpha::State_::State_(const State_& s)
  : I_(0.0),
		I_AMPA_(s.I_AMPA_),
    I_NMDA_(s.I_NMDA_),
    I_GABAA_(s.I_GABAA_)
{
  for ( size_t i = 0 ; i < STATE_VEC_SIZE ; ++i )
    y[i] = s.y[i];
}

mynest::izhik_cond_alpha::State_& mynest::izhik_cond_alpha::State_::operator=(const State_& s)
{
  if ( this == &s )  // avoid assignment to self
    return *this;

  for ( size_t i = 0 ; i < STATE_VEC_SIZE ; ++i )
    y[i] = s.y[i];

  I_ = s.I_;
  I_AMPA_ = s.I_AMPA_;
  I_NMDA_ = s.I_NMDA_;
  I_GABAA_ = s.I_GABAA_;
  //r = s.r;
  return *this;
}

mynest::izhik_cond_alpha::State_::~State_()
{
}

mynest::izhik_cond_alpha::Buffers_::Buffers_(izhik_cond_alpha& n)
  : logger_(n),
    s_(0),
    c_(0),
    e_(0)
{
  // The other member variables are left uninitialised or are
  // automatically initialised by their default constructor.
}

mynest::izhik_cond_alpha::Buffers_::Buffers_(const Buffers_&, izhik_cond_alpha& n)
  : logger_(n),
    s_(0),
    c_(0),
    e_(0)
{
  // The other member variables are left uninitialised or are
  // automatically initialised by their default constructor.
}

/* ---------------------------------------------------------------- 
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void mynest::izhik_cond_alpha::Parameters_::get(DictionaryDatum &dd) const
{
  def<double>(dd,names::V_m,          V_m);
  def<nest::double_t>(dd, "V_r",      V_r);
  def<nest::double_t>(dd, "V_t",      V_t);
  def<double>(dd,names::C_m,          C_m);
  def<nest::double_t>(dd, "c",        c);
  def<nest::double_t>(dd, "k",        k);
  def<nest::double_t>(dd, "V_peak",   V_peak);
  def<double>(dd,names::I_e,          I_e);

  def<nest::double_t>(dd, "u",      	u);
  def<nest::double_t>(dd, "V_b",      V_b);
  def<nest::double_t>(dd, "a",        a);
  def<nest::double_t>(dd, "b_1",      b_1);
  def<nest::double_t>(dd, "b_2",      b_2);
  def<nest::double_t>(dd, "p_1",      p_1);
  def<nest::double_t>(dd, "p_2",      p_2);
  def<nest::double_t>(dd, "u_const",  u_const);
  def<nest::double_t>(dd, "d",        d);

  def<nest::double_t>(dd, "E_AMPA",          E_AMPA);
  def<nest::double_t>(dd, "E_NMDA",          E_NMDA);
  def<nest::double_t>(dd, "E_GABAA",         E_GABAA);
  def<nest::double_t>(dd, "tau_synAMPA",          tau_synAMPA);
  def<nest::double_t>(dd, "tau_synNMDA",          tau_synNMDA);
  def<nest::double_t>(dd, "tau_synGABAA",         tau_synGABAA);

  def<nest::double_t>(dd, "NMDA_Vact",     NMDA_Vact);
  def<nest::double_t>(dd, "NMDA_Sact",     NMDA_Sact);

}

void mynest::izhik_cond_alpha::Parameters_::set(const DictionaryDatum& dd)
{
  // allow setting the membrane potential
	updateValue<double>(dd,names::V_m,          V_m);
  updateValue<nest::double_t>(dd, "V_r",      V_r);
  updateValue<nest::double_t>(dd, "V_t",      V_t);
  updateValue<double>(dd,names::C_m,          C_m);
  updateValue<nest::double_t>(dd, "c",        c);
  updateValue<nest::double_t>(dd, "k",        k);
  updateValue<nest::double_t>(dd, "V_peak",   V_peak);
  updateValue<double>(dd,names::I_e,          I_e);

  updateValue<nest::double_t>(dd, "u",        u);
  updateValue<nest::double_t>(dd, "V_b",      V_b);
  updateValue<nest::double_t>(dd, "a",        a);
  updateValue<nest::double_t>(dd, "b_1",      b_1);
  updateValue<nest::double_t>(dd, "b_2",      b_2);
  updateValue<nest::double_t>(dd, "p_1",      p_1);
  updateValue<nest::double_t>(dd, "p_2",      p_2);
  updateValue<nest::double_t>(dd, "u_const",  u_const);
  updateValue<nest::double_t>(dd, "d",        d);

  updateValue<nest::double_t>(dd, "E_AMPA",          E_AMPA);
  updateValue<nest::double_t>(dd, "E_NMDA",          E_NMDA);
  updateValue<nest::double_t>(dd, "E_GABAA",         E_GABAA);
  updateValue<nest::double_t>(dd, "tau_synAMPA",          tau_synAMPA);
  updateValue<nest::double_t>(dd, "tau_synNMDA",          tau_synNMDA);
  updateValue<nest::double_t>(dd, "tau_synGABAA",         tau_synGABAA);
  updateValue<nest::double_t>(dd, "NMDA_Vact",     NMDA_Vact);
  updateValue<nest::double_t>(dd, "NMDA_Sact",     NMDA_Sact);


  if ( V_r >= V_t )
    throw BadProperty("Reset potential must be smaller than instantaneous threshold.");
    
  if ( C_m <= 0 )
    throw BadProperty("Capacitance must be strictly positive.");
    
  //if ( t_ref < 0 )
  //  throw BadProperty("Refractory time cannot be negative.");
      
  if ( tau_synAMPA <= 0 || tau_synNMDA <= 0 || tau_synGABAA <= 0 )
    throw BadProperty("All time constants must be strictly positive.");
}

void mynest::izhik_cond_alpha::State_::get(DictionaryDatum &d) const
{
  def<double>(d, names::V_m, y[V_M]); // Membrane potential
  def<nest::double_t>(d, "u", y[u]); // Recovery variable
}

void mynest::izhik_cond_alpha::State_::set(const DictionaryDatum& d, const Parameters_&)
{
  updateValue<double>(d, names::V_m, y[V_M]);
  updateValue<nest::double_t>(d, "u", y[u]); // Recovery variable
}


/* ---------------------------------------------------------------- 
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

mynest::izhik_cond_alpha::izhik_cond_alpha()
  : Archiving_Node(), 
    P_(), 
    S_(P_),
    B_(*this)
{
}

mynest::izhik_cond_alpha::izhik_cond_alpha(const izhik_cond_alpha& n)
  : Archiving_Node(n), 
    P_(n.P_), 
    S_(n.S_),
    B_(n.B_, *this)
{
}

mynest::izhik_cond_alpha::~izhik_cond_alpha()
{
  // GSL structs only allocated by init_nodes_(), so we need to protect destruction
  if ( B_.s_ ) gsl_odeiv_step_free(B_.s_);
  if ( B_.c_ ) gsl_odeiv_control_free(B_.c_);
  if ( B_.e_ ) gsl_odeiv_evolve_free(B_.e_);
}

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void mynest::izhik_cond_alpha::init_node_(const Node& proto)
{
  const izhik_cond_alpha& pr = downcast<izhik_cond_alpha>(proto);
  P_ = pr.P_;
  S_ = pr.S_;
}

void mynest::izhik_cond_alpha::init_state_(const Node& proto)
{
  const izhik_cond_alpha& pr = downcast<izhik_cond_alpha>(proto);
  S_ = pr.S_;
}

void mynest::izhik_cond_alpha::init_buffers_()
{
  Archiving_Node::clear_history();

  B_.spikes_AMPA_.clear();       // includes resize
  B_.spikes_NMDA_.clear();       // includes resize
  B_.spikes_GABAA_.clear();       // includes resize
  B_.currents_.clear();        // includes resize

  B_.logger_.init();

  nest::Archiving_Node::clear_history();

  B_.step_ = Time::get_resolution().get_ms();
  B_.IntegrationStep_ = B_.step_;

  static const gsl_odeiv_step_type* T1 = gsl_odeiv_step_rkf45;
  
  if ( B_.s_ == 0 )
    B_.s_ = gsl_odeiv_step_alloc (T1, State_::STATE_VEC_SIZE);
  else 
    gsl_odeiv_step_reset(B_.s_);
    
  if ( B_.c_ == 0 )  
    B_.c_ = gsl_odeiv_control_y_new (1e-3, 0.0);
  else
    gsl_odeiv_control_init(B_.c_, 1e-3, 0.0, 1.0, 0.0);
    
  if ( B_.e_ == 0 )  
    B_.e_ = gsl_odeiv_evolve_alloc(State_::STATE_VEC_SIZE);
  else 
    gsl_odeiv_evolve_reset(B_.e_);
  
  B_.sys_.function  = izhik_cond_alpha_dynamics;
  B_.sys_.jacobian  = NULL;
  B_.sys_.dimension = State_::STATE_VEC_SIZE;
  B_.sys_.params    = reinterpret_cast<void*>(this);
}

void mynest::izhik_cond_alpha::calibrate()
{
  V_.PSConInit_AMPA  = 1.0 * numerics::e / P_.tau_synAMPA;
  V_.PSConInit_NMDA  = 1.0 * numerics::e / P_.tau_synNMDA;
  V_.PSConInit_GABAA = 1.0 * numerics::e / P_.tau_synGABAA;
  //V_.RefractoryCounts = Time(Time::ms(P_.t_ref)).get_steps();
  
  //assert(V_.RefractoryCounts >= 0);  // since t_ref >= 0, this can only fail in error
}

/* ---------------------------------------------------------------- 
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void mynest::izhik_cond_alpha::update(Time const & origin, const nest::long_t from, const nest::long_t to)
{
   
  assert(to >= 0 && (delay) from < Scheduler::get_min_delay());
  assert(from < to);

  for ( nest::long_t lag = from ; lag < to ; ++lag )
  {
    
    double t = 0.0;

    // numerical integration with adaptive step size control:
    // ------------------------------------------------------
    // gsl_odeiv_evolve_apply performs only a single numerical
    // integration step, starting from t and bounded by step;
    // the while-loop ensures integration over the whole simulation
    // step (0, step] if more than one integration step is needed due
    // to a small integration step size;
    // note that (t+IntegrationStep > step) leads to integration over
    // (t, step] and afterwards setting t to step, but it does not
    // enforce setting IntegrationStep to step-t; this is of advantage
    // for a consistent and efficient integration across subsequent
    // simulation intervals
    while ( t < B_.step_ )
    {
      const int status = gsl_odeiv_evolve_apply(B_.e_, B_.c_, B_.s_, 
			   &B_.sys_,             // system of ODE
			   &t,                   // from t
			    B_.step_,            // to t <= step
			   &B_.IntegrationStep_, // integration step size
			    S_.y); 	         // neuronal state

      if ( status != GSL_SUCCESS )
        throw GSLSolverFailure(get_name(), status);
    }

    // refractoriness and spike generation 
    //if ( S_.r )
    //{// neuron is absolute refractory
	  //  --S_.r;
	  //  S_.y[State_::V_M] = P_.V_reset;  // clamp potential
    //}
    //else
      // neuron is not absolute refractory
    if ( S_.y[State_::V_M] >= P_.V_peak )
      {
    	// v is reset to c (v=c)
    	S_.y[State_::V_M] = P_.c;

		  //S_.r              = V_.RefractoryCounts;

    	// d is added to u (u=u + d)
    	S_.y[State_::u] += P_.d;

    	// log spike with Archiving_Node
    	set_spiketime(Time::step(origin.get_steps()+lag+1));
	
    	SpikeEvent se;
    	network()->send(*this, se, lag);
      }

    // Here incomming spikes are added. It is setup such that AMPA, NMDA
    // and GABA recieves from one receptor each.
    S_.y[State_::DG_AMPA] += B_.spikes_AMPA_.get_value(lag)* V_.PSConInit_AMPA;
    S_.y[State_::DG_NMDA] += B_.spikes_NMDA_.get_value(lag)* V_.PSConInit_NMDA;
    S_.y[State_::DG_GABAA] += B_.spikes_GABAA_.get_value(lag) * V_.PSConInit_GABAA;

    // set new input current
    V_.I_CURR = B_.currents_.get_value(lag);
   
    // log state data
    B_.logger_.record_data(origin.get_steps() + lag);

  }
}

void mynest::izhik_cond_alpha::handle(SpikeEvent & e)
{
  assert(e.get_delay() > 0);
  // Assert that port is 0 or 1 (SUP_SPIKE_RECEPTOR (3)- MIN_SPIKE_RECEPTOR (1)
  // = 2). AMPA =1, NMDA =2
  assert(0 <= e.get_rport() && e.get_rport() < SUP_SPIKE_RECEPTOR - MIN_SPIKE_RECEPTOR);

  if (e.get_rport() == AMPA - MIN_SPIKE_RECEPTOR)
  	B_.spikes_AMPA_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
  			e.get_weight() * e.get_multiplicity() );
  else if (e.get_rport() == NMDA - MIN_SPIKE_RECEPTOR)
    B_.spikes_NMDA_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
    		e.get_weight() * e.get_multiplicity() );
  else
  	B_.spikes_GABAA_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
  			e.get_weight() * e.get_multiplicity() );
}

void mynest::izhik_cond_alpha::handle(CurrentEvent& e)
{
  assert(e.get_delay() > 0);

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), 
			 e.get_weight() * e.get_current());

  assert(e.get_delay() > 0);
  // Assert that port is 0 (SUP_SPIKE_RECEPTOR (4) - MIN_SPIKE_RECEPTOR (3) = 1)
  assert(0 <= e.get_rport() && e.get_rport() < SUP_CURR_RECEPTOR - MIN_CURR_RECEPTOR);

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
					e.get_weight() * e.get_current());
}

void mynest::izhik_cond_alpha::handle(DataLoggingRequest& e)
{
  B_.logger_.handle(e);
}

#endif //HAVE_GSL
