/*
 *  izhik_cond_beta.cpp
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


#include "izhik_cond_beta.h"

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

nest::RecordablesMap<mynest::izhik_cond_beta> mynest::izhik_cond_beta::recordablesMap_;

namespace nest   // template specialization must be placed in namespace
{
  /*
   * Override the create() method with one call to RecordablesMap::insert_() 
   * for each quantity to be recorded.
   */
  template <>
  void RecordablesMap<mynest::izhik_cond_beta>::create()
  {
  	// use standard names whereever you can for consistency!
  	// Recording current seems
  	insert_(names::V_m,
  			&mynest::izhik_cond_beta::get_y_elem_<mynest::izhik_cond_beta::State_::V_M>);
  	insert_(Name("u"),
  			&mynest::izhik_cond_beta::get_y_elem_<mynest::izhik_cond_beta::State_::u>);
  	insert_(Name("g_AMPA"),
  			&mynest::izhik_cond_beta::get_y_elem_<mynest::izhik_cond_beta::State_::G_AMPA>);
  	insert_(Name("g_NMDA"),
  			&mynest::izhik_cond_beta::get_y_elem_<mynest::izhik_cond_beta::State_::G_NMDA>);
  	insert_(Name("g_GABAA_1"),
  			&mynest::izhik_cond_beta::get_y_elem_<mynest::izhik_cond_beta::State_::G_GABAA_1>);
  	insert_(Name("g_GABAA_2"),
  			&mynest::izhik_cond_beta::get_y_elem_<mynest::izhik_cond_beta::State_::G_GABAA_2>);

	insert_(Name("I"), &mynest::izhik_cond_beta::get_I_);
	insert_(Name("I_AMPA"), &mynest::izhik_cond_beta::get_I_AMPA_);
	insert_(Name("I_NMDA"), &mynest::izhik_cond_beta::get_I_NMDA_);
	insert_(Name("I_GABAA_1"), &mynest::izhik_cond_beta::get_I_GABAA_1_);
  insert_(Name("I_GABAA_2"), &mynest::izhik_cond_beta::get_I_GABAA_2_);
  insert_(Name("I_V_clamp"), &mynest::izhik_cond_beta::get_I_V_clamp_);

}
}

/* ---------------------------------------------------------------- 
 * Iteration function
 * ---------------------------------------------------------------- */

extern "C"
inline int mynest::izhik_cond_beta_dynamics(double, const double y[], double f[], void* pnode)
{ 
	// some shorthands
	//typedef mynest::izhik_cond_beta         N;
	typedef mynest::izhik_cond_beta::State_ S;

	// get access to node so we can almost work as in a member class
	assert(pnode);
	mynest::izhik_cond_beta& node =  *(reinterpret_cast<mynest::izhik_cond_beta*>(pnode));

	// easier access to membrane potential and recovery variable
	const nest::double_t& V = y[S::V_M];
	const nest::double_t& u = y[S::u];

	// y[] here is---and must be---the state vector supplied by the integrator,
	// not the state vector in the node, node.S_.y[].

	// The following code is verbose for the sake of clarity. We assume that a
	// good compiler will optimize the verbosity away.
	const nest::double_t I_AMPA = - y[S::G_AMPA] * ( V - node.P_.AMPA_E_rev );
	const nest::double_t I_NMDA = - y[S::G_NMDA] * ( V - node.P_.NMDA_E_rev )
      		/ ( 1 + std::exp( (node.P_.NMDA_Vact - V)/node.P_.NMDA_Sact ) );
  const nest::double_t I_GABAA_1 = - y[S::G_GABAA_1] * ( V - node.P_.GABAA_1_E_rev );
  const nest::double_t I_GABAA_2 = - y[S::G_GABAA_2] * ( V - node.P_.GABAA_2_E_rev );

	// Set state variable used for recording AMPA, NMDA and GABAA current
	// contributions
	node.S_.I_AMPA_    = I_AMPA;
	node.S_.I_NMDA_    = I_NMDA;
  node.S_.I_GABAA_1_ = I_GABAA_1;
  node.S_.I_GABAA_2_ = I_GABAA_2;

	// Total input current from synapses and external input
	node.S_.I_ =  I_AMPA + I_NMDA + I_GABAA_1 + I_GABAA_2 + node.B_.I_stim_ ;

	// If voltage clamp. I_ is feed back if in v clamp mode.
	if ( node.P_.V_clamp == 1 )
		node.S_.I_V_clamp_ = node.P_.k*( node.P_.V_clamp_at - node.P_.V_r )*( node.P_.V_clamp_at - node.P_.V_t ) - u  + node.S_.I_ + node.P_.I_e;
	else
		node.S_.I_V_clamp_ = 0;

  // Neuron dynamics
	// dV_m/dt
	f[ S::V_M ] = ( node.P_.k*( V - node.P_.V_r )*(  V - node.P_.V_t ) - u + node.S_.I_ + node.P_.I_e - node.S_.I_V_clamp_ )/node.P_.C_m ;

	// du/dt
	// If V is less than Vb then b=b1 and p=p1 else b=b2 and p=p2
	if ( V < node.P_.V_b )
		f[ S::u ] = node.P_.a*( node.P_.b_1*std::pow( V - node.P_.V_b, node.P_.p_1 ) - u + node.P_.u_const);
	else
		f[ S::u ] = node.P_.a*( node.P_.b_2*std::pow( V - node.P_.V_b, node.P_.p_2 ) - u + node.P_.u_const);


	// Synapse dynamics
	// dg_AMPA/dt
	f[ S::G_AMPA ] 		= -y[ S::G_AMPA ] / node.P_.AMPA_Tau_decay;
	
	// Synapse dynamics
	// d dg_AMPA/dt, dg_AMPA/dt
	f[ S::DG_AMPA ] = -y[S::DG_AMPA] / node.P_.AMPA_Tau_rise;
	f[ S::G_AMPA ]  = y[S::DG_AMPA] - y[S::G_AMPA]/node.P_.AMPA_Tau_decay;

	// d dg_NMDA/dt, dg_NMDA/dt
	f[ S::DG_NMDA ] = -y[S::DG_NMDA] / node.P_.NMDA_Tau_rise;
	f[ S::G_NMDA ]  =  y[S::DG_NMDA] - y[S::G_NMDA]/node.P_.NMDA_Tau_decay;

	// d dg_GABAA_1/dt, dg_GABAA_1/dt
	f[ S::DG_GABAA_1 ] = -y[S::DG_GABAA_1] / node.P_.GABAA_1_Tau_rise;
	f[ S::G_GABAA_1 ]  =  y[S::DG_GABAA_1] - y[S::G_GABAA_1]/node.P_.GABAA_1_Tau_decay;

// d dg_GABAA_1/dt, dg_GABAA_1/dt
	f[ S::DG_GABAA_2 ] = -y[S::DG_GABAA_2] / node.P_.GABAA_2_Tau_rise;
	f[ S::G_GABAA_2 ]  =  y[S::DG_GABAA_2] - y[S::G_GABAA_2]/node.P_.GABAA_2_Tau_decay;

	return GSL_SUCCESS;
}

/* ---------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

mynest::izhik_cond_beta::Parameters_::Parameters_()
: V_clamp   (   0    ),    		 //!< 0 or 1, no/yes for voltage clamp, i.e. when 1 then membrane voltage will be clamped at V_clamp_at
  V_clamp_at(  -70   ), 			 //!< Membrane potential in mV

	V_m			  ( -50.0  ), 		   //!< Membrane potential in mV
  V_r			  ( -50.0  ),        //!< Membrane resting potential
  V_t			  ( -40.0  ),			   //!< Instantaneous activation threshold
  C_m			  (  100.0 ),        //!< Capacity of the membrane in pF
  k				  (  0.7   ),				 //!< Factor determining instantaneous I-V relation
  c				  ( -50.0  ),				 //!< Reset value membrane potential
  kc				(  0.0   ),				 //!< Slope of voltage reset point
  V_peak	  (  30.0  ),        //!< spike cut of value
  I_e       (  0.0   ),        //!< pA

  u			    (  0.0   ), 			 //!< Recovery variable pA
  V_b       ( -50.0  ),			   //!< Recovery variable voltage threshold mV
  a				  (  0.03  ),			   //!< Time constant slow dynamics
  b_1			  ( -2.0   ),				 //!< Slope factor 1 slow dynamics
  b_2			  ( -2.0   ),				 //!< Slope factor 2 slow dynamics
  p_1			  (  1.0   ),				 //!< Polynomial voltage dependency factor 1
  p_2			  (  1.0   ),				 //!< Polynomial voltage dependency factor 2
  d_1			  (  100.0 ),        //!< Slow variable change at spike when u > du
  d_2			  (  100.0 ),        //!< slow variable change at spike when u < du
  d_u       (  0.0   ),        //!< if u <du then d = d1 else d=d2
  d_1_bound (  0.0   ),        //!< If 1 then d=d1 if abs( u - du ) <d1 else  d=-(u -du)
  u_const   (  0.0   ),				 //!< Constant current recovery variable

  AMPA_E_rev         (  0.0    ),  // mV
  AMPA_g_peak        (  1.0    ),  // mV
  AMPA_Tau_rise      (  1.0    ),  // ms
  AMPA_Tau_decay     (  3.0    ),  // ms
    
  NMDA_E_rev         (  0.0    ),  // mV
  NMDA_g_peak        (  1.0    ),  // mV
  NMDA_Tau_rise      (  10.0   ),  // ms
  NMDA_Tau_decay     (  100.0  ),  // ms
  NMDA_Vact          ( -58.0   ), // mV
  NMDA_Sact          (  2.5    ), // mV

  GABAA_1_E_rev        (-70.0    ),  // mV
  GABAA_1_g_peak       (  1.0    ),  // mV
  GABAA_1_Tau_rise     (  1.0    ),  // ms
  GABAA_1_Tau_decay    (  4.0    ),  // ms
  
  GABAA_2_E_rev        (-70.0    ),  // mV
  GABAA_2_g_peak       (  1.0    ),  // mV
  GABAA_2_Tau_rise     (  1.0    ),  // ms
  GABAA_2_Tau_decay    (  4.0    )  // ms

{
	recordablesMap_.create();
}

mynest::izhik_cond_beta::State_::State_(const Parameters_& p)
: I_(0.0),
  I_AMPA_(0.0),
  I_NMDA_(0.0),
  I_GABAA_1_(0.0),
	I_GABAA_2_(0.0),
	I_V_clamp_(0.0)
{
	y[V_M] = p.V_r;  // initialize to resting potential
	y[u] = p.u;   // initialize to u
	for ( size_t i = 2 ; i < STATE_VEC_SIZE ; ++i )
		y[i] = 0;
}

mynest::izhik_cond_beta::State_::State_(const State_& s)
: I_(         s.I_         ),
  I_AMPA_(    s.I_AMPA_    ),
  I_NMDA_(    s.I_NMDA_    ),
  I_GABAA_1_( s.I_GABAA_1_ ),
  I_GABAA_2_( s.I_GABAA_2_ ),
  I_V_clamp_( s.I_V_clamp_ )

{
	for ( size_t i = 0 ; i < STATE_VEC_SIZE ; ++i )
		y[i] = s.y[i];
}

mynest::izhik_cond_beta::State_& mynest::izhik_cond_beta::State_::operator=(const State_& s)
{
	if ( this == &s )  // avoid assignment to self
		return *this;

	for ( size_t i = 0 ; i < STATE_VEC_SIZE ; ++i )
		y[i] = s.y[i];

	I_         = s.I_;
	I_AMPA_    = s.I_AMPA_;
	I_NMDA_    = s.I_NMDA_;
  I_GABAA_1_ = s.I_GABAA_1_;
  I_GABAA_2_ = s.I_GABAA_2_;
  I_V_clamp_ = s.I_V_clamp_;

  //r = s.r;
	return *this;
}

mynest::izhik_cond_beta::State_::~State_()
{
}

mynest::izhik_cond_beta::Buffers_::Buffers_(izhik_cond_beta& n)
: logger_(n),
  s_(0),
  c_(0),
  e_(0)
{
	// The other member variables are left uninitialised or are
	// automatically initialised by their default constructor.
}

mynest::izhik_cond_beta::Buffers_::Buffers_(const Buffers_&, izhik_cond_beta& n)
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

void mynest::izhik_cond_beta::Parameters_::get(DictionaryDatum &dd) const
{
	def<nest::double_t>(dd, "V_clamp",      V_clamp);
	def<nest::double_t>(dd, "V_clamp_at",   V_clamp_at);

	def<double>(dd,names::V_m,          V_m);
	def<nest::double_t>(dd, "V_r",      V_r);
	def<nest::double_t>(dd, "V_t",      V_t);
	def<double>(dd,names::C_m,          C_m);
	def<nest::double_t>(dd, "k",        k);
	def<nest::double_t>(dd, "c",        c);
	def<nest::double_t>(dd, "kc",       kc);
	def<nest::double_t>(dd, "V_peak",   V_peak);
	def<double>(dd,names::I_e,          I_e);

	def<nest::double_t>(dd, "u",      	 u);
	def<nest::double_t>(dd, "V_b",       V_b);
	def<nest::double_t>(dd, "a",         a);
	def<nest::double_t>(dd, "b_1",       b_1);
	def<nest::double_t>(dd, "b_2",       b_2);
	def<nest::double_t>(dd, "p_1",       p_1);
	def<nest::double_t>(dd, "p_2",       p_2);
	def<nest::double_t>(dd, "d_1",       d_1);	      //!< Slow variable change at spike when u > du
	def<nest::double_t>(dd, "d_2",       d_2);        //!< slow variable change at spike when u < du
	def<nest::double_t>(dd, "d_u",       d_u);        //!< if u <du then d = d1 else d=d2
	def<nest::double_t>(dd, "d_1_bound", d_1_bound);  //!< If 1 then d=d1 if abs( u - du ) <d1 elseif 0 d=-(u -du)
	def<nest::double_t>(dd, "u_const",   u_const);

	def<nest::double_t>(dd, "AMPA_E_rev",      AMPA_E_rev);
	def<nest::double_t>(dd, "AMPA_g_peak",     AMPA_g_peak);
	def<nest::double_t>(dd, "AMPA_Tau_rise",   AMPA_Tau_rise);
	def<nest::double_t>(dd, "AMPA_Tau_decay",  AMPA_Tau_decay);
	
	def<nest::double_t>(dd, "NMDA_E_rev",      NMDA_E_rev);
	def<nest::double_t>(dd, "NMDA_g_peak",     NMDA_g_peak);
	def<nest::double_t>(dd, "NMDA_Tau_rise",   NMDA_Tau_rise);
	def<nest::double_t>(dd, "NMDA_Tau_decay",  NMDA_Tau_decay);
	def<nest::double_t>(dd, "NMDA_Vact",       NMDA_Vact);
	def<nest::double_t>(dd, "NMDA_Sact",       NMDA_Sact);

  def<nest::double_t>(dd, "GABAA_1_E_rev",     GABAA_1_E_rev);
	def<nest::double_t>(dd, "GABAA_1_g_peak",    GABAA_1_g_peak);
	def<nest::double_t>(dd, "GABAA_1_Tau_rise",  GABAA_1_Tau_rise);
	def<nest::double_t>(dd, "GABAA_1_Tau_decay", GABAA_1_Tau_decay);

  def<nest::double_t>(dd, "GABAA_2_E_rev",     GABAA_2_E_rev);
	def<nest::double_t>(dd, "GABAA_2_g_peak",    GABAA_2_g_peak);
	def<nest::double_t>(dd, "GABAA_2_Tau_rise",  GABAA_2_Tau_rise);
	def<nest::double_t>(dd, "GABAA_2_Tau_decay", GABAA_2_Tau_decay);
}

void mynest::izhik_cond_beta::Parameters_::set(const DictionaryDatum& dd)
{
	updateValue<nest::double_t>(dd, "V_clamp",      V_clamp);
	updateValue<nest::double_t>(dd, "V_clamp_at",   V_clamp_at);

	// allow setting the membrane potential
	updateValue<double>(dd,names::V_m,          V_m);
	updateValue<nest::double_t>(dd, "V_r",      V_r);
	updateValue<nest::double_t>(dd, "V_t",      V_t);
	updateValue<double>(dd,names::C_m,          C_m);
	updateValue<nest::double_t>(dd, "c",        c);
	updateValue<nest::double_t>(dd, "k",        k);
	updateValue<nest::double_t>(dd, "kc",       kc);
	updateValue<nest::double_t>(dd, "V_peak",   V_peak);
	updateValue<double>(dd,names::I_e,          I_e);

	updateValue<nest::double_t>(dd, "u",         u);
	updateValue<nest::double_t>(dd, "V_b",       V_b);
	updateValue<nest::double_t>(dd, "a",         a);
	updateValue<nest::double_t>(dd, "b_1",       b_1);
	updateValue<nest::double_t>(dd, "b_2",       b_2);
	updateValue<nest::double_t>(dd, "p_1",       p_1);
	updateValue<nest::double_t>(dd, "p_2",       p_2);
	updateValue<nest::double_t>(dd, "d_1",       d_1);	      //!< Slow variable change at spike when u > du
	updateValue<nest::double_t>(dd, "d_2",       d_2);        //!< slow variable change at spike when u < du
	updateValue<nest::double_t>(dd, "d_u",       d_u);        //!< if u <du then d = d1 else d=d2
	updateValue<nest::double_t>(dd, "d_1_bound", d_1_bound);  //!< If 1 then d=d1 if abs( u - du ) <d1 elseif 0 d=-(u -du)
	updateValue<nest::double_t>(dd, "u_const",   u_const);

	updateValue<nest::double_t>(dd, "AMPA_E_rev",      AMPA_E_rev);
	updateValue<nest::double_t>(dd, "AMPA_g_peak",     AMPA_g_peak);
	updateValue<nest::double_t>(dd, "AMPA_Tau_rise",   AMPA_Tau_rise);
	updateValue<nest::double_t>(dd, "AMPA_Tau_decay",  AMPA_Tau_decay);

	updateValue<nest::double_t>(dd, "NMDA_E_rev",      NMDA_E_rev);
	updateValue<nest::double_t>(dd, "NMDA_g_peak",     NMDA_g_peak);
	updateValue<nest::double_t>(dd, "NMDA_Tau_rise",   NMDA_Tau_rise);
	updateValue<nest::double_t>(dd, "NMDA_Tau_decay",  NMDA_Tau_decay);
	updateValue<nest::double_t>(dd, "NMDA_Vact",       NMDA_Vact);
	updateValue<nest::double_t>(dd, "NMDA_Sact",       NMDA_Sact);

	updateValue<nest::double_t>(dd, "GABAA_1_E_rev",     GABAA_1_E_rev);
	updateValue<nest::double_t>(dd, "GABAA_1_g_peak",    GABAA_1_g_peak);
	updateValue<nest::double_t>(dd, "GABAA_1_Tau_rise",  GABAA_1_Tau_rise);
	updateValue<nest::double_t>(dd, "GABAA_1_Tau_decay", GABAA_1_Tau_decay);

	updateValue<nest::double_t>(dd, "GABAA_2_E_rev",     GABAA_2_E_rev);
	updateValue<nest::double_t>(dd, "GABAA_2_g_peak",    GABAA_2_g_peak);
	updateValue<nest::double_t>(dd, "GABAA_2_Tau_rise",  GABAA_2_Tau_rise);
	updateValue<nest::double_t>(dd, "GABAA_2_Tau_decay", GABAA_2_Tau_decay);


	if ( V_r > V_t )
		throw BadProperty("Reset potential must be smaller than instantaneous threshold.");

	if ( C_m <= 0 )
		throw BadProperty("Capacitance must be strictly positive.");

	//if ( t_ref < 0 )
	//  throw BadProperty("Refractory time cannot be negative.");

  if ( AMPA_Tau_decay    <= 0 ||
  		 NMDA_Tau_decay    <= 0 ||
  		 GABAA_1_Tau_decay <= 0 ||
  		 GABAA_2_Tau_decay <= 0 )
		throw BadProperty("All time constants must be strictly positive.");


}

void mynest::izhik_cond_beta::State_::get(DictionaryDatum &d) const
{
	def<double>(d, names::V_m, y[V_M]); // Membrane potential
	def<nest::double_t>(d, "u", y[u]); // Recovery variable
}

void mynest::izhik_cond_beta::State_::set(const DictionaryDatum& d, const Parameters_&)
{
	updateValue<double>(d, names::V_m, y[V_M]);
	updateValue<nest::double_t>(d, "u", y[u]); // Recovery variable
}


/* ---------------------------------------------------------------- 
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

mynest::izhik_cond_beta::izhik_cond_beta()
: Archiving_Node(),
  P_(),
  S_(P_),
  B_(*this)
{
}

mynest::izhik_cond_beta::izhik_cond_beta(const izhik_cond_beta& n)
: Archiving_Node(n),
  P_(n.P_),
  S_(n.S_),
  B_(n.B_, *this)
{
}

mynest::izhik_cond_beta::~izhik_cond_beta()
{
	// GSL structs only allocated by init_nodes_(), so we need to protect destruction
	if ( B_.s_ ) gsl_odeiv_step_free(B_.s_);
	if ( B_.c_ ) gsl_odeiv_control_free(B_.c_);
	if ( B_.e_ ) gsl_odeiv_evolve_free(B_.e_);
}

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void mynest::izhik_cond_beta::init_node_(const Node& proto)
{
	const izhik_cond_beta& pr = downcast<izhik_cond_beta>(proto);
	P_ = pr.P_;
	S_ = pr.S_;
}

void mynest::izhik_cond_beta::init_state_(const Node& proto)
{
	const izhik_cond_beta& pr = downcast<izhik_cond_beta>(proto);
	S_ = pr.S_;
}

void mynest::izhik_cond_beta::init_buffers_()
{

	B_.spikes_AMPA_.clear();       // includes resize
	B_.spikes_NMDA_.clear();       // includes resize
  B_.spikes_GABAA_1_.clear();    // includes resize
  B_.spikes_GABAA_2_.clear();    // includes resize
	B_.currents_.clear();          // includes resize

	B_.logger_.reset();


	nest::Archiving_Node::clear_history();

	B_.step_ = Time::get_resolution().get_ms();

	// We must integrate this model with high-precision to obtain decent results
	B_.IntegrationStep_ = std::min(0.01, B_.step_);
	//B_.IntegrationStep_ = B_.step_;

	static const gsl_odeiv_step_type* T1 = gsl_odeiv_step_rkf45;

	if ( B_.s_ == 0 )
		B_.s_ = gsl_odeiv_step_alloc (T1, State_::STATE_VEC_SIZE);
	else
		gsl_odeiv_step_reset(B_.s_);

	// Lower tolerance
	if ( B_.c_ == 0 )
		B_.c_ = gsl_odeiv_control_y_new (1e-6, 1e-6); // Changed from (1e-3, 0)
	else
		gsl_odeiv_control_init(B_.c_, 1e-6, 1e-6, 0.0, 1.0); // Changed from ( 1e-3, 0.0, 1.0, 0.0)

	/*
	if ( B_.c_ == 0 )
			B_.c_ = gsl_odeiv_control_y_new (1e-3, 0.0);
	else
			gsl_odeiv_control_init(B_.c_,  1e-3, 0.0, 1.0, 0.0);
	*/

	if ( B_.e_ == 0 )
		B_.e_ = gsl_odeiv_evolve_alloc(State_::STATE_VEC_SIZE);
	else
		gsl_odeiv_evolve_reset(B_.e_);

	B_.sys_.function  = izhik_cond_beta_dynamics;
	B_.sys_.jacobian  = NULL;
	B_.sys_.dimension = State_::STATE_VEC_SIZE;
	B_.sys_.params    = reinterpret_cast<void*>(this);

	B_.I_stim_ = 0.0;
}

nest::double_t mynest::izhik_cond_beta::get_synapse_constant(nest::double_t Tau_rise,
		nest::double_t Tau_decay,
		nest::double_t g_peak)
{
	// Factor used to account for the missing 1/((1/Tau_2)-(1/Tau_1)) term
	// in the ht_neuron_dynamics integration of the synapse terms.
	// See: Exact digital simulation of time-invariant linear systems
	// with applications to neuronal modeling, Rotter and Diesmann,
	// section 3.1.2.
	nest::double_t exact_integration_adjustment = (1/Tau_decay)-(1/Tau_rise);

	nest::double_t t_peak = (Tau_decay*Tau_rise)*std::log(Tau_decay/Tau_rise)/(Tau_decay-Tau_rise);
	nest::double_t normalisation_factor = 1/(std::exp(-t_peak/Tau_rise) - std::exp(-t_peak/Tau_decay));

	return g_peak*normalisation_factor*exact_integration_adjustment;
}


void mynest::izhik_cond_beta::calibrate()
{
	B_.logger_.init();  // ensures initialization in case mm connected after Simulate

// NOTE: code below initializes conductance step size for incoming pulses.
	// Variable and function names need to be changed!
	V_.cond_steps_.resize(SUP_SPIKE_RECEPTOR-1);

	V_.cond_steps_[AMPA-1] =
			get_synapse_constant(P_.AMPA_Tau_rise, P_.AMPA_Tau_decay, P_.AMPA_g_peak);

	V_.cond_steps_[NMDA-1] =
			get_synapse_constant(P_.NMDA_Tau_rise, P_.NMDA_Tau_decay, P_.NMDA_g_peak);

	V_.cond_steps_[GABAA_1-1] =
			get_synapse_constant(P_.GABAA_1_Tau_rise, P_.GABAA_1_Tau_decay, P_.GABAA_1_g_peak);

	V_.cond_steps_[GABAA_2-1] =
			get_synapse_constant(P_.GABAA_2_Tau_rise, P_.GABAA_2_Tau_decay, P_.GABAA_2_g_peak);
}

/* ---------------------------------------------------------------- 
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void mynest::izhik_cond_beta::update(Time const & origin, const nest::long_t from, const nest::long_t to)
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



			// spikes are handled inside the while-loop
			// due to spike-driven adaptation
			if ( S_.y[State_::V_M] >= P_.V_peak )
			{
				// v is reset to c + u*kc
				S_.y[State_::V_M] = P_.c + S_.y[State_::u]*P_.kc;


				// Set d, if u < d_u then d=d1 else d=d2
				if ( S_.y[State_::u] < P_.d_u )
				{
          // If d1 < -(u-du) d->d1 else d=-(u-du)
					if ( P_.d_1_bound == 1 )
							S_.y[State_::u] += ( P_.d_1 < -( S_.y[State_::u] - P_.d_u ) ) ? P_.d_1 : -( S_.y[State_::u] - P_.d_u ) ;
					else
						S_.y[State_::u] += P_.d_1;
				}
				else
					S_.y[State_::u] += P_.d_2;


				// log spike with Archiving_Node
				set_spiketime(Time::step(origin.get_steps()+lag+1));
				SpikeEvent se;
				network()->send(*this, se, lag);
			}
		}

		// Here incomming spikes are added. It is setup such that AMPA, NMDA
		// and GABA receives from one receptor each.
		S_.y[State_::DG_AMPA]  += B_.spikes_AMPA_.get_value(lag)* V_.cond_steps_[AMPA-1];
		S_.y[State_::DG_NMDA]  += B_.spikes_NMDA_.get_value(lag)* V_.cond_steps_[NMDA-1];
		S_.y[State_::DG_GABAA_1] += B_.spikes_GABAA_1_.get_value(lag) * V_.cond_steps_[GABAA_1-1];
		S_.y[State_::DG_GABAA_2] += B_.spikes_GABAA_2_.get_value(lag) * V_.cond_steps_[GABAA_2-1];

		// set new input current
    B_.I_stim_ = B_.currents_.get_value(lag);

		// log state data
		B_.logger_.record_data(origin.get_steps() + lag);

	}
}

void mynest::izhik_cond_beta::handle(SpikeEvent & e)
{
	assert(e.get_delay() > 0);
	// Assert that port is 0 or 1 (SUP_SPIKE_RECEPTOR (3)- MIN_SPIKE_RECEPTOR (1)
	// = 2). AMPA =1, NMDA =2
	assert(0 <= e.get_rport() && e.get_rport() < SUP_SPIKE_RECEPTOR - MIN_SPIKE_RECEPTOR);

  // If AMPA
	if (e.get_rport() == AMPA - MIN_SPIKE_RECEPTOR)
		B_.spikes_AMPA_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
				e.get_weight() * e.get_multiplicity() );

  // If NMDA
	else if (e.get_rport() == NMDA - MIN_SPIKE_RECEPTOR)
		B_.spikes_NMDA_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
				e.get_weight() * e.get_multiplicity() );

	// If GABAA_1
  else if (e.get_rport() == GABAA_1 - MIN_SPIKE_RECEPTOR)
  	B_.spikes_GABAA_1_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
  			e.get_weight() * e.get_multiplicity() );

  // If GABAA_2
  else
  	B_.spikes_GABAA_2_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
  			e.get_weight() * e.get_multiplicity() );

}

void mynest::izhik_cond_beta::handle(CurrentEvent& e)
{
	assert(e.get_delay() > 0);

	// add weighted current; HEP 2002-10-04
	B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
			e.get_weight() * e.get_current());

	assert(e.get_delay() > 0);
	// Assert that port is 0 (SUP_SPIKE_RECEPTOR (4) - MIN_SPIKE_RECEPTOR (3) = 1)
	assert(0 <= e.get_rport() && e.get_rport() < SUP_CURR_RECEPTOR - MIN_CURR_RECEPTOR);

}

void mynest::izhik_cond_beta::handle(DataLoggingRequest& e)
{
	B_.logger_.handle(e);
}

#endif //HAVE_GSL
