/*
 *  my_iaf_cond_exp.cpp
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


#include "my_iaf_cond_exp.h"

#ifdef HAVE_GSL

#include "exceptions.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
#include "universal_data_logger_impl.h"
#include <limits>

#include <iomanip>
#include <iostream>
#include <cstdio>

using namespace nest; //added

/* ---------------------------------------------------------------- 
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap<mynest::my_iaf_cond_exp> mynest::my_iaf_cond_exp::recordablesMap_;

namespace nest   // template specialization must be placed in namespace
{
/*
 * Override the create() method with one call to RecordablesMap::insert_()
 * for each quantity to be recorded.
 */
template <>
void RecordablesMap<mynest::my_iaf_cond_exp>::create()
{
	// use standard names whereever you can for consistency!
	// Recording current seems
	insert_(names::V_m,
			&mynest::my_iaf_cond_exp::get_y_elem_<mynest::my_iaf_cond_exp::State_::V_M>);
	insert_(Name("g_AMPA"),
			&mynest::my_iaf_cond_exp::get_y_elem_<mynest::my_iaf_cond_exp::State_::G_AMPA>);
	insert_(Name("g_NMDA"),
			&mynest::my_iaf_cond_exp::get_y_elem_<mynest::my_iaf_cond_exp::State_::G_NMDA>);
	insert_(Name("g_GABAA_1"),
			&mynest::my_iaf_cond_exp::get_y_elem_<mynest::my_iaf_cond_exp::State_::G_GABAA_1>);
	insert_(Name("g_GABAA_2"),
			&mynest::my_iaf_cond_exp::get_y_elem_<mynest::my_iaf_cond_exp::State_::G_GABAA_2>);

	insert_(Name("I"        ), &mynest::my_iaf_cond_exp::get_I_);
	insert_(Name("I_AMPA"   ), &mynest::my_iaf_cond_exp::get_I_AMPA_);
	insert_(Name("I_NMDA"   ), &mynest::my_iaf_cond_exp::get_I_NMDA_);
	insert_(Name("I_GABAA_1"), &mynest::my_iaf_cond_exp::get_I_GABAA_1_);
	insert_(Name("I_GABAA_2"), &mynest::my_iaf_cond_exp::get_I_GABAA_2_);
	insert_(Name("I_V_clamp"), &mynest::my_iaf_cond_exp::get_I_V_clamp_);

	insert_(names::t_ref_remaining,
			&mynest::my_iaf_cond_exp::get_r_);
}
}

/* ---------------------------------------------------------------- 
 * Iteration function
 * ---------------------------------------------------------------- */

extern "C"
inline int mynest::my_iaf_cond_exp_dynamics(double, const double y[], double f[], void* pnode)
{ 
	// some shorthands
	typedef mynest::my_iaf_cond_exp         N;
	typedef mynest::my_iaf_cond_exp::State_ S;

	// get access to node so we can almost work as in a member class
	assert(pnode);
	mynest::my_iaf_cond_exp& node =  *(reinterpret_cast<mynest::my_iaf_cond_exp*>(pnode));

	// easier access to membrane potential
	const nest::double_t& V = y[S::V_M];

	// y[] here is---and must be---the state vector supplied by the integrator,
	// not the state vector in the node, node.S_.y[].

	// The following code is verbose for the sake of clarity. We assume that a
	// good compiler will optimize the verbosity away.
	const nest::double_t I_AMPA = - y[S::G_AMPA] * ( V - node.P_.AMPA_E_rev );
	const nest::double_t I_NMDA = - y[S::G_NMDA] * ( V - node.P_.NMDA_E_rev )
      						/ ( 1 + std::exp( (node.P_.NMDA_Vact - V)/node.P_.NMDA_Sact ) );
	const nest::double_t I_GABAA_1 = - y[S::G_GABAA_1] * ( V - node.P_.GABAA_1_E_rev );
	const nest::double_t I_GABAA_2 = - y[S::G_GABAA_2] * ( V - node.P_.GABAA_2_E_rev );

	const nest::double_t I_leak    = - node.P_.g_L * ( V - node.P_.E_L  );

	// Set state variable used for recording AMPA, NMDA and GABAA current
	// contributions
	node.S_.I_AMPA_    = I_AMPA;
	node.S_.I_NMDA_    = I_NMDA;
	node.S_.I_GABAA_1_ = I_GABAA_1;
	node.S_.I_GABAA_2_ = I_GABAA_2;

	node.S_.I_ =  I_AMPA + I_NMDA + I_GABAA_1 + I_GABAA_2 + node.B_.I_stim_ ;


	// dV_m/dt
	f[0]= ( I_leak + node.S_.I_ + node.P_.I_e ) / node.P_.C_m;

	// Synapse dynamics
	// dg_AMPA/dt
	f[ S::G_AMPA ] 		= -y[ S::G_AMPA ] / node.P_.AMPA_Tau_decay;

	// dg_NMDA/dt
	f[ S::G_NMDA ] 		= -y[ S::G_NMDA ] / node.P_.NMDA_Tau_decay;

	// dg_GABAA_1/dt
	f[ S::G_GABAA_1 ] = -y[ S::G_GABAA_1 ] / node.P_.GABAA_1_Tau_decay;

	// dg_GABAA_2/dt
	f[ S::G_GABAA_2 ] = -y[ S::G_GABAA_2 ] / node.P_.GABAA_2_Tau_decay;




	return GSL_SUCCESS;
}

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

mynest::my_iaf_cond_exp::Parameters_::Parameters_()
: V_th           (-55.0    ),  // mV
  V_reset        (-60.0    ),  // mV
  t_ref          (  2.0    ),  // ms
  g_L            ( 16.6667 ),  // nS
  C_m            (250.0    ),  // pF
  E_L            (-75.0    ),  // mV
  I_e            (  0.0    ),   // pA
 
  AMPA_E_rev     			(  0.0   ),  	// mV
  AMPA_Tau_decay 			(  3.0   ),  	// ms

  NMDA_E_rev     			(  0.0   ), 	// mV
  NMDA_Tau_decay 			(  100.0 ), 	// ms
  NMDA_Vact      			( -58.0  ),  	// mV
  NMDA_Sact           (  2.5   ),  	// mV

  GABAA_1_E_rev       (-70.0    ),  // mV
  GABAA_1_Tau_decay   (  4.0    ),  // ms

  GABAA_2_E_rev       (-70.0    ),  // mV
  GABAA_2_Tau_decay   (  4.0    )  	// ms


{
	recordablesMap_.create();
}

mynest::my_iaf_cond_exp::State_::State_(const Parameters_& p)
:I_AMPA_(0.0),
 I_NMDA_(0.0),
 I_GABAA_1_(0.0),
 I_GABAA_2_(0.0),
 I_V_clamp_(0.0),
 r(0)
{
	y[V_M] = p.E_L;  // initialize to reversal potential
	for ( size_t i = 1 ; i < STATE_VEC_SIZE ; ++i )
		y[i] = 0;
}

mynest::my_iaf_cond_exp::State_::State_(const State_& s)
:I_(s.I_),
 I_AMPA_(  s.I_AMPA_  ),
 I_NMDA_(  s.I_NMDA_  ),
 I_GABAA_1_(s.I_GABAA_1_),
 I_GABAA_2_(s.I_GABAA_2_),
 I_V_clamp_(s.I_V_clamp_),
 r(s.r)
{
	for ( size_t i = 0 ; i < STATE_VEC_SIZE ; ++i )
		y[i] = s.y[i];
}

mynest::my_iaf_cond_exp::State_& mynest::my_iaf_cond_exp::State_::operator=(const State_& s)
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
	r = s.r;
	return *this;
}

mynest::my_iaf_cond_exp::State_::~State_()
{
}

mynest::my_iaf_cond_exp::Buffers_::Buffers_(my_iaf_cond_exp& n)
: logger_(n),
  s_(0),
  c_(0),
  e_(0)
{
	// The other member variables are left uninitialised or are
	// automatically initialised by their default constructor.
}

mynest::my_iaf_cond_exp::Buffers_::Buffers_(const Buffers_&, my_iaf_cond_exp& n)
: logger_(n),
  s_(0),
  c_(0),
  e_(0)
{
	// The other member variables are left uninitialised or are
	// automatically initialised by their default constructor.
}

/* ---------------------------------------------------------------- 
 * Paramater and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void mynest::my_iaf_cond_exp::Parameters_::get(DictionaryDatum &dd) const
{
	def<double>(dd,names::V_th,         V_th);
	def<double>(dd,names::V_reset,      V_reset);
	def<double>(dd,names::t_ref,        t_ref);
	def<double>(dd,names::g_L,          g_L);
	def<double>(dd,names::E_L,          E_L);
	def<double>(dd,names::C_m,          C_m);
	def<double>(dd,names::I_e,          I_e);

	def<nest::double_t>(dd, "AMPA_E_rev",         AMPA_E_rev);
	def<nest::double_t>(dd, "AMPA_Tau_decay",     AMPA_Tau_decay);

	def<nest::double_t>(dd, "NMDA_E_rev",         NMDA_E_rev);
	def<nest::double_t>(dd, "NMDA_Tau_decay",     NMDA_Tau_decay);
	def<nest::double_t>(dd, "NMDA_Vact",          NMDA_Vact);
	def<nest::double_t>(dd, "NMDA_Sact",       		NMDA_Sact);

	def<nest::double_t>(dd, "GABAA_1_E_rev",     	GABAA_1_E_rev);
	def<nest::double_t>(dd, "GABAA_1_Tau_decay", 	GABAA_1_Tau_decay);

	def<nest::double_t>(dd, "GABAA_2_E_rev",     	GABAA_2_E_rev);
	def<nest::double_t>(dd, "GABAA_2_Tau_decay", 	GABAA_2_Tau_decay);


}

void mynest::my_iaf_cond_exp::Parameters_::set(const DictionaryDatum& dd)
{
	// allow setting the membrane potential
	updateValue<double>(dd,names::V_th,    V_th);
	updateValue<double>(dd,names::V_reset, V_reset);
	updateValue<double>(dd,names::t_ref,   t_ref);
	updateValue<double>(dd,names::g_L,     g_L);
	updateValue<double>(dd,names::E_L,     E_L);
	updateValue<double>(dd,names::C_m,     C_m);
	updateValue<double>(dd,names::I_e,     I_e);


	updateValue<nest::double_t>(dd, "AMPA_E_rev",        AMPA_E_rev);
	updateValue<nest::double_t>(dd, "AMPA_Tau_decay",    AMPA_Tau_decay);

	updateValue<nest::double_t>(dd, "NMDA_E_rev",        NMDA_E_rev);
	updateValue<nest::double_t>(dd, "NMDA_Tau_decay",    NMDA_Tau_decay);
	updateValue<nest::double_t>(dd, "NMDA_Vact",         NMDA_Vact);
	updateValue<nest::double_t>(dd, "NMDA_Sact",         NMDA_Sact);

	updateValue<nest::double_t>(dd, "GABAA_1_E_rev",     GABAA_1_E_rev);
	updateValue<nest::double_t>(dd, "GABAA_1_Tau_decay", GABAA_1_Tau_decay);

	updateValue<nest::double_t>(dd, "GABAA_2_E_rev",     GABAA_2_E_rev);
	updateValue<nest::double_t>(dd, "GABAA_2_Tau_decay", GABAA_2_Tau_decay);



	if ( V_reset >= V_th )
		throw BadProperty("Reset potential must be smaller than threshold.");

	if ( C_m <= 0 )
		throw BadProperty("Capacitance must be strictly positive.");

	if ( t_ref < 0 )
		throw BadProperty("Refractory time cannot be negative.");

	if ( AMPA_Tau_decay    <= 0 ||
			NMDA_Tau_decay    <= 0 ||
			GABAA_1_Tau_decay <= 0 ||
			GABAA_2_Tau_decay <= 0 )
		throw BadProperty("All time constants must be strictly positive.");

}

void mynest::my_iaf_cond_exp::State_::get(DictionaryDatum &d) const
{
	def<double>(d, names::V_m, y[V_M]); // Membrane potential
}

void mynest::my_iaf_cond_exp::State_::set(const DictionaryDatum& d, const Parameters_&)
{
	updateValue<double>(d, names::V_m, y[V_M]);
}


/* ---------------------------------------------------------------- 
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

mynest::my_iaf_cond_exp::my_iaf_cond_exp()
: Archiving_Node(),
  P_(),
  S_(P_),
  B_(*this)
{
}

mynest::my_iaf_cond_exp::my_iaf_cond_exp(const my_iaf_cond_exp& n)
: Archiving_Node(n),
  P_(n.P_),
  S_(n.S_),
  B_(n.B_, *this)
{
}

mynest::my_iaf_cond_exp::~my_iaf_cond_exp()
{
	// GSL structs only allocated by init_nodes_(), so we need to protect destruction
	if ( B_.s_ ) gsl_odeiv_step_free(B_.s_);
	if ( B_.c_ ) gsl_odeiv_control_free(B_.c_);
	if ( B_.e_ ) gsl_odeiv_evolve_free(B_.e_);
}

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void mynest::my_iaf_cond_exp::init_node_(const Node& proto)
{
	const my_iaf_cond_exp& pr = downcast<my_iaf_cond_exp>(proto);
	P_ = pr.P_;
	S_ = pr.S_;
}

void mynest::my_iaf_cond_exp::init_state_(const Node& proto)
{
	const my_iaf_cond_exp& pr = downcast<my_iaf_cond_exp>(proto);
	S_ = pr.S_;
}

void mynest::my_iaf_cond_exp::init_buffers_()
{
	Archiving_Node::clear_history();

	B_.spikes_AMPA_.clear();       // includes resize
	B_.spikes_NMDA_.clear();       // includes resize
	B_.spikes_GABAA_1_.clear();    // includes resize
	B_.spikes_GABAA_2_.clear();    // includes resize
	B_.currents_.clear();          // includes resize

	B_.logger_.reset();

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

	B_.sys_.function  = my_iaf_cond_exp_dynamics;
	B_.sys_.jacobian  = NULL;
	B_.sys_.dimension = State_::STATE_VEC_SIZE;
	B_.sys_.params    = reinterpret_cast<void*>(this);
}

void mynest::my_iaf_cond_exp::calibrate()
{
	B_.logger_.init();  // ensures initialization in case mm connected after Simulate

 V_.RefractoryCounts = Time(Time::ms(P_.t_ref)).get_steps();

	assert(V_.RefractoryCounts >= 0);  // since t_ref >= 0, this can only fail in error
}

/* ---------------------------------------------------------------- 
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void mynest::my_iaf_cond_exp::update(Time const & origin, const nest::long_t from, const nest::long_t to)
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
		if ( S_.r )
		{// neuron is absolute refractory
			--S_.r;
			S_.y[State_::V_M] = P_.V_reset;  // clamp potential
		}
		else
			// neuron is not absolute refractory
			if ( S_.y[State_::V_M] >= P_.V_th )
			{
				S_.r              = V_.RefractoryCounts;
				S_.y[State_::V_M] = P_.V_reset;

				// log spike with Archiving_Node
				set_spiketime(Time::step(origin.get_steps()+lag+1));

				SpikeEvent se;
				network()->send(*this, se, lag);
			}

		// Here incomming spikes are added. It is setup such that AMPA, NMDA
		// and GABA recieves from one receptor each.
		S_.y[State_::G_AMPA]    += B_.spikes_AMPA_.get_value(lag);
		S_.y[State_::G_NMDA]    += B_.spikes_NMDA_.get_value(lag);
		S_.y[State_::G_GABAA_1] += B_.spikes_GABAA_1_.get_value(lag);
		S_.y[State_::G_GABAA_2] += B_.spikes_GABAA_2_.get_value(lag);

		// set new input current
		B_.I_stim_ = B_.currents_.get_value(lag);

		// log state data
		B_.logger_.record_data(origin.get_steps() + lag);

	}
}

void mynest::my_iaf_cond_exp::handle(SpikeEvent & e)
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

void mynest::my_iaf_cond_exp::handle(CurrentEvent& e)
{
	assert(e.get_delay() > 0);

	// add weighted current; HEP 2002-10-04
	B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
			e.get_weight() * e.get_current());

	assert(e.get_delay() > 0);
	// Assert that port is 0 (SUP_SPIKE_RECEPTOR (4) - MIN_SPIKE_RECEPTOR (3) = 1)
	assert(0 <= e.get_rport() && e.get_rport() < SUP_CURR_RECEPTOR - MIN_CURR_RECEPTOR);

}

void mynest::my_iaf_cond_exp::handle(DataLoggingRequest& e)
{
	B_.logger_.handle(e);
}

#endif //HAVE_GSL
