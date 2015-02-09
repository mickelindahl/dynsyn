/*
 *  pif_psc_alpha.cpp
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2008 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */

#include "exceptions.h"
#include "pif_psc_alpha.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
#include "analog_data_logger_impl.h"
#include "lockptrdatum.h"

#include <limits>

using namespace nest;

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */
    
mynest::pif_psc_alpha::Parameters_::Parameters_()
  : C_m    (250.0),  // pF
    I_e    (  0.0),  // nA
    tau_syn(  2.0),  // ms
    V_th   (-55.0),  // mV
    V_reset(-70.0),  // mV
    t_ref  (  2.0)   // ms   
  {}

mynest::pif_psc_alpha::State_::State_(const Parameters_& p)
  : V_m       (p.V_reset),
    dI_syn    (0.0),
    I_syn     (0.0),
    I_ext     (0.0),
    refr_count(0  )
{}

/* ---------------------------------------------------------------- 
 * Paramater and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void mynest::pif_psc_alpha::Parameters_::get(DictionaryDatum &d) const
{
  (*d)[names::C_m    ] = C_m;
  (*d)[names::I_e    ] = I_e;
  (*d)[names::tau_syn] = tau_syn;
  (*d)[names::V_th   ] = V_th;
  (*d)[names::V_reset] = V_reset;
  (*d)[names::t_ref  ] = t_ref;
}

void mynest::pif_psc_alpha::Parameters_::set(const DictionaryDatum& d)
{
  updateValue<double>(d, names::C_m, C_m);
  updateValue<double>(d, names::I_e, I_e);
  updateValue<double>(d, names::tau_syn, tau_syn);
  updateValue<double>(d, names::V_th   , V_th);
  updateValue<double>(d, names::V_reset, V_reset);
  updateValue<double>(d, names::t_ref, t_ref);

  if ( C_m <= 0 )
    throw nest::BadProperty("The membrane capacitance must be strictly positive.");

  if ( tau_syn <= 0 )
    throw nest::BadProperty("The synaptic time constant must be strictly positive.");

  if ( V_reset >= V_th )
    throw nest::BadProperty("The reset potential must be below threshold.");
  
  if ( t_ref < 0 )
    throw nest::BadProperty("The refractory time must be at least one simulation step.");  
}

void mynest::pif_psc_alpha::State_::get(DictionaryDatum &d) const
{
  // Only the membrane potential is shown in the status; one could show also the other
  // state variables
  (*d)[names::V_m] = V_m;   
}

void mynest::pif_psc_alpha::State_::set(const DictionaryDatum& d, const Parameters_& p)
{
  // Only the membrane potential can be set; one could also make other state variables
  // settable.
  updateValue<double>(d, names::V_m, V_m);
}

/* ---------------------------------------------------------------- 
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

mynest::pif_psc_alpha::pif_psc_alpha()
  : Node(), 
    P_(), 
    S_(P_)
{}

mynest::pif_psc_alpha::pif_psc_alpha(const pif_psc_alpha& n)
  : Node(n), 
    P_(n.P_), 
    S_(n.S_)
{}

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void mynest::pif_psc_alpha::init_node_(const Node& proto)
{
  const pif_psc_alpha& pr = downcast<pif_psc_alpha>(proto);
  P_ = pr.P_;
  S_ = pr.S_;
}

void mynest::pif_psc_alpha::init_state_(const Node& proto)
{
  const pif_psc_alpha& pr = downcast<pif_psc_alpha>(proto);
  S_ = pr.S_;
}

void mynest::pif_psc_alpha::init_buffers_()
{
  B_.spikes.clear();    // includes resize
  B_.currents.clear();  // include resize
  B_.potentials.clear_data(); // includes resize
}

void mynest::pif_psc_alpha::calibrate()
{
  const double h  = Time::get_resolution().get_ms(); 
  const double eh = std::exp(-h/P_.tau_syn);
  const double tc = P_.tau_syn / P_.C_m;

  // compute matrix elements, all other elements 0
  V_.P11 = eh;
  V_.P22 = eh;
  V_.P21 = h * eh;
  V_.P30 = h / P_.C_m;
  V_.P31 = tc * ( P_.tau_syn - (h+P_.tau_syn) * eh );
  V_.P32 = tc * ( 1 - eh );
  // P33_ is 1
    
  // initial value ensure normalization to max amplitude 1.0
  V_.pscInitialValue = 1.0 * numerics::e / P_.tau_syn;

  // refractory time in steps
  V_.t_ref_steps = Time(Time::ms(P_.t_ref)).get_steps();
  assert(V_.t_ref_steps >= 0);  // since t_ref_ >= 0, this can only fail in error
}

/* ---------------------------------------------------------------- 
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void mynest::pif_psc_alpha::update(Time const& slice_origin, 
                                   const nest::long_t from_step, 
                                   const nest::long_t to_step)
{
  for ( long lag = from_step ; lag < to_step ; ++lag )
  {
    // order is important in this loop, since we have to use the old values
    // (those upon entry to the loop) on right hand sides everywhere
    
    // update membrane potential
    if ( S_.refr_count == 0 )  // neuron absolute not refractory
      S_.V_m +=   V_.P30 * ( S_.I_ext + P_.I_e )
                  + V_.P31 * S_.dI_syn
                  + V_.P32 * S_.I_syn;
    else 
     --S_.refr_count;  // count down refractory time

    // update synaptic currents
    S_.I_syn   = V_.P21 * S_.dI_syn + V_.P22 * S_.I_syn;
    S_.dI_syn *= V_.P11;

    // log membrane potential
    B_.potentials.record_data(slice_origin.get_steps()+lag, S_.V_m);
    
    // check for threshold crossing
    if ( S_.V_m >= P_.V_th )
    {
      // reset neuron
      S_.refr_count = V_.t_ref_steps;
      S_.V_m        = P_.V_reset;
      
      // send spike
      SpikeEvent se;
      network()->send(*this, se, lag);
    }

    // add synaptic input currents for this step 
    S_.dI_syn += V_.pscInitialValue * B_.spikes.get_value(lag);

    // set new input current
    S_.I_ext = B_.currents.get_value(lag);
    
    // log membrane potential
    B_.potentials.record_data(slice_origin.get_steps()+lag, S_.V_m);
  }  
}                           

void mynest::pif_psc_alpha::handle(SpikeEvent & e)
{
  assert(e.get_delay() > 0);

  B_.spikes.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
                      e.get_weight());
}

void mynest::pif_psc_alpha::handle(CurrentEvent& e)
{
  assert(e.get_delay() > 0);

  B_.currents.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), 
		                    e.get_weight() * e.get_current());
}

// Do not move this function as inline to h-file. It depends on 
// analog_data_logger_impl.h being included here.
void mynest::pif_psc_alpha::handle(PotentialRequest& e)
{
  B_.potentials.handle(*this, e);  // the logger does this for us
}