/*
 *  my_poisson_generator.cpp
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2004 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */

#include "my_poisson_generator.h"
#include "network.h"
#include "dict.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "exceptions.h"

using namespace nest;
/* ----------------------------------------------------------------
 * Default constructors defining default parameter
 * ---------------------------------------------------------------- */

nest::my_poisson_generator::Parameters_::Parameters_()
  : rate_(0.0    )  // pA
{}


/* ----------------------------------------------------------------
 * Parameter extraction and manipulation functions
 * ---------------------------------------------------------------- */

void nest::my_poisson_generator::Parameters_::get(DictionaryDatum &d) const
{
  def<double>(d, names::rate, rate_);
}

void nest::my_poisson_generator::Parameters_::set(const DictionaryDatum& d)
{
  updateValue<double>(d, names::rate, rate_);
  if ( rate_ < 0 )
    throw BadProperty("The rate cannot be negative.");
}


/* ----------------------------------------------------------------
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

nest::my_poisson_generator::my_poisson_generator()
  : Node(),
    device_(),
    P_()
{}

nest::my_poisson_generator::my_poisson_generator(const my_poisson_generator& n)
  : Node(n),
    device_(n.device_),
    P_(n.P_)
{}


/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void nest::my_poisson_generator::init_node_(const Node& proto)
{
  const my_poisson_generator& pr = downcast<my_poisson_generator>(proto);

  device_.init_parameters(pr.device_);

  P_ = pr.P_;
}

void nest::my_poisson_generator::init_state_(const Node& proto)
{
  const my_poisson_generator& pr = downcast<my_poisson_generator>(proto);

  device_.init_state(pr.device_);
}

void nest::my_poisson_generator::init_buffers_()
{
  device_.init_buffers();
}

void nest::my_poisson_generator::calibrate()
{
  device_.calibrate();

  // rate_ is in Hz, dt in ms, so we have to convert from s to ms
  V_.poisson_dev_.set_lambda(Time::get_resolution().get_ms() * P_.rate_ * 1e-3);
}


/* ----------------------------------------------------------------
 * Update function and event hook
 * ---------------------------------------------------------------- */

void nest::my_poisson_generator::update(Time const & T, const long_t from, const long_t to)
{
  assert(to >= 0 && (delay) from < Scheduler::get_min_delay());
  assert(from < to);

  if ( P_.rate_ <= 0 )
    return;

  for ( long_t lag = from ; lag < to ; ++lag )
  {
    if ( !device_.is_active( T + Time::step(lag) ) )
      continue;  // no spike at this lag

    DSSpikeEvent se;
    network()->send(*this, se, lag);
  }
}

void nest::my_poisson_generator::event_hook(DSSpikeEvent& e)
{
  librandom::RngPtr rng = net_->get_rng(get_thread());
  ulong_t n_spikes = V_.poisson_dev_.uldev(rng);

  if ( n_spikes > 0 ) // we must not send events with multiplicity 0
  {
    e.set_multiplicity(n_spikes);
    e.get_receiver().handle(e);
  }
}
