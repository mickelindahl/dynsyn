/*
 *  stdp_dop_connection.h
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

#ifndef STDP_DOP_CONNECTION_H
#define STDP_DOP_CONNECTION_H

/* BeginDocumentation
  Name: stdp_dop_synapse - Synapse type for spike-timing dependent
   plasticity with dopamine modulation.

  Description:
   stdp_synapse is a connector to create synapses with spike time 
   dependent plasticity (as defined in [1]). Here the weight dependence
   exponent can be set separately for potentiation and depression.

  Examples:
   multiplicative STDP mu_plus = mu_minus = 1.0
   additive STDP       mu_plus = mu_minus = 0.0
   Guetig STDP         mu_plus = mu_minus = [0,1]
   van Rossum STDP     mu_plus = 0.0 mu_minus = 1.0 

  Parameters:
   tau_plus     time constant of STDP window, potentiation 
                (tau_minus defined in post-synaptic neuron)
   lambda       step size
   alpha        asymmetry parameter (scales depressing increments as alpha*lambda)
   mu_plus      weight dependence exponent, potentiation
   mu_minus     weight dependence exponent, depression
   Wmax         maximum allowed weight

  Transmits: SpikeEvent

  References:
   [1] Guetig et al. (2003) Learning Input Correlations through Nonlinear
       Temporally Asymmetric Hebbian Plasticity. Journal of Neuroscience

	Based on stdp_synapse by Moritz Helias, Abigail Morrison March 2006

  FirstVersion: June 2010
  Author: Mikael Lindahl, Wiebke Potjan
  SeeAlso: synapsedict, tsodyks_synapse, static_synapse
 */

#include "connection_het_wd.h"
#include "archiving_node.h"
#include "generic_connector.h"
#include <cmath>

namespace mynest
{
class stdp_dop_connection : public nest::ConnectionHetWD
{

public:
	/**
	 * Default Constructor.
	 * Sets default values for all parameters. Needed by GenericConnectorModel.
	 */
	stdp_dop_connection();

	/**
	 * Copy constructor.
	 * Needs to be defined properly in order for GenericConnector to work.
	 */
	stdp_dop_connection(const stdp_dop_connection &);

	/**
	 * Default Destructor.
	 */
	~stdp_dop_connection() {}

	void check_connection(nest::Node & s, nest::Node & r, nest::port receptor_type, nest::double_t t_lastspike);

	/**
	 * Get all properties of this connection and put them into a dictionary.
	 */
	void get_status(DictionaryDatum & d) const;

	/**
	 * Set properties of this connection from the values given in dictionary.
	 */
	void set_status(const DictionaryDatum & d, nest::ConnectorModel &cm);

	/**
	 * Set properties of this connection from position p in the properties
	 * array given in dictionary.
	 */
	void set_status(const DictionaryDatum & d, nest::index p, nest::ConnectorModel &cm);

	/**
	 * Create new empty arrays for the properties of this connection in the given
	 * dictionary. It is assumed that they are not existing before.
	 */
	void initialize_property_arrays(DictionaryDatum & d) const;

	/**
	 * Append properties of this connection to the given dictionary. If the
	 * dictionary is empty, new arrays are created first.
	 */
	void append_properties(DictionaryDatum & d) const;

	/**
	 * Send an event to the receiver of this connection.
	 * \param e The event to send
	 * \param t_lastspike Point in time of last spike sent.
	 * \param cp common properties of all synapses (empty).
	 */
	void send(nest::Event& e, nest::double_t t_lastspike, const nest::CommonSynapseProperties &cp);

	// overloaded for all supported event types
	using nest::Connection::check_event;
	void check_event(nest::SpikeEvent&) {}

private:

	nest::double_t facilitate_(nest::double_t w, nest::double_t kplus);
	nest::double_t depress_(nest::double_t w, nest::double_t kminus);

	// data members of each connection
	nest::double_t tau_plus_;
	nest::double_t lambda_;
	nest::double_t alpha_;
	nest::double_t mu_plus_;
	nest::double_t mu_minus_;
	nest::double_t Wmax_;
	nest::double_t Kplus_;

};


inline
nest::double_t stdp_dop_connection::facilitate_(nest::double_t w, nest::double_t kplus)
{
	nest::double_t norm_w = (w / Wmax_) + (lambda_ * std::pow(1.0 - (w/Wmax_), mu_plus_) * kplus);
	return norm_w < 1.0 ? norm_w * Wmax_ : Wmax_;
}

inline 
nest::double_t stdp_dop_connection::depress_(nest::double_t w, nest::double_t kminus)
{
	nest::double_t norm_w = (w / Wmax_) - (alpha_ * lambda_ * std::pow(w/Wmax_, mu_minus_) * kminus);
	return norm_w > 0.0 ? norm_w * Wmax_ : 0.0;
}


inline 
void stdp_dop_connection::check_connection(nest::Node & s, nest::Node & r, nest::port receptor_type, nest::double_t t_lastspike)
{
	nest::ConnectionHetWD::check_connection(s, r, receptor_type, t_lastspike);

	// For a new synapse, t_lastspike contains the point in time of the last spike.
	// So we initially read the history(t_last_spike - dendritic_delay, ...,  T_spike-dendritic_delay]
	// which increases the access counter for these entries.
	// At registration, all entries' access counters of history[0, ..., t_last_spike - dendritic_delay] will be 
	// incremented by the following call to Archiving_Node::register_stdp_connection().
	// See bug #218 for details.
	// Here do not change name to register_stdp_dop_connection
	r.register_stdp_connection(t_lastspike - nest::Time(nest::Time::step(delay_)).get_ms());
}

/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param p The port under which this connection is stored in the Connector.
 * \param t_lastspike Time point of last spike emitted
 */
inline
void stdp_dop_connection::send(nest::Event& e, nest::double_t t_lastspike, const nest::CommonSynapseProperties &)
{
	// synapse STDP depressing/facilitation dynamics

	nest::double_t t_spike = e.get_stamp().get_ms();
	// t_lastspike_ = 0 initially
	nest::double_t dendritic_delay = nest::Time(nest::Time::step(delay_)).get_ms();

	//get spike history in relevant range (t1, t2] from post-synaptic neuron
	std::deque<nest::histentry>::iterator start;
	std::deque<nest::histentry>::iterator finish;

	// For a new synapse, t_lastspike contains the point in time of the last spike.
	// So we initially read the history(t_last_spike - dendritic_delay, ...,  T_spike-dendritic_delay]
	// which increases the access counter for these entries.
	// At registration, all entries' access counters of history[0, ..., t_last_spike - dendritic_delay] have been 
	// incremented by Archiving_Node::register_stdp_dop_connection(). See bug #218 for details.
	target_->get_history(t_lastspike - dendritic_delay, t_spike - dendritic_delay,
			&start, &finish);
	//facilitation due to post-synaptic spikes since last pre-synaptic spike
	nest::double_t minus_dt;
	while (start != finish)
	{
		minus_dt = t_lastspike - (start->t_ + dendritic_delay);
		start++;
		if (minus_dt == 0)
			continue;
		weight_ = facilitate_(weight_, Kplus_ * std::exp(minus_dt / tau_plus_));
	}

	//depression due to new pre-synaptic spike
	weight_ = depress_(weight_, target_->get_K_value(t_spike - dendritic_delay));

	e.set_receiver(*target_);
	e.set_weight(weight_);
	e.set_delay(delay_);
	e.set_rport(rport_);
	e();


	Kplus_ = Kplus_ * std::exp((t_lastspike - t_spike) / tau_plus_) + 1.0;
}

} // of namespace nest

#endif // of #ifndef STDP_DOP_CONNECTION_H
