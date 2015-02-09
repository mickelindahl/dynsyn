/*
 *  tsodyks_connection.h
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

#ifndef TSODYKS_STOCASTIC_CONNECTION_H
#define TSODYKS_STOCASTIC_CONNECTION_H



/* BeginDocumentation
  Name: tsodyks_synapse - Synapse type with short term plasticity.

  Description:
   This synapse model implements synaptic short-term depression and short-term facilitation
   according to [1]. In particular it solves Eqs (3) and (4) from this paper in an
   exact manner.

   Synaptic depression is motivated by depletion of vesicles in the readily releasable pool
   of synaptic vesicles (variable x in equation (3)). Synaptic facilitation comes about by
   a presynaptic increase of release probability, which is modeled by variable U in Eq (4).
   The original interpretation of variable y is the amount of glutamate concentration in
   the synaptic cleft. In [1] this variable is taken to be directly proportional to the
   synaptic current caused in the postsynaptic neuron (with the synaptic weight w as a
   proportionality constant). In order to reproduce the results of [1] and to use this
   model of synaptic plasticity in its original sense, the user therefore has to ensure
   the following conditions:

   1.) The postsynaptic neuron must be of type iaf_psc_exp or iaf_tum_2000, because
   these neuron models have a postsynaptic current which decays exponentially.

   2.) The time constant of each tsodyks_synapse targeting a particular neuron
   must be chosen equal to that neuron's synaptic time constant. In particular that means
   that all synapses targeting a particular neuron have the same parameter tau_psc.

   However, there are no technical restrictions using this model of synaptic plasticity
   also in conjunction with neuron models that have a different dynamics for their synaptic
   current or conductance. The effective synaptic weight, which will be transmitted
   to the postsynaptic neuron upon occurrence of a spike at time t is u(t)*x(t)*w, where
   u(t) and x(t) are defined in Eq (3) and (4), w is the synaptic weight specified upon
   connection.
   The interpretation is as follows: The quantity u(t)*x(t) is the release probability
   times the amount of releasable synaptic vesicles at time t of the presynaptic neuron's
   spike, so this equals the amount of transmitter expelled into the synaptic cleft.
   The amount of transmitter than relaxes back to 0 with time constant tau_psc of the
   synapse's variable y.
   Since the dynamics of y(t) is linear, the postsynaptic neuron can reconstruct from the
   amplitude of the synaptic impulse u(t)*x(t)*w the full shape of y(t).
   The postsynaptic neuron, however, might choose to have a synaptic current that is not
   necessarily identical to the concentration of transmitter y(t) in the synaptic cleft.
   It may realize an arbitrary postsynaptic effect depending on y(t).

   Parameters: 
     The following parameters can be set in the status dictionary:
     Us         double - maximum probability of realease [0,1]
     tau_pscs   double - time constants of synaptic current in ms
     tau_facs   double - time constant for facilitation in ms
     tau_recs   double - time constant for depression in ms
     xs         double - initial fraction of synaptic vesicles in the readily releasable pool [0,1]
     ys         double - initial fraction of synaptic vesicles in the synaptic cleft [0,1]

  References:
   [1] Tsodyks, Uziel, Markram (2000) Synchrony Generation in Recurrent Networks
       with Frequency-Dependent Synapses. Journal of Neuroscience, vol 20 RC50

  Transmits: SpikeEvent

  FirstVersion: March 2006
  Author: Moritz Helias
  SeeAlso: synapsedict, stdp_synapse, static_synapse, iaf_psc_exp, iaf_tum_2000
 */


/**
 * Class representing a synapse with Tsodyks short term plasticity.
 * A suitale Connector containing these connections can be obtained from the template GenericConnector.
 */

#include "connection_het_wd.h"
#include "archiving_node.h"
#include "generic_connector.h"
#include <cmath>

namespace mynest 
{
class tsodyks_stocastic_connection : public nest::ConnectionHetWD
{
public:

	/**
	 * Default Constructor.
	 * Sets default values for all parameters. Needed by GenericConnectorModel.
	 */
	tsodyks_stocastic_connection();

	/**
	 * Default Destructor.
	 */
	~tsodyks_stocastic_connection() {}

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
	 * \param cp Common properties to all synapses (empty).
	 */
	void send(nest::Event& e, nest::double_t t_lastspike, const nest::CommonSynapseProperties &cp);

	// overloaded for all supported event types
	using nest::Connection::check_event;
	void check_event(nest::SpikeEvent&) {}

private:
	nest::double_t tau_psc_;   //!< [ms] time constant of postsyn current
	nest::double_t tau_fac_;   //!< [ms] time constant for fascilitation
	nest::double_t tau_rec_;   //!< [ms] time constant for recovery
	nest::double_t U_;         //!< asymptotic value of probability of release
	nest::double_t x_;         //!< amount of resources in recovered state
	nest::double_t y_;         //!< amount of resources in active state
	nest::double_t u_;         //!< actual probability of release
	nest::double_t p_;         //!< probability
	nest::double_t w_;         //!< weight
	nest::double_t m_;         //!< number of synaptic contacts
};


/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param p The port under which this connection is stored in the Connector.
 * \param t_lastspike Time point of last spike emitted
 */
inline
void tsodyks_stocastic_connection::send(nest::Event& e, nest::double_t t_lastspike, const nest::CommonSynapseProperties &)
{
	nest::double_t h = e.get_stamp().get_ms() - t_lastspike;

	// t_lastspike_ = 0 initially
	// this has no influence on the dynamics, IF y = z = 0 initially
	// !!! x != 1.0 -> z != 0.0 -> t_lastspike_=0 has influence on dynamics

	// propagator
	// TODO: use expm1 here instead, where applicable
	nest::double_t Puu = (tau_fac_ == 0.0) ? 0.0 : std::exp(-h/tau_fac_);
	nest::double_t Pyy = std::exp(-h/tau_psc_);
	nest::double_t Pzz = std::exp(-h/tau_rec_);

	//double_t Pzy = (Pyy - Pzz) * tau_rec_ / (tau_psc_ - tau_rec_);
	nest::double_t Pxy = ((Pzz - 1.0)*tau_rec_ - (Pyy - 1.0)*tau_psc_) / (tau_psc_ - tau_rec_);
	nest::double_t Pxz = 1.0 - Pzz;

	nest::double_t z = 1.0 - x_ - y_;

	// propagation t_lastspike -> t_spike
	// don't change the order !

	u_ *= Puu;
	x_ += Pxy * y_ + Pxz * z;
	//z = Pzz * z_ + Pzy * y_;
	y_ *= Pyy;

	// delta function u
	u_ += U_*(1.0-u_);

  p_ =u_;

  w_=0;
  for ( nest::long_t i = 1 ; i <= m_ ; ++i )
  {
  	// Generate random number
  	nest::double_t r = std::rand();
  	nest::double_t r_max=(RAND_MAX+1);
  	r=r/r_max;
//r=r/(RAND_MAX+1);
  	// Assign weight depending on r_ and p_
  	if (p_>r)
  		w_ +=U_/m_;

  }

	// postsynaptic current step caused by incoming spike
	//nest::double_t delta_y_tsp = u_*x_;
  nest::double_t delta_y_tsp = w_*x_;
	// delta function x, y
	x_ -= delta_y_tsp;
	y_ += delta_y_tsp;

	// send the spike to the target
	e.set_receiver(*target_);
	e.set_weight( weight_ * delta_y_tsp );
	e.set_delay( delay_ );
	e.set_rport( rport_ );
	e();
}

} // namespace

#endif // TSODYKS_STOCASTIC_CONNECTION_H
