#ifndef DROP_ODD_SPIKE_CONNECTION_H
#define DROP_ODD_SPIKE_CONNECTION_H

/*
 *  drop_odd_spike_connection.h
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

#include "connection_het_wd.h"

/* BeginDocumentation
  Name: drop_odd_spike - Synapse dropping spikes with odd time stamps.

  Description:
  This synapse will not deliver any spikes with odd time stamps, while spikes with even
  time stamps go through unchanged. 
  
  Remarks:
  This synapse type is provided only for illustration purposes in MyModule.

  SeeAlso: synapsedict
*/

namespace mynest {

/**
 * Connection class for illustration purposes.
 */
class DropOddSpikeConnection : public nest::ConnectionHetWD
{
 public:

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  DropOddSpikeConnection() : ConnectionHetWD() {}

  /**
   * Default Destructor.
   */
  ~DropOddSpikeConnection() {}

  /**
   * Send an event to the receiver of this connection.
   * @param e The event to send
   * @param t_lastspike Point in time of last spike sent.
   * @param cp Common properties to all synapses (empty).
   */
  void send(nest::Event& e, nest::double_t t_lastspike, 
            const nest::CommonSynapseProperties &cp);

  //! Defining this as empty means we can handle spike events
  using Connection::check_event;  // see http://www.gotw.ca/gotw/005.htm
  void check_event(nest::SpikeEvent&) {}

};


inline
void DropOddSpikeConnection::send(nest::Event& e, nest::double_t last, 
                                  const nest::CommonSynapseProperties &props)
{
  if ( e.get_stamp().get_steps() % 2 ) // stamp is odd, drop it
    return;
    
  // Even time stamp, we send the spike using the normal sending mechnism  
  // send the spike to the target
  nest::ConnectionHetWD::send(e, last, props);
}
 
} // namespace

#endif // drop_odd_spike_connection.h
