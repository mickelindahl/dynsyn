/*
 *  tsodyks_stocastic_connection.cpp
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

#include "network.h"
#include "dictdatum.h"
#include "connector_model.h"
#include "common_synapse_properties.h"
#include "tsodyks_stocastic_connection.h"
#include "event.h"

namespace mynest
{

  tsodyks_stocastic_connection::tsodyks_stocastic_connection() :
		nest::ConnectionHetWD(),
    tau_psc_(3.0),
    tau_fac_(0.0),
    tau_rec_(800.0),
    U_(0.5),
    x_(1.0),
    y_(0.0),
    u_(0.0),

    p_(1.0),
    w_(1.0),
    m_(1.0)
  { }

  void tsodyks_stocastic_connection::get_status(DictionaryDatum & d) const
  {
  	nest::ConnectionHetWD::get_status(d);

    def<nest::double_t>(d, "U", U_);
    def<nest::double_t>(d, "tau_psc", tau_psc_);
    def<nest::double_t>(d, "tau_rec", tau_rec_);
    def<nest::double_t>(d, "tau_fac", tau_fac_);
    def<nest::double_t>(d, "x", x_);
    def<nest::double_t>(d, "y", y_);
    def<nest::double_t>(d, "u", u_);

    def<nest::double_t>(d, "p", p_);
    def<nest::double_t>(d, "w", w_);
    def<nest::double_t>(d, "m", m_);
  }
  
  void tsodyks_stocastic_connection::set_status(const DictionaryDatum & d, nest::ConnectorModel &cm)
  {
  	nest::ConnectionHetWD::set_status(d, cm);

    updateValue<nest::double_t>(d, "U", U_);
    updateValue<nest::double_t>(d, "tau_psc", tau_psc_);
    updateValue<nest::double_t>(d, "tau_rec", tau_rec_);
    updateValue<nest::double_t>(d, "tau_fac", tau_fac_);
    updateValue<nest::double_t>(d, "p", p_);
    updateValue<nest::double_t>(d, "w", w_);
    updateValue<nest::double_t>(d, "m", m_);


    nest::double_t x = x_;
    nest::double_t y = y_;
    updateValue<double_t>(d, "x", x);
    updateValue<double_t>(d, "y", y);

    if (x + y > 1.0)
    {
      cm.network().message(SLIInterpreter::M_ERROR,
			   "tsodyks_stocastic_connection::set_status()", "x + y must be <= 1.0.");
      throw nest::BadProperty();
    }
    else
    {
      x_ = x;
      y_ = y;
    }

    updateValue<nest::double_t>(d, "u", u_);
  }

  /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */  
  void tsodyks_stocastic_connection::set_status(const DictionaryDatum & d, nest::index p, nest::ConnectorModel &cm)
  {
  	nest::ConnectionHetWD::set_status(d, p, cm);

  	nest::set_property<nest::double_t>(d, "Us", p, U_);
  	nest::set_property<nest::double_t>(d, "tau_pscs", p, tau_psc_);
  	nest::set_property<nest::double_t>(d, "tau_facs", p, tau_fac_);
  	nest::set_property<nest::double_t>(d, "tau_recs", p, tau_rec_);

  	nest::set_property<nest::double_t>(d, "ps", p, p_);
  	nest::set_property<nest::double_t>(d, "ws", p, w_);
  	nest::set_property<nest::double_t>(d, "ms", p, m_);


  	nest::double_t x = x_;
    nest::double_t y = y_;
    nest::set_property<nest::double_t>(d, "xs", p, x);
    nest::set_property<nest::double_t>(d, "ys", p, y);

    if (x + y > 1.0)
    {
      cm.network().message(SLIInterpreter::M_ERROR, 
			   "tsodyks_stocastic_connection::set_status()",
			   "x + y must be <= 1.0.");
      throw nest::BadProperty();
    }
    else
    {
      x_ = x;
      y_ = y;
    }

    nest::set_property<nest::double_t>(d, "us", p, u_);
  }

  void tsodyks_stocastic_connection::initialize_property_arrays(DictionaryDatum & d) const
  {
  	nest::ConnectionHetWD::initialize_property_arrays(d);

    initialize_property_array(d, "Us"); 
    initialize_property_array(d, "tau_pscs");
    initialize_property_array(d, "tau_facs");
    initialize_property_array(d, "tau_recs");  
    initialize_property_array(d, "xs"); 
    initialize_property_array(d, "ys");
    initialize_property_array(d, "us");

    initialize_property_array(d, "ps");
    initialize_property_array(d, "ws");
    initialize_property_array(d, "ms");
  }

  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void tsodyks_stocastic_connection::append_properties(DictionaryDatum & d) const
  {
  	nest::ConnectionHetWD::append_properties(d);

    append_property<double_t>(d, "Us", U_); 
    append_property<double_t>(d, "tau_pscs", tau_psc_);
    append_property<double_t>(d, "tau_facs", tau_fac_);
    append_property<double_t>(d, "tau_recs", tau_rec_);  
    append_property<double_t>(d, "xs", x_); 
    append_property<double_t>(d, "ys", y_);
    append_property<double_t>(d, "us", u_);

    append_property<double_t>(d, "ps", x_);
    append_property<double_t>(d, "ws", y_);
    append_property<double_t>(d, "ms", u_);
  }

} // of namespace nest
