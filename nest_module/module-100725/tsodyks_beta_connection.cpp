/*
 *  tsodyks_beta_connection.cpp
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

#include "tsodyks_beta_connection.h"
#include "network.h"
#include "connector_model.h"

namespace mynest
{

  TsodyksConnection::TsodyksConnection() :
    nest::ConnectionHetWD(),
    tau_psc_rise_(3.0),
    tau_psc_decay_(3.0),
    tau_fac_(0.0),
    tau_rec_(800.0),
    U_(0.5),
    x_(1.0),
    y_(0.0),
    u_(0.0)
  { }

  void TsodyksConnection::get_status(DictionaryDatum & d) const
  {
    nest::ConnectionHetWD::get_status(d);

    def<nest::double_t>(d, "U",             U_);
    def<nest::double_t>(d, "tau_psc_rise",  tau_psc_rise_);
    def<nest::double_t>(d, "tau_psc_decay", tau_psc_decay_);
    def<nest::double_t>(d, "tau_rec",       tau_rec_);
    def<nest::double_t>(d, "tau_fac",       tau_fac_);
    def<nest::double_t>(d, "x",             x_);
    def<nest::double_t>(d, "y",             y_);
    def<nest::double_t>(d, "u",             u_);
  }
  
  void TsodyksConnection::set_status(const DictionaryDatum & d, nest::ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, cm);

    updateValue<nest::double_t>(d, "U", U_);
    updateValue<nest::double_t>(d, "tau_psc_rise", tau_psc_rise_);
    updateValue<nest::double_t>(d, "tau_psc_decay", tau_psc_decay_);
    updateValue<nest::double_t>(d, "tau_rec", tau_rec_);
    updateValue<nest::double_t>(d, "tau_fac", tau_fac_);

    nest::double_t x = x_;
    nest::double_t y = y_;
    updateValue<nest::double_t>(d, "x", x);
    updateValue<nest::double_t>(d, "y", y);

    if (x + y > 1.0)
    {
      cm.network().message(SLIInterpreter::M_ERROR,
			   "TsodyksConnection::set_status()", "x + y must be <= 1.0.");
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
  void TsodyksConnection::set_status(const DictionaryDatum & d, nest::index p, nest::ConnectorModel &cm)
  {
  	nest::ConnectionHetWD::set_status(d, p, cm);

    nest::set_property<nest::double_t>(d, "Us", p, U_);
    nest::set_property<nest::double_t>(d, "tau_psc_rises", p, tau_psc_rise_);
    nest::set_property<nest::double_t>(d, "tau_psc_decays", p, tau_psc_decay_);
    nest::set_property<nest::double_t>(d, "tau_facs", p, tau_fac_);
    nest::set_property<nest::double_t>(d, "tau_recs", p, tau_rec_);


    nest::double_t x = x_;
    nest::double_t y = y_;
    nest::set_property<nest::double_t>(d, "xs", p, x);
    nest::set_property<nest::double_t>(d, "ys", p, y);

    if (x + y > 1.0)
    {
      cm.network().message(SLIInterpreter::M_ERROR, 
			   "TsodyksConnection::set_status()", 
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

  void TsodyksConnection::initialize_property_arrays(DictionaryDatum & d) const
  {
  	nest::ConnectionHetWD::initialize_property_arrays(d);

  	initialize_property_array(d, "Us"             );
  	initialize_property_array(d, "tau_psc_rises"  );
  	initialize_property_array(d, "tau_psc_decays" );
  	initialize_property_array(d, "tau_facs"       );
  	initialize_property_array(d, "tau_recs"       );
  	initialize_property_array(d, "xs"             );
  	initialize_property_array(d, "ys"             );
  	initialize_property_array(d, "us"             );
  }

  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void TsodyksConnection::append_properties(DictionaryDatum & d) const
  {
  	nest::ConnectionHetWD::append_properties(d);

  	append_property<nest::double_t>(d, "Us",             U_             );
  	append_property<nest::double_t>(d, "tau_psc_rises",  tau_psc_rise_  );
  	append_property<nest::double_t>(d, "tau_psc_decays", tau_psc_decay_ );
  	append_property<nest::double_t>(d, "tau_facs",       tau_fac_       );
  	append_property<nest::double_t>(d, "tau_recs",       tau_rec_       );
  	append_property<nest::double_t>(d, "xs",             x_             );
  	append_property<nest::double_t>(d, "ys",             y_             );
  	append_property<nest::double_t>(d, "us",             u_             );
  }

} // of namespace mynest
