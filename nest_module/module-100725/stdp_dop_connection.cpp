
/*
 *  stdp_dop_connection.cpp
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
#include "stdp_dop_connection.h"
#include "event.h"


namespace mynest
{

  stdp_dop_connection::stdp_dop_connection() :
    nest::ConnectionHetWD(),
    tau_plus_(20.0),
    lambda_(0.01),
    alpha_(1.0),
    mu_plus_(1.0),
    mu_minus_(1.0),
    Wmax_(100.0),
    Kplus_(0.0)
  { }


  stdp_dop_connection::stdp_dop_connection(const stdp_dop_connection &rhs) :
    nest::ConnectionHetWD(rhs)
  {
    tau_plus_ = rhs.tau_plus_;
    lambda_ = rhs.lambda_;
    alpha_ = rhs.alpha_;
    mu_plus_ = rhs.mu_plus_;
    mu_minus_ = rhs.mu_minus_;
    Wmax_ = rhs.Wmax_;
    Kplus_ = rhs.Kplus_;
  }

  void stdp_dop_connection::get_status(DictionaryDatum & d) const
  {
    nest::ConnectionHetWD::get_status(d);
    def<nest::double_t>(d, "tau_plus", tau_plus_);
    def<nest::double_t>(d, "lambda", lambda_);
    def<nest::double_t>(d, "alpha", alpha_);
    def<nest::double_t>(d, "mu_plus", mu_plus_);
    def<nest::double_t>(d, "mu_minus", mu_minus_);
    def<nest::double_t>(d, "Wmax", Wmax_);
    def<nest::double_t>(d, "Kplus", Kplus_);
  }

  void stdp_dop_connection::set_status(const DictionaryDatum & d, nest::ConnectorModel &cm)
  {
    nest::ConnectionHetWD::set_status(d, cm);
    updateValue<nest::double_t>(d, "tau_plus", tau_plus_);
    updateValue<nest::double_t>(d, "lambda", lambda_);
    updateValue<nest::double_t>(d, "alpha", alpha_);
    updateValue<nest::double_t>(d, "mu_plus", mu_plus_);
    updateValue<nest::double_t>(d, "mu_minus", mu_minus_);
    updateValue<nest::double_t>(d, "Wmax", Wmax_);
    updateValue<nest::double_t>(d, "Kplus", Kplus_);
  }

   /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */
  void stdp_dop_connection::set_status(const DictionaryDatum & d, nest::index p, nest::ConnectorModel &cm)
  {
    nest::ConnectionHetWD::set_status(d, p, cm);

    nest::set_property<nest::double_t>(d, "tau_pluss", p, tau_plus_);
    nest::set_property<nest::double_t>(d, "lambdas", p, lambda_);
    nest::set_property<nest::double_t>(d, "alphas", p, alpha_);
    nest::set_property<nest::double_t>(d, "mu_pluss", p, mu_plus_);
    nest::set_property<nest::double_t>(d, "mu_minuss", p, mu_minus_);
    nest::set_property<nest::double_t>(d, "Wmaxs", p, Wmax_);
    nest::set_property<nest::double_t>(d, "Kpluss", p, Kplus_);
  }

  void stdp_dop_connection::initialize_property_arrays(DictionaryDatum & d) const
  {
  	nest::ConnectionHetWD::initialize_property_arrays(d);

  	nest::initialize_property(d, "tau_pluss");
  	nest::initialize_property(d, "lambdas");
  	nest::initialize_property(d, "alphas");
  	nest::initialize_property(d, "mu_pluss");
  	nest::initialize_property(d, "mu_minuss");
  	nest::initialize_property(d, "Wmaxs");
  	nest::initialize_property(d, "Kpluss");
  }

  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void stdp_dop_connection::append_properties(DictionaryDatum & d) const
  {
  	nest::ConnectionHetWD::append_properties(d);

  	nest::append_property<nest::double_t>(d, "tau_pluss", tau_plus_);
  	nest::append_property<nest::double_t>(d, "lambdas", lambda_);
  	nest::append_property<nest::double_t>(d, "alphas", alpha_);
  	nest::append_property<nest::double_t>(d, "mu_pluss", mu_plus_);
  	nest::append_property<nest::double_t>(d, "mu_minuss", mu_minus_);
  	nest::append_property<nest::double_t>(d, "Wmaxs", Wmax_);
  	nest::append_property<nest::double_t>(d, "Kpluss", Kplus_);
  }

} // of namespace mynest
