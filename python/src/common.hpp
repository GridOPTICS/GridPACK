// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 * Copyright (c) 2013 Battelle Memorial Institute
 * Licensed under modified BSD License. A copy of this license can be found
 * in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created February 10, 2025 by Perkins
// -------------------------------------------------------------

#ifndef _common_hpp_
#define _common_hpp_

#include <mpi4py/mpi4py.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <boost/shared_ptr.hpp>

#include <gridpack/environment/environment.hpp>
#include <gridpack/configuration/configuration.hpp>
#include <gridpack/parallel/communicator.hpp>

namespace gp = gridpack;
namespace gpp = gridpack::parallel;
namespace gpu = gridpack::utility;


/// A functor to keep smart pointers from deleting their pointer
struct null_deleter
{
  void operator()(void const *) const { }
};


// GridPACK uses Boost smart pointers, so let's use those here
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>, false);

// -------------------------------------------------------------
// gridpack::utility::Configuration wrapping
//
// pybind11 gets really confused with the Configuration and
// Configuration::Cursor being the same type.  
// -------------------------------------------------------------

// -------------------------------------------------------------
//  class ConfigurationCursorWrapper
// -------------------------------------------------------------
struct ConfigurationCursorWrapper {
  gpu::Configuration::CursorPtr the_cursor;
};

#endif
