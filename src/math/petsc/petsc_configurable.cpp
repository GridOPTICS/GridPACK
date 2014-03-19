// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
/**
 * @file   petsc_configurable.cpp
 * @author William A. Perkins
 * @date   2014-03-19 09:50:01 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <cctype>
#include <iostream>
#include <boost/tokenizer.hpp>
#include <boost/random/random_device.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/mpi/collectives.hpp>
#include <petscsys.h>

#include "petsc_exception.hpp"
#include "petsc_configurable.hpp"

namespace gridpack {
namespace math {

static const bool verbose(false);

// -------------------------------------------------------------
//  class PETScConfigurable
// -------------------------------------------------------------

// -------------------------------------------------------------
// PETScConfigurable static members
// -------------------------------------------------------------

/// The PETSc option prefix key
const std::string PETScConfigurable::p_prefixKey("PETScPrefix");

/// The PETSc options string key
const std::string PETScConfigurable::p_optionsKey("PETScOptions");

// -------------------------------------------------------------
// PETScConfigurable::p_generatePrefix
// -------------------------------------------------------------
std::string
PETScConfigurable::p_generatePrefix(void)
{
  static const int prefix_len(4);
  static std::string chars("abcdefghijklmnopqrstuvwxyz");
  boost::random::random_device rng;
  boost::random::uniform_int_distribution<> index_dist(0, chars.size() - 1);
  std::string result("");
  for(int i = 0; i < prefix_len; ++i) {
    result += chars[index_dist(rng)];
  }
  return result;
}


std::string
PETScConfigurable::p_generatePrefix(const parallel::Communicator& comm)
{
  std::string result;
  if (comm.rank() == 0) {
    result = p_generatePrefix();
  }
  boost::mpi::broadcast(comm, result, 0);
  return result;
}

// -------------------------------------------------------------
// PETScConfigurable constructor / destructor
// -------------------------------------------------------------
PETScConfigurable::PETScConfigurable(const parallel::Communicator& comm)
  : p_comm(comm)
{}

PETScConfigurable::~PETScConfigurable(void)
{
  try {
    this->p_unprocessOptions();
  } catch (...) {
    // just eat it
  }
}

// -------------------------------------------------------------
// PETScConfigurable::build
// -------------------------------------------------------------
void
PETScConfigurable::build(utility::Configuration::CursorPtr props)
{
  this->p_processOptions(props);
  this->p_build(p_prefix);
}

// -------------------------------------------------------------
// prefixOptionMaybe
// -------------------------------------------------------------
/// Prepend the prefix to the specified option, if it looks like an option
/** 
 * 
 * 
 * @param prefix option prefix
 * @param opt original option string
 * 
 * @return modified option string
 */
static std::string 
prefixOptionMaybe(const std::string& prefix, const std::string& opt)
{
  std::string result(opt);

  if (opt[0] == '-' && islower(opt[1])) {
    result = "-";
    result.append(prefix);
    result.append(opt.substr(1));
  }      

  return result;
}

// -------------------------------------------------------------
// PETScConfigurable::p_processOptions
// -------------------------------------------------------------
/** 
 * 
 * 
 * @param comm 
 * @param props 
 * 
 * @return PETSc options prefix to use
 */
void
PETScConfigurable::p_processOptions(utility::Configuration::CursorPtr props)
{
  if (!props) return;

  p_prefix = props->get(p_prefixKey, p_generatePrefix(p_comm));
  if (*p_prefix.rbegin() != '_') {
    p_prefix.append("_");
  }

  std::string optsorig, optsmod, optsfmt;
  optsorig = props->get(p_optionsKey, "");

  boost::char_separator<char> sep(" \t\f\n\r\v", "");
  boost::tokenizer<boost::char_separator<char> > 
    opttok(optsorig, sep);
  boost::tokenizer<boost::char_separator<char> >::iterator o;
  for (o = opttok.begin(); o != opttok.end(); ++o) {
    optsfmt.append(*o);
    optsfmt.append(" ");
    optsmod.append(prefixOptionMaybe(p_prefix, *o));
    optsmod.append(" ");
  }

  if (verbose) {
    std::cout << "p_processOptions:  in: " << optsorig << std::endl;
    std::cout << "p_processOptions: fmt: " << optsfmt << std::endl;
    std::cout << "p_processOptions: out: " << optsmod << std::endl;
  }
  p_loadedOptions = optsmod;

  PetscErrorCode ierr(0);
  try {
    ierr = PetscOptionsInsertString(p_loadedOptions.c_str()); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
  return;
}

// -------------------------------------------------------------
// PETScConfigurable::p_unprocessOptions
// -------------------------------------------------------------
void
PETScConfigurable::p_unprocessOptions(void)
{
  boost::char_separator<char> sep(" \t\f\n\r\v", "");
  boost::tokenizer<boost::char_separator<char> > 
    opttok(p_loadedOptions, sep);
  boost::tokenizer<boost::char_separator<char> >::iterator o;
  for (o = opttok.begin(); o != opttok.end(); ++o) {
    if ((*o)[0] == '-' && islower((*o)[1])) {
      PetscErrorCode ierr(0);
      try {
        if (verbose) {
          std::cout << "p_unprocessOptions: removing \"" 
                    << *o << "\"" << std::endl;
        }
        ierr = PetscOptionsClearValue(o->c_str());
      } catch (const PETSc::Exception& e) {
        throw PETScException(ierr, e);
      }
    }
  } 
  p_loadedOptions.clear();
}



} // namespace math
} // namespace gridpack
