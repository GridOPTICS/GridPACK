// -------------------------------------------------------------
/**
 * @file   petsc_configuration.cpp
 * @author William A. Perkins
 * @date   2013-10-08 08:11:18 d3g096
 * 
 * @brief Implementation of routines for handling PETSc options
 * through Configuration
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

#include "petsc_configuration.hpp"
#include "petsc_exception.hpp"

namespace gridpack {
namespace math {

/// The PETSc option prefix key
static const std::string petscPrefixKey("PETScPrefix");

/// The PETSc options string key
static const std::string petscOptionsKey("PETScOptions");

// -------------------------------------------------------------
// generatePrefix
// -------------------------------------------------------------
static std::string
generatePrefix(void)
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
generatePrefix(const parallel::Communicator& comm)
{
  std::string result;
  if (comm.rank() == 0) {
    result = generatePrefix();
  }
  broadcast(comm, result, 0);
  return result;
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
// petscProcessOptions
// -------------------------------------------------------------
/** 
 * 
 * 
 * @param comm 
 * @param props 
 * 
 * @return PETSc options prefix to use
 */
std::string
petscProcessOptions(const parallel::Communicator& comm,
                    utility::Configuration::Cursor *props)
{
  std::string prefix("");

  if (props == NULL) return prefix;

  prefix = props->get(petscPrefixKey, generatePrefix(comm));
  if (*prefix.rbegin() != '_') {
    prefix.append("_");
  }

  std::string optsorig, optsmod, optsfmt;
  optsorig = props->get(petscOptionsKey, "");

  boost::char_separator<char> sep(" \t\f\n\r\v", "");
  boost::tokenizer<boost::char_separator<char> > 
    opttok(optsorig, sep);
  boost::tokenizer<boost::char_separator<char> >::iterator o;
  for (o = opttok.begin(); o != opttok.end(); ++o) {
    optsfmt.append(*o);
    optsfmt.append(" ");
    optsmod.append(prefixOptionMaybe(prefix, *o));
    optsmod.append(" ");
  }

  if (false) {
    std::cout << "petscProcessOptions:  in: " << optsorig << std::endl;
    std::cout << "petscProcessOptions: fmt: " << optsfmt << std::endl;
    std::cout << "petscProcessOptions: out: " << optsmod << std::endl;
  }

  PetscErrorCode ierr(0);
  try {
    ierr = PetscOptionsInsertString(optsmod.c_str()); CHKERRXX(ierr);
  } catch (const PETSc::Exception& e) {
    throw PETScException(ierr, e);
  }
  
  return prefix;
}

} // namespace math
} // namespace gridpack
