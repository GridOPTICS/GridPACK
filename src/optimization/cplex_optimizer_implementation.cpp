// -------------------------------------------------------------
// file: cplex_optimizer_implementation.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created September 28, 2015 by Yilin Fang
// Last Change: 2015-11-06 07:19:46 d3g096
// -------------------------------------------------------------


#include <algorithm>
#include <fstream>
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <cstring>
#include <iostream>
#include <string>
#include <ilcplex/ilocplex.h>
#include "cplex_optimizer_implementation.hpp"


namespace gridpack {
namespace optimization {

// -------------------------------------------------------------
//  class CPlexOptimizerImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// CPlexOptimizerImplementation:: constructors / destructor
// -------------------------------------------------------------
CPlexOptimizerImplementation::CPlexOptimizerImplementation(const parallel::Communicator& comm)
  : LPFileOptimizerImplementation(comm)
{
}

CPlexOptimizerImplementation::~CPlexOptimizerImplementation(void)
{
}

// -------------------------------------------------------------
// CPlexOptimizerImplementation::p_solve
// -------------------------------------------------------------
void
CPlexOptimizerImplementation::p_solve(const p_optimizeMethod& m)
{
  std::string tmpname(p_temporaryFileName());
  std::ofstream tmp;

  tmp.open(tmpname.c_str());
  if (!tmp) {
    std::string msg("Cannot open temporary file: ");
    msg += tmpname.c_str();
    throw gridpack::Exception(msg);
  }
  p_write(m, tmp);
  tmp.close();
  
  std::cout << tmpname << std::endl;
  IloEnv env;
  IloModel model(env);
  IloCplex cplex(env);
  IloObjective obj;
  IloNumVarArray var(env);
  IloRangeArray rng(env);
  cplex.importModel(model, tmpname.c_str(), obj,var,rng);
  cplex.extract(model);
  if ( !cplex.solve() ) {
     env.error() << "Failed to optimize LP" << std::endl;
     throw(-1);
  }

  IloNumArray vals(env);
  cplex.getValues(vals,var);
//  env.out() << "solution vector = " << vals << std::endl;

  IloInt n(vals.getSize());
  unsigned sz; 
  for (IloInt i = 0; i < n; ++i) {
    std::string vname(var[i].getName());
    std::string buff;
    buff = vname;
    sz = buff.size();
    buff.resize(sz+20-sz,' ');
    VariablePtr v(p_allVariables[vname]);
    SetVariableInitial vset(vals[i]);
    v->accept(vset);
    std::cout << buff << " = " << vals[i] << std::endl;
  }


}

} // namespace optimization
} // namespace gridpack
