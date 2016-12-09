// -------------------------------------------------------------
// file: cplex_optimizer_implementation.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created September 28, 2015 by Yilin Fang
// Last Change: 2016-12-08 14:59:37 d3g096
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
  LPFileOptimizerImplementation::p_solve(m);

  if (p_runMaybe) {
    IloEnv env;
    IloModel model(env);
    IloCplex cplex(env);
    IloObjective obj;
    IloNumVarArray var(env);
    IloRangeArray rng(env);
    cplex.importModel(model, p_outputName.c_str(), obj,var,rng);
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
}

} // namespace optimization
} // namespace gridpack
