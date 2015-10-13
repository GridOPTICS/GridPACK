// -------------------------------------------------------------
// file: glpk_optimizer_implementation.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created September 16, 2015 by William A. Perkins
// Last Change: 2015-10-13 15:34:10 d3g096
// -------------------------------------------------------------


#include <algorithm>
#include <fstream>
#include <glpk.h>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include "glpk_optimizer_implementation.hpp"


namespace gridpack {
namespace optimization {

// -------------------------------------------------------------
//  class GLPKOptimizerImplementation
// -------------------------------------------------------------

// -------------------------------------------------------------
// GLPKOptimizerImplementation:: constructors / destructor
// -------------------------------------------------------------
GLPKOptimizerImplementation::GLPKOptimizerImplementation(const parallel::Communicator& comm)
  : LPFileOptimizerImplementation(comm)
{
}

GLPKOptimizerImplementation::~GLPKOptimizerImplementation(void)
{
}

// -------------------------------------------------------------
// GLPKOptimizerImplementation::p_solve
// -------------------------------------------------------------
void
GLPKOptimizerImplementation::p_solve(const p_optimizeMethod& m)
{
  parallel::Communicator comm(this->communicator());
  int nproc(comm.size());
  int me(comm.rank());
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
  
  int ierr;
  glp_prob *lp = glp_create_prob();
  std::cout << tmpname << std::endl;
  ierr = glp_read_lp(lp, NULL, tmpname.c_str());
  if (ierr != 0) {
    std::string msg = 
      boost::str(boost::format("GLPK LP parse failure, code = %d") % ierr);
    glp_delete_prob(lp);
    throw gridpack::Exception(msg);
  }

  ierr = glp_simplex(lp, NULL);

  if (ierr != 0) {
    std::string msg = 
      boost::str(boost::format("GLPK optimizer failure, code = %d") % ierr);
    glp_delete_prob(lp);
    throw gridpack::Exception(msg);
  }

  comm.barrier();
  for (int p = 0; p < nproc; ++p) {
    if (p == me) {
      std::cout << "Optimimal variable values (process " << me << "):" << std::endl;

      int varnum(glp_get_num_cols(lp));
      for (int idx = 1; idx <= varnum; ++idx) {
        std::string gname(glp_get_col_name(lp, idx));
        VariablePtr v(p_allVariables[gname]);
        std::string vname(v->name());
        
        std::cout << gname << " " << vname << " " 
                  << glp_get_col_prim(lp, idx) << " " 
                  << glp_get_col_dual(lp, idx) << " "
                  << std::endl;
        
        SetVariableInitial vset(glp_get_col_prim(lp, idx));
        v->accept(vset);
      } 
      VariableTable vtab(std::cout);
      BOOST_FOREACH(VarMap::value_type& i, p_allVariables) {
        i.second->accept(vtab);
      }
    }
    comm.barrier();
  }
  glp_delete_prob(lp);
}

} // namespace optimization
} // namespace gridpack
