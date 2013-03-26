// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: petsc_vector_implementation.h
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 26, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

// SCCS ID: $Id$ Battelle PNL

#ifndef _petsc_vector_implementation_h_
#define _petsc_vector_implementation_h_

#include "vector_implementation.h"


// -------------------------------------------------------------
//  class PETScVectorImplementation
// -------------------------------------------------------------
class PETScVectorImplementation 
  : public VectorImplementation {
public:

  /// Default constructor.
  PETScVectorImplementation(void);

  /// Destructor
  ~PETScVectorImplementation(void);

protected:

  /// Set an individual element (specialized)
  void set_element_(const int& i, const double& x);

  /// Set an several elements (specialized)
  void set_elements_(cont int& n, const int *i, const double *x);

  /// Add to an individual element (specialized)
  void add_element_(const int& i, const double& x);

  /// Add to an several elements (specialized)
  void add_elements_(const int& n, const int *i, const double *x);

  /// Make all the elements zero (specialized)
  void zero_(void);

};


#endif
