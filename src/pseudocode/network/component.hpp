// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   component.hpp
 * @author William A. Perkins
 * @date   Fri Mar 22 12:04:21 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 22, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#ifndef _component_hpp_
#define _component_hpp_

#include "gridpack/utility/indexed.h"

namespace gridpack {
namespace network {

// -------------------------------------------------------------
//  class GeneratorInterface
// -------------------------------------------------------------
class GeneratorInterface : public utility::UnCopyable {
public:

  /// Default constructor.
  GeneratorInterface(void);

  /// Destructor
  ~GeneratorInterface(void);

};



// -------------------------------------------------------------
//  class BusInterface
// -------------------------------------------------------------
/**
 * 
 * 
 */
class BusInterface : public utility::Indexed, public utility::UnCopyable {
public:

  /// Default constructor.
  BusInterface(void);

  /// Destructor
  virtual ~BusInterface(void);

  /// Way to get the size of admittance matrix contribution
  /** 
   * This returns the size (of a square sub matrix) that is the
   * contribution of this bus to the admittance matrix.
   * 
   * 
   * @return size of the square region this bus needs in admittance matrix
   */
  int admittance_size(void) const
  {
    return this->admittance_size_();
  }

  /// Way to get admittance matrix contribution

  void admittance_contribution(double *submatrix) const
  {
    return this->admittance_contribution_(submatrix);
  }

protected:
  
  /// Specialized way to get size of admittance matrix contribution
  virtual int admittance_size_(void) const = 0;

  /// Specialized way to get admittance matrix contribution
  virtual void admittance_contribution_(double *submatrix) const = 0;

private:

  /// Protected copy constructor to avoid unwanted copies.
  BusInterface(const BusInterface& old);
};

// -------------------------------------------------------------
//  class BranchInterface
// -------------------------------------------------------------
class BranchInterface : public utility::Indexed {
public:

  /// Default constructor.
  BranchInterface(void);

  /// Destructor
  virtual ~BranchInterface(void);

private:

  /// Protected copy constructor to avoid unwanted copies.
  BranchInterface(const BranchInterface& old);

};





} // network namespace
} // gridpack namespace
#endif
