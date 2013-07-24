// -------------------------------------------------------------
/**
 * @file   component_visitor.hpp
 * @author Kevin A. Glass
 * @date   Fri Apr  19 13:36:28 2013
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Created April  19, 2013 by Kevin A. Glass
// Last Change:
// -------------------------------------------------------------

#ifndef _bus_visitor_hpp_
#define _bus_visitor_hpp_

namespace gridpack {
namespace network {

class MatrixInterface;


// -------------------------------------------------------------
//  class ComponenVisitor
// -------------------------------------------------------------
/**
 * 
 * To be safe, these should be used simultaneously on all processes.  
 */

class BusVisitor {
public:

  /// Default constructor.
  BusVisitor() : matrixImpl_(NULL), nRows(0), nCols(0) {};
  /// Destructor
  virtual ~BusVisitor(void){};

  virtual void math::MatrixImplementation();

  virtual void mapData(MatrixInterface & interface) {

  }
protected:
  virtual void getMatrixSize(network::ComponentNetwork & network) {
      BusCountVisitor        visitor;

  }
private:
  math::MatrixImplementation      *  matrixImpl_;
  int                                nRows;
  int                                nCols;
};

} // namespace math
} // namespace gridpack

#endif
