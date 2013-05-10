// -------------------------------------------------------------
/**
 * @file   base_mapper.hpp
 * @author Bruce Palmer
 * @date   May 3, 2013
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _base_mapper_h_
#define _base_mapper_h_

// -------------------------------------------------------------
//  class BaseMapper:
//  This class is the base class for creating matrices and
//  vectors from a collection of network components. It only
//  requires that each of the components implements the
//  MatVecInterface
// -------------------------------------------------------------

#include "smart_ptr.hpp"
class BaseMapper {
public:

  /**
   * Constructor
   */
  BaseMapper(void);

  /**
   * Destructor
   */
  ~BaseMapper(void);

  /**
   * Construct a matrix from a single field implementing a MatVecInterface
   * @param field: field of components implementing the MatVecInterface
   */
  createMatrix(stlplus::smart_ptr<BaseField<stlplus::smart_ptr<MatVecInterface> > > field);

  /**
   * Construct a matrix from a two fields implementing a MatVecInterface
   * @param field1: first field of components implementing the MatVecInterface
   * @param field2: second field of components implementing the MatVecInterface
   */
  createMatrix(stlplus::smart_ptr<BaseField<stlplus::smart_ptr<MatVecInterface> > > field1,
               stlplus::smart_ptr<BaseField<stlplus::smart_ptr<MatVecInterface> > > field2);

  /**
   * Construct a vector from a single field implementing a MatVecInterface
   * @param field: field of components implementing the MatVecInterface
   */
  createVector(stlplus::smart_ptr<BaseField<stlplus::smart_ptr<MatVecInterface> > > field);
};

#endif
