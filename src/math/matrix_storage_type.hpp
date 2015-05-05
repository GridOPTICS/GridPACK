// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: matrix_storage_type.hpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created January 28, 2015 by William A. Perkins
// Last Change: 2015-02-09 11:37:36 d3g096
// -------------------------------------------------------------


#ifndef _matrix_storage_type_hpp_
#define _matrix_storage_type_hpp_

namespace gridpack {
namespace math {

/// The types of matrices that can be created
/**
 * The gridpack::math library provides two storage schemes for
 * matrices. This is used by Matrix and MatrixImplementation
 * subclasses.
 *
 * The actual storage scheme and memory used is dependent upon the
 * underlying math library implementation.
 * 
 */
enum MatrixStorageType { 
  Dense,                      /**< dense matrix storage scheme */
  Sparse                      /**< sparse matrix storage scheme */
};

} // namespace math
} // namespace gridpack

#endif
