// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: petsc_types.hpp
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created February 18, 2015 by William A. Perkins
// Last Change: 2015-02-18 08:03:29 d3g096
// -------------------------------------------------------------


#ifndef _petsc_types_hpp_
#define _petsc_types_hpp_


#include <petsc.h>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <gridpack/utilities/complex.hpp>

namespace gridpack {
namespace math {

/// A flag (type) to denote whether the library can be used
/**
 * Some operations can be passed directly to the underlying library
 * if the TheType is the same as the PETSc type @e or the PETSc type
 * is complex.  This type computes and stores that flag. 
 * 
 */
template <typename T>
struct UsePetscLibrary
  : public boost::mpl::bool_<
            boost::is_same<T, PetscScalar>::value ||
            boost::is_same<ComplexType, PetscScalar>::value >::type
{
};

/// The number of library elements used to represent a single vector element
template <typename T>
struct PetscElementSize
  : public boost::mpl::if_< UsePetscLibrary<T>, 
                            boost::mpl::int_<1>, 
                            boost::mpl::int_<2> >::type
{
};


} // namespace math
} // namespace gridpack

#endif
