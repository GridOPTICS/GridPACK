// -------------------------------------------------------------
/**
 * @file   implementation_visitor.cpp
 * @author William A. Perkins
 * @date   2013-06-11 14:16:24 d3g096
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------


#include <boost/assert.hpp>
#include "implementation_visitor.hpp"

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class ImplementationVisitor
// -------------------------------------------------------------

// -------------------------------------------------------------
// ImplementationVisitor:: constructors / destructor
// -------------------------------------------------------------
ImplementationVisitor::ImplementationVisitor()
  : utility::Uncopyable()
{
}

ImplementationVisitor::~ImplementationVisitor(void)
{
}

// -------------------------------------------------------------
// ImplementationVisitor::visit
// -------------------------------------------------------------
/** 
 * 
 * 
 */
void 
ImplementationVisitor::visit(VectorImplementation&)
{
  BOOST_ASSERT(false);
}

void 
ImplementationVisitor::visit(PETScVectorImplementation&)
{
  BOOST_ASSERT(false);
}

void 
ImplementationVisitor::visit(MatrixImplementation&)
{
  BOOST_ASSERT(false);
}

void 
ImplementationVisitor::visit(PETScMatrixImplementation&)
{
  BOOST_ASSERT(false);
}

// -------------------------------------------------------------
//  class ConstImplementationVisitor
// -------------------------------------------------------------

// -------------------------------------------------------------
// ConstImplementationVisitor:: constructors / destructor
// -------------------------------------------------------------
ConstImplementationVisitor::ConstImplementationVisitor()
  : utility::Uncopyable()
{
}

ConstImplementationVisitor::~ConstImplementationVisitor(void)
{
}

// -------------------------------------------------------------
// ConstImplementationVisitor::visit
// -------------------------------------------------------------
void 
ConstImplementationVisitor::visit(const VectorImplementation&)
{
  BOOST_ASSERT(false);
}

void 
ConstImplementationVisitor::visit(const PETScVectorImplementation&)
{
  BOOST_ASSERT(false);
}

void 
ConstImplementationVisitor::visit(const MatrixImplementation&)
{
  BOOST_ASSERT(false);
}

void 
ConstImplementationVisitor::visit(const PETScMatrixImplementation&)
{
  BOOST_ASSERT(false);
}


} // namespace math
} // namespace gridpack

