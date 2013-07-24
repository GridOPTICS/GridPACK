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

#ifndef matrix_interface_hpp_
#define matrix_interface_hpp_

namespace gridpack {
namespace network {


// -------------------------------------------------------------
//  class ComponenVisitor
// -------------------------------------------------------------
/**
 * 
 * To be safe, these should be used simultaneously on all processes.  
 */

class MatrixInterface {
public:

  /// Default constructor.
    MatrixInterface(int ni, int nj, int mi, int mj) :
        ni_(ni),
        nj_(nj),
        mi_(mi),
        mj_(mj_),
        size_((nj - ni)*(mj - ni)),
        x_(new double[size_])
    {}

    MatrixInterface() :
        ni_(0),
        nj_(0),
        mi_(0),
        mj_(0),
        size_(0),
        x_(NULL)
    {}

    /// Destructor
  virtual ~MatrixInterface(void);

  /// Allow visits by component count visitors
  virtual void accept(BusCountVisitor& visitor){accept_(visitor);}
  virtual void accept(BranchCountVisitor& visitor){accept_(visitor);}
  virtual void accept(ComponentCountVisitor& visitor){accept_(visitor);}

  virtual void accept(BusSizeVisitor& visitor){accept_(visitor);}
  virtual void accept(BranchSizeVisitor& visitor){accept_(visitor);}
  virtual void accept(ComponentSizeVisitor& visitor){accept_(visitor);}

  virtual void accept(BusDataVisitor& visitor){accept_(visitor);}
  virtual void accept(BranchDataVisitor& visitor){accept_(visitor);}
  virtual void accept(ComponentDataVisitor& visitor){accept_(visitor);}

protected:
  // the matrix interface has to deal with three types of visitors
  //     the number of each component type (count)
  //     the size of each component's contribution to matrix/vector (size)
  //     a copy of each component's data contribution (data)
  virtual void accept_(BusCountVisitor & visitor){};
  virtual void accept_(ComponentSizeVisitor & visitor){};
  virtual void accept_(ComponentCountVisitor & visitor){};

  virtual void accept_(BusSizeVisitor & visitor){};
  virtual void accept_(BranchSizeVisitor & visitor){};
  virtual void accept_(ComponentSizeVisitor & visitor){};

  virtual void accept_(BusDataVisitor & visitor){};
  virtual void accept_(BranchDataVisitor & visitor){};
  virtual void accept_(ComponentDataVisitor & visitor){};

  virtual void setSize(ComponentSizeVisitor & visitor){
      visitor->setSize(size);
  }
private:
  int                ni_;
  int                nj_;
  int                mi_;
  int                mj_
  int                size_;
  double           * x_;

};  /* END MATRIX_INTERFACE */

} // namespace math
} // namespace gridpack

#endif
