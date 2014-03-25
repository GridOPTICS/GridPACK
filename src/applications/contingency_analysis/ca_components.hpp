/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ca_components.hpp
 * @author Bruce Palmer
 * @date   2014-03-05 14:50:50 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _ca_components_h_
#define _ca_components_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/utilities/complex.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/component/data_collection.hpp"
#include "gridpack/applications/components/y_matrix/ymatrix_components.hpp"
#include "gridpack/applications/components/pf_matrix/pf_components.hpp"

namespace gridpack {
namespace  contingency_analysis {

class CABus
  : public gridpack::powerflow::PFBus
{
  public:
    /**
     *  Simple constructor
     */
    CABus(void);

    /**
     *  Simple destructor
     */
    ~CABus(void);

    /**
     * Set voltage limits on bus. Voltage magnitudes that are outside limits
     * represent a contingency violation
     * @param vMin minimum voltage
     * @param vMax maximum voltage
     */
    void setVoltageLimits(double vMin, double vMax);

    /**
     * Write output from buses to standard out
     * @param string (output) string with information to be printed out
     * @param bufsize size of string buffer in bytes
     * @param signal an optional character string to signal to this
     * routine what about kind of information to write
     * @return true if bus is contributing string to output, false otherwise
     */
    bool serialWrite(char *string, const int bufsize, const char *signal = NULL);

    /**
     * Return the size of the buffer used in data exchanges on the network.
     * For this problem, the power flow data need to be exchanged plus data
     * for bool that keeps track of the isolated status of the bus
     * @return size of buffer
     */
    int getXCBufSize(void);

    /**
     * Assign pointers for powerflow data and isolated status variable
     */
    void setXCBuf(void *buf);

    /**
     * Voltage min and max limits
     */
    double p_vMin;
    double p_vMax;
    bool *p_isolated;

private:


  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar  & boost::serialization::base_object<gridpack::powerflow::PFBus>(*this)
      & p_vMin
      & p_vMax;
  }  

};

class CABranch
  : public gridpack::powerflow::PFBranch {
  public:
    /**
     *  Simple constructor
     */
    CABranch(void);

    /**
     *  Simple destructor
     */
    ~CABranch(void);

private:


  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar  & boost::serialization::base_object<gridpack::powerflow::PFBranch>(*this);
  }  

};


}     // contingency_analysis
}     // gridpack

BOOST_CLASS_EXPORT_KEY(gridpack::contingency_analysis::CABus);
BOOST_CLASS_EXPORT_KEY(gridpack::contingency_analysis::CABranch);


#endif
