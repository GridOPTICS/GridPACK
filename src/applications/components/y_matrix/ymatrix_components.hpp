/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ymatrix_components.hpp
 * @author Bruce Palmer
 * @date   2016-07-14 13:27:00 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#define USE_ACOPF

#ifndef _ymatrix_components_h_
#define _ymatrix_components_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/utilities/complex.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/component/data_collection.hpp"

namespace gridpack {
namespace ymatrix {

enum YMatrixMode{YBus};

class YMBus
  : public gridpack::component::BaseBusComponent {
  public:
    /**
     *  Simple constructor
     */
    YMBus(void);

    /**
     *  Simple destructor
     */
    ~YMBus(void);

    /**
     * Return size of matrix block on the diagonal contributed by
     * component
     * @param isize, jsize: number of rows and columns of matrix block
     * @return: false if network component does not contribute matrix
     *        element
     */
    bool matrixDiagSize(int *isize, int *jsize) const;

    /**
     * Return the values of for a diagonal matrix block. The values are
     * returned in row-major order
     * @param values: pointer to matrix block values
     * @return: false if network component does not contribute
     *        matrix element
     */
    bool matrixDiagValues(ComplexType *values);

    /**
     * Set values of YBus matrix. These can then be used in subsequent
     * calculations
     */
    void setYBus(void);

    /**
     * Modify diagonal values of matrix.
     * @param rval real part of diagonal matrix element
     * @param ival imaginary part of diagonal matrix element
     */
    void setYBusDiag(double rval, double ival);

    /**
     * Get values of YBus matrix. These can then be used in subsequent
     * calculations
     */
    gridpack::ComplexType getYBus(void);

    /**
     * Load values stored in DataCollection object into YMBus object. The
     * DataCollection object will have been filled when the network was created
     * from an external configuration file
     * @param data: DataCollection object contain parameters relevant to this
     *       bus that were read in when network was initialized
     */
    void load(const boost::shared_ptr<gridpack::component::DataCollection> &data);

    /**
     * Set the mode to control what matrices and vectors are built when using
     * the mapper
     * @param mode: enumerated constant for different modes
     */
    void setMode(int mode);

    /**
     * Return whether or not a bus is isolated
     * @return true if bus is isolated
     */
    bool isIsolated(void) const;

    /**
     * Change isolated status of bus
     * @param flag true if bus is isolated
     */
    void setIsolated(const bool flag);

    /**
     * Get shunt values
     * @param gl shunt GL value
     * @param bl shunt BL value
     */
    void getShuntValues(double *bl, double *gl) const;

    /**
     * Set internal parameters inside the Y-bus component
     * @param name character string describing component to be modified
     * @param value of parameter to be modified
     * @param idx index (if necessary) of variable to be modified
     */
    void setParam(std::string name, double value, int idx);

  private:
    double p_shunt_gs;
    double p_shunt_bs;
    bool p_shunt;
    int p_mode;
    bool p_isolated;

    // p_v and p_a are initialized to p_voltage and p_angle respectively,
    // but may be subject to change during the NR iterations
    double p_ybusr, p_ybusi;

private:


  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<gridpack::component::BaseBusComponent>(*this)
      & p_shunt_gs
      & p_shunt_bs
      & p_shunt
      & p_mode
      & p_isolated
      & p_ybusr & p_ybusi;
  }  

};

class YMBranch
  : public gridpack::component::BaseBranchComponent {
  public:
    /**
     *  Simple constructor
     */
    YMBranch(void);

    /**
     *  Simple destructor
     */
    ~YMBranch(void);

    /**
     * Return size of off-diagonal matrix block contributed by the component
     * for the forward/reverse directions
     * @param isize, jsize: number of rows and columns of matrix block
     * @return: false if network component does not contribute matrix element
     */
    bool matrixForwardSize(int *isize, int *jsize) const;
    bool matrixReverseSize(int *isize, int *jsize) const;

    /**
     * Return the values of the forward/reverse matrix block. The values are
     * returned in row-major order
     * @param values: pointer to matrix block values
     * @return: false if network component does not contribute matrix element
     */
    bool matrixForwardValues(ComplexType *values);
    bool matrixReverseValues(ComplexType *values);

    /**
     * Set values of YBus matrix. These can then be used in subsequent
     * calculations
     */
    void setYBus(void);

    /**
     * Get values of YBus matrix. These can then be used in subsequent
     * calculations
     */
    gridpack::ComplexType getForwardYBus(void);
    gridpack::ComplexType getReverseYBus(void);

    /**
     * Load values stored in DataCollection object into YMBranch object. The
     * DataCollection object will have been filled when the network was created
     * from an external configuration file
     * @param data: DataCollection object contain parameters relevant to this
     *       branch that were read in when network was initialized
     */
    void load(const boost::shared_ptr<gridpack::component::DataCollection> &data);

    /**
     * Return the complex admittance of the branch
     * @return: complex addmittance of branch
     */
    gridpack::ComplexType getAdmittance(void);

    /**
     * Return transformer contribution from the branch to the calling
     * bus
     * @param bus: pointer to the bus making the call
     * @return: contribution from transformers to Y matrix
     */
    gridpack::ComplexType getTransformer(YMBus *bus);

    /**
     * Return the contribution to a bus from shunts
     * @param bus: pointer to the bus making the call
     * @return: contribution to Y matrix from shunts associated with branches
     */
    gridpack::ComplexType getShunt(YMBus *bus);

    /**
     * Return contributions to Y-matrix from a specific transmission element
     * @param tag character string for transmission element
     * @param Yii contribution from "from" bus
     * @param Yij contribution from line element
     */
    void getLineElements(const std::string tag,
        gridpack::ComplexType *Yii, gridpack::ComplexType *Yij);

    /**
     * Return contributions to Y-matrix from a specific transmission element at to end
     * @param tag character string for transmission element
     * @param Yii contribution from "from" bus
     * @param Yij contribution from line element
     */
    void getRvrsLineElements(const std::string tag,
        gridpack::ComplexType *Yii, gridpack::ComplexType *Yij);

    /**
     * Return status of all transmission elements
     * @return vector containing status of transmission elements
     */
    std::vector<bool> getLineStatus();

    /**
     * Return tags of all transmission elements
     * @return vector containging tag of transmission elements
     */
    std::vector<std::string> getLineTags();

    /**
     * Set the status of a transmission element based on its tag name
     * @param tag name of transmission element
     * @param status that transmission element should be set to
     * @return false if no transmission element with that name exists
     */
    bool setLineStatus(std::string tag, bool status);

    /**
     * Set the mode to control what matrices and vectors are built when using
     * the mapper
     * @param mode: enumerated constant for different modes
     */
    void setMode(int mode);

    /**
     * Get branch susceptance (charging)
     * @param tag string identifier for transmission element
     * @return value of susceptance
     */
    double getSusceptance(std::string tag);

    /**
     * Set internal parameters inside the Y-branch component
     * @param name character string describing component to be modified
     * @param value of parameter to be modified
     * @param idx index (if necessary) of variable to be modified
     */
    void setParam(std::string name, double value, int idx);

#ifdef USE_ACOPF
    /**
     * Return components from individual transmission elements
     * @param yffr list of real parts of Yff
     * @param yffr list of imaginary parts of Yff
     * @param yttr list of real parts of Ytt
     * @param yttr list of imaginary parts of Ytt
     * @param yftr list of real parts of Yft
     * @param yftr list of imaginary parts of Yft
     * @param ytfr list of real parts of Ytf
     * @param ytfr list of imaginary parts of Ytf
     * @param switched flag on whether line is switched or not
     */
    void getYElements(std::vector<double> &yffr, std::vector<double> &yffi,
                      std::vector<double> &yttr, std::vector<double> &ytti,
                      std::vector<double> &yftr, std::vector<double> &yfti,
                      std::vector<double> &ytfr, std::vector<double> &ytfi,
                      std::vector<bool> &switched);
#endif

  private:
    std::vector<double> p_reactance;
    std::vector<double> p_resistance;
    std::vector<double> p_tap_ratio;
    std::vector<double> p_phase_shift;
    std::vector<double> p_charging;
    std::vector<double> p_shunt_admt_g1;
    std::vector<double> p_shunt_admt_b1;
    std::vector<double> p_shunt_admt_g2;
    std::vector<double> p_shunt_admt_b2;
    std::vector<bool> p_xform, p_shunt;
    int p_mode;
    double p_ybusr_frwd, p_ybusi_frwd;
    double p_ybusr_rvrs, p_ybusi_rvrs;
    std::vector<bool> p_branch_status;
    std::vector<bool> p_switched;
    std::vector<std::string> p_tag;
    int p_elems;
    bool p_isolated;
    bool p_active;
#ifdef USE_ACOPF
    std::vector<double> p_yffr, p_yffi;
    std::vector<double> p_yttr, p_ytti;
    std::vector<double> p_yftr, p_yfti;
    std::vector<double> p_ytfr, p_ytfi;
#endif

private:


  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<gridpack::component::BaseBranchComponent>(*this)
#ifdef USE_ACOPF
      & p_yffr & p_yffi
      & p_yttr & p_ytti
      & p_yftr & p_yfti
      & p_ytfr & p_ytfi
#endif
      & p_reactance
      & p_resistance
      & p_tap_ratio
      & p_phase_shift
      & p_charging
      & p_shunt_admt_g1
      & p_shunt_admt_b1
      & p_shunt_admt_g2
      & p_shunt_admt_b2
      & p_xform & p_shunt
      & p_mode
      & p_ybusr_frwd & p_ybusi_frwd
      & p_ybusr_rvrs & p_ybusi_rvrs
      & p_branch_status
      & p_switched
      & p_tag
      & p_elems
      & p_isolated
      & p_active;
  }  

};


}     // ymatrix
}     // gridpack

BOOST_CLASS_EXPORT_KEY(gridpack::ymatrix::YMBus)
BOOST_CLASS_EXPORT_KEY(gridpack::ymatrix::YMBranch)


#endif
