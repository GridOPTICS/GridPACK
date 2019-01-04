/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   ds_components.hpp
 * @author Shrirang Abhyankar
 * @date   2016-07-14 14:17:34 d3g096
 * 
 * @brief  
 * Network component defitions
 * 
 */
// -------------------------------------------------------------

#ifndef _ds_components_h_
#define _ds_components_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/component/data_collection.hpp"
#include "gridpack/network/base_network.hpp"

namespace gridpack {
namespace dsimplicit {

  enum DSMode{INIT_X,RESIDUAL_EVAL,XVECTOBUS,XDOTVECTOBUS,FAULT_EVAL};

class DSBus
  : public gridpack::component::BaseBusComponent {
  public:
    /**
     *  Simple constructor
     */
  DSBus(void);

    /**
     *  Simple destructor
     */
  ~DSBus(void);

/**
  *  Check if the bus is isolated. Returns true if the bus is isolated

*/
  bool isIsolated(void) const;

    /**
     * Load values stored in DataCollection object into DSBus object. The
     * DataCollection object will have been filled when the network was created
     * from an external configuration file
     * @param data: DataCollection object contain parameters relevant to this
     *       bus that were read in when network was initialized
     */
    void load(const boost::shared_ptr<gridpack::component::DataCollection> &data);

  /**
   * Set the model to control what matrices and vectors and built when using the
   * mapper
   * @param mode: enumerated constant for different modes
   */
  void setMode(int mode);

    /**
     * Return size of matrix block on the diagonal contributed by component
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
     * Return size of vector block contributed by component
     * @param isize: number of vector elements
     * @return: false if network component does not contribute
     *        vector element
     */
    bool vectorSize(int *isize) const;

    /**
     * Return the values of the vector block
     * @param values: pointer to vector values
     * @return: false if network component does not contribute
     *        vector element
     */
    bool vectorValues(ComplexType *values);

    /**
     * Set the internal values of the voltage magnitude and phase angle. Need this
     * function to push values from vectors back onto buses 
     * @param values array containing voltage magnitude and angle
     */
    void setValues(gridpack::ComplexType *values);

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
      * Get voltages in the rectangular form VD, VQ
      * @param double VD - real part of complex voltage at this bus
      * @param double VQ - imaginary part of complex voltage at this bus
      */
  void getVoltagesRectangular(double*,double*) const;

  /**
     Get number of variables for this bus
  */
  void getNvar(int*) const;

  /**
    Set the shift value provided by TS
  */
  void setTSshift(double);

  /**
     Set buffer size for exchange
  */
  int getXCBufSize(void);

  void setXCBuf(void*);

  /*
    Add Bus shunts
  */
  void addBusShunt(double Gs,double Bs)
  {
    p_gl += Gs;
    p_bl += Bs;
  }
  private:
  // Anything declared here should be set in the Archive class in exactly the same order!!
  // Data needed for calculations
  double p_gl,p_bl; // Shunt conductance and susceptance (p.u.)
  double p_pl,p_ql; // Active and reactive load (p.u)
  double p_Vm0,p_Va0;     // Voltage magnitude and angle at t=0 used for building constant impedance load
  int    p_ngen;    // Number of generators incident on this bus
  int    p_nactivegen; // Number of active generators (status=1) on this bus
  bool   p_isolated;   // flag for isolated bus
  int    p_mode; // factory mode
  double p_TSshift;  // shift value provided by TSIJacobian. 
  double p_ws;   // synchronous speed
  std::vector<int> p_gstatus; // Generator status
  std::vector<double> p_pg, p_qg; // Real and reactive generator output
  std::vector<double> p_mbase; // Generator machine base
  std::vector<double> p_Rs;     // Machine stator resistance (converted to MVAbase)
  std::vector<double> p_Xdp;   // Machine transient reactance (converted to MVAbase)
  std::vector<double> p_H;     // Machine Inertia constant (converted to MVAbase)
  std::vector<double> p_D;     // Machine damping coeffient (converted to MVAbase)
  std::vector<double> p_Pm;    // Mechanical power input
  std::vector<double> p_Ep;    // Machine internal emf

  // Variables
  double p_VD,p_VQ; // Real and imaginary part of bus voltage
  std::vector<double> p_delta,p_dw; // Machine angle and speed deviation */

  // xdot from PETSc
  std::vector<double> p_deltadot,p_dwdot;
  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<gridpack::component::BaseBusComponent>(*this)
      & p_gl & p_bl
      & p_pl & p_ql
      & p_Vm0 & p_Va0
      & p_ngen
      & p_nactivegen
      & p_isolated
      & p_mode
      & p_TSshift
      & p_ws
      & p_gstatus
      & p_pg & p_qg
      & p_mbase
      & p_Rs
      & p_Xdp
      & p_H
      & p_D
      & p_Pm
      & p_Ep
      & p_VD & p_VQ
      & p_delta & p_dw
      & p_deltadot & p_dwdot;
  }  

}; // End of DSBus

class DSBranch
  : public gridpack::component::BaseBranchComponent {
  public:
    /**
     *  Simple constructor
     */
    DSBranch(void);

    /**
     *  Simple destructor
     */
    ~DSBranch(void);

    /**
     * Load values stored in DataCollection object into DSBranch object. The
     * DataCollection object will have been filled when the network was created
     * from an external configuration file
     * @param data: DataCollection object contain parameters relevant to this
     *       branch that were read in when network was initialized
     */
    void load(const boost::shared_ptr<gridpack::component::DataCollection> &data);

  /**
   * Set the model to control what matrices and vectors and built when using the
   * mapper
   * @param mode: enumerated constant for different modes
   */
  void setMode(int mode);

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
     * Write output from branches to standard out
     * @param string (output) string with information to be printed out
     * @param bufsize size of string buffer in bytes
     * @param signal an optional character string to signal to this
     * routine what about kind of information to write
     * @return true if branch is contributing string to output, false otherwise
     */
    bool serialWrite(char *string, const int bufsize, const char *signal = NULL);

  /**
   * Get the forward self admittance of the branch
   * @param double - The forward self conductance Gff
   * @param double - The forward self susceptance Bff
   */
  bool getForwardSelfAdmittance(double*,double*);

  /**
   * Get the reserve self admittance of the branch
   * @param double - The reserve self conductance Gtt
   * @param double - The reverse self susceptance Btt
   */
  bool getReverseSelfAdmittance(double*,double*);

  /**
   * Get the forward transfer admittance of the branch
   * @param double - The forward transfer conductance Gft
   * @param double - The forward transfer susceptance Bft
   */
  bool getForwardTransferAdmittance(double*,double*);

  /**
   * Get the reserve transfer admittance of the branch
   * @param double - The reserve transfer conductance Gtf
   * @param double - The reverse transfer susceptance Btf
   */
  bool getReverseTransferAdmittance(double*,double*);

  private:
  int p_nparlines; // Number of parallel lines
  std::vector<int> p_status; // Status of the lines
  std::vector<std::string> p_cktid; // circuit id
  // Each line can be defined by the complex two-port network

  // |If| =  |  Yff Yft | |Vf|
  // |It|    |  Ytf Ytt | |Vt|
  // where Yxx = Gxx + sqrt(-1)*Bxx
  std::vector<double> p_Gff, p_Bff;
  std::vector<double> p_Gft, p_Bft;
  std::vector<double> p_Gtf, p_Btf;
  std::vector<double> p_Gtt, p_Btt;

  int p_mode;

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<gridpack::component::BaseBranchComponent>(*this)
      & p_nparlines
      & p_status
      & p_cktid
      & p_Gff & p_Bff
      & p_Gft & p_Bft
      & p_Gtf & p_Btf
      & p_Gtt & p_Btt
      & p_mode;
    
  }  

};


/// The type of network used in the examples application
typedef network::BaseNetwork<DSBus, DSBranch > DSNetwork;


}     // dsimplicit
}     // gridpack

BOOST_CLASS_EXPORT_KEY(gridpack::dsimplicit::DSBus)
BOOST_CLASS_EXPORT_KEY(gridpack::dsimplicit::DSBranch)


#endif
