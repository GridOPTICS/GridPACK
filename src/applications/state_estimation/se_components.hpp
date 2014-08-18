/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   se_components.hpp
 * @author Yousu Chen 
 * @date   2/24/2014 
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _se_components_h_
#define _se_components_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/include/gridpack.hpp"
#include "gridpack/applications/components/y_matrix/ymatrix_components.hpp"

namespace gridpack {
namespace state_estimation{

enum SEMode{YBus,Jacobian_H, R_inv};

struct Measurement
{
  char p_type[4];
  //std::string p_type;
  int p_busid;
  int p_fbusid;
  int p_tbusid;
  //std::string p_ckt;
  char p_ckt[3];
  double p_value;
  double p_deviation;
};

class SEBus
  : public gridpack::ymatrix::YMBus
{
  public:
    /**
     *  Simple constructor
     */
    SEBus(void);

    /**
     *  Simple destructor
     */
    ~SEBus(void);

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
     * Return the size of the buffer used in data exchanges on the network.
     * For this problem, the voltage magnitude and phase angle need to be exchanged
     * @return size of buffer
     */
    int getXCBufSize(void);

    /**
     * Assign pointers for voltage magnitude and phase angle
     */
    void setXCBuf(void *buf);

    /**
     * Set values of YBus matrix. These can then be used in subsequent
     * calculations
     */
    void setYBus(void);

    /**
     * Get values of YBus matrix. These can then be used in subsequent
     * calculations
     */
    gridpack::ComplexType getYBus(void);

    /**
     * Load values stored in DataCollection object into SEBus object. The
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
     * Return the value of the voltage magnitude on this bus
     * @return voltage magnitude
     */
    double getVoltage(void);

    /**
     * Return the complex voltage on this bus
     * @return the complex voltage
     */
    ComplexType getComplexVoltage(void);

    /**
     * Return the value of the phase angle on this bus
     * @return: phase angle
     */
    double getPhase(void);

    /**
     * Return whether or not the bus is a PV bus (V held fixed in powerflow
     * equations)
     * @return true if bus is PV bus
     */
    bool isPV(void);

    /**
     * Return whether or not a bus is isolated
     * @return true if bus is isolated
     */
    bool isIsolated(void) const;

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
     * Add a measurement to the bus
     * @param measurement a measurement struct that will be used to assign
     * internal paramters
     */
    void addMeasurement(Measurement measurement);

    /**
     * Return number of rows in matrix from component
     * @return number of rows from component
     */
    int matrixNumRows() const;

    /**
     * Return number of columns in matrix from component
     * @return number of columnsows from component
     */
    int matrixNumCols() const;

    /**
     * Set row indices corresponding to the rows contributed by this
     * component
     * @param irow index of row contributed by this component (e.g. if component
     * contributes 3 rows then irow is between 0 and 2)
     * @param idx matrix index of row irow
     */
    void matrixSetRowIndex(int irow, int idx);

    /**
     * Set column indices corresponding to the columns contributed by this
     * component
     * @param icol index of column contributed by this component (e.g. if component
     * contributes 3 columns then icol is between 0 and 2)
     * @param idx matrix index of column icol
     */
    void matrixSetColIndex(int icol, int idx);

    /**
     * Get the row index corresponding to the rows contributed by this component
     * @param irow index of row contributed by this component (e.g. if component
     * contributes 3 rows then irow is between 0 and 2)
     * @return matrix index of row irow
     */
    int matrixGetRowIndex(int idx);

    /**
     * Get the column index corresponding to the columns contributed by this component
     * @param icol index of column contributed by this component (e.g. if component
     * contributes 3 columns then icol is between 0 and 2)
     * @return matrix index of column icol
     */
    int matrixGetColIndex(int idx);

    /**
     * Return the number of matrix values contributed by this component
     * @return number of matrix values
     */
    int matrixNumValues() const;

    /**
     * Return values from a matrix block
     * @param values: pointer to matrix block values
     * @param rows: pointer to matrix block rows
     * @param cols: pointer to matrix block cols
    */
    void matrixGetValues(ComplexType *values, int *rows, int *cols);

    /**
     * Return values from a vector
     * @param values: pointer to vector values (z-h(x))
     * @param idx: pointer to vector index 
     */
    void VectorGetElementValues(ComplexType *values, int *idx);

    /**
     * Configure buses with state estimation parameters. These can be used in
     * other methods
     */
    void configureSE(void);

  private:
    double p_shunt_gs;
    double p_shunt_bs;
    bool p_shunt;
    bool p_load;
    int p_mode;

    // p_v and p_a are initialized to p_voltage and p_angle respectively,
    // but may be subject to change during the NR iterations
    double p_v, p_a;
    double p_theta; //phase angle difference
    double p_ybusr, p_ybusi;
    double p_P0, p_Q0; //double p_sbusr, p_sbusi;
    double p_angle;   // initial bus angle read from parser
    double p_voltage; // initial bus voltage read from parser
    // newly added priavate variables:
    std::vector<double> p_pg, p_qg;
    std::vector<int> p_gstatus;
    std::vector<double> p_vs;
    std::vector<int> p_gid;
    double p_pl, p_ql;
    double p_sbase;
    double p_Pinj, p_Qinj;
    bool p_isPV;
    int p_numElements;
    std::vector<int> p_colJidx;
    std::vector<int> p_rowJidx;
    std::vector<int> p_colRidx;
    std::vector<int> p_rowRidx;
    std::vector<Measurement> p_meas;

    /**
     * Variables that are exchanged between buses
     */
    double* p_vMag_ptr;
    double* p_vAng_ptr;

private:


  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
//    ar & boost::serialization::base_object<gridpack::component::BaseBusComponent>(*this)
    ar  & boost::serialization::base_object<gridpack::ymatrix::YMBus>(*this)
      & p_shunt_gs
      & p_shunt_bs
      & p_shunt
      & p_load
      & p_mode
      & p_v & p_a & p_theta
      & p_ybusr & p_ybusi
      & p_P0 & p_Q0
      & p_angle & p_voltage
      & p_pg & p_qg
      & p_gstatus
      & p_vs & p_gid
      & p_pl & p_ql
      & p_sbase
      & p_Pinj & p_Qinj
      & p_isPV
      & p_numElements
      & p_colJidx
      & p_rowJidx
      & p_colRidx
      & p_rowRidx;
  }  

};

class SEBranch
  : public gridpack::ymatrix::YMBranch {
  public:
    /**
     *  Simple constructor
     */
    SEBranch(void);

    /**
     *  Simple destructor
     */
    ~SEBranch(void);

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
    gridpack::ComplexType getYBus(void);

    /**
     * Load values stored in DataCollection object into SEBranch object. The
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
    gridpack::ComplexType getTransformer(SEBus *bus);

    /**
     * Return the contribution to a bus from shunts
     * @param bus: pointer to the bus making the call
     * @return: contribution to Y matrix from shunts associated with branches
     */
    gridpack::ComplexType getShunt(SEBus *bus);

    /**
     * Set the mode to control what matrices and vectors are built when using
     * the mapper
     * @param mode: enumerated constant for different modes
     */
    void setMode(int mode);

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
     * Add a measurement to the branch
     * @param measurement a measurement struct that will be used to assign
     * internal paramters
     */
    void addMeasurement(Measurement measurement);

    /**
     * Return contribution to constraints
     * @param v: voltage at the other bus
     * @param theta: angle difference between two buses
     */
    void getVTheta(gridpack::state_estimation::SEBus *bus, double *v, double *theta);

    /**
     * Return contribution to constraints
     * @param v1, v2: voltages at buses
     * @param theta: angle difference between two buses
     */
    void getV1V2Theta(gridpack::state_estimation::SEBranch *branch, double *v1, double *v2, double *theta);

    /**
     * Return number of rows in matrix from component
     * @return number of rows from component
     */
    int matrixNumRows() const;

    /**
     * Return number of columns in matrix from component
     * @return number of columnsows from component
     */
    int matrixNumCols() const;

    /**
     * Set row indices corresponding to the rows contributed by this
     * component
     * @param irow index of row contributed by this component (e.g. if component
     * contributes 3 rows then irow is between 0 and 2)
     * @param idx matrix index of row irow
     */
    void matrixSetRowIndex(int irow, int idx);

    /**
     * Set column indices corresponding to the columns contributed by this
     * component
     * @param icol index of column contributed by this component (e.g. if component
     * contributes 3 columns then icol is between 0 and 2)
     * @param idx matrix index of column icol
     */
    void matrixSetColIndex(int icol, int idx);

    /**
     * Get the row index corresponding to the rows contributed by this component
     * @param irow index of row contributed by this component (e.g. if component
     * contributes 3 rows then irow is between 0 and 2)
     * @return matrix index of row irow
     */
    int matrixGetRowIndex(int idx);

    /**
     * Get the column index corresponding to the columns contributed by this component
     * @param icol index of column contributed by this component (e.g. if component
     * contributes 3 columns then icol is between 0 and 2)
     * @return matrix index of column icol
     */
    int matrixGetColIndex(int idx);

    /**
     * Return the number of matrix values contributed by this component
     * @return number of matrix values
     */
    virtual int matrixNumValues() const;

    /**
     * Return values from a matrix block
     * @param values: pointer to matrix block values
     * @param rows: pointer to matrix block rows
     * @param cols: pointer to matrix block cols
     */
    void matrixGetValues(ComplexType *values, int *rows, int *cols);

    /**
     * Return values from a vector
     * @param values: pointer to vector values (z-h(x))
     * @param idx: pointer to vector index 
     */
    void VectorGetElementValues(ComplexType *values, int *idx);

    /**
     * Configure branches with state estimation parameters. These can be used in
     * other methods
     */
    void configureSE(void);

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
    double p_theta;
    double p_sbase;
    std::vector<bool> p_branch_status;
    std::vector<std::string> p_tag;
    int p_elems;
    bool p_active;
    int p_numElements;
    std::vector<int> p_colJidx;
    std::vector<int> p_rowJidx;
    std::vector<int> p_colRidx;
    std::vector<int> p_rowRidx;
    std::vector<Measurement> p_meas;

private:


  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
//    ar & boost::serialization::base_object<gridpack::component::BaseBranchComponent>(*this)
    ar  & boost::serialization::base_object<gridpack::ymatrix::YMBranch>(*this)
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
      & p_theta
      & p_sbase
      & p_branch_status
      & p_tag
      & p_elems
      & p_active
      & p_numElements
      & p_colJidx
      & p_rowJidx
      & p_colRidx
      & p_rowRidx;
  }  

};


/// The type of network used in the contingency analysis application
typedef network::BaseNetwork<SEBus, SEBranch > SENetwork;


}     // state_estimation
}     // gridpack

BOOST_CLASS_EXPORT_KEY(gridpack::state_estimation::SEBus);
BOOST_CLASS_EXPORT_KEY(gridpack::state_estimation::SEBranch);


#endif
