/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   kds_components.hpp
 * @author Da Meng and Yousu Chen
 * @date   2019-12-03 07:24:03 d3g096
 * 
 * @brief  
 * 
 * @ Modified by Xinya Li 6/30/2015 
 */
// -------------------------------------------------------------

#ifndef _kalman_components_h_
#define _kalman_components_h_

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/utilities/complex.hpp"
#include "gridpack/component/base_component.hpp"
#include "gridpack/component/data_collection.hpp"
#include "gridpack/applications/components/y_matrix/ymatrix_components.hpp"

namespace gridpack {
namespace kalman_filter{

enum KalmanMode{YBus, EnsembleX, Perturbation, YC, RefreshY,
                E_Ensemble1, E_Ensemble2, E_Ensemble3,
                V1, V2, V3, HX, HA, 
                X_INC,X_Update,Measurements,
                onFY, posFY};

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

class KalmanBus
  : public gridpack::ymatrix::YMBus
{
  public:
    /**
     *  Simple constructor
     */
    KalmanBus(void);

    /**
     *  Simple destructor
     */
    ~KalmanBus(void);

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
     * Set the number of ensemble samples on the buses
     * @param nsize number of ensembles
     */
    void setEnsembleSize(int nsize);

    /**
     * Set the distribution width
     * @param sigma width of guassian distribution
     * @param noise width of guassian distribution
     */
    void setGaussianWidth(double sigma, double noise);

    /**
     * Return number of rows and columns in matrix from component
     * Number of columns must be the same for all components
     * @return size of block contributed by component
     */
    void slabSize(int *rows, int *cols) const;

    /**
     * Set indices corresponding to the rows contributed by this
     * component
     * @param irow index of row contributed by this component (e.g. if component
     * contributes 3 rows then irow is between 0 and 2)
     * @param idx row index of row irow
     */
    void slabSetRowIndex(int irow, int idx);

    /**
     * Get list of row indices from component
     * @param idx list of row indices that component maps onto
     */
    void slabGetRowIndices(int *idx);

    /**
     * Get a list of row values contributed by this component and their
     * indices
     * @param values list of values for rows
     * @param idx indices for the matrix rows
     */
    void slabGetValues(std::vector<ComplexType*> &values, int *idx);

    /**
     * Transfer slab values to component
     * @param values list of slab values
     */
    void slabSetValues(ComplexType **values);

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
     * Load values stored in DataCollection object into KalmanBus object. The
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
     * Check to see if a fault event applies to this bus and set an internal
     * flag marking the bus as the "from" or "to" bus for the event
     * @param from_idx index of "from" bus for fault event
     * @param to_idx index of "to" bus for fault event
     * @param branch_ptr pointer to branch on which fault occurs
     */
    void setEvent(int from_idx, int to_idx,
        gridpack::component::BaseBranchComponent* branch_ptr);

    /**
     * Clear fault event from bus
     */
    void clearEvent();

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
     * Set voltage angle time series data
     * @param ang array containing time series
     */
    void setVAngSeries(double *ang);

    /**
     * Set voltage magnitude time series data
     * @param mag array containing time series
     */
    void setVMagSeries(double *mag);

    /**
     * Set time increment and number of timesteps
     * @param delta_t time step increment
     * @param nsteps number of steps in time series
     */
    void setTimeSteps(double delta_t, int nsteps);
    void printT();

    /**
     * Set current step
     * @param istep current step
     */
    void setCurrentTimeStep(int istep);

    /**
     * Create the ensemble of rotor phases and rotor speeds 
     */
    void createEnsemble();

    /**
     * Evaluate dX_dt_1 and then get X2
     */
    void evaluateX2();

    /**
     * Evaluate dX_dt_2 and then get X3
     */
    void evaluateX3();

    /**
     * The number of generators on this bus
     * @return number of generators
     */
    int numGenerators();

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

    /**
     * Additional generator parameters
     */
    int p_ngen;
    std::vector<double> p_omega_0;
    std::vector<double> p_delta_0;
    std::vector<double> p_omega_n;
    std::vector<double> p_delta_n;
    std::vector<double> p_Pm;
    std::vector<double> p_Dm;
    std::vector<double> p_mva, p_r, p_dstr, p_dtr;
    std::vector<double> p_hmach;
    std::vector<double> p_Ebase;
    std::vector<std::string> p_genID;

    /**
     * Arrays used in Kalman Filter analysis
     */
    double* p_mag_series;
    double* p_ang_series;
    double p_delta_t;
    int p_t_len;
    int p_nEnsemble;
    std::vector<int> p_X_idx;
    std::vector<int> p_E_idx;
    std::vector<int> p_V_idx;
    std::vector<int> p_V3_idx;
    std::vector<int> p_HX_idx;
    double p_sigma;
    double p_noise;
    int p_currentStep;
    double **p_omega1;
    double **p_delta1;
    double **p_omega_diff;
    double **p_delta_diff;
    ComplexType **p_V1;
    double **p_delta1_dt;
    double **p_omega1_dt;
    double **p_omega2;
    double **p_delta2;
    ComplexType **p_V2;
    double **p_delta2_dt;
    double **p_omega2_dt;
    double **p_omega3;
    double **p_delta3;
    ComplexType *p_V3;

    bool p_from_flag, p_to_flag;
    gridpack::component::BaseBranchComponent* p_branch;

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
      & p_omega_0
      & p_delta_0
      & p_omega_n
      & p_delta_n
      & p_Pm & p_Dm
      & p_ngen
      & p_mva & p_r
      & p_dstr & p_dtr
      & p_from_flag & p_to_flag;
  }  

};

class KalmanBranch
  : public gridpack::ymatrix::YMBranch {
  public:
    // Small utility structure to encapsulate information about fault events
    struct Event{
      std::string name;
      double start,end;
      int from_idx, to_idx;
      double step;
    };

    /**
     *  Simple constructor
     */
    KalmanBranch(void);

    /**
     *  Simple destructor
     */
    ~KalmanBranch(void);

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
     * Load values stored in DataCollection object into KalmanBranch object. The
     * DataCollection object will have been filled when the network was created
     * from an external configuration file
     * @param data: DataCollection object contain parameters relevant to this
     *       branch that were read in when network was initialized
     */
    void load(const boost::shared_ptr<gridpack::component::DataCollection> &data);

    /**
     * Set the mode to control what matrices and vectors are built when using
     * the mapper
     * @param mode: enumerated constant for different modes
     */
    void setMode(int mode);

    /**
     * Return the updating factor that will be applied to the ybus matrix at
     * the clear fault phase
     * @return: value of update factor
     */
    gridpack::ComplexType getPosfy11YbusUpdateFactor(int sw2_2, int sw3_2);
    gridpack::ComplexType getUpdateFactor();

    /**
     * Check to see if an event applies to this branch and set appropriate
     * internal parameters
     * @param event a struct containing parameters that describe a fault
     * event in a dyanamic simulation
     */
    void setEvent(const Event &event);

    /**
     * Write output from branches to standard out
     * @param string (output) string with information to be printed out
     * @param bufsize size of string buffer in bytes
     * @param signal an optional character string to signal to this
     * routine what about kind of information to write
     * @return true if branch is contributing string to output, false otherwise
     */
    bool serialWrite(char *string, const int bufsize, const char *signal = NULL);

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
    bool p_event;
    int p_numElements;
    std::vector<int> p_colJidx;
    std::vector<int> p_rowJidx;
    std::vector<int> p_colRidx;
    std::vector<int> p_rowRidx;
    std::vector<int> p_vecZidx;
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
      & p_event
      & p_numElements
      & p_colJidx
      & p_rowJidx
      & p_colRidx
      & p_rowRidx
      & p_vecZidx;
  }  

};


}     // kalman_filter
}     // gridpack

BOOST_CLASS_EXPORT_KEY(gridpack::kalman_filter::KalmanBus)
BOOST_CLASS_EXPORT_KEY(gridpack::kalman_filter::KalmanBranch)


#endif
