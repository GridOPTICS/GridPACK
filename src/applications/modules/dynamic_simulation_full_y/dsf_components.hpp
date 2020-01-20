/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   dsf_components.hpp
 * @author Shuangshuang Jin 
 * @Last modified:   May 13, 2015
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _dsf_components_h_
#define _dsf_components_h_

//#define USE_FNCS
/**
 * Some preprocessor string declarations. These will need to be put in an
 * include file someplace else. Just declare them here for the time being.
 */

#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/applications/components/y_matrix/ymatrix_components.hpp"
#include "generator_factory.hpp"
#include "relay_factory.hpp"
#include "load_factory.hpp"

namespace gridpack {
namespace dynamic_simulation {

enum DSMode{YBUS, YL, YDYNLOAD, PG, onFY, posFY, jxd, make_INorton_full, bus_relay, branch_relay};

// Small utility structure to encapsulate information about fault events
struct Event{
  double start,end; // start and end times of fault
  double step;      // time increment of fault (not used?)
  char tag[3];      // 2-character identifier of line or generator
  bool isGenerator; // fault is a generator failure
  int bus_idx;      // index of bus hosting generator
  bool isLine;      // fault is a line failure
  int from_idx;     // "from" bus of line
  int to_idx;       // "to" bus of line
};

class DSFullBranch;
class DSFullBus;

class DSFullBus
  : public gridpack::ymatrix::YMBus
{
  public:

#ifdef USE_FNCS
    /**
     * Data package for use with FNCS framework
     */
    typedef struct{int busID;
      gridpack::ComplexType voltage;
    } voltage_data;
#endif

    /**
     *  Simple constructor
     */
    DSFullBus(void);

    /**
     *  Simple destructor
     */
    ~DSFullBus(void);

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

    void setValues(ComplexType *values);

    /**
     * Set values of YBus matrix. These can then be used in subsequent
     * calculations
     */
    void setYBus(void);

    /**
     * Set initial values of vectors for integration. 
     * These can then be used in subsequent calculations
     * @param ts time step
     */
    void initDSVect(double ts);

    /**
     * Update values for vectors in each integration time step (Predictor)
     * @param flag initial step if true
     */
    void predictor_currentInjection(bool flag);

    /**
     * Update values for vectors in each integration time step (Corrector)
     * @param flag initial step if true
     */
    void corrector_currentInjection(bool flag);

    /**
     * Update values for vectors in each integration time step (Predictor)
     * @param t_inc time step increment
     * @param flag initial step if true
     */
    void predictor(double t_inc, bool flag);

    /**
     * Update values for vectors in each integration time step (Corrector)
     * @param t_inc time step increment
     * @param flag initial step if true
     */
    void corrector(double t_inc, bool flag);
	
	/**
     * Update dynamic load internal relays action
     */
	void dynamicload_post_process(double t_inc, bool flag);

    /**
     * Get rotor angle of generators
     */
    double getAngle();

    /**
     * Set volt from volt_full
     */
    void setVolt(bool flag);
	
	/**
     * compute bus frequency and set it to dynamic load models
     */
	void updateFreq (double delta_t);

    /**
     * Get values of YBus matrix. These can then be used in subsequent
     * calculations
     */
    gridpack::ComplexType getYBus(void);

    /**
     * Load values stored in DataCollection object into DSFullBus object. The
     * DataCollection object will have been filled when the network was created
     * from an external configuration file
     * @param data: DataCollection object contain parameters relevant to this
     *       bus that were read in when network was initialized
     */
    void load(const boost::shared_ptr<gridpack::component::DataCollection> &data);
	
 	/**
     * load parameters for the extended buses from composite load model
     */
	void LoadExtendedCmplBus(const boost::shared_ptr<gridpack::component::DataCollection> &data);
	
	/**
     * set voltage for the extended buses from composite load model
     */
	void setExtendedCmplBusVoltage(const boost::shared_ptr<gridpack::component::DataCollection> &data);

    /**
     * Set the mode to control what matrices and vectors are built when using
     * the mapper
     * @param mode: enumerated constant for different modes
     */
    void setMode(int mode);

    /**
     * Return the value of the voltage magnitude on this bus
     * @return: voltage magnitude
     */
    double getVoltage(void);

    /**
     * Return the value of the phase angle on this bus
     * @return: phase angle
     */
    double getPhase(void);
	
	/**
     * Return the value of whether the bus is an extended bus due to compositeload
	 * return true if this is an extended bus due to compositeload
     */
	int checkExtendedLoadBus(void);
	
	/**
     * Set the value of the voltage magnitude on this bus
     */
    void setVoltage(double mag);

    /**
     * Set the value of the phase angle on this bus
     */
    void setPhase(double ang);
	
	/**
     * Set the point of the related extended transformer branch of this bus, due to composite load model
     */
    void setCmplXfmrPt(gridpack::dynamic_simulation::DSFullBranch* p_CmplXfmr);
	
	/**
     * Set the point of the related extended feeder branch of this bus, due to composite load model
     */
    void setCmplXfeederPt(gridpack::dynamic_simulation::DSFullBranch* p_CmplFeeder);
	
	/**
     * Return the complex value of the voltage on this bus
     * @return: complex value of the voltage
     */
    gridpack::ComplexType getComplexVoltage(void); //renke add
	
	/**
     * compute the value of the voltage frequency on this bus
     * @return: voltage frequency
     */
    void computeBusVolFrequency(double timestep); //renke add
	
	/**
     * return the value of the voltage frequency on this bus
     * @return: voltage frequency
     */
    double getBusVolFrequency(void); //renke add
	
	void setBusVolFrequencyFlag(bool flag); //renke add
	
	void printbusvoltage (void); //renke add
	
	void setWideAreaFreqforPSS(double freq); //renke hard coded;
	
	
	/**
     * update the relay status associate with this bus
     */
    bool updateRelay(bool flag, double delta_t); //renke add
	
	/**
     * update the old bus voltage with this bus
     */
    void updateoldbusvoltage (); //renke add

    /**
     * Return the number of generators on this bus
     * @return number of generators on bus
     */
    int getNumGen(void);

    /**
     * Return whether or not a bus is isolated
     * @return true if bus is isolated
     */
    bool isIsolated(void) const;

    /**
     * Set values of the IFunction on this bus (gen)
     */
    void setIFunc(void);

    /**
     * Set values of the IJaco on this bus (gen)
     */
    void setIJaco(void);
	
	void setBranchRelayFromBusStatus(bool sta);
	void setBranchRelayToBusStatus(bool sta);
	void setRelayTrippedbranch(gridpack::component::BaseBranchComponent* branch_ptr);
	void clearRelayTrippedbranch();

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
    bool serialWrite(char *string, const int bufsize, const char *signal);

    /**
     * Set an internal parameter that specifies that the rotor speed and angle
     * for the generator corresponding to the string tag are to be printed to
     * output
     * @param tag 2-character identifier of generator
     * @param flag set to true to monitor generator
     */
    void setWatch(std::string tag, bool flag);

    /**
     * Add constant impedance load admittance to diagonal elements of Y-matrix
     */
    void addLoadAdmittance();

    /**
     * Set load on bus
     * @param pl real load
     * @param ql imaginary load
     */
    void setLoad(double pl, double ql);

    /**
     * Get load on bus
     * @param pl real load
     * @param ql imaginary load
     */
    void getLoad(double *pl, double *ql);
	
	bool checkisolated();

    /**
     * Set value of real power on individual generators
     * @param tag generator ID
     * @param value new value of real power
     * @param data data collection object associated with bus
     */
    void setGeneratorRealPower(std::string tag, double value,
        gridpack::component::DataCollection *data);

    /**
     * Set value of real power on individual loads
     * @param tag load ID
     * @param value new value of real power
     * @param data data collection object associated with bus
     */
    void setLoadRealPower(std::string tag, double value,
        gridpack::component::DataCollection *data);

    /**
     * Return a list of watched generators
     * @return list of generator tags
     */
    std::vector<std::string> getWatchedGenerators();

    /**
     * Return a vector of watched values
     * @return rotor angle and speed for all watched generators on bus
     */
    std::vector<double> getWatchedValues();

    /**
     * Check generators for frequency violations
     * @param start time at which monitoring begins
     * @param time current time
     * @return true if no violation has occured
     */
    bool checkFrequency(double start, double time);

    /**
     * Check generators for frequency violations
     * @param limit maximum allowable frequency excursion
     * @return true if no violation has occured
     */
    bool checkFrequency(double limit);

    /**
     * Scale value of real power on all generators
     * @param character ID for generator
     * @param value scale factor for real power
     * @param data data collection object for bus holding generators
     */
    void scaleGeneratorRealPower(std::string tag, double value,
        boost::shared_ptr<gridpack::component::DataCollection> data);

    /**
     * Scale value of real power on loads
     * @param character ID for load
     * @param value scale factor for real power
     */
    void scaleLoadPower(std::string tag, double value);

    /**
     * Reset power for generators and loads back to original values
     * @param data data collection object for bus
     */
    void resetPower(boost::shared_ptr<gridpack::component::DataCollection> data);

    /**
     * Get available margin for generator
     * @param tag character ID for generator
     * @param current initial generation
     * @param slack amount generation can be reduced
     * @param excess amount generation can be increased
     * @param status current status of generator
     */
    void getGeneratorMargins(std::vector<std::string> &tag,
        std::vector<double> &current, std::vector<double> &slack,
        std::vector<double> &excess,std::vector<int> &status);

    /**
     * Get current value of loads
     * @param tag character ID for load
     * @param pl initial value of load real power
     * @param ql initial value of load reactive power
     * @param status current status of load
     */
    void getLoadPower(std::vector<std::string> &tag, std::vector<double> &pl,
        std::vector<double> &ql, std::vector<int> &status);

    /**
     * Get list of generator IDs
     * @return vector of generator IDs
     */
    std::vector<std::string> getGenerators();

    /**
     * Get list of load IDs
     * @return vector of load IDs
     */
    std::vector<std::string> getLoads();

    /**
     * Get area parameter for bus
     * @return bus area index
     */
    int getArea();

    /**
     * Get zone parameter for bus
     * @return bus zone index
     */
    int getZone();

    /**
     * Label bus as a source for real time path rating
     * @param flag identify bus as source
     */
    void setSource(bool flag);

    /**
     * Label bus as a sink for real time path rating
     * @param flag identify bus as sink
     */
    void setSink(bool flag);

    /**
     * Store scale factor
     * @param scale factor for scaling generation or
     loads
     */
    void setScale(double scale);

#ifdef USE_FNCS
    /**
     * Retrieve an opaque data item from component.
     * @param data item to retrieve from component
     * @param signal string to control behavior of routine (currently ignored)
     * @return true if component is returning data item, false otherwise
     */
    bool getDataItem(void *data, const char *signal = NULL);
#endif

  private:
    double p_shunt_gs;
    double p_shunt_bs;
    bool p_shunt;
    int p_mode;
    double p_theta; // phase angle difference
    double p_ybusr, p_ybusi;
    double p_angle, p_voltage;
	double p_busvolfreq, pbusvolfreq_old; //renke add, bus voltage frequency at current and previous timesteps
    bool p_load;
    double p_pl, p_ql;
	double p_loadimpedancer, p_loadimpedancei;
    double p_sbase;
    bool p_isGen;
    int p_area;
    int p_zone;
    bool p_source;
    bool p_sink;
    double p_rtpr_scale;
    std::vector<double> p_pg, p_qg, p_savePg, p_negpg, p_negqg;
    std::vector<double> p_mva, p_r, p_dstr, p_dtr, p_gpmin, p_gpmax;
    int p_ngen, p_negngen;
    int p_ndyn_load, p_npowerflow_load;
    int p_type;
    gridpack::ComplexType p_permYmod;
    bool p_from_flag, p_to_flag;
	bool p_branchrelay_from_flag, p_branchrelay_to_flag;
	bool p_busrelaytripflag;
	int p_bextendedloadbus; // whether it is an extended load bus with composite load model
	                        // -1: normal bus
							//  1: LOW_SIDE_BUS
							//  2: LOAD_BUS
    std::vector<std::string> p_genid;
    std::vector<std::string> p_loadid;
    std::vector<int> p_gstatus;
    std::vector<bool> p_watch;

    // DAE related variables
    //double user_eqprime, user_pmech, user_gen_d0, user_gen_h; // User app context variables
    //int user_ngen; // User app context variables
    std::vector<double> p_h, p_d0;
    //std::vector<double> x, xdot; // DAE variables
    std::vector<gridpack::ComplexType> p_pelect, p_volt;

    std::vector<gridpack::ComplexType> p_elect_final, p_mech_final;
    std::vector<gridpack::ComplexType> p_mac_ang_final, p_mac_spd_final;

    std::vector<gridpack::ComplexType> p_mac_ang_s0, p_mac_spd_s0;
    std::vector<gridpack::ComplexType> p_mac_ang_s1, p_mac_spd_s1;
    std::vector<gridpack::ComplexType> p_dmac_ang_s0, p_dmac_spd_s0;
    std::vector<gridpack::ComplexType> p_dmac_ang_s1, p_dmac_spd_s1;
    std::vector<gridpack::ComplexType> p_eqprime, p_pmech;

    std::vector<gridpack::ComplexType> p_eprime_s0, p_eprime_s1;

    std::vector<gridpack::ComplexType> p_INorton;
    gridpack::ComplexType p_volt_full;
	gridpack::ComplexType p_volt_full_old; //renke add
	double p_volt_full_old_real, p_volt_full_old_imag; //renke add
	bool bcomputefreq; // renke add

    gridpack::component::BaseBranchComponent* p_branch;
	gridpack::component::BaseBranchComponent* p_relaytrippedbranch; //renke add

    //std::vector<boost::shared_ptr<gridpack::dynamic_simulation::BaseGenerator> >
      //p_generators;
    std::vector<boost::shared_ptr<gridpack::dynamic_simulation::BaseGeneratorModel> >
      p_generators;
	  
    std::vector<boost::shared_ptr<gridpack::dynamic_simulation::BaseRelayModel> >
      p_loadrelays;   // renke add  

    std::vector<boost::shared_ptr<gridpack::dynamic_simulation::BaseLoadModel> >
      p_loadmodels;
	
    std::vector<double> p_powerflowload_p;	
    std::vector<double> p_powerflowload_p_save;	
	std::vector<double> p_powerflowload_q;	
    std::vector<double> p_powerflowload_q_save;	
   std::vector<int> p_powerflowload_status;
	
	gridpack::dynamic_simulation::DSFullBranch* p_CmplXfmrBranch, *p_CmplFeederBranch;  // the point to the transformer and feeder branch due to the added composite load model
	gridpack::dynamic_simulation::DSFullBus* p_CmplXfmrBus, *p_CmplFeederBus; // the point to the transformer and feeder bus due to the added composite load model
	gridpack::ComplexType p_CmplXfmrBusVolt_cplx; // value of the initialized complex voltage at the transformer bus due to the added composite load model
	double p_CmplXfmr_xxf, p_CmplXfmr_tap; // parameter of the extended transformer due to the added composite load model
	double loadMVABase; // load MVA Base is the load MVA base if the bus is an extend load bus (LOW_SIDE_BUS or LOAD_BUS) due to composite load model 

    bool p_isolated;

    // variables for monitoring the frequency of generators to find out
    // if any are going out of bounds
    std::vector<double> p_previousFrequency;
    std::vector<double> p_upIntervalStart;
    std::vector<bool> p_upStartedMonitoring;
    std::vector<double> p_downIntervalStart;
    std::vector<bool> p_downStartedMonitoring;


    friend class boost::serialization::access;

    template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
      {
        ar &
          boost::serialization::base_object<gridpack::ymatrix::YMBus>(*this)
          & p_shunt_gs
          & p_shunt_bs
          & p_shunt
          & p_mode
          & p_theta
          & p_ybusr & p_ybusi
          & p_angle & p_voltage
          & p_load
          & p_pl & p_ql
          & p_sbase
          & p_isGen
          & p_area & p_zone
          & p_source & p_sink
          & p_rtpr_scale
          & p_pg & p_qg & p_savePg
          & p_mva & p_r & p_dstr & p_dtr & p_gpmin & p_gpmax
          & p_ngen & p_type & p_permYmod
          & p_from_flag & p_to_flag
          & p_h & p_d0
          & p_pelect & p_eprime_s0 & p_eprime_s1 & p_eqprime & p_pmech 
          & p_mac_ang_s0 & p_mac_spd_s0
          & p_mac_ang_s1 & p_mac_spd_s1
          & p_dmac_ang_s0 & p_dmac_spd_s0
          & p_elect_final & p_mech_final
          & p_mac_ang_final & p_mac_spd_final
          & p_genid
          & p_previousFrequency
          & p_upIntervalStart & p_upStartedMonitoring
          & p_downIntervalStart & p_downStartedMonitoring;
      }

};

class DSFullBranch
  : public gridpack::ymatrix::YMBranch {
  public:

    /**
     *  Simple constructor
     */
    DSFullBranch(void);

    /**
     *  Simple destructor
     */
    ~DSFullBranch(void);

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
     * Load values stored in DataCollection object into DSFullBranch object. The
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
    gridpack::ComplexType getTransformer(DSFullBus *bus);

    /**
     * Return the contribution to a bus from shunts
     * @param bus: pointer to the bus making the call
     * @return: contribution to Y matrix from shunts associated with branches
     */
    gridpack::ComplexType getShunt(DSFullBus *bus);

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
     * Return the updating factor that will be applied to the ybus matrix at
     * the branch relay trip phase
     * @return: value of update factor
     */
	gridpack::ComplexType getBranchRelayTripUpdateFactor(); //renke add

    /**
     * Check to see if an event applies to this branch and set appropriate
     * internal parameters
     * @param event a struct containing parameters that describe a fault
     * event in a dyanamic simulation
     */
    void setEvent(const Event &event);
	
	/**
     * update branch current
     */
	void updateBranchCurrent(); //RENKE ADD
	bool updateRelay(bool flag, double delta_t); //renke add
	
	/**
     * Set parameters of the transformer branch due to composite load model
     */
	void SetCmplXfmrBranch(double dx, double dtap); //RENKE ADD
	
	/**
     * Set parameters of the feeder branch due to composite load model
     */
	void SetCmplFeederBranch(double dr, double dx); //RENKE ADD
	
	/*
	* print the content of the DSFullBranch
	*/
	void printDSFullBranch();
	
	/**
     * check the type of the extended load branch type variable: p_bextendedloadbranch
     */
	int checkExtendedLoadBranchType(void);

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
	std::vector<int> p_newtripbranchcktidx;
	
    int p_mode;
    double p_ybusr_frwd, p_ybusi_frwd;
    double p_ybusr_rvrs, p_ybusi_rvrs;
    double p_theta;
    double p_sbase;
    std::vector<int> p_branch_status;
    int p_elems;
    bool p_active;
    bool p_event;
	bool p_branchrelaytripflag;
	int  p_bextendedloadbranch; //whether this branch is added by composite load model
								// -1: normal branch
								//  1: transformer branch added by composite load model
								//  2: feeder branch added by composite load model
	
	std::vector<boost::shared_ptr<gridpack::dynamic_simulation::BaseRelayModel> >
      p_linerelays;   // renke add 
	std::vector<int> p_relaybranchidx; //renke add 
    std::vector<std::string> p_ckt; //renke_add
	std::vector<gridpack::ComplexType> p_branchcurrent; //renke add
	gridpack::ComplexType p_branchfrombusvolt; //renke add
	gridpack::ComplexType p_branchtobusvolt; //renke add

    friend class boost::serialization::access;

    template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
      {
        ar &
          boost::serialization::base_object<gridpack::ymatrix::YMBranch>(*this)
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
          & p_theta & p_sbase
          & p_branch_status
          & p_elems & p_active & p_event;
      }
};

/// The type of network used in the dynamic_simulation application
typedef network::BaseNetwork<DSFullBus, DSFullBranch > DSFullNetwork;

}     // dynamic_simulation
}     // gridpack

BOOST_CLASS_EXPORT_KEY(gridpack::dynamic_simulation::DSFullBus)
BOOST_CLASS_EXPORT_KEY(gridpack::dynamic_simulation::DSFullBranch)

#endif
