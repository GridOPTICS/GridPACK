/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */

#ifndef _dsf_components_h_
#define _dsf_components_h_

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

  enum DSMode{YBUS, YL, YDYNLOAD, PG, onFY, posFY, jxd, make_INorton_full, bus_relay, branch_relay, branch_trip_action, bus_Yload_change_P, bus_Yload_change_Q, BUSFAULTON,BUSFAULTOFF, LINESTATUSCHANGE, GENSTATUSCHANGE,INIT_V};

// Utility structure to encapsulate information about fault events
struct Event{
  bool isGenerator;     // event is a generator status change
  bool isBus;           // event is on the bus
  bool isLine;          // event is a line status change
  bool isBusFault;      // event is a bus fault
  bool isLineStatus;    // event is a line status change
  bool isGenStatus;     // event is a gen status change

  // fault information
  double start;         // start time of fault
  double end;           // end time of fault
  int bus_idx;          // fault or generator bus number
  double Gfault;       // Fault conductance 
  double Bfault;       // Fault susceptance

  // Line status change
  double time;          // line or generator status change time
  std::string tag;      // 2-character identifier of line or generator
  int from_idx;      // "from" bus of line
  int to_idx;        // "to" bus of line

  int status;        // Line or generator status

  double step;
  Event(void)
    : start(0.0), end(0.0), time(0.0),
      tag(3, '\0'), isGenerator(false), isBus(false), isLine(false), isBusFault(false), isLineStatus(false), isGenStatus(false) ,
      bus_idx(-1), from_idx(-1), to_idx(-1), Gfault(0.0),
      Bfault(99999), status(1)
  {
    tag = "1";
  }
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
     * Add shunt
     * @param gs shunt Gshunt value
     * @param bs shunt Bshunt value
     */
  void addShunt(double gs,double bs);

  /**
     * set fault
     * @param flag => true = fault on, false otherwise
     * @param gfault => fault conductance (pu)
     * @param bfault => fault susceptance (pu)
     */
  void setFault(double gfault, double bfault);

  /**
     getGenNum - Get the generator number given generator id
     
     @param:   gen_id - generator id
     @return:  idx - the internal index for the generator that can be used to access its parameters (-1 if not found)
  **/
  int getGenNum(std::string id);

  /**
     getGenStatus - Get the generator status 
   
     @param:  ckt_id - generator id
     @return: status - generator status
   
  **/
  int getGenStatus(std::string id);

    /**
     * set generator status
     * @param id => generator id
     * @param status => generator status
     */
  void setGenStatus(std::string id, int status);

  /**
   * clear fault
   */
  void clearFault();

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
	
	void setGeneratorObPowerBaseFlag(bool generator_observationpower_systembase);

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
	
	void setBranchTripAction(int from_idx, int to_idx,
		gridpack::component::BaseBranchComponent* branch_ptr);
		
	void clearBranchTripAction();
	
	void applyConstYLoad_Change_P(double loadPChangeMW);
	void applyConstYLoad_Change_Q(double loadPChangeMVAR);
	void clearConstYLoad_Change_P();
	void clearConstYLoad_Change_Q();
	
	bool setConstYLoadtoZero_P( );
	bool setConstYLoadtoZero_Q( );
	
	//bool setConstYLoadtoValue_P( double impedancer);
	//bool setConstYLoadtoValue_Q( double impedancei);
	bool setConstYLoadtoValue( double impedancer, double impedancei);

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
     * Set status variable on individual generators
     * @param tag generator ID
     * @param value new value of status
     * @param data data collection object associated with bus
     */
    void setGeneratorStatus(std::string tag, int status,
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
     * Return a vector of watched values
     * @return rotor angle and speed for generator with given tag (id)
     */
   std::vector<double> getWatchedValues(std::string tag);


    /**
     * Return rotor speed and angle for a specific generator
     * @param idx index of generator
     * @param speed generator rotor speed
     * @param angle generator rotor angle
	 * @param speed generator real power
     * @param angle generator reactive power
     */
    void getWatchedValues(int idx, double *speed, double *angle, double *genP, double *genQ);

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
     * Return pointer to generator corresponding to two-character id
     * @param id two character id labeling generator
     * @return pointer to generator device
     */
    gridpack::dynamic_simulation::BaseGeneratorModel *getGenerator(std::string id);

    /**
     * Return pointer to load corresponding to two-character id
     * @param id two character id labeling load
     * @return pointer to load device
     */
    gridpack::dynamic_simulation::BaseLoadModel *getLoad(std::string id);

    /**
     * Return pointer to relay corresponding to two-character id
     * @param id two character id labeling relay
     * @return pointer to relay device
     */
    gridpack::dynamic_simulation::BaseRelayModel *getRelay(std::string id);

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
     * Get list of IDs for dynamic loads
     * @return vector of dynamic load IDs
     */
    std::vector<std::string> getDynamicLoads();

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
	
	/**
	 * execute load scattering, the P and Q values of the STATIC load at certain buses will be changed to the values of 
	 * the loadP and loadQ
	*/
	void scatterInjectionLoad(double loadP, double loadQ);
	
	/**
	 * execute load scattering, the P and Q values of the STATIC load at certain buses will be changed to the values of 
	 * the loadP and loadQ, this function keeps the Y load component of the bus still at the bus, while only compenstates the difference
	*/
	void scatterInjectionLoad_compensateY(double loadP, double loadQ);
	
	/**
	* execute load scattering, the current values of the STATIC load at certain buses will be changed to the values of 
	* the curR and curI
	*/
	void scatterInjectionLoadConstCurrent(double curR, double curI);
	
	/**
     * apply load shedding for the loads in this bus
     */
	void applyLoadShedding(std::string loadid, double percentage);
	
	/**
	 * execute constant Y load shedding 	 
	 * percentage: float load shed percentage, for example -0.2 means shed 20%
	 */
	void applyConstYLoadShedding(double percentage );
	
	/**
	 * set the wide area control signals of the PSS of a certain generator
	 * input bus_number: generator bus number
	 * input bus_number: generator gen ID
	 * input wideAreaControlSignal:  wide area control signal for the PSS of the generator
	 */
	void setWideAreaControlSignal(std::string genid, double wideAreaControlSignal);
	
	/**
     * execute Grid Forming Inverter control parameters adjustment at this bus
	 * input controlTyp: 0: GFI mp adjust; 1: GFI mq adjust; 2: GFI Pset adjust; 3: GFI Qset adjust; others: invalid
	 * input bus_number: GFI gen ID
	 * input newParValScaletoOrg: GFI new parameter scale factor to the very initial parameter value at the begining of dynamic simulation
     */
    void applyGFIAdjustment(int controlType, std::string genid, double newParValScaletoOrg);
	
     /**
     * apply generator tripping for the specific generator with genid in this bus
     */
	void applyGeneratorTripping(std::string genid);


   /**
    * return the fraction online from dynamic load model
    * @param idx index of dynamic load model
    * @return fraction of dynamic load that is online
    */
   double getOnlineLoadFraction(int idx);

   /**
    * return the total load on the bus
    * @param total_p total real power load on bus
    * @param total_q total reactive power load on bus
    */
    void getTotalLoadPower(double &total_p, double &total_q) const;

    /**
     * return the real and reactive power produced by the generator indicated by
     * the tag variable
     * @param tag 2-character identifier for the generator
     * @param pg real power produced by generator
     * @param qg reactive power produced by generator
     * @return false if no generator corresponds to tag value.
     */
    bool getGeneratorPower(std::string tag, double &pg, double &qg) const;

   /**
    * return the power generated on the bus
    * @param total_p total active power generated on bus
    * @param total_q total reactive power generated on bus
    */
    void getTotalGeneratorPower(double &total_p, double &total_q) const;

  /**
     update the diag value contributions on line status change
     @param : ybr_self - contribution for line status change
  **/
  void diagValuesInsertForLineStatusChange(gridpack::ComplexType Ybr_self);
  
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
    // Used for faults only
    double p_gfault; // Fault conductance
    double p_bfault; // Fault susceptance
    // Used for line status change only
    bool p_line_status_change; // A line connected to this bus is changing its status
  gridpack::ComplexType p_yii; // This value must be inserted in the Ybus matrix for the line status change
    // Used for generator status change only
    bool p_gen_status_change; // A generaator connected to this bus is changing its status
  gridpack::ComplexType p_ygen; // This value must be inserted in the Ybus matrix for the gen status change

    double p_shunt_gs;
    double p_shunt_bs;
    bool p_shunt;
    int p_mode;
    double p_theta; // phase angle difference
    double p_ybusr, p_ybusi;
    double p_angle, p_voltage;
	double p_busvolfreq, pbusvolfreq_old; //renke add, bus voltage frequency at current and previous timesteps
    bool p_load;
    double p_pl, p_ql; //static loads for Y-bus formation, check the ::load() function for more info.
	double p_loadimpedancer, p_loadimpedancei;
	double p_Yload_change_r, p_Yload_change_i; // renke add, to enable constant-Y load changes during the dynamic simulation
	                                           // per unit value based on system MVA 100, increase 500 MW load, p_Yload_change_r = 5.0
	bool p_bconstYLoadSheddingFlag; //renke add, whether the static constant Y load of the bus has been shedding, true: has been shed
    double 	remainConstYLoadPerc; // renke add, remaining percent of the static constant Y load of the bus
	bool p_bscatterinjload_flag; //renke add, whether the static load of the bus could be modified at each time step 
	bool p_bscatterinjload_flag_compensateY; //renke add, whether the static load of the bus could be modified at each time step, this flag is corresponding to the
											 // scatterInjectionLoad_compensateY function, 
											 // which keep the Y load component still at the bus, while only compenstate the difference
	double p_scatterinjload_p, p_scatterinjload_q; //renke add, the value of the static load of the bus modified at each time step
	bool p_bscatterinjloadconstcur_flag; //renke add, whether the static load of the bus could be modified at each time step as constant current load
	double p_scatterinjload_constcur_r, p_scatterinjload_constcur_i; //renke add, the value of the static load of the bus modified at each time step as constant current load
    double p_sbase;
    bool p_isGen;
    int p_area;
    int p_zone;
    bool p_source;
    bool p_sink;
    double p_rtpr_scale;
    std::vector<double> p_pg, p_qg, p_savePg, p_negpg, p_negqg, p_genpg_nodynmodel, p_genqg_nodynmodel;
    std::vector<bool> p_gen_nodynmodel; // Generator does not have a dynamic model specificed in the dyr file
    std::vector<double> p_mva, p_r, p_dstr, p_dtr, p_gpmin, p_gpmax;
    int p_ngen, p_negngen, p_ngen_nodynmodel;
    int p_ndyn_load, p_npowerflow_load;
    int p_type;
    gridpack::ComplexType p_permYmod;
    bool p_from_flag, p_to_flag;
	bool p_Yload_change_P_flag, p_Yload_change_Q_flag; // renke add, indicating whether this bus has constant-Y load P or Q change,to enable constant-Y load change during dyn. simu.
	bool p_bConstYLoadSettoZero_P, p_bConstYLoadSettoZero_Q; // renke add, indicating whether this bus's constant-Y load P or Q has been already set to zero, check the function
															 // setConstYLoadtoZero_P to see the usage
															 
	bool p_bConstYLoadSettoValue; // renke add, indicating whether this bus's constant-Y load P or Q has been already set to some value, check the function
															 // setConstYLoadtoValue to see the usage
	bool p_branchrelay_from_flag, p_branchrelay_to_flag;
	bool p_branchtripaction_from_flag, p_branchtripaction_to_flag;
	bool p_busrelaytripflag;
	int p_bextendedloadbus; // whether it is an extended load bus with composite load model
	                        // -1: normal bus
							//  1: LOW_SIDE_BUS
							//  2: LOAD_BUS
    std::vector<std::string> p_genid;
    std::vector<std::string> p_loadid;
    std::vector<std::string> p_dynamic_loadid;
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
	std::vector<gridpack::component::BaseBranchComponent*> p_vec_tripactionbranch; // the branch connect to this bus that will be tripped for the branch trip action
	

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

    bool p_busfault; // Flag for bus fault

    friend class boost::serialization::access;

    template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
      {
        ar &
          boost::serialization::base_object<gridpack::ymatrix::YMBus>(*this)
	  & p_gfault
	  & p_bfault
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
     * Return the updating factor that will be applied to the ybus matrix at
     * the branch trip action
     * @return: value of update factor
     */
	gridpack::ComplexType getBranchTripActionUpdateFactorForBus(int busNo); //renke add
	gridpack::ComplexType getBranchTripActionUpdateFactorForward(); //renke add
	
	gridpack::ComplexType getBranchTripActionUpdateFactorReverse(); //renke add
	
	/**
	 * Clear fault event from branch
	*/
	void clearEvent();
	
    /**
     * Check to see if an event applies to this branch and set appropriate
     * internal parameters
     * @param event a struct containing parameters that describe a fault
     * event in a dyanamic simulation
     */
    void setEvent(const Event &event);
	
	bool setBranchTripAction(const std::string ckt_tag);
	bool setBranchTripAction();
	void clearBranchTripAction();
		
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
	
     /**
      * Return contributions to Y-matrix from a specific transmission element
      * @param tag character string for transmission element
      * @param Yjj contribution at "to" bus
      * @param Yji contribution for ji_th Y-matrix element
      */
  void getRvrsLineElements(const std::string tag,gridpack::ComplexType *Yjj, gridpack::ComplexType *Yji);

  /**
   * Return contributions to Y-matrix from a specific transmission element
   * @param tag character string for transmission element
   * @param Yii contribution at "from" bus
   * @param Yij contribution for ij_th Y-matrix element
   */
  void getLineElements(const std::string tag,
		       gridpack::ComplexType *Yii, gridpack::ComplexType *Yij);

  /**
     setLineStatus - Sets the line status and updates the associated
     branch and bus objects. 
   
     @param: ckt_id - circuit id
     @param: status - new line status
   
     Note: This method is used to
     update the branch status and update the bus/branch
     objects. It sets up values in the bus and branch objects
     so that incrementMatrix method called on the network Ybus
     uses these values to remove the branch contributions from
     the Y-bus matrix
  **/
  void setLineStatus(std::string ckt_id, int status);
  
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
	std::vector<bool> p_switched;
	std::vector<int> p_newtripbranchcktidx;
	std::vector<int> p_vec_tripaction_branchcktidx;
	
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
	bool p_branchactiontripflag;  // whether this branch will be tripped due to a branch trip action
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

  bool p_line_status_change; // Flag to indicate line status change
  gridpack::ComplexType p_yft,p_ytf; // Ybus off-diagonal contributions from this line
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
