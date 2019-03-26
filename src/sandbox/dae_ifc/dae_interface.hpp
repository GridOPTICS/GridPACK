namespace gridpack{
namespace component{

class DAEDeviceInterface {
  public:

    /**
     * Constructor
     */
    DAEBaseInterface(void);

    /**
     * Destructor
     */
    ~DAEBaseInterface(void);

    /**
     * functions for time dependent values
     */

    /**
     * Return number of time-dependent variables
     */
    int numTimeDependentVariables(void);

    /**
     * Return vector of current values for time-dependent variables
     */
    void currentValues(std::vector<double> &values);

    /**
     * Return vector of time derivatives of current values for
     * time-dependent variables
     */
    void currentTimeDerivatives(std::vector<double> &values);

    /**
     * Set current values
     */
    void setCurrentValues(std::vector<double> &values);

    /**
     * Functions for algebraic variables
     */

    /**
     * Return the number of algebraic variables
     */
    int numAlgebraicVariables(void);
    
    /**
     * Set the current values of algebraic variables
     */
    void setAlgebraicValues(std::vector<double> &values);

    /**
     * Get the current values of algebraic variables
     */
    void getAlgebraicValues(std::vector<double> &values);

    /**
     * Get Jacobian block for algebraic values. Assume column or row major form
     * (TBD)
     */
    void getJacobian(std::vector<double> &values);
};

class DAEBusInterface: DAEDeviceInterface {
  public:

    /**
     * Constructor
     */
    DAEBaseInterface(void);

    /**
     * Destructor
     */
    ~DAEBaseInterface(void);

    /**
     * functions for time dependent values
     */

    /**
     * Return number of time-dependent variables
     */
    int numTimeDependentVariables(void);

    /**
     * Return vector of current values for time-dependent variables
     */
    void currentValues(std::vector<double> &values);

    /**
     * Return vector of time derivatives of current values for
     * time-dependent variables
     */
    void currentTimeDerivatives(std::vector<double> &values);

    /**
     * Set current values
     */
    void setCurrentValues(std::vector<double> &values);

    /**
     * Functions for algebraic variables
     */

    /**
     * Return the number of algebraic variables
     */
    int numAlgebraicVariables(void);
    
    /**
     * Set the current values of algebraic variables
     */
    void setAlgebraicValues(std::vector<double> &values);

    /**
     * Get the current values of algebraic variables
     */
    void getAlgebraicValues(std::vector<double> &values);

    /**
     * Get Jacobian block for algebraic values. Assume column or row major form
     * (TBD)
     */
    bool getJacobian(std::vector<double> &values);

    /**
     * Do we need to have something on the branches to account for off-diagonal
     * elements in the Jacobian?
     */

};

/**
 * This class is currently assuming that there are potentially some variables
 * to the DAE solver that are contributed directly by the bus. If this is not
 * true, then this class could be simplified considerably.
 */

class DAEBusInterface : DAEDeviceInterface {
  public:

    /**
     * Add another DAEDeviceInterface object to internal list of devices
     */
    void addDevice(boost::shared_ptr<DAEDeviceInterface> device);

    /**
     * Return total number of time dependent variables on all devices.
     */
    int totalTimeDependentVariables();

    /**
     * get the current values in the device list in the values vector
     */
    void totalCurrentValues(std::vector<double> &values);

    /**
     * set values in the device list using the values vector
     */
    void setTotalCurrentValues(std::vector<double> &values)

    /**
     * append the current time derivatives in the device list to the values
     * vector
     */
    void totalCurrentTimeDerivatives(std::vector<double> &values); 

    /**
     * return total number of algebraic variables
     */
    int totalAlgebraicVariables();

    /**
     * get algebraic variables in the device list
     */
    void getTotalAlgebraicValues(std::values<double> &values);

    /**
     * set algebraic variables in the device list
     */
    void setTotalAlgebraicValues(std::values<double> &values);

    /*
     * Get Jacobian block for algebraic values. Assume column or row major form
     * (TBD)
     */
    bool getTotalJacobian(std::vector<double> &values);

  private:

    std::vector<boost::shared_ptr<DAEDeviceInterface> p_devices;
};
} // component
} // gridpack
