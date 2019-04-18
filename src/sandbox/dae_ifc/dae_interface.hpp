namespace gridpack{
namespace component{

class DAEBaseInterface {
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

class DAEBusInterface: DAEBaseInterface {
  public:

    /**
     * Constructor
     */
    DAEBusInterface(void);

    /**
     * Destructor
     */
    ~DAEBusInterface(void);

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

    /**
     * Do we need to have something on the branches to account for off-diagonal
     * elements in the Jacobian?
     */

    /**
     * Add a new device to the bus
     */
    void addDevice(boost::shared_ptr<DAEBaseInterface> device);

  private:

    std::vector<boost::shared_ptr<DAEBaseInterface> p_devices;
};


} // component
} // gridpack
