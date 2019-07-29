namespace gridpack{
namespace component{

template <_Data_type>
class BaseDevice {
  public:
    /**
     * All variables are treated on an equal footing, which matches the
     * interfaces in both PETSc and SUNDIALS. This should be the only thing that
     * is exposed to the user, other details can be buried within the DAE solver
     * class.
     */

    /**
     * Constructor
     */
    BaseDevice(void);

    /**
     * Destructor
     */
    ~BaseDevice(void);

    /**
     * Initialize variables. Assign initial values, which may be determined
     * elsewhere, to all variables
     * @param values vector of initial values
     */
    void initializeVariables(std::vector<_Data_type> &values);

    /**
     * Return number of variables from the device. This includes both
     * time-dependent and algebraic variables.
     */
    int numVariables(void);

    /**
     * Return vector of current values. The ordering of variables needs to
     * remain consistent across all calls and it is up to the user to ensure
     * that this happens
     * @param values vector of current values of the DAE variables
     */
    void getDeviceValues(std::vector<_Data_type> &values);

    /**
     * Return vector of time derivatives of current values for all DAE
     * variables. Algebraic variables have a time derivative of zero.
     * @param values time derivative of all variables at the current state value
     */
    void getDeviceTimeDerivatives(std::vector<_Data_type> &values);

    /**
     * Set current values of all DAE values
     * @param values vector of current values
     * @param indices index values of variables
     */
    void setDeviceVariableValues(std::vector<int> &indices,
        std::vector<_Data_type> &values);

    /**
     * Get Jacobian block for variables for non-zero matrix elements.
     * Variable number is based on the number of variables in the buses
     * and must be consistent with ordering used in other calls. Only indices
     * within the block need to considered, global positioning of the block is
     * handled by the DAE solver. The diagonal blocks are created on buses and
     * the forward and reverse blocks are created on branches.
     * @param idx vector of row indices of non-zero elements
     * @param jdx vector of column indices of non-zero elements
     * @param values vector on non-zero elements of Jacobian matrix block
     * @param time_const variable provided by time integrator
     */
    void getDeviceJacobianValues(std::vector<int> &idx, std::vector<int> &jdx,
        std::vector<_Data_type> &values, double time_const);
};
} // component
} // gridpack
