/**
 * Constructor
 */
gridpack::component::DAEBaseInterface::DAEBaseInterface(void)
{
}

/**
 * Destructor
 */
gridpack::component::DAEBaseInterface::~DAEBaseInterface(void)
{
}

/**
 * initialize DAE variables. Assign initial values, which may be determined
 * elsewhere, to all DAE variables
 * @param values vector of initial values
 */
void gridpack::component::DAEBaseInterface::initializeDAEVariables(
    std::vector<_Data_type> &values)
{
}

/**
 * Return number of DAE variables. This includes both time-dependent and
 * algebraic variables.
 */
int gridpack::component::DAEBaseInterface::numVariables(void)
{
  return 0;
}

/**
 * Return vector of current values. The ordering of variables needs to
 * remain consistent across all calls and it is up to the user to ensure
 * that this happens
 * @param values vector of current values of the DAE variables
 */
void gridpack::component::DAEBaseInterface::getVariableValues(
    std::vector<_Data_type> &values)
{
  values.clear();
}

/**
 * Return vector of time derivatives of current values for all DAE
 * variables. Algebraic variables have a time derivative of zero.
 * @param values time derivative of all variables at the current state value
 */
void gridpack::component::DAEBaseInterface::getTimeDerivatives(
    std::vector<_Data_type> &values)
{
  values.clear();
}

/**
 * Set current values of all DAE values
 * @param values vector of current values
 */
void gridpack::component::DAEBaseInterface::setVariableValues(
    std::vector<_Data_type> &values)
{
}

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
 */
void gridpack::component::DAEBaseInterface::getDiagJacobian(
    std::vector<int> &idx, std::vector<int> &jdx,
    std::vector<_Data_type> &values)
{
  idx.clear();
  jdx.clear();
  values.clear();
}
void gridpack::component::DAEBaseInterface::getForwardJacobian(
    std::vector<int> &idx, std::vector<int> &jdx,
    std::vector<_Data_type> &values)
{
  idx.clear();
  jdx.clear();
  values.clear();
}
void gridpack::component::DAEBaseInterface::getReverseJacobian(
    std::vector<int> &idx, std::vector<int> &jdx,
    std::vector<_Data_type> &values)
{
  idx.clear();
  jdx.clear();
  values.clear();
}
