// Simple outline of data collection object
class DataCollection {
public:
  /**
   * Simple constructor
   */
  DataCollection(void);

  /**
   * Simple destructor
   */
  ~DataCollection(void);

  /**
   *  Add variables to DataCollection object
   *  @param name: name given to data element
   *  @param value: value of data element
   */
  void AddValue(char *name, int value);
  void AddValue(char *name, long value);
  void AddValue(char *name, bool value);
  void AddValue(char *name, char *value);
  void AddValue(char *name, float value);
  void AddValue(char *name, double value);
  void AddValue(char *name, SingleComplex value);
  void AddValue(char *name, DoubleComplex value);

  /**
   *  Modify current value of existing data element in
   *  DataCollection object
   *  @param name: name of data element
   *  @param value: new value of data element
   *  @return: false if no element of the correct name and type exists in
   *  DataCollection object
   */
  bool SetValue(char *name, int value);
  bool SetValue(char *name, long value);
  bool SetValue(char *name, bool value);
  bool SetValue(char *name, char *value);
  bool SetValue(char *name, float value);
  bool SetValue(char *name, double value);
  bool SetValue(char *name, SingleComplex value);
  bool SetValue(char *name, DoubleComplex value);

  /**
   *  Retrieve current value of existing data element in
   *  DataCollection object
   *  @param name: name of data element
   *  @param value: current value of data element
   *  @return: false if no element of the correct name and type exists in
   *  DataCollection object
   */
  bool GetValue(char *name, int *value);
  bool GetValue(char *name, long *value);
  bool GetValue(char *name, bool *value);
  bool GetValue(char *name, char *std::string);
  bool GetValue(char *name, float *value);
  bool GetValue(char *name, double *value);
  bool GetValue(char *name, SingleComplex *value);
  bool GetValue(char *name, DoubleComplex *value);
private:
  std::map<std::string, int> p_ints; 
  std::map<std::string, long> p_longs; 
  std::map<std::string, bool> p_bools; 
  std::map<std::string, std::string> p_strings; 
  std::map<std::string, float> p_floats; 
  std::map<std::string, double> p_doubles; 
  std::map<std::string, SingleComplex> p_singleComplex; 
  std::map<std::string, DoubleComplex> p_doubleComplex; 
}

main(int argc, char **argv)
{
  // Read configuration parameters from external file
  Configure configure;
  if (!configure.ImportConfiguration("powerflow.in")) {
     //Throw error
  }
  // Read network from external file
  ReaderPTI reader;
  BaseNetwork *network = reader->GetNetwork("network.pti");
  if (!network) {
    //Throw error
  }

  // Partition network amongst processors and add ghost buses and branches
  Partioner partition;
  partition.PartitionNetwork(network);

  PowerflowFactory *factory = new PowerflowFactory(network);
  factory->Load();

  factory->SetState();

  PowerflowSolver *solver = new PowerflowSolver(network, factory);

  solver->solve();

  // Do something to export results from application
}

PowerflowFactory::Load(void) 
{
  // a code fragment illustrating how you would move raw data from a data
  // collection object to actual bus and branch model objects
  BusField<BusModel> *buses = network->GetBusField("buses");
  BranchField<BranchModel> *branches = network->GetBranchField("branches");

  // busData and branchData are fields containing objects of type DataCollection
  BusField<DataCollection> *busData = network->GetBusField("raw_bus_data");
  BranchField<DataCollection> *branchData = network->GetBranchField("raw_branch_data");

  int numBus = buses->Size();
  int numBranch = branches->Size();

  // load method is something that accesses data elements in a DataCollection
  // object and moves it into the corresponding bus or branch object. Load
  // method is defined as BusModel::load(DataCollection data) and similarly for
  // a BranchModel
  for (i=0; i<numBus; i++) {
    (*buses)[i]->load((*busData)[i]);
  }
  for (i=0; i<numBranch; i++) {
    (*branches)[i]->load((*branchData)[i]);
  }
}

PowerflowFactory::SetState(void) 
{
  for (i=0; i<numBus; i++) {
    (*buses)[i]->evaluateDiag(network);
  }
  for (i=0; i<numBranch; i++) {
    (*branches)[i]->evaluateOffDiag(network);
  }
}
