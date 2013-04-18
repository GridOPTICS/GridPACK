// Simple outline of data collection object
class DataCollection {
public:
  DataCollection(void);
  ~DataCollection(void);
  void AddValue(char* name, int value);
  void AddValue(char* name, double value);

  bool GetValue(char* name, int *value);
  bool GetValue(char* name, double *value);
private:
  std::map<std::string, int> p_ints; 
  std::map<std::string, double> p_doubles; 
}

  // a code fragment illustrating how you would move raw data from a data
  // collection object to actual bus and branch model objects

  BusField *buses = network->GetBusField("buses");
  BranchField *branches = network->GetBranchField("branches");

  // busData and branchData are fields containing objects of type DataCollection
  BusField *busData = network->GetBusField("raw_bus_data");
  BranchField *branchData = network->GetBranchField("raw_branch_data");

  int numBus = buses->Size();
  int numBranch = branches->Size();

  // load method is something that access data elements in a DataCollection
  // object and moves it into the corresponding bus or branch object. Load
  // method is defined as BusModel::load(DataCollection data) and similarly for
  // a BranchModel
  for (i=0; i<numBus; i++) {
    buses[i]->load(busData[i]);
  }
  for (i=0; i<numBranch; i++) {
    branches[i]->load(branchData[i]);
  }
