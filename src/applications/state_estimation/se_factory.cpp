#define BUSES "buses"
#define BRANCHES "branches"
#define MEASUREMENTS "measurements"
#define RAW_BUS_DATA "raw_bus_data"
#define RAW_BRANCH_DATA "raw_branch_data"

SEFactory::SEFactory(BaseNetwork *network)
{
  p_network = network;
  BusField<BusModel> *buses;
  if (!p_network->getBusField(BUSES)) {
    buses = new BusField<BusModel>;
    p_network->addBusField(BUSES,buses);
  }
  BranchField<BranchModel> *branches;
  if (!p_network->getBrancheField(BRANCHES)) {
    branches = new BusField<BranchModel>;
    p_network->addBranchField(BRANCHES,branches);
  }
  BusField<MeasurementModel> *measurements;
  if (!p_network->getBusField(MEASUREMENTS)) {
    measurements = new BusField<MeasurementModel>;
    p_network->addBusField(MEASUREMENTS,measurements);
  }
}

SEFactory::~SEFactory(void)
{
}

SEFactory::load(void) 
{
  // a code fragment illustrating how you would move raw data from a data
  // collection object to actual bus and branch model objects
  BusField<BusModel> *buses = p_network->getBusField(BUSES);
  BranchField<BranchModel> *branches = p_network->getBranchField(BRANCHES);
  BusField<MeasurementModel> *measurements = p_network->getBusField(MEASUREMENTS)

  // busData and branchData are fields containing objects of type DataCollection
  BusField<DataCollection> *busData = p_network->GetBusField(RAW_BUS_DATA);
  BranchField<DataCollection> *branchData =
           p_network->GetBranchField(RAW_BRANCH_DATA);

  int numBus = buses->Size();
  int numBranch = branches->Size();

  // load method is something that accesses data elements in a DataCollection
  // object and moves it into the corresponding bus or branch object. Load
  // method is defined as BusModel::load(DataCollection *data) and similarly for
  // a BranchModel or MeasurementModel
  for (i=0; i<numBus; i++) {
    (*buses)[i]->load((*busData)[i]);
    (*measurements)[i]->load((*busData)[i])
  }
  for (i=0; i<numBranch; i++) {
    (*branches)[i]->load((*branchData)[i]);
  }
}

SEFactory::setState(void) 
{
  for (i=0; i<numBus; i++) {
    (*buses)[i]->evaluateDiag(p_network);
  }
  for (i=0; i<numBranch; i++) {
    (*branches)[i]->evaluateOffDiag(p_network);
  }
  for (i=0; i<numBus; i++) {
    (*measurements)[i]->evaluateMeasurement(p_network);
  }
}
