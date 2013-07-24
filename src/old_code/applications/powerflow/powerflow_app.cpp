main(int argc, char **argv)
{
  // Read configuration parameters from external file
  Configure configure;
  if (!configure.importConfiguration("powerflow.in")) {
     //Throw error
  }
  // Read network from external file
  ReaderPTI reader;
  BaseNetwork *network = reader->getNetwork("network.pti");
  if (!network) {
    //Throw error
  }

  // Partition network amongst processors and add ghost buses and branches
  Partioner partition;
  partition.partitionNetwork(network);

  PowerflowFactory *factory = new PowerflowFactory(network);
  factory->Load();

  factory->SetState();

  // Create Y and H matrices for state estimation problem
  Mapper map;
  BusField<BusModel> *buses = network->getBusField(BUSES);
  BranchField<BranchModel> *branches = network->getBranchField(BRANCHES);
  Matrix *YMatrix = map.createMatrix(buses, branches);

  PowerflowSolver *solver = new PowerflowSolver(network, factory);

  solver->solve();

  // Do something to export results from application
}
