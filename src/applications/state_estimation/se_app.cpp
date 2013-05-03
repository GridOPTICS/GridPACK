main(int argc, char **argv)
{
  // Read configuration parameters from external file
  Configure configure;
  if (!configure.importConfiguration("state_est.in")) {
     //Throw error
  }
  // Read network from external file
  ReaderPTI reader;
  BaseNetwork *network = reader->getNetwork("se_network.pti");
  if (!network) {
    //Throw error
  }

  // Partition network amongst processors and add ghost buses and branches
  Partioner partition;
  partition.partitionNetwork(network);

  SEFactory *factory = new SEFactory(network);
  factory->load();

  factory->setState();

  // Create Y and H matrices for state estimation problem
  Mapper map;
  BusField<BusModel> *buses = network->getBusField(BUSES);
  BranchField<BranchModel> *branches = network->getBranchField(BRANCHES);
  BusField<MeasurementModel> *measurements = network->getBusField(MEASUREMENTS);
  Matrix *YMatrix = map.createMatrix(buses, branches);
  Matrix *AnotherYMatrix = map.createMatrix(branches, buses);
  // Note: YMatrix and AnotherYMatrix should be the same matrix
  Matrix *HMatrix = map.createMatrix(measurements);
  // YMatrix and HMatrix are created using overloaded versions of createMatrix

  // Not yet sure how SESolver would interact with the YMatrix and HMatrix
  // but the solver will depend on them.
  SESolver *solver = new SESolver(network, factory);

  solver->solve();

  // Do something to export results from application
}
