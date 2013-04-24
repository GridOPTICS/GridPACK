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
