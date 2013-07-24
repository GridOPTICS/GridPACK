// -------------------------------------------------------------
// Powerflow::solve()
// -------------------------------------------------------------
void
Powerflow::solve()
{

  // Before we get here, the admittance matrix and power/current
  // injection vector have been correctly formed, including the
  // selection of the slack bus.  There is also an initial estimate. 

  math::Matrix *Ybus;
  math::Vector *Cinj;
  math::Vector *Vest;

  // Or, the smart (but not shared pointer) of your choice
  boost::scoped_ptr<PowerflowSolver> solver;

  std::string solver_option(this->get_param<std::string>("solver"));
  
  // This should actually be coded as a factory:
  // solver.reset(PowerflowSolverFactory(solver_option));

  if (solver_option == "linear") {
    solver.reset(new LinearPowerflowSolver(*Ybus, *Cinj, *Vest));
  } else if (solver_option == "nonlinear") {
    solver.reset(new NonlinearPowerflowSolver(*Ybus, *Cinj, *Vest));
  } else if (solver_option == "newtonraphson") {
    solver.reset(new NetwonRaphsonPowerflowSolver(*Ybus, *Cinj, *Vest));
  } else if (solver_option == "...") {
    // whatever
  } else {
    solver.reset(new ReliablePowerflowSolver(*Ybus, *Cinj, *Vest));
  }

  solver->solve();

  math::Vector *V;
  solver->solution(*V);

  // compute voltage magnitudes and phase angle

  // compute branch currents / power

  // put values back on network
  
}
