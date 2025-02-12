// -------------------------------------------------------------
// file: gridpack.cpp
// -------------------------------------------------------------
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created January 24, 2020 by Perkins
// Last Change: 2025-05-06 13:10:00 d3g096
// -------------------------------------------------------------

#include "common.hpp"

#include <gridpack/configuration/no_print.hpp>
#include <gridpack/parallel/task_manager.hpp>
#include <gridpack/timer/coarse_timer.hpp>

// Some temporary hacks

// #define RHEL_OPENMPI_HACK 1
#ifdef RHEL_OPENMPI_HACK

// This stupidity is needed on RHEL7 with stock OpenMPI packages
// installed.

#include <dlfcn.h>


// -------------------------------------------------------------
// stupid_openmpi_hack
// from https://github.com/baidu-research/tensorflow-allreduce/issues/4
// -------------------------------------------------------------
static void 
stupid_openmpi_hack(void)
{
  void *handle = NULL;
  int mode = RTLD_NOW | RTLD_GLOBAL;

  // GNU/Linux and others 
#ifdef RTLD_NOLOAD
      mode |= RTLD_NOLOAD;
#endif
  if (!handle) handle = dlopen("libmpi.so.20", mode);
  if (!handle) handle = dlopen("libmpi.so.12", mode);
  if (!handle) handle = dlopen("libmpi.so.1", mode);
  if (!handle) handle = dlopen("libmpi.so.0", mode);
  if (!handle) handle = dlopen("libmpi.so", mode);
}

#endif

// Hack to return a value from nextTask functions
class TaskCounter
{
public:
  TaskCounter(){};
  ~TaskCounter(){};
  int task_id;
};

// Wrapper class to deal with method pointer arguments (nextTask)
class TaskManagerWrapper
{
public:
  TaskManagerWrapper(gpp::Communicator &comm)
  {
    p_tskmgr.reset(new gpp::TaskManager(comm));
  }
  ~TaskManagerWrapper()
  { }
  void set(int ntask)
  {
    p_tskmgr->set(ntask);
  }
  bool nextTask(TaskCounter &next)
  {
    return p_tskmgr->nextTask(&next.task_id);
  }
  bool nextTask(gpp::Communicator &comm, TaskCounter &next)
  {
    return p_tskmgr->nextTask(comm,&next.task_id);
  }
  void cancel()
  {
    p_tskmgr->cancel();
  }
  void printStats()
  {
    p_tskmgr->printStats();
  }
private:
  boost::shared_ptr<gpp::TaskManager> p_tskmgr;
};

// -------------------------------------------------------------
// get_mpi_comm
// Return a MPI communicator from mpi4py communicator object.
// From https://gitlab.com/robertodr/pybind11-mpi4py/-/blob/main/src/pb11mpi.cpp
// -------------------------------------------------------------
MPI_Comm *get_mpi_comm(py::object py_comm) {
  auto comm_ptr = PyMPIComm_Get(py_comm.ptr());
  
  if (!comm_ptr)
    throw py::error_already_set();
  
  return comm_ptr;
}

// -------------------------------------------------------------
// gridpack module
// -------------------------------------------------------------

// routines to initialize submodules
extern void init_gridpack_ds(py::module& gpm);
extern void init_gridpack_hadrec(py::module& gpm);
extern void init_gridpack_emt(py::module& gpm);
extern void init_gridpack_pf(py::module& gpm);

PYBIND11_MODULE(gridpack, gpm) {
  gpm.doc() = "GridPACK Module";

#ifdef RHEL_OPENMPI_HACK
  stupid_openmpi_hack();
#endif

  // initialize mpi4py's C-API
  if (import_mpi4py() < 0) {
    // mpi4py calls the Python C API
    // we let pybind11 give us the detailed traceback
    throw py::error_already_set();
  }
  
  // -------------------------------------------------------------
  // gridpack.Envronment
  // -------------------------------------------------------------

  py::class_<gp::Environment, boost::shared_ptr<gp::Environment> >
    env(gpm, "Environment");

  env.doc() = "GridPACK parallel environment";

  env
    .def(py::init<>([]()
                    { return boost::shared_ptr<gp::Environment>
                        (new gp::Environment(0, NULL)); }),
      "Create a parallel environment with default MPI communicator")
    
    .def(py::init<>([](py::object py_comm)
                    { MPI_Comm *c = get_mpi_comm(py_comm);
                      return boost::shared_ptr<gp::Environment>
                        (new gp::Environment(0, NULL, *c)); }),
      "Create a parallel environment with specified MPI communicator")
    ;
  

  // -------------------------------------------------------------
  // gridpack.NoPrint
  // -------------------------------------------------------------
  py::class_<gp::NoPrint, std::unique_ptr<gp::NoPrint, py::nodelete> >
    np(gpm, "NoPrint");

  np.doc() = "Device to control output volume";

  np
    .def(py::init([](){
                    return std::unique_ptr<gp::NoPrint, py::nodelete>
                      (gp::NoPrint::instance());
    }), "Create, or get current, NoPrint instance")
    .def("status", &gp::NoPrint::status,
         "Get current status")
    .def("setStatus", &gp::NoPrint::setStatus,
         "Set status")
    ;

  // -------------------------------------------------------------
  // gridpack.Communicator
  // -------------------------------------------------------------
  py::class_<gpp::Communicator> comm(gpm, "Communicator");
  comm.doc() = "GridPACK parallel communicator";

  comm
    .def(py::init<>())
    .def("size", &gpp::Communicator::size)
    .def("rank", &gpp::Communicator::rank)
    .def("worldRank", &gpp::Communicator::worldRank)
    .def("barrier", &gpp::Communicator::barrier)
    .def("sync", &gpp::Communicator::sync)
    .def("divide", &gpp::Communicator::divide)
    .def("split", &gpp::Communicator::split)
    ;

  // -------------------------------------------------------------
  // gridpack.CoarseTimer
  // -------------------------------------------------------------
  py::class_<gpu::CoarseTimer, std::unique_ptr<gpu::CoarseTimer, py::nodelete> >
    timer(gpm, "CoarseTimer");
  timer.doc() = "A general purpose execution timer";

  timer
    .def(py::init([]() {
                    return std::unique_ptr<gpu::CoarseTimer, py::nodelete>
                      (gpu::CoarseTimer::instance());
                  }))
    .def("createCategory", &gpu::CoarseTimer::createCategory)
    .def("start", &gpu::CoarseTimer::start)
    .def("stop", &gpu::CoarseTimer::stop)
    .def("dump", &gpu::CoarseTimer::dump)
    .def("currentTime", &gpu::CoarseTimer::currentTime)
    .def("configTimer", &gpu::CoarseTimer::configTimer)
    .def("dumpProfile",
         static_cast<void ( gpu::CoarseTimer::*)(int) const >(&gpu::CoarseTimer::dumpProfile) )
    .def("dumpProfile",
         static_cast<void ( gpu::CoarseTimer::*)(std::string) const >(&gpu::CoarseTimer::dumpProfile) )
    ;

  // -------------------------------------------------------------
  // gridpack.Configuration
  // -------------------------------------------------------------
  py::class_< ConfigurationCursorWrapper>
    confwrap(gpm, "ConfigurationCursor");

  confwrap.doc() = "Cursor into a specific Configuration path";
  confwrap
    .def("get",
         [](ConfigurationCursorWrapper& self,
            const gpu::Configuration::KeyType& key) -> py::object {
           std::string s;
           if (self.the_cursor->get(key, &s)) {
             py::handle py_s =
               PyUnicode_DecodeLatin1(s.data(), s.length(), nullptr);
             return py::reinterpret_steal<py::str>(py_s);
           }
           return py::cast<py::none>(Py_None);
         })
    .def("getCursor",
         [](ConfigurationCursorWrapper& self,
            const gpu::Configuration::KeyType& key) -> py::object {
           gpu::Configuration::CursorPtr s = self.the_cursor->getCursor(key);
           if (s) {
             ConfigurationCursorWrapper result;
             result.the_cursor = s;
             return py::cast(result);
           } 
           return py::cast<py::none>(Py_None);
         })
    ;
  
  py::class_<gpu::Configuration, std::unique_ptr<gpu::Configuration, py::nodelete>>
    conf (gpm, "Configuration");

  conf.doc() = "Configuration database";

  conf
    .def(py::init([]() {
                    return std::unique_ptr<gpu::Configuration, py::nodelete>
                      (gpu::Configuration::configuration());
                  }))
    .def_property_readonly_static("KeySep",
                                  [](py::object) {return gpu::Configuration::KeySep; })
    .def("open",
         py::overload_cast<const std::string&, gpp::Communicator>
         (&gpu::Configuration::open))
    .def("getCursor",
         [] (gpu::Configuration& self, const std::string& path) {
           ConfigurationCursorWrapper result;
           result.the_cursor = self.getCursor(path);
           return(result);
         })
    .def("get",
         [](gpu::Configuration& self,
            const gpu::Configuration::KeyType& key) -> py::object {
           std::string s;
           if (self.get(key, &s)) {
             py::handle py_s =
               PyUnicode_DecodeLatin1(s.data(), s.length(), nullptr);
             return py::reinterpret_steal<py::str>(py_s);
           }
           return py::cast<py::none>(Py_None);
         })
    ;
    
  // -------------------------------------------------------------
  // gridpack.TaskCounter
  // -------------------------------------------------------------
  py::class_<TaskCounter> tc(gpm, "TaskCounter");
  tc.doc() = "";
  tc
    .def(py::init<>())
    .def_readwrite("task_id", &TaskCounter::task_id)
    ;

  // -------------------------------------------------------------
  // gridpack.TaskManager
  // -------------------------------------------------------------
  py::class_<TaskManagerWrapper>
    tm(gpm, "TaskManager");
  tm.doc() = "";
  tm
    .def(py::init<gpp::Communicator&>(),
         "Create an instance using the specified communicator")
    .def("set", &TaskManagerWrapper::set,
         "Specify total number of tasks and set task manager to zero")
    .def("nextTask",
         (bool (TaskManagerWrapper::*)(TaskCounter&))
         &TaskManagerWrapper::nextTask,
         "Get the next task from the task manager.")
    .def("nextTask",
         (bool (TaskManagerWrapper::*)(gpp::Communicator&, TaskCounter&))
         &TaskManagerWrapper::nextTask,
         "Get the next task from the task manager on a (sub)communicator.")
    .def("cancel", &TaskManagerWrapper::cancel)
    .def("printStats", &TaskManagerWrapper::printStats,
         "Print out statistics on how tasks are distributed on processors")
    ;

  // -------------------------------------------------------------
  // dynamic simulation application module
  // -------------------------------------------------------------
  init_gridpack_ds(gpm);

  // -------------------------------------------------------------
  // hadrec application module
  // -------------------------------------------------------------
  init_gridpack_hadrec(gpm);

  // -------------------------------------------------------------
  // emt application module
  // -------------------------------------------------------------
  init_gridpack_emt(gpm);

  // -------------------------------------------------------------
  // powerflow application module
  // -------------------------------------------------------------
  init_gridpack_pf(gpm);

}
