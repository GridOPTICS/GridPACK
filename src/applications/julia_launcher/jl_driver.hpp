/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/**
 * @file   jl_driver.hpp
 * @author Yousu Chen
 * @date   2026-01-09
 *
 * @brief  Driver for launching Julia workers to solve OPF tasks
 *
 * This application uses GridPACK's TaskManager to distribute
 * (scenario, contingency) tasks across MPI ranks, with each rank
 * spawning a Julia worker to solve the AC-OPF problem.
 */

#ifndef _jl_driver_hpp_
#define _jl_driver_hpp_

#include <string>
#include <vector>
#include "gridpack/parallel/communicator.hpp"
#include "gridpack/parallel/task_manager.hpp"
#include "gridpack/parallel/global_vector.hpp"

namespace gridpack {
namespace julia_launcher {

/**
 * Structure representing a single task from the task queue
 */
struct Task {
    int task_id;
    int scenario_id;
    std::string contingency_type;
    std::string contingency_id;
    std::string contingency_name;
};

/**
 * Driver class for Julia-based OPF task execution
 */
class JuliaLauncherDriver {
public:
    /**
     * Default constructor using world communicator
     */
    JuliaLauncherDriver(void);

    /**
     * Constructor with custom communicator
     * @param comm MPI communicator to use
     */
    JuliaLauncherDriver(gridpack::parallel::Communicator comm);

    /**
     * Destructor
     */
    ~JuliaLauncherDriver(void);

    /**
     * Parse command line arguments and configure the driver
     * @param argc argument count
     * @param argv argument values
     */
    void configure(int argc, char **argv);

    /**
     * Execute the task distribution and Julia worker spawning
     */
    void execute(void);

private:
    /**
     * Read task queue from CSV file
     * @param filepath path to task_queue.csv
     */
    void readTaskQueue(const std::string &filepath);

    /**
     * Spawn a Julia worker to process a single task
     * @param task_id ID of the task to process
     * @return true if task completed successfully
     */
    bool runJuliaWorker(int task_id);

    /**
     * Print usage information
     */
    void printUsage(void);

    // Configuration parameters
    std::string p_taskQueueFile;   ///< Path to task_queue.csv
    std::string p_workdir;         ///< Working directory
    std::string p_juliaProject;    ///< Julia project path
    std::string p_juliaScript;     ///< Path to task_worker.jl
    std::string p_instance;        ///< PGLib instance name
    std::string p_pRange;          ///< P perturbation range
    std::string p_qRange;          ///< Q perturbation range

    // Task list
    std::vector<Task> p_tasks;     ///< All tasks from queue

    // MPI communicator
    gridpack::parallel::Communicator p_comm;
};

} // namespace julia_launcher
} // namespace gridpack

#endif
