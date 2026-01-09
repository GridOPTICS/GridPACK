/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/**
 * @file   jl_driver.cpp
 * @author Yousu Chen
 * @date   2026-01-09
 *
 * @brief  Implementation of Julia launcher driver
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include "jl_driver.hpp"

namespace gridpack {
namespace julia_launcher {

/**
 * Default constructor
 */
JuliaLauncherDriver::JuliaLauncherDriver(void)
{
    p_pRange = "0.9,1.1";
    p_qRange = "0.9,1.1";
}

/**
 * Constructor with communicator
 */
JuliaLauncherDriver::JuliaLauncherDriver(gridpack::parallel::Communicator comm)
    : p_comm(comm)
{
    p_pRange = "0.9,1.1";
    p_qRange = "0.9,1.1";
}

/**
 * Destructor
 */
JuliaLauncherDriver::~JuliaLauncherDriver(void)
{
}

/**
 * Print usage information
 */
void JuliaLauncherDriver::printUsage(void)
{
    if (p_comm.rank() == 0) {
        std::cout << "Usage: julia_launcher.x [options]" << std::endl;
        std::cout << "Options:" << std::endl;
        std::cout << "  --task_queue=FILE    Path to task_queue.csv (required)" << std::endl;
        std::cout << "  --workdir=DIR        Working directory (required)" << std::endl;
        std::cout << "  --julia_project=DIR  Julia project path (required)" << std::endl;
        std::cout << "  --julia_script=FILE  Path to task_worker.jl (optional)" << std::endl;
        std::cout << "  --instance=NAME      PGLib instance name (optional)" << std::endl;
        std::cout << "  --p_range=A,B        P perturbation range (default: 0.9,1.1)" << std::endl;
        std::cout << "  --q_range=A,B        Q perturbation range (default: 0.9,1.1)" << std::endl;
    }
}

/**
 * Parse command line arguments
 */
void JuliaLauncherDriver::configure(int argc, char **argv)
{
    for (int i = 1; i < argc; i++) {
        std::string arg(argv[i]);

        if (arg.find("--task_queue=") == 0) {
            p_taskQueueFile = arg.substr(13);
        } else if (arg.find("--workdir=") == 0) {
            p_workdir = arg.substr(10);
        } else if (arg.find("--julia_project=") == 0) {
            p_juliaProject = arg.substr(16);
        } else if (arg.find("--julia_script=") == 0) {
            p_juliaScript = arg.substr(15);
        } else if (arg.find("--instance=") == 0) {
            p_instance = arg.substr(11);
        } else if (arg.find("--p_range=") == 0) {
            p_pRange = arg.substr(10);
        } else if (arg.find("--q_range=") == 0) {
            p_qRange = arg.substr(10);
        } else if (arg == "--help" || arg == "-h") {
            printUsage();
            exit(0);
        }
    }

    // Validate required arguments
    if (p_taskQueueFile.empty() || p_workdir.empty() || p_juliaProject.empty()) {
        if (p_comm.rank() == 0) {
            std::cerr << "Error: Missing required arguments" << std::endl;
            printUsage();
        }
        exit(1);
    }

    // Set default julia script path if not specified
    if (p_juliaScript.empty()) {
        p_juliaScript = p_workdir + "/task_worker.jl";
    }
}

/**
 * Read task queue from CSV file
 */
void JuliaLauncherDriver::readTaskQueue(const std::string &filepath)
{
    std::ifstream file(filepath.c_str());
    if (!file.is_open()) {
        if (p_comm.rank() == 0) {
            std::cerr << "Error: Cannot open task queue file: " << filepath << std::endl;
        }
        exit(1);
    }

    std::string line;

    // Skip header line
    std::getline(file, line);

    // Parse task lines
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        Task task;
        std::stringstream ss(line);
        std::string token;

        // Parse CSV: task_id,scenario_id,contingency_type,contingency_id,contingency_name
        std::getline(ss, token, ',');
        task.task_id = std::atoi(token.c_str());

        std::getline(ss, token, ',');
        task.scenario_id = std::atoi(token.c_str());

        std::getline(ss, task.contingency_type, ',');
        std::getline(ss, task.contingency_id, ',');
        std::getline(ss, task.contingency_name, ',');

        // Remove trailing newline/carriage return if present
        if (!task.contingency_name.empty() &&
            (task.contingency_name.back() == '\n' || task.contingency_name.back() == '\r')) {
            task.contingency_name.pop_back();
        }

        p_tasks.push_back(task);
    }

    file.close();
}

/**
 * Spawn a Julia worker to process a single task
 */
bool JuliaLauncherDriver::runJuliaWorker(int task_id)
{
    // Build command
    std::stringstream cmd;
    cmd << "julia --project=" << p_juliaProject << " "
        << p_juliaScript << " "
        << "--task_id=" << task_id << " "
        << "--workdir=" << p_workdir;

    if (!p_instance.empty()) {
        cmd << " --instance=" << p_instance;
    }

    cmd << " --p_range=" << p_pRange
        << " --q_range=" << p_qRange;

    // Execute command
    int ret = std::system(cmd.str().c_str());

    return (ret == 0);
}

/**
 * Execute task distribution and Julia worker spawning
 */
void JuliaLauncherDriver::execute(void)
{
    int rank = p_comm.rank();
    int nprocs = p_comm.size();

    // Read task queue (all ranks read it)
    readTaskQueue(p_taskQueueFile);
    int ntasks = static_cast<int>(p_tasks.size());

    if (rank == 0) {
        std::cout << "========================================" << std::endl;
        std::cout << "GridPACK Julia Launcher" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "Tasks: " << ntasks << std::endl;
        std::cout << "MPI ranks: " << nprocs << std::endl;
        std::cout << "Working directory: " << p_workdir << std::endl;
        std::cout << "Julia project: " << p_juliaProject << std::endl;
        std::cout << "========================================" << std::endl;
    }

    // Synchronize before starting
    p_comm.barrier();

    // Set up task manager for dynamic load balancing
    gridpack::parallel::TaskManager taskmgr(p_comm);
    taskmgr.set(ntasks);

    // Track local results
    std::vector<int> local_task_ids;
    std::vector<int> local_results;

    // Process tasks dynamically
    int task_id;
    int tasks_completed = 0;

    while (taskmgr.nextTask(&task_id)) {
        bool success = runJuliaWorker(task_id);

        local_task_ids.push_back(task_id);
        local_results.push_back(success ? 1 : 0);

        tasks_completed++;

        // Print progress (only from rank that completed the task)
        std::cout << "Rank " << rank << ": Task " << task_id
                  << " (" << p_tasks[task_id].contingency_name << ") - "
                  << (success ? "SUCCESS" : "FAILED") << std::endl;
    }

    // Synchronize after all tasks complete
    p_comm.barrier();

    // Collect results using GlobalVector
    gridpack::parallel::GlobalVector<int> results(p_comm);
    if (!local_task_ids.empty()) {
        results.addElements(local_task_ids, local_results);
    }
    results.upload();

    // Summary (rank 0 only)
    if (rank == 0) {
        std::vector<int> all_ids;
        std::vector<int> all_results;

        for (int i = 0; i < ntasks; i++) {
            all_ids.push_back(i);
        }
        results.getData(all_ids, all_results);

        int n_success = 0;
        int n_failed = 0;
        for (int i = 0; i < ntasks; i++) {
            if (all_results[i] == 1) {
                n_success++;
            } else {
                n_failed++;
            }
        }

        std::cout << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "SUMMARY" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "Total tasks: " << ntasks << std::endl;
        std::cout << "Successful: " << n_success << std::endl;
        std::cout << "Failed: " << n_failed << std::endl;
        std::cout << "Success rate: "
                  << (100.0 * n_success / ntasks) << "%" << std::endl;
        std::cout << "========================================" << std::endl;
    }

    // Print per-processor task statistics
    taskmgr.printStats();
}

} // namespace julia_launcher
} // namespace gridpack
