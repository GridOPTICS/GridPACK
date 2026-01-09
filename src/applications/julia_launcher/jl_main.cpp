/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
/**
 * @file   jl_main.cpp
 * @author Yousu Chen
 * @date   2026-01-09
 *
 * @brief  Main entry point for Julia launcher application
 */

#include "mpi.h"
#include "gridpack/environment/environment.hpp"
#include "jl_driver.hpp"

/**
 * Main entry point for julia_launcher.x
 *
 * Usage:
 *   mpirun -np N julia_launcher.x \
 *       --task_queue=workdir/task_queue.csv \
 *       --workdir=workdir \
 *       --julia_project=/path/to/ExaGrid-DataGen \
 *       --instance=pglib_opf_case24_ieee_rts
 */
int main(int argc, char **argv)
{
    // Initialize GridPACK environment (handles MPI init)
    gridpack::Environment env(argc, argv);

    // Create and run the Julia launcher driver
    gridpack::julia_launcher::JuliaLauncherDriver driver;
    driver.configure(argc, argv);
    driver.execute();

    return 0;
}
