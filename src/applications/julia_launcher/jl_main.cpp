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
 *
 * Uses same initialization pattern as contingency_analysis/ca_main.cpp
 * to ensure proper MPI/GA/PETSc setup.
 */

#include "mpi.h"
#include <ga.h>
#include <macdecls.h>
#include "gridpack/math/math.hpp"
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
    // Initialize MPI libraries (same pattern as CA app)
    int ierr = MPI_Init(&argc, &argv);

    // Initialize Global Arrays
    GA_Initialize();
    int stack = 200000, heap = 200000;
    MA_init(C_DBL, stack, heap);

    // Initialize Math libraries (PETSc)
    gridpack::math::Initialize(&argc, &argv);

    // Create and run the Julia launcher driver
    gridpack::julia_launcher::JuliaLauncherDriver driver;
    driver.configure(argc, argv);
    driver.execute();

    // Cleanup
    GA_Terminate();
    gridpack::math::Finalize();
    ierr = MPI_Finalize();

    return 0;
}
