# GridPACK Julia Launcher Application

**Author:** Yousu Chen
**Date:** January 2026

> **Note:** This README covers local development setup. Perlmutter/HPC deployment instructions will be added later.

## Overview

This application uses GridPACK's TaskManager to distribute OPF tasks across MPI ranks, with each rank spawning a Julia worker to solve the AC-OPF problem using PowerModels.jl. 

Features:
- TaskManager provides dynamic load balancing across ranks
- Each task spawns a new Julia process (~2-5s overhead)

## Dependencies

- GridPACK parallel components (TaskManager, GlobalVector): https://github.com/GridOPTICS/GridPACK
- Julia with PowerModels.jl: https://github.com/lanl-ansi/PowerModels.jl
- ExaGrid-DataGen: https://github.com/exanauts/ExaGrid-DataGen

## Architecture

```
GridPACK (MPI)                    Julia Workers
  Rank 0                           task_worker.jl 
    TaskManager    ──spawns───       PowerModels 
    GlobalVector                     Ipopt/other optimizers      
  ...                              ...         
  Rank N                           task_worker.jl 
    TaskManager    ──spawns────      PowerModels 
    GlobalVector                     Ipopt/other optimizers
  ...                              ...          
```

## Building

```bash
cd /path/to/GridPACK/build
cmake ..
make julia_launcher.x
```

## Usage, using pglib_opf_case24_ieee_rts as an example

### Step 1: Generate Task Queue (Julia)

```bash
julia --project=/path/to/ExaGrid-DataGen \
    /path/to/ExaGrid-DataGen/gridpack_integration/task_queue.jl \
    --instance=pglib_opf_case24_ieee_rts \
    --n_scenarios=10 \
    --output_dir=workdir
```

### Step 2: Run Julia Launcher (GridPACK)

```bash
mpirun -np 8 julia_launcher.x \
    --task_queue=workdir/task_queue.csv \
    --workdir=workdir \
    --julia_project=/path/to/ExaGrid-DataGen \
    --instance=pglib_opf_case24_ieee_rts
```

### Step 3: Aggregate Results (Julia)

```bash
julia --project=/path/to/ExaGrid-DataGen \
    /path/to/ExaGrid-DataGen/gridpack_integration/result_aggregator.jl \
    --workdir=workdir \
    --output_dir=final_output
```



