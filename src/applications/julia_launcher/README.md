# Julia Launcher Application

## Overview

This application uses GridPACK's TaskManager to distribute OPF tasks across MPI ranks, with each rank spawning a Julia worker to solve the AC-OPF problem using PowerModels.jl.

## Architecture

```
GridPACK (MPI)                    Julia Workers
┌─────────────────┐              ┌─────────────────┐
│ Rank 0          │              │ task_worker.jl  │
│   TaskManager   │──spawns──────│   PowerModels   │
│   GlobalVector  │              │   Ipopt         │
├─────────────────┤              └─────────────────┘
│ Rank 1          │              ┌─────────────────┐
│   TaskManager   │──spawns──────│ task_worker.jl  │
│   GlobalVector  │              │   PowerModels   │
├─────────────────┤              └─────────────────┘
│ ...             │              │ ...             │
└─────────────────┘              └─────────────────┘
```

## Building

```bash
cd /path/to/GridPACK/build
cmake ..
make julia_launcher.x
```

## Usage

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

## Command Line Options

| Option | Required | Description |
|--------|----------|-------------|
| `--task_queue` | Yes | Path to task_queue.csv |
| `--workdir` | Yes | Working directory with results/ |
| `--julia_project` | Yes | Path to Julia project (ExaGrid-DataGen) |
| `--julia_script` | No | Path to task_worker.jl (default: workdir/task_worker.jl) |
| `--instance` | No | PGLib instance name |
| `--p_range` | No | P perturbation range (default: 0.9,1.1) |
| `--q_range` | No | Q perturbation range (default: 0.9,1.1) |

## Files

| File | Description |
|------|-------------|
| `jl_main.cpp` | Main entry point |
| `jl_driver.hpp` | Driver class header |
| `jl_driver.cpp` | Driver implementation |
| `CMakeLists.txt` | Build configuration |

## Dependencies

- GridPACK parallel components (TaskManager, GlobalVector)
- Julia with PowerModels.jl, Ipopt.jl
- MPI

## Performance Notes

- Each task spawns a new Julia process (~2-5s overhead)
- For better performance, consider persistent Julia workers
- TaskManager provides dynamic load balancing across ranks
