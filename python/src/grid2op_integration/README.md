# GridPACK ↔ Grid2Op Integration (Python)

This package exposes a Grid2Op-compatible GridPACK backend to train and evaluate agents against detailed grid dynamics.

## 1. Prerequisites

- **OS/Toolchain**:  
  - Linux (tested on CentOS 8 / RHEL 8)  
  - **GCC 11.2.0**  
  - **OpenMPI 4.1.4**  
  - **CMake 3.29.0**  
  - Make (from GNU coreutils, ≥4.3)

- **Python**: 3.9–3.11 recommended (use a dedicated virtual environment or conda environment).  

- [**GridPACK**](https://grid2op.readthedocs.io/):
  1. Load required modules (see below).  
  2. Install dependencies:
     ```bash
     source install_gridpack_deps.sh
     ```
  3. Build GridPACK core libraries (no Python):
     ```bash
     # Set install_gridpack=true and install_gridpack_python=false
     source install_gridpack.sh
     ```
  4. Build GridPACK Python bindings:
     ```bash
     # Set install_gridpack=false and install_gridpack_python=true
     source install_gridpack.sh
     ```
  5. Verify that the `gridpack` Python module can be imported:
     ```bash
     python -c "import gridpack; print('GridPACK OK')"
     ```  

- [**Grid2Op**](https://github.com/rte-france/Grid2Op) (tested with v1.12.1):
  1. Install into your virtual environment:
  ```bash
  pip install grid2op
  ```
  2. Verify installation:
  ```bash
  python -c "import grid2op; print('Grid2Op OK, version:', grid2op.__version__)"
  ```
## 2. File reference

- **`fidvr_study.py`**  
  Primary entry point for running experiments on Fault-Induced Delayed Voltage Recovery (FIDVR).  
  Uses `config.py` and the GridPACK–Grid2Op integration to set up environments, train/evaluate  
  agents, and log results in different modes (`simulate`, `train`, `deploy`).

- **`generate_results.py`**  
  Post-processing and analysis tool. Loads logs and model outputs from training/deployment runs,  
  computes metrics, and generates summary results and plots.

- **`build_prods_charac.py`**  
  Utility for building and characterizing production datasets. Often used to prepare case-specific  
  data (e.g., load characteristics or baseline trajectories) needed for experiments.


## 3. Running FIDVR Experiments

Fault-Induced Delayed Voltage Recovery (FIDVR) occurs when post-fault dynamics (especially from induction motors and load recovery)  
cause bus voltages to remain depressed for several seconds even after a fault is cleared.  
Understanding and mitigating this phenomenon is a critical application of learning-based  
control in power systems. With this integration:

- **GridPACK** provides the physics-based dynamic simulation of faults and load recovery.  
- **Grid2Op** provides the agent interface (observations, actions, rewards) for testing  
  different control strategies.

The main driver script, `fidvr_study.py`, connects these  
two components together through the configuration settings given in `config.py`. The script sets up experiments where agents attempt to stabilize the system following simulated faults.

### 3.1 Command-line arguments

- `--grid-xml` (str)  
  GridPACK XML case to use. Choices:  
  `input_9bus.xml`, `input_9bus_fault.xml`,  
  `input_39bus_IBR.xml`, `input_39bus_IBR_fault.xml`,  
  `input_145.xml`, `input_145_fault.xml`  
  **Default:** `input_9bus.xml`

- `--mode` (str)  
  Execution mode: `simulate`, `train`, or `deploy`.  
  **Default:** `simulate`

- `--outdir` (str)  
  Directory where logs, models, and CSV outputs are written.  
  **Default:** `runs/dqn_grid2op`

- `--episodes` (int)  
  Number of episodes to run in `simulate` or `deploy` mode.  
  **Default:** `1`

- `--max-steps` (int)  
  Maximum number of steps per episode (for `simulate`/`deploy`).  
  **Default:** `5`

- `--total-steps` (int)  
  Total number of environment steps to use for training.  
  **Default:** `TOTAL_STEPS` (value defined in code)

- `--model-path` (str)  
  Path to a trained model file (used in `deploy` mode).  
  **Default:** `./models/dqn_grid2op.pt`

- `--log-sim-during-train` (flag)  
  If set, logs per-step simulation data during training for analysis.  
  **Default:** `False`

- `--seed` (int)  
  Random seed for reproducibility.  
  **Default:** `123`

### 3.2 Usage

Run `fidvr_study.py` in one of the following **modes** (set with `--mode`):

- **simulate**: run the environment with a no-op (do-nothing) policy  
- **train**: learn a DQN agent from interactions with the environment  
- **deploy**: load a trained model and evaluate its performance  

#### Example runs

**Simulate** (environment only, no learning):
```bash
python fidvr_study.py --mode simulate --max-steps 10 --grid-xml input_9bus_fault.xml
```

**Train**  (learn a DQN agent):
```bash
python fidvr_study.py --mode train --total-steps 60000 --log-sim-during-train --max-steps 60 --grid-xml input_9bus_fault.xml
```

**Deploy**  (evaluate a trained DQN policy):
```bash
python fidvr_study.py --mode deploy --max-steps 10 --grid-xml input_9bus_fault.xml --episodes 10 --model-path ./models/dqn_grid2op.pt
```

### 3.3 What happens under the hood?

1. `config.py` loads integration defaults.  
2. The XML is parsed to locate its `.dyr` and `.raw` files.  
3. The matching `grid2op_<FILENAME>/` folder is located automatically.  
4. A closed-loop experiment is launched where the agent must respond to faults and accelerate voltage recovery."

