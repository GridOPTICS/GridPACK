# pygridpack examples

Five self-contained demos, the first four one application each.  They are meant to
be read: each takes the shortest honest path from an input file to a number
you can check, with no framework of its own.

| Demo | Application | Case | Shows |
|---|---|---|---|
| `01_power_flow.py` | power flow | IEEE 14 (`IEEE14_PTIv33_rated.raw`: every branch rated) | solve, read the solved network, report limit violations |
| `02_contingency_screening.py` | contingency analysis | IEEE 118 | 179 N-1 line outages ranked by what they break |
| `03_dynamic_simulation.py` | dynamic simulation | 9-bus, 3-machine | XML-declared fault, whole run in C++, series from the watch file |
| `04_stepped_event.py` | dynamic simulation | Kundur two-area | Python drives the step loop and trips a generator at run time |
| `05_composed_study.py` | all three above | IEEE 14 | base case, N-1 screen, worst outage opened in a stepped run |

## Running them

```sh
pip install -e python/[all]          # needs GRIDPACK_DIR set to a GridPACK install
cd python/examples
python 01_power_flow.py
```

Each demo copies its inputs into `$TMPDIR/gridpack-demo-<case>/`, works there,
and prints the path of anything it writes.  Inputs come from the GridPACK
source tree; set `GRIDPACK_TEST_DATA` to that tree if you installed without it.

## Running them in parallel

```sh
mpiexec -np 4 python 02_contingency_screening.py
```

How many ranks a demo can take is set by the case, not by your machine:

| Demo | Verified | Why |
|---|---|---|
| `01_power_flow.py` | `-np 4` | |
| `02_contingency_screening.py` | `-np 4` | contingencies are handed out across ranks, so more ranks finishes the screen sooner |
| `03_dynamic_simulation.py` | `-np 3` | a 9-bus network will not partition across four ranks |
| `04_stepped_event.py` | serial | see below |
| `05_composed_study.py` | `-np 2`, screen only | as `04`; in parallel it stops after the screen |

**A network too small for the rank count hangs.**  `03` at `-np 4` never
finishes -- both ranks spin at 100% CPU and the run has to be killed.  The
stock `dsf.x` driver does the same thing on the same input at `-np 4`, and
both finish normally at `-np 3`, so this is a GridPACK property rather than a
pygridpack one.  Keep a few buses of work per rank and it does not arise; the
118-bus screen in `02` is comfortable at four.

## Which inputs carry observations

`04` reads its time series from `<observations>` in the XML.  Only the four
Kundur two-area inputs under `data_sets/input/ds/` declare that block; without
it the stepper's channel list comes back empty and `to_csv()` writes only a
time column.  `03` sidesteps this by reading the generator watch file, which
`<generatorWatch>` controls instead.
