# GridPACK

GridPACK is a software framework consisting of a set of modules
designed to simplify the development of programs that model the power
grid and run on parallel, high performance computing platforms. The
modules are available as a library and consist of components for
setting up and distributing power grid networks, support for modeling
the behavior of individual buses and branches in the network,
converting the network models to the corresponding algebraic
equations, and parallel routines for manipulating and solving large
algebraic systems. Additional modules support input and output as well
as basic profiling and error management.  

<!---See the [GridPACK home page](https://www.gridpack.org) for more information.-->

## Installing

GridPACK can be built in a number of different ways. The `docs/notes` subdirectory has instructions for building GridPACK on different clusters. Also, it includes a [step-by-step guide](docs/notes/install.md) for building, installing GridPACK and its dependencies.

## Usage
See [User manual](https://www.gridpack.org/wiki/images/9/9a/GridPACK_User_Manual_3.4.pdf) for a deep dive on GridPACK internals. One can also use GridPACK through its Python interface.

## Applications
GridPACK includes a number of different power system applications. The two most commonly used and well-developed are:
- AC Power Flow
- Dynamics Simulation
- Contingency Analysis

Other applications, that are in development or not full featured are
- Dynamic security assessment
- State estimation

## Authors
- Bruce Palmer
- William Perkins
- Yousu Chen
- Renke Huang
- Yuan Liu
- Shuangshuang Jin
- Shrirang Abhyankar

## Acknowledgement
GridPACK has been developed through funding from various sources over the years.
- PNNL LDRD Future Grid Initiative
- DOE OE [Advanced Grid Modeling (AGM)](https://www.energy.gov/oe/advanced-grid-modeling) program
- [Grid Modernization Laboratory Consortium](https://www.energy.gov/gmi/grid-modernization-lab-consortium)
- DOE EERE [Solar Energy Technologies Office](https://www.energy.gov/eere/solar/solar-energy-technologies-office)
- DOE EERE [Wind Energy Technologies Office](https://www.energy.gov/eere/wind/wind-energy-technologies-office)

## Copyright
Copyright &copy; 2017, Battelle Memorial Institute.

GridPACK<sup>TM</sup> is a free software distributed under a BSD 2-clause license. You may reuse, modify, and redistribute the software. See the [license](LICENSE) file for details.


## Disclaimer
This material was prepared as an account of work sponsored by an agency of the United States Government.  Neither the United States Government nor the United States Department of Energy, nor Battelle, nor any of their employees, nor any jurisdiction or organization that has cooperated in the development of these materials, makes any warranty, express or implied, or assumes any legal liability or responsibility for the accuracy, completeness, or usefulness or any information, apparatus, product, software, or process disclosed, or represents that its use would not infringe privately owned rights.
Reference herein to any specific commercial product, process, or service by trade name, trademark, manufacturer, or otherwise does not necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or any agency thereof, or Battelle Memorial Institute. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or any agency thereof.