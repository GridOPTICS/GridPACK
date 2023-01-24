## Ubuntu Personal Package Archives (PPA) Installation
GridPACK is relatively easy to install on some older versions of the
[Ubuntu Linux](https://www.ubuntu.com/)
system. Working packages for most [GridPACK
prerequisites](../required/PREREQUISITES.md)
are available from standard repositories. If you do not
need to modify the GridPACK source, e.g. just use the GridPACK example
applications or build your own application, you can install the binary package
from the
[GridPACK PPA](https://launchpad.net/~wperkins/+archive/ubuntu/gridpack-ppa).

GridPACK packages are available for the long term support (LTS) 64-bit (AMD64)
Ubuntu series 16.04 (xenial) and 18.04 (bionic).  If you are installing Ubuntu
on a system or virtual machine, download the
[16.04](http://releases.ubuntu.com/16.04)
or [18.04](http://releases.ubuntu.com/18.04) LTS Desktop AMD64
distribution and follow the
[installation
instructions](https://tutorials.ubuntu.com/tutorial/tutorial-install-ubuntu-desktop#0).
GridPACK packages install a version of Global
Arrays (GA) specifically for GridPACK use.  This will conflict with applications
that depend on the Ubuntu GA packages, but this should be a rarity for GridPACK
users. 

*You will need super user or sudo privileges for this installation.*

## GridPACK Installation on Ubuntu Linux 16.04 or 18.04 LTS

Add the PPA to your system, and install GridPACK with

```
sudo add-apt-repository ppa:wperkins/gridpack-ppa
sudo apt-get update
sudo apt-get install gridpack-dev
```

Make sure your installation works by building and running one of the several
example applications, e.g. power flow:

```
mkdir tmpbuild
cd tmpbuild
cmake /usr/share/gridpack/example/powerflow
make
mpiexec -np 2 ./powerflow.x 118.xml
```

Also, you can run the installed powerflow application on this same input:

```
mpiexec -np 2 /usr/bin/pf.x 118.xml
```

Other GridPACK applications that can be found in `/usr/bin` are
contingency analysis (`ca.x`), dynamic simulation (`dsf.x`),
kalman filter (`kds.x`), and state estimation (`stes.x`).
Additional example inputs for the GridPACK applications can be found under
`/usr/share/gridpack/example`. The GridPACK libraries and include files
are all found under `/usr`, so specify this directory as your GridPACK
installation directory when compiling you own applications.

GridPACK can be removed with

```
sudo apt-get purge gridpack-dev
sudo apt autoremove
```

