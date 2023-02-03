## Frequently asked questions

#### Where can I download GridPACK?
Release tarballs are available
[here](https://github.com/GridOPTICS/GridPACK/releases). You can also download
gridpack directly from the Github repository. The latest relatively stable work
is located in the [develop branch](https://github.com/GridOPTICS/GridPACK/tree/develop).

If you obtain GridPACK directly from the Github repository, it is currently
necessary to run the submodule command to download
some files used in the CMake configuration process. After cloning GridPACK, `cd`
into the top level GridPACK directory and type

```
git submodule update --init
```

#### What kind of platforms can I use to run GridPACK?
GridPACK is designed to run on computers using the Linux operating system. It
can also be built fairly easily on Apple Macintosh computers running some
version of Mac OS. Use a terminal on a Mac to access a Linux-like environment.

#### How do I build GridPACK?
Instructions for building and running GridPACK can be found
[here](BASIC_INSTALL.md). This page describes a generic basic install and
contains links to pages providing more information. There are also links
describing builds on specific platforms.

#### Where are examples of applications of GridPACK software
Applications built using the GridPACK frameworks are included in the GridPACK
distributon and can be found in the
[$GRIDPACK/src/applications](../../src/applications) directory. Most
of the applications in this directory are built using GridPACK
[modules](../../src/applications/modules), which
represent fundamental power grid calculations such as power flow and dynamic
simulation. Additional examples can be found in the
[$GRIDPACK/src/applications/examples](../../src/applications/examples)
directory and represent more self-contained examples of applications built with
GridPACK.
