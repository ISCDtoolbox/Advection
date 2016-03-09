# advect [![Build Status](https://travis-ci.org/ICStoolbox/Advection.svg?branch=master)](https://travis-ci.org/ICStoolbox/Advection)
Advect is a program for solving linear advection problems in two and three dimensions.

#### Installation

1. Install the [ICS Commons Library](https://github.com/ICStoolbox/Commons) on your system. 
Please refer to the instructions provided on the ICS Commons Library page in order to install this library.

2. download the zip archive of Advection or clone this repository:

   ` git clone https://github.com/ICStoolbox/Advection.git `

   navigate to the downloaded directory: 

   ` cd Advection `

   create a build directory and compile with cmake
   ```
   mkdir build
   cd build
   cmake ..
   make
   make install
   ```

#### Usage
After compiling advect as described above, you should have an executable file in your $HOME/bin directory. If your PATH variable is correctly set to this directory, advect can be called with the following syntax:

    usage: advect [+/-v | -h] [-dt step] source[.mesh] [-s data[.sol]] [-o output[.sol]]
    
The square braces indicate optional arguments. Some commands have flags, some others do not.

The options and flags are:
```
  --help       show the syntax and exit.
  --version    show the version and date of release and exit.

  -dt step     time step (time units)
  -v           suppress any message (for use with function call).
  +v           increase the verbosity level for output.

  source.mesh    name of the mesh file
  data.sol       name of file containing the initial solution or boundary conditions
  output.sol     name of the output file
```

A full description of all parameters and options that can be specified in the command line or in a parameter file [file.advect] can be found in the project [wiki](https://github.com/ICStoolbox/Advection/wiki) (Coming soon...).

#### Quickstart (Coming soon...)
You can test the installation and look at examples by entering the [demos](demos) directory and running the program:

    cd demos/2d
    advection test.mesh -dt 0.001 -s test.sol -c test.chi.sol -o test.chi.sol

that will produce an output that will look like:
```
user:~/code/Advection/demos/2d$ advect test.mesh -dt 0.001 -s test.sol -c test.chi.sol -o test.chi.sol
 - ADVECT, Release 3.0a, Feb. 19, 2016
   (C) Copyright 2007- , ICS-SU

 - LOADING DATA
    test.mesh: 647 vertices, 1213 triangles
    test.sol : 647 data vectors
    test.chi.sol : 647 data scalar
    Adjacency table:  3639 updated
 - COMPLETED: 0.003s

 ** MODULE ADVECT: 3.0a
    Time stepping: 0.005
    Solving: 0 characteristics
 ** COMPLETED: 0.000s

 - WRITING DATA
    test.chi.sol: 647 data vectors
 - COMPLETED: 0.001s

 ** Cumulative time: 0.004s.
```

#### Authors & contributors
* advect has been initiated by Thi Thu Cuc Bui, Charles Dapogny and Pascal Frey (Université Pierre et Marie Curie).
* Contributors to this project are warmly welcomed. 

#### License
advect is given under the [terms of the GNU Lesser General Public License] (LICENSE.md).
