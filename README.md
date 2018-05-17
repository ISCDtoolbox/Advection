# advect [![Build Status](https://travis-ci.org/ISCDtoolbox/Advection.svg?branch=master)](https://travis-ci.org/ISCDtoolbox/Advection)
Advect is a program for solving linear advection problems in two and three dimensions.

#### Installation

1. install the [ISCD Commons Library](https://github.com/ISCDtoolbox/Commons) on your system. 
Please refer to the instructions provided on the ISCD Commons Library page in order to install this library.

2. download the zip archive of Advection or clone this repository:

   ` git clone https://github.com/ISCDtoolbox/Advection.git `

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

    usage: advect [+/-v | -h] [-dt step] source[.mesh] [-c function[.sol]] [-s data[.sol]] [-o output[.sol]]
    
The square braces indicate optional arguments. Some commands have flags, some others do not.

The options and flags are:
```
  --help       show the syntax and exit.
  --version    show the version and date of release and exit.

  -dt step     time step (time units)
  -nocfl       avoid truncation of the time period for advection due to cfl condition.
  -noex        deactivates the default feature whereby characteristic lines are extrapolated outside the domain when the input velocity field causes them to do so.  
  -v           suppress any message (for use with function call).
  +v           increase the verbosity level for output.

  source.mesh    name of the mesh file
  function.sol   name of file containing the (scalar) values to be advected
  data.sol       name of file containing the velocity field
  output.sol     name of the output file containing the function values
```

A full description of all parameters and options that can be specified in the command line or in a parameter file [file.advect] can be found in the project [wiki](https://github.com/ISCDtoolbox/Advection/wiki) (Coming soon...).

#### Quickstart (Coming soon...)
You can test the installation and look at examples by entering the [demos](demos) directory and running the program:

    cd demos
    advect test.mesh -dt 0.01 -s test.sol -c test.chi.sol -o test.chi.sol

that will produce an output that will look like:
```
user:~/code/Advection/demos/2d$ advect test.mesh -dt 0.01 -s test.sol -c test.chi.sol -o test.chi.sol
 - ADVECT, Release 3.0a, Feb. 19, 2016
   (C) Copyright 2007- , ICS-SU

 - LOADING DATA
    test.mesh: 2225 vertices, 4222 triangles
    test.sol : 2225 data vectors
    test.chi.sol : 2225 data scalar
    Adjacency table:  12666 updated
 - COMPLETED: 0.013s

 ** MODULE ADVECT: 3.0a
    Time stepping: 0.001
    Solving: 2225 characteristics
 ** COMPLETED: 0.009s

 - WRITING DATA
    test.chi.sol: 2225 data vectors
 - COMPLETED: 0.010s

 ** Cumulative time: 0.032s.
```

#### Authors & contributors
* advect has been initiated by Thi Thu Cuc Bui, Charles Dapogny and Pascal Frey (Université Pierre et Marie Curie).
* Contributors to this project are warmly welcomed. 

#### License
advect is given under the [terms of the GNU Lesser General Public License] (LICENSE.md).
