# advect
Advect is a program for solving linear advection problems in two and three dimensions.

#### Installation

##### Automatic installation

1. download the zip archive of Advection or clone this repository:

   ` git clone https://github.com/ICStoolbox/Advection.git `

   navigate to the downloaded directory: 

   ` cd Advection `

2. execute the installation script, which will install the [ICS Commons Library](https://github.com/ICStoolbox/Commons) on your system, along with the Advection library and executable, located in ~/lib/ and ~/bin/:

   ` sh install.sh `

##### Manual installation

1. you will need to install the [ICS Commons Library](https://github.com/ICStoolbox/Commons) on your system. 
Please refer to the instructions provided on the ICS Commons Library page in order to install this library.

2. download the zip archive of Advection or clone this repository:

   ` git clone https://github.com/ICStoolbox/Advection.git `

   navigate to the downloaded directory: 

   ` cd Advection `

   then create build directory and create Makefile
   ```
   mkdir build
   cd build
   cmake ..
   make
   ```

   if no errors are produced, install the binary and library

   ` make install ` 

#### Usage
After compiling nstokes as described above, you should have an executable file in your $HOME/bin directory. If your PATH variable is correctly set to this directory, advect can be called with the following syntax:

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

A full description of all parameters and options that can be specified in the command line or in a parameter file [file.nstokes] can be found in the project [wiki](https://github.com/ICStoolbox/NavierStokes/wiki).

#### Quickstart
You can test the installation and look at examples by entering the [demos](demos) directory and running the program:

    cd demos/2d
    advection test.mesh

that will produce an output that will look like:
```
user:~/code/Advection/demos/2d$ advect test.mesh
 - ADVECT, Release 5.2a, Jan. 29, 2016
   (C) Copyright 2006- , ICS-SU

 - LOADING DATA
    cavity.mesh: 555 vertices, 74 edges, 1034 triangles
    cavity.nstokes: 2 parameters
 - COMPLETED: 0.002s

 ** MODULE NSTOKES: 5.2a
    Matrix and right-hand side assembly
    Solving linear system:
     pressure: res=6.862905e-07, nit=36
     velocity: res=7.920541e-07, nit=26
 ** COMPLETED: 0.120s

 - WRITING DATA
    cavity.solb: 555 data vectors
 - COMPLETED: 0.008s

 ** Cumulative time: 0.130s.
```

#### Authors & contributors
* advect has been initiated by Thi Thu Cuc Bui, Charles Dapogny and Pascal Frey (Université Pierre et Marie Curie).
* Contributors to this project are warmly welcomed. 

#### License
advect is given under the [terms of the GNU Lesser General Public License] (LICENSE.md).
