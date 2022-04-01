# Advection
Advection is a program for solving linear advection problems in two and three dimensions.

This repo is a fork of the [original
Advection](https://github.com/ISCDtoolbox/Advection) to be used as a dependency in
[Rodin](https://github.com/cbritopacheco/rodin).

# Building

```bash
git clone --recursive https://github.com/cbritopacheco/Advection
cd Advection
mkdir build && cd build
cmake ..
make -j4
```

# License
Advection is given under the [terms of the GNU Lesser General Public License] (LICENSE.md).
