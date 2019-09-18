# Spherical Harmonic Map

This package includes the prototype code for implementing least squares conformal maps.

![alt text](data/brain.jpg?raw=true "Spherical harmonic map")

## Build

[![Build status](https://ci.appveyor.com/api/projects/status/6nyv0sobm4k0ey2j?svg=true)](https://ci.appveyor.com/project/icemiliang/spherical-harmonic)


In the root directory, run:
```
$ rm -r build
$ mkdir build
$ cd build
$ cmake ..
$ make
```

The program has been tested on Ubuntu 16.04 with g++ 5.4.0.

## Usage
```
./main ../data/brain.obj output
```

## References
Gu, Xianfeng David. Computational conformal geometry. Edited by Shing-Tung Yau. Somerville, Mass, USA: International Press, 2008.

## Contact
Please contact Liang Mi icemiliang@gmail.com for any issues. 
