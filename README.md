# Spherical Harmonic Map

This package includes the prototype code for implementing spherical harmonic maps.

Harmonic maps between genus-0 surfaces are conformal. Thus, this code can be used for computing conformal maps from any genus-0 surface to a unit sphere.

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

The output is three files: star.obj tuette.obj and harmonic.obj. star is the initial map that projects all vertices onto the unit sphere. Tuette map is an intermidiate step. Harmonic map is the final result.

## Reference
Gu, Xianfeng David. Computational conformal geometry. Edited by Shing-Tung Yau. Somerville, Mass, USA: International Press, 2008.

## Contact
Please contact Liang Mi icemiliang@gmail.com for any issues. 
