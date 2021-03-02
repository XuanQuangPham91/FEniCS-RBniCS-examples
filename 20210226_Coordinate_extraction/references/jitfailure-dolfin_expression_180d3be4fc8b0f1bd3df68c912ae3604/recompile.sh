#!/bin/bash
# Execute this file to recompile locally
c++ -Wall -shared -fPIC -std=c++11 -O3 -fno-math-errno -fno-trapping-math -ffinite-math-only -I/usr/lib/slepcdir/slepc3.12/x86_64-linux-gnu-real/include -I/usr/lib/petscdir/petsc3.12/x86_64-linux-gnu-real/include -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/include/hdf5/openmpi -I/usr/include/eigen3 -I/usr/lib/python3/dist-packages/ffc/backends/ufc -I/home/xuanquang/.cache/dijitso/include dolfin_expression_180d3be4fc8b0f1bd3df68c912ae3604.cpp -L/usr/lib/x86_64-linux-gnu/openmpi/lib -L/usr/lib/petscdir/petsc3.12/x86_64-linux-gnu-real/lib -L/usr/lib/slepcdir/slepc3.12/x86_64-linux-gnu-real/lib -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -L/home/xuanquang/.cache/dijitso/lib -Wl,-rpath,/home/xuanquang/.cache/dijitso/lib -lmpi -lmpi_cxx -lpetsc_real -lslepc_real -lm -ldl -lz -lsz -lhdf5 -lboost_timer -ldolfin -olibdijitso-dolfin_expression_180d3be4fc8b0f1bd3df68c912ae3604.so