# FEnics-RBniCS-examples

Provide personal standard code for FEniCS-RBniCS cases. 


# 🛠️ Installation

1. Install [FEniCS](https://fenics.readthedocs.io/en/latest/installation.html#debian-ubuntu-packages) on Ubuntu via Ubuntu Personal Package Archives (PPA)

```bash
sudo apt-get install --no-install-recommends software-properties-common
sudo add-apt-repository ppa:fenics-packages/fenics
sudo apt-get update
sudo apt-get install --no-install-recommends fenics
```

2. Install [PETSc](https://www.mcs.anl.gov/petsc/download/index.html) and [SLEPc](https://slepc.upv.es/download/)

* Prerequisites

```bash
sudo apt-get install valgrind
sudo apt-get install gfortran
sudo apt-get install python3-distutils
```

* Download **PETSc** and **SLEPc** (should be in /home/xuanquang/Downloads)
* Install **petsc-3.14.1** locally on a user-defined folder

```bash
mkdir /home/gaumap/Packages
cd Packages
tar -xzf /home/gaumap/Downloads/petsc-3.14.1.tar.gz
cd petsc-3.14.1
./configure --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --download-mpich --download-fblaslapack
make all check
```
  * alway aware the version of petsc and slepc for the coincide directories

* Install **slepc-3.14.0** locally on a user-defined folder

```bash
cd /home/gaumap/Packages
tar -xzf /home/gaumap/Downloads/slepc-3.14.0.tar.gz
cd slepc-3.14.0
export SLEPC_DIR=/home/gaumap/Packages/slepc-3.18.0
export PETSC_DIR=/home/gaumap/Packages/petsc
export PETSC_ARCH=arch-linux-c-debug
./configure
make all check
```
clear the evironment, use below
```
SLEPC_DIR
```

 * arch-linux2-gnu-c-debug: double check if exist this folder in petsc folder
	 * It can be  arch-linux-c-debug 
 * xuanquang: this is ubuntu username, you can change to your personal username
 * ```make all check``` in this case will show error. Because it cannot find the installed folder in root. Just ignore this warning.

3. Install **mpi** and **pip3**

```bash
sudo apt install mpi
sudo apt install python3-pip
```

4. Install **petsc4py** and **slepc4py**

```bash
pip3 install petsc4py
pip3 install slepc4py
```

5. Install [RBniCS](https://www.rbnicsproject.org/)

* Go back to **Packages** directory

```bash
git clone https://github.com/RBniCS/RBniCS.git
cd RBniCS
sudo python3 setup.py install
```


<!--stackedit_data:
eyJoaXN0b3J5IjpbLTc0MDY1MTU0MiwtMTE0NTU3MTUyMSwtMT
UyNDI3NDU0MiwxNzczODA5MjU3LC0xNjIxNjc1ODMsLTYwMDAw
MDY5OCwxMTkzMTUxMDE3LC0xNjEyMjM4Njc4LDExOTMxNTEwMT
csNzU1NTUxLC0xODgyNzE4MjMzLC01OTQ5MDAxODddfQ==
-->
