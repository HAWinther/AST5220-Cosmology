# AST5220-Cosmology

This repository contains C++ templates for making an Einstein-Boltzmann solver (a CAMB like code). This is used for the course AST5220 "Cosmology II" at ITA Univeristy of Oslo. The aim of this course is for the students to learn how to do cosmology in both theory and practice by making their own Einstein-Boltzmann solver.

---

For an introduction to the project, C++ and the tools that we provide see http://folk.uio.no/hansw/AST5220/notes/about.html

For the first milestone see http://folk.uio.no/hansw/AST5220/notes/milestone1.html

---

Compile the code running [ make ]. If you want to compile this on your computer you need to install the GSL library first. See below for instructions.

If you get it compiled then run it as ( ./cmb ) It will crash with "Spline eta has not been created". That is fine, its one of your task to implement this.

See Examples.cpp, run the examples as ( make examples; ./examples ; ), for examples on how to make splines, solve ODEs, etc.
and the functionality of the stuff supplied with this template

The code runs from Main.cpp. See this file for the order of things to implement.

In the last module you will also need a library for Bessel functions (the GSL one often fails for very large inputs), if 
COMPLEX\_BESSEL is defined in the makefile then you will have to have this library installed:

https://github.com/joeydumont/complex\_bessel

If you don't have this then just comment out the two lines below ( "Add bessel function library" ) in the Makefile - its only needed in the end.

---
INSTALL GSL LOCALLY:
---

Run the following commands in order:

- Go the the home directory:

cd $HOME

- Make a local folder to hold libraries:

mkdir local

- Enter this directory:

cd local

- Download the code (if you don't have wget you need to get the file to this dir by other means):

wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.6.tar.gz

- Unzip the code:

tar -xvf gsl-2.6.tar.gz

- You should now have the gsl-2.6 folder. Enter it:

cd gsl-2.6

- Run the configure script:

./configure --prefix=$HOME/local

- Compile and install it:

make ; make install

- In the CMB code Makefile change the include and lib paths to point to the library:

INC  = -I$(HOME)/local/include
LIBS = -L$(HOME)/local/lib -lgsl -lgslcblas

- If this fails with "libgsl.so not found" then run the command:

export LD\_LIBRARY\_PATH="$LD\_LIBRARY\_PATH:$HOME/local/lib"

and try to run ./cmb again and it should work. To avoid having
to run this command every time you open a new terminal open
the $HOME/.bashrc file and add this line to the end of the file
and it will load everytime you open a new window

---
INSTALL COMPLEX BESSEL LOCALLY:
---

- Download the source

cd $HOME

mkdir local

cd local

git clone https://github.com/joeydumont/complex\_bessel

cd complex\_bessel

See the README.md file for how to proceed

