# Cosmology II
## The Cosmic Microwave Background and the Large Scale Structure of our Universe in theory and practice

This repository contains C++ templates for making an Einstein-Boltzmann solver (a CAMB like code). This is used for the course AST5220 "Cosmology II" at ITA University of Oslo. The aim of this course is for students to learn both the theory, the physics and the numerics. We derive all the equations and discuss the physics in the lectures and then the students have to implement and solve them in a numerical code that will ultimately lead to matter and CMB power spectra. 

For the current version of the course C++ is the main language and you are strongly reccomended to use this as I can much more easily help you out and most of the course website only contains information for this C++ template. However you are free to do this project in any language you want (and you don't even have to use the template we provide if you don't want to), just be aware that I'm only good at C, C++, Fortran and Python so if you do it any other language you are on your own. See below for templates in other languages: 

# Fortran

Templates for how to do this in Fortran (written by Hans Kristian Eriksen) can be found in the Fortran\_template directory together with PDFs outlining what to do for each milestone.

# Python?

Every year there is someone who absolutely wants to do it in Python and very often it ends up not going so well (though some people manage to do it extremely well). I have therefore also included some simple templates for how you could do this in Python and you can find these in python\_template. However its strongly reccomended that you do not do this in Python even if that is what you know best. First of all its extremely useful to learn a more low level language which will make you a better human being and secondly the last two milestones are going to be super slow unless you really know what you are doing and how to speed up slow parts in Python (and is willing to spend time doing this). If you just naively implement things its likely going to take 10-20 minutes for a proper run of the full code (which will be very challenging to debug). Having said all that, this is not kindergarden so make your own choice.

# Website
All relevant information about the project and the different milestones can be found on this [website](https://cmb.wintherscoming.no/).

# Compiling

Compile the code running [ make ]. If you want to compile this on your computer you need to install the [GSL library](ftp://ftp.gnu.org/gnu/gsl/) first. See below for instructions if you haven't installed a library lik this before.

The code runs from Main.cpp and then proceeds to go through the different milestones one by one untill we have the CMB power spectra in the end. If you get it compiled then run it as [ ./cmb ]. It will crash with "Error: Spline eta has not been created". That is fine, it's one of your task to implement this. 

See Examples.cpp - and run the examples as [ make examples; ./examples ; ] - for examples on how to make splines, solve ODEs and the functionality of the stuff supplied with this template.

For the last milestone you need to compute spherical bessel functions (see the function j\_ell(n, x) in Utils.cpp). There are several options for this: if you have a C++17 compiler (use -std=c++17 instead of c++11 in the Makefile) then you can use provided by the C++ standard std::sph\_bessel(n, x). Otherwise the GSL library provides a function gsl\_sf\_bessel\_jl(n, x) for this that we use in the template. This implementation has problems with underflow for small x and large n, but we correct this using known asymptotical expressions. The last option is to use another library like for example the [Complex Bessel](https://github.com/joeydumont/complex_bessel) library.

# How to install GSL

See [this](https://solarianprogrammer.com/) for how to install it on a Windows machine. On Linux or a Mac you can either use a package manager or install it directly as follows:

- Go the the home directory:

cd $HOME

- Make a local folder to hold libraries:

mkdir local

- Enter this directory:

cd local

- Download the code (if you don't have wget you need to get the file to this dir by other means):

wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.6.tar.gz

- Untar the code:

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
and it will load everytime you open a new window.

