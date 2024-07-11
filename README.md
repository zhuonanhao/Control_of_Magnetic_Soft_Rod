# Control_of_Magnetic_Soft_Rod


Compile and build:
------------------

Instructions for Ubuntu:
(1) To run this code you need Eigen, OpenGL, MKL, and Lapack. Lapack is usually preinstalled on your computer.
Eigen can be found at http://eigen.tuxfamily.org/index.php?title=Main_Page

(2) Create a file named "Makefile". The content of the "Makefile" should be the same "Makefile_Sample" except that you will need to change the path to eigen from "/usr/local/include/eigen3/" to the location of eigen in your system.

(3) Open a terminal, "cd" to this folder and run the command "make" (without the quotes).

(4) To start the simulation, run the command "./simDER option.txt" (without the quotes). More on option.txt later.
