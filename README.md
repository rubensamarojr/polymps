# EMPS_MESH
Explicit and Weakly Compressible Moving Particle Semi-implicit method with Polygon wall

## Requirements
- [Eigen](http://eigen.tuxfamily.org/)
- [libigl](https://github.com/libigl/libigl)
- [JSON for Modern C++](https://github.com/nlohmann/json)

Eigen, libigl and JSON for Modern C++ are header-only libraries already located in folder **include**.

## MPS input files
1. It is necessary to create a file (extension **.grid**) with:
- First line: **0**
- Second line: **number of particles**
- Next lines in the columns: **material ID** coordinates of particles **X** **Y** **Z** the initial fluid velocities (generally 0.0) **VX** **VY** **VZ** the initial particle pressure (generally 0.0) **P** and presure average (generally 0.0) **PAV**

2. It is necessary to create a file (extension **.json**) and set foldernames and filenames, as well as physical and numerical parameters.

There are some examples in the folder **input**.

## Usage

Code compiled and tested on Linux CentOS 7.

You can build the project in GNU/Linux using the makefile. Follow these steps (CPU version):

Clone this repository into your system `git clone https://github.com/rubensaaj/EMPS_MESH.git`

In a terminal, go to the folder **EMPS_MESH**.

Edit the `Makefile` file (if necessary) with a text editor.

Make sure the environment is clean and ready to compile
```bash
make clean
```
Execute make
```bash
make all
```
Run the code as
```bash
./bin/main
```
Type the name of the json input file (located in input directory)

This code writes pvd (header file) and corresponding vtu files as output.
You can visualize them by open the pvd file with [Paraview](https://www.paraview.org).
