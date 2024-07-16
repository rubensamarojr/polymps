![logo_PolyMPS](logo_polymps.png)

A C++ code for numerical modelling of free-surface flow :ocean: using Explicit, Weakly Compressible or Incompressible Moving Particle Simulation/Semi-implicit [MPS](https://doi.org/10.13182/NSE96-A24205) method. Boundary walls can be modeled by polygonal mesh (triangles) or discrete layers of wall and dummy (ghost) particles.

<img src='/output/3D_dam_1610_h300_lo0p0100_pnd2_adj2_LNJ_SITE.gif' width="96%">

[**Quick start**](#contents) | [**Examples**](https://github.com/rubensamarojr/polymps/tree/master/input) | [**Paper**](https://doi.org/10.1016/j.simpa.2022.100376) | [**Citation**](#how-to-cite-polymps)

## Forum

On the [PolyMPS Discussions](https://github.com/rubensamarojr/polymps/discussions) page you can ask questions, discuss about simulation issues, share ideas, and interact with other members.

## Contents
- [Dependencies](#dependencies)
- [PolyMPS workflow](#polymps-workflow)
- [Input files](#input-files)
    - [Foldernames, Filenames, Physical and Numerical parameters](#foldernames-filenames-physical-and-numerical-parameters)
    - [Solid domain](#solid-domain)
    - [Fluid domain](#fluid-domain)
- [Compile](#compile)
- [Run](#run)
    - [on LINUX](#on-linux)
    - [on WINDOWS](#on-windows)
- [Json input filename](#json-input-filename)
- [Output](#output)
- [Directories](#directories)
- [License](#license)

## Dependencies

- [GCC (GNU Compiler Collection)](https://gcc.gnu.org)
- [Eigen](http://eigen.tuxfamily.org)
- [libigl](https://github.com/libigl/libigl)
- [JSON for Modern C++](https://github.com/nlohmann/json)

In order to install C++ compiler (GCC) on **windows**, we recommend to install [Cygwin](https://cygwin.com). You can find [here](https://www3.ntu.edu.sg/home/ehchua/programming/howto/Cygwin_HowTo.html) how to install Cygwin. 
The following Packages should be selected during the Cygwin installation:
- automake
- gcc-core
- gcc-fortran
- gcc-g++
- gdb
- libstdc++
- make

Eigen, libigl and JSON for Modern C++ are third party [header-only](https://en.wikipedia.org/wiki/Header-only) libraries, i.e., they do not need to be separately compiled, packaged and installed to be used :heart_eyes:.

## PolyMPS workflow

<img src='/output/workflow.jpg' width="100%">

## Input files

Have a look at some examples in the folder [**input**](https://github.com/rubensamarojr/polymps/tree/master/input).

In order to run these examples, please, go to the directory [**input/grid**](https://github.com/rubensamarojr/polymps/tree/master/input/grid) and extract the compressed folder [**grid.zip**](https://github.com/rubensamarojr/polymps/blob/master/input/grid/grid.zip) in the grid directory itself.

### Foldernames, Filenames, Physical and Numerical parameters
It is necessary to create a file (extension **.json**) and set all parameters.

### Solid domain
1. ... using **triangular meshes**. It is necessary to create a file (extension **.stl**) with informations about the initial geometry.
2. ... using **particles**. Necessary to add one layer of wall particles (material ID = 2) and two layers of dummy particles (material ID = 3) in the **.grid** file.

### Fluid domain
It is necessary to create a file (extension **.grid**) with informations about the initial geometry and some numerical and physical parameters:
- First line: **0**
- Second line: **number of particles**
- Next lines in the columns:
    - material **ID**
    - initial coordinates of particles **X** **Y** **Z** 
    - initial fluid velocities **VX** **VY** **VZ** (generally 0.0 0.0 0.0)
    - initial particle pressure **P** (generally 0.0)
    - initial presure average **PAV** (generally 0.0)

## Compile

Code compiled and tested on Windows 7, and Linux CentOS 7 and Ubuntu64.

You can build the project in GNU/Linux using the makefile. Follow these steps (CPU version):

Clone this repository into your system using `terminal in Linux`, and [Git BASH](https://gitforwindows.org) `or command prompt (cmd) in Windows`
```bash
git clone https://github.com/rubensaaj/polymps.git
```
Go to the folder **polymps**
```bash
cd polymps
```
Now, clone the third party libraries. You must run two commands
```bash
git submodule init
```
to initialize your local configuration file, and 
```bash
git submodule update
```
to fetch all the data from the third party libraries.

Edit the `Makefile` file (if necessary) with a text editor.

Make sure the environment is clean and ready to compile
```bash
make clean
```
Execute make
```bash
make all
```

This should create a binary `main` in folder **bin**

## Run

:warning:FIRST OF ALL, go to the directory [**input/grid**](https://github.com/rubensamarojr/polymps/tree/master/input/grid) and extract the compressed folder [**grid.zip**](https://github.com/rubensamarojr/polymps/blob/master/input/grid/grid.zip) in the grid directory itself. Check if **dam1610_3D_fluid_lo0p010_mps.grid** contains the data mentioned before in [Fluid domain](#fluid-domain).

### on LINUX

In the terminal, type
```bash
./bin/main
```
### on WINDOWS

You can do this in two ways:

1st way - In the command prompt, type
```bash
bin\main.exe
```
2nd way - Move the *main.exe* from the folder **bin** to the root folder **polymps**. After that, double click on *main.exe*.

## Json input filename
Type the name of the json input file (located in input directory), e.g.
```bash
MpsInputExample
```

You can specify a different case by changing numerical and physical parameters in the input json file. Also, input grid and stl files should be updated according your problem. We recommend that you rename the json file, e.g. **case_02.json**, and set a new name to output folder.
After that, you can [run](#run) PolyMPS at any time, and type the name of the new json input file, e.g.
```bash
case_02
```

## Output
This code writes pvd (header file) and corresponding vtu files as output. Look in the [**output**](https://github.com/rubensamarojr/polymps/tree/master/output) directory.
You can visualize them by open the pvd file with [Paraview](https://www.paraview.org) :eyes:.

## Directories

The PolyMPS contains several files and directories:

| File/Folder | Description |
| --- | --- |
| eigen | library for linear algebra: matrices, vectors, numerical solvers, and related algorithms |
| [include](https://github.com/rubensamarojr/polymps/tree/master/include) | header files |
| [input](https://github.com/rubensamarojr/polymps/tree/master/input) |	simple input examples (json, grid and stl files). Grid files compressed in a folder|
| libigl | geometry processing library |
| json | file that uses human-readable text to store and transmit data objects |
| [output](https://github.com/rubensamarojr/polymps/tree/master/output) |ouput files (pvd, vtu and txt files) |
| [src](https://github.com/rubensamarojr/polymps/tree/master/src) |	source files |
| [LICENSE](https://github.com/rubensamarojr/polymps/blob/master/LICENSE) |	MIT License |
| [Makefile](https://github.com/rubensamarojr/polymps/blob/master/Makefile) | set of tasks to compile the program |
| README |text file |

## How to cite PolyMPS

Please star :star: this project if it helps you.

If you publish work using PolyMPS, please cite the following reference:
```latex
@Article{Amaro2022,
  author    = {Amaro Jr., Rubens Augusto and Cheng, Liang-Yee},
  title     = {PolyMPS - An open source {CFD} solver based on Polygon walls in Moving Particle Semi-implicit ({MPS}) method},
  year      = {2022},
  journal   = {Software Impacts},
  volume    = {14},
  doi       = {10.1016/j.simpa.2022.100376},
  publisher = {Elsevier {BV}},
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
