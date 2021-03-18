# EMPS_MESH

Explicit and Weakly Compressible Moving Particle Semi-implicit method with Polygon wall

<img src='/output/dambreak01_movie.gif' width="500">

## Requirements

- [GCC (GNU Compiler Collection)](https://gcc.gnu.org/)
- [Eigen](http://eigen.tuxfamily.org/)
- [libigl](https://github.com/libigl/libigl)
- [JSON for Modern C++](https://github.com/nlohmann/json)

In order to install C++ compiler (GCC) on **windows**, we recommend to install [Cygwin](https://cygwin.com/). You can find [here](https://www3.ntu.edu.sg/home/ehchua/programming/howto/Cygwin_HowTo.html) how to install Cygwin. 
The following Packages should be installed:
- automake
- gcc
- gcc-c++
- gdb
- libstdc++
- make

Eigen, libigl and JSON for Modern C++ are third party header-only libraries.

## MPS input files

1. It is necessary to create a file (extension **.grid**) with informations about the initial geometry and some numerical and physical parameters:
- First line: **0**
- Second line: **number of particles**
- Next lines in the columns: **material ID** coordinates of particles **X** **Y** **Z** the initial fluid velocities (generally 0.0) **VX** **VY** **VZ** the initial particle pressure (generally 0.0) **P** and presure average (generally 0.0) **PAV**

2. It is necessary to create a file (extension **.json**) and set foldernames and filenames, as well as physical and numerical parameters.

There are some examples in the folder **input**.

## Compile

Code compiled and tested on Windows 7, and Linux CentOS 7 and Ubuntu64.

You can build the project in GNU/Linux using the makefile. Follow these steps (CPU version):

Clone this repository into your system using `terminal in Linux`, and [Git BASH](https://gitforwindows.org/) `or command prompt (cmd) in Windows`
```bash
git clone https://github.com/rubensaaj/EMPS_MESH.git
```
Go to the folder **EMPS_MESH**
```bash
cd EMPS_MESH
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
Run the code as
- LINUX
```bash
./bin/main
```
- WINDOWS
```bash
bin\main.exe
```

Type the name of the json input file (located in input directory), e.g.

```bash
MpsInputExample
```

This code writes pvd (header file) and corresponding vtu files as output.
You can visualize them by open the pvd file with [Paraview](https://www.paraview.org).

## Directories

The EMPS_MESH contains several files and directories:

| File/Folder | Description |
| --- | --- |
| eigen | library for linear algebra: matrices, vectors, numerical solvers, and related algorithms |
| include | header files |
| input |	simple input examples (json, grid and stl files) |
| libigl | geometry processing library |
| json | file that uses human-readable text to store and transmit data objects |
| output |ouput files (pvd, vtu and txt files) |
| src |	source files |
| LICENSE |	MIT License |
| Makefile | set of tasks to compile the program |
| README |text file |

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
