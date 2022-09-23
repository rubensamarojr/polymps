/**
 * @defgroup   MPSINPUTOUTPUT Input/Output
 *
 * @brief      This file implements Input/Output file managers.
 * @details    Read input data from .json, .stl and .grid files. Write particles data in .vtu files and bucktes in vtk files.
 *
 * @author     Rubens Amaro
 * @date       2022
 */

#ifndef MPS_INCLUDE_INPUTOUTPUT_H_
#define MPS_INCLUDE_INPUTOUTPUT_H_

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "MpsParticleSystem.h"
#include "MpsParticle.h"
#include "MpsPressure.h"

/**
 * @brief      This class manage input/output files.
 * @details    Read input data from .json, .stl and .grid files. Write particles data in .vtu files and bucktes in vtk files.
 */
class MpsInputOutput {
public:
	/**
	 * @brief      Constructs a new instance.
	 */
	MpsInputOutput();
	
	/**
	 * @brief      Destroys the object.
	 */
	virtual ~MpsInputOutput();

	double getTime();

	void displayInfo(MpsParticleSystem *PSystem, MpsParticle *Particles, MpsPressure* PPress, const int intervalIter);

	/**
	 * @brief      Reads input data from input file .json and set variables.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 */
	void readInputFile(MpsParticleSystem *PSystem, MpsParticle *Particles);

	/**
	 * @brief      Reads data from input file .grid to class MpsParticle and allocate memory to its attributes.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 */
	void readMpsParticleFile(MpsParticleSystem *PSystem, MpsParticle *Particles);

	/**
	 * @brief      Calls the functions to write output files.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 */
	void writeOutputFiles(MpsParticleSystem *PSystem, MpsParticle *Particles);
	
	/**
	 * @brief      Writes output files .prof in ascii format.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @note       Implemented but not used.
	 */
	void writeProfAscii(MpsParticleSystem *PSystem, MpsParticle *Particles);

	/**
	 * @brief      SwapEnd.
	 * @param      var   The variable
	 * @tparam     T     { description }
	 * @see        <a href="https://stackoverflow.com/questions/105252" target="_blank">How do I convert between big-endian and little-endian values in C++?</a>
	 * @see        <a href="https://stackoverflow.com/questions/10913666" target="_blank">Error writing binary VTK files</a>
	 * @see        <a href="https://stackoverflow.com/questions/55829282" target="_blank">Write vtk file in binary format</a>
	 */
	template <typename T> void SwapEnd(T& var);

	/**
	 * @brief      Writes output files .vtu (Paraview) in binary format.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @warning    This function is not working.
	 */
	void writeVtuBinary(MpsParticleSystem *PSystem, MpsParticle *Particles);

	/**
	 * @brief      Writes output files .vtu (Paraview) in binary format.
	 * @details    Only free-surface particles are considered to be shown.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @warning    This function is not working.
	 */
	void writeVtuBinaryFreeSurface(MpsParticleSystem *PSystem, MpsParticle *Particles);
	
	/**
	 * @brief      Writes output files .vtu (Paraview) in ascii format.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 */
	void writeVtuAscii(MpsParticleSystem *PSystem, MpsParticle *Particles);
	
	/**
	 * @brief      Writes output files .vtu (Paraview) in ascii format.
	 * @details    Only free-surface particles are considered to be shown.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 */
	void writeVtuAsciiFreeSurface(MpsParticleSystem *PSystem, MpsParticle *Particles);

	/**
	 * @brief      Write header for vtu files.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 */
	void writePvd(MpsParticleSystem *PSystem, MpsParticle *Particles);

	/**
	 * @brief      Writes buckets in .vtk file (Paraview).
	 * @details    The bounding box of the entire analysis domain is filled with uniform buckets of cubes. 
	 * Particles are assigned to corresponding buckets. Therefore, each bucket has a sequence of particle 
	 * data assigned to it as a linked-list.
	 * @see        <a href="https://doi.org/10.1007/s40571-015-0059-2" target="_blank">Performance improvements of differential operators code for MPS method on GPU</a>
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 */
	void writeBuckets(MpsParticleSystem *PSystem, MpsParticle *Particles);

	// Pressure sensors

	/**
	 * @brief      Writes pressure at sensors in .txt file..
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @warning    Not working.
	 */
	void writePressSensors(MpsParticleSystem *PSystem, MpsParticle *Particles);
	
	/**
	 * @brief      Writes forces on walls in .txt file.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @warning    Not implemented.
	 */
	//void writeForce(MpsParticle *Particles);

	/**
	 * @brief      Writes header in .txt files.
	 * @param      PSystem    The particle system
	 * @param      Particles  The particles data
	 * @warning    Not working.
	 */
	void writeHeaderTxtFiles(MpsParticleSystem *PSystem, MpsParticle *Particles);

	/**
	 * @brief      Convert String to Char Array
	 * @param      out_char  The Char Array
	 */
	void stringToChar(char *out_char);

	/**
	 * @brief      Delete all files inside the output folder.
	 */
	void deleteDirectoryFiles();
	
	std::string meshRigidFilename;
	std::string meshDeformableFilename;
	std::string meshForcedFilename;

protected:
	
private:
	// Files
	// MpsParticle data json file
	FILE* js;
	// MpsParticle data vtu files
	FILE* fp;
	// Pressure txt file
	FILE* pressTxtFile;
	// Force txt file
	FILE* forceTxtFile;

	std::string gridFilename;
	std::string forceTxtFilename;
	std::string pressTxtFilename;
	std::string vtuOutputFoldername;
	
	int vtuType;
	bool freeSurfWall;
	bool txtPress;
	bool txtForce;
	bool outputPnd;
	bool outputNeigh;
	bool outputDeviation;
	bool outputConcentration;
	bool outputAuxiliar;
	bool outputNonNewtonian;
};

#endif // MPS_INCLUDE_INPUTOUTPUT_H_