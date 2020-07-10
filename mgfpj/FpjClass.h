// Forward projection class

#pragma once

#include <string>
#include <vector>

namespace mango
{
	struct Config
	{
		/*********************************************************
		* input and output directory and files
		*********************************************************/
		std::string		inputDir;
		std::string		outputDir;
		std::vector<std::string>	inputFiles;
		std::vector<std::string>	outputFiles;

		/*********************************************************
		* image parameters
		*********************************************************/
		int		imgDim;					// number of rows/cols of reconstructed images
		float	pixelSize;				// image pixel size [mm]
		int		sliceCount;				// number of slice in each image file

		/*********************************************************
		* geometry and detector parameters
		*********************************************************/
		float	sid;					// source to isocenter distance [mm]
		bool	nonuniformSID;			// whether the sids are nonuniform across views
		std::string		sidFile;		// name of the jsonc file of explicit sids
		float	sdd;					// source to detector distance [mm]
		bool	nonuniformSDD;			// whether the sdds are nonuniform across views
		std::string		sddFile;		// name of the jsonc file of explicit sdds

		float	startAngle = 0;			// angle position of source for the first view [degree]
		int		detEltCount;			// number of detector elements
		int		views;					// number of views
		float	totalScanAngle;			// total scan angle for short scan [degree]
		bool		nonuniformScanAngle;	// whether the scan angles are nonuniform across views
		std::string		scanAngleFile;		// name of the jsonc file of explicit scan angle values

		float	detEltSize;				// physical size of detector element [mm]
		float	detOffCenter;			// the position (coordinate) of center of detector

		int		oversampleSize = 1;		// oversample size

	};


	class FpjClass
	{
	public:
		static Config config;

	private:
		// array of sdds across views [mm]
		static float* sdd_array;
		// array of sids across views [mm]
		static float* sid_array;
		// array of detector element coordinate
		static float* u;
		// array of each view angle [radius]
		static float* beta;


	private:
		float* image = nullptr;
		float* sinogram_large = nullptr;	// sinogram for oversampling
		float* sinogram = nullptr;			// sinogram for final output

	public:
		FpjClass();
		~FpjClass();

		// Read config file
		void ReadConfigFile(const char* filename);

		// Initialize parameters
		void InitParam();

		// Read image file (raw format)
		void ReadImageFile(const char* filename);

		// Save sinogram data to file (raw format)
		void SaveSinogram(const char* filename);

		// Forward projection, using bilinear interpolation
		void ForwardProjectionBilinear();

	};

	

}