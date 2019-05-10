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
		unsigned	imgDim;					// number of rows/cols of reconstructed images
		float		pixelSize;				// image pixel size [mm]
		unsigned	sliceCount;				// number of slice in each image file

		/*********************************************************
		* geometry and detector parameters
		*********************************************************/
		float		sid;					// source to isocenter distance [mm]
		float		sdd;					// source to detector distance [mm]

		float		startAngle = 0;			// angle position of source for the first view [degree]
		unsigned	detEltCount;			// number of detector elements
		unsigned	views;					// number of views

		float		detEltSize;				// physical size of detector element [mm]
		float		detOffCenter;			// the position (coordinate) of center of detector

	};


	class FpjClass
	{
	public:
		static Config config;

	private:
		// array of detector element coordinate
		static float* u;
		// array of each view angle [radius]
		static float* beta;


	private:
		float* image = nullptr;
		float* sinogram = nullptr;

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