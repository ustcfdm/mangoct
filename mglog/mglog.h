#pragma once

#include "stdafx.h"

struct Config
{
	/*********************************************************
	* input and out directory and files
	*********************************************************/

	std::string		inputDir;
	std::string		outputDir;

	std::string		bkgFile;
	std::vector<std::string>	inputFiles;
	std::vector<std::string>	outputFiles;

	// For dual energy subtraction, we need more variables. Below are for high energy and subtraction.
	// For subtraction, e.g. bkgFile - bkgFile2. The result is in outputFiles3.
	bool			dualEnergySubtraction;
	std::string		bkgFile2;
	std::vector<std::string>	inputFiles2;
	std::vector<std::string>	outputFiles2;
	std::vector<std::string>	outputFiles3;
	std::vector<std::string>	outputTypes;

	/*********************************************************
	* sinogram and slice parameters
	*********************************************************/

	unsigned		bkgViews;
	unsigned		objViews;

	unsigned		rebinSize;
	unsigned		sliceCount;
	unsigned		sliceStartIdx;
	unsigned		sliceThickness;

	/*********************************************************
	* correct sinogram artifacts
	*********************************************************/

	bool		interpolateWhiteLines;
	bool		smoothoutMarginalData;
	unsigned	smoothoutLeftIdx;
	unsigned	smoothoutRightIdx;

	/*********************************************************
	*  EVI file parameters, usually you do not need to change
	*********************************************************/

	unsigned	detectorWidth;
	unsigned	detectorHeight;
	unsigned	offsetToFirstImage;
	unsigned	gap;

}config;


