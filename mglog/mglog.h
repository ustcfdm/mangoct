#pragma once

#include "stdafx.h"

struct Config
{
	/*********************************************************
	* input and out director and files
	*********************************************************/

	std::string		inputDir;
	std::string		outputDir;

	std::string		bkgFile;
	std::vector<std::string>	inputFiles;

	std::string		outputFilePrefix;
	std::string		outputFilePostfix;
	std::vector<std::string>	outputFileReplace;

	/*********************************************************
	* sinogram and slice parameters
	*********************************************************/

	unsigned		bkgViews;
	unsigned		objViews;

	unsigned		rebinSize;
	unsigned		sliceNumber;
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


