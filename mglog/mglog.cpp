// mglog.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "mglog.h"

namespace fs = boost::filesystem;
namespace mg = mango;
namespace js = rapidjson;

void LoadConfigFile(js::Document& d, Config& config);


int main(int argc, char* argv[])
{

	if (argc == 1)
	{
		fprintf(stderr, "No config files!\n");
		return -1;
	}


	for (int configIdx = 1; configIdx < argc; configIdx++)
	{

		//////////////////////////////////////////////////////////////////////////////
		// Step 1: load the config file
		//////////////////////////////////////////////////////////////////////////////
		std::ifstream ifs(argv[configIdx]);
		js::IStreamWrapper isw(ifs);
		js::Document doc;
		doc.ParseStream<js::kParseCommentsFlag | js::kParseTrailingCommasFlag>(isw);

		LoadConfigFile(doc, config);

	}


	return 0;
}


// acquire the list of file names in the dir that matches filter
std::vector<std::string> GetInputFileNames(const std::string& dir, const std::string& filter)
{
	std::vector<std::string> filenames;

	if (!fs::is_directory(dir))
	{
		return filenames;
	}

	fs::directory_iterator end_iter;
	std::regex e(filter);

	for (auto&& fe : fs::directory_iterator(dir))
	{
		std::string file = fe.path().filename().string();

		if (std::regex_match(file, e))
		{
			filenames.push_back(file);
		}
	}
	return filenames;
}

// acquire the list of output file names, replace substring in input file names, add prefix and postfix
std::vector<std::string> GetOutputFileNames(const std::vector<std::string>& inputFileNames, const std::vector<std::string>& replace, const std::string& prefix)
{
	std::vector<std::string> outputFiles;

	for (size_t fileIdx = 0; fileIdx < inputFileNames.size(); fileIdx++)
	{
		std::string outputFile = inputFileNames[fileIdx];
		for (size_t i = 0; i < replace.size() / 2; i++)
		{
			auto pos = outputFile.find(replace[2 * i]);
			if (pos == std::string::npos)
			{
				fprintf(stderr, "Did not find substring \"%s\" to be replaced!\n", replace[2 * i].c_str());
				exit(2);
			}
			outputFile.replace(pos, replace[2 * i].size(), replace[2 * i + 1]);
		}
		outputFiles.push_back(prefix + outputFile);
	}

	return outputFiles;
}

// load the config parameters in json document d into struct config
void LoadConfigFile(js::Document& d, Config& config)
{
	/*********************************************************
	* input and output director and files
	*********************************************************/

	config.inputDir = d["InputDir"].GetString();
	config.outputDir = d["OutputDir"].GetString();

	// check inputDir and outputDir
	fs::path inDir(config.inputDir);
	fs::path outDir(config.outputDir);
	if (!fs::exists(inDir))
	{
		fprintf(stderr, "Input directory %s does not exist!\n", config.inputDir.c_str());
		exit(1);
	}
	if (!fs::exists(outDir))
	{
		fprintf(stderr, "Output directory %s does not exist!\n", config.inputDir.c_str());
		exit(1);
	}

	config.bkgFile = d["BackgroundFile"].GetString();
	// check background file
	if (!fs::exists(inDir / config.bkgFile))
	{
		fprintf(stderr, "Cannot find backgorund file %s!\n", config.bkgFile.c_str());
		exit(1);
	}

	// get input file names
	config.inputFiles = GetInputFileNames(config.inputDir, d["InputFiles"].GetString());
	if (config.inputFiles.size() == 0)
	{
		fprintf(stderr, "No file name like %s!\n", d["InputFiles"].GetString());
		exit(1);
	}

	// get output file names
	const js::Value& replaceJs = d["OutputFileReplace"];
	std::vector<std::string> replace;
	for (js::SizeType i = 0; i < replaceJs.Size(); i++)
	{
		replace.push_back(replaceJs[i].GetString());
	}

	config.outputFiles = GetOutputFileNames(config.inputFiles, replace, d["OutputFilePrefix"].GetString());


	/*********************************************************
	* sinogram and slice parameters
	*********************************************************/

	config.bkgViews = d["BackgroundViews"].GetUint();
	config.objViews = d["ObjectViews"].GetUint();
	config.rebinSize = d["RebinSize"].GetUint();
	config.sliceCount = d["SliceCount"].GetUint();
	config.sliceStartIdx = d["SliceStartIdx"].GetUint();
	config.sliceThickness = d["SliceThickness"].GetUint();

	/*********************************************************
	* correct sinogram artifacts
	*********************************************************/

	config.interpolateWhiteLines = d["InterpolateWhiteLines"].GetBool();
	config.smoothoutMarginalData = d["SmoothoutMarginalData"].GetBool();
	config.smoothoutLeftIdx = d["SmoothoutLeftIdx"].GetUint();
	config.smoothoutRightIdx = d["SmoothoutRightIdx"].GetUint();

	/*********************************************************
	*  EVI file parameters
	*********************************************************/

	config.detectorWidth = d["DetectorWidth"].GetUint();
	config.detectorHeight = d["DetectorHeight"].GetUint();
	config.offsetToFirstImage = d["OffsetToFirstImage"].GetUint();
	config.gap = d["Gap"].GetUint();

}