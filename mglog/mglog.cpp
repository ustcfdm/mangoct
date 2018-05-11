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
std::vector<std::string> GetFileNames(const std::string& dir, const std::string& filter)
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

// load the config parameters in json document d into struct config
void LoadConfigFile(js::Document& d, Config& config)
{
	/*********************************************************
	* input and out director and files
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
	config.inputFiles = GetFileNames(config.inputDir, d["InputFiles"].GetString());
	if (config.inputFiles.size()==0)
	{
		fprintf(stderr, "No file name like %s!\n", d["InputFiles"].GetString());
		exit(1);
	}

}