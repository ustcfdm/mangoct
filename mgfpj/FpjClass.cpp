#include "FpjClass.h"
#include "FpjClass_Agent.cuh"
#include "stdafx.h"

mango::Config mango::FpjClass::config;
float* mango::FpjClass::u = nullptr;
float* mango::FpjClass::beta = nullptr;

mango::FpjClass::FpjClass()
{
}

mango::FpjClass::~FpjClass()
{
	FreeMemory_Agent(u);
	FreeMemory_Agent(beta);

	FreeMemory_Agent(image);
	FreeMemory_Agent(sinogram);
}


// acquire the list of file names in the dir that matches filter
std::vector<std::string> GetInputFileNames(const std::string& dir, const std::string& filter)
{
	namespace fs = std::experimental::filesystem;

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

// read and parse config file
void mango::FpjClass::ReadConfigFile(const char * filename)
{
	namespace fs = std::experimental::filesystem;
	namespace js = rapidjson;

	// load the config file
	std::ifstream ifs(filename);
	if (!ifs)
	{
		printf("Cannot open config file '%s'!\n", filename);
		exit(-2);
	}
	rapidjson::IStreamWrapper isw(ifs);
	rapidjson::Document doc;
	doc.ParseStream<js::kParseCommentsFlag | js::kParseTrailingCommasFlag>(isw);

#pragma region input and output directory and files

	config.inputDir = doc["InputDir"].GetString();
	config.outputDir = doc["OutputDir"].GetString();

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
		fprintf(stderr, "Output directory %s does not exist!\n", config.outputDir.c_str());
		exit(1);
	}

	// get input file names
	config.inputFiles = GetInputFileNames(config.inputDir, doc["InputFiles"].GetString());
	if (config.inputFiles.size() == 0)
	{
		fprintf(stderr, "No file name like %s!\n", doc["InputFiles"].GetString());
		exit(1);
	}

	// get output file names
	const js::Value& replaceJs = doc["OutputFileReplace"];
	std::vector<std::string> replace;
	for (js::SizeType i = 0; i < replaceJs.Size(); i++)
	{
		replace.push_back(replaceJs[i].GetString());
	}

	config.outputFiles = GetOutputFileNames(config.inputFiles, replace, doc["OutputFilePrefix"].GetString());

#pragma endregion


#pragma region image parameters
	config.imgDim = doc["ImageDimension"].GetInt();

	// get pixel size
	if (doc.HasMember("PixelSize"))
	{
		config.pixelSize = doc["PixelSize"].GetFloat();
	}
	else if (doc.HasMember("ImageSize"))
	{
		config.pixelSize = doc["ImageSize"].GetFloat() / config.imgDim;
	}
	else
	{
		fprintf(stderr, "Did not find PixelSize or ImageSize! Please check your config file: %s.\n", filename);
		exit(1);
	}

	config.sliceCount = doc["SliceCount"].GetInt();

#pragma endregion


#pragma region  geometry and detector parameters
	config.sid = doc["SourceIsocenterDistance"].GetFloat();
	config.sdd = doc["SourceDetectorDistance"].GetFloat();

	if (doc.HasMember("StartAngle"))
	{
		config.startAngle = doc["StartAngle"].GetFloat();
	}

	if (doc.HasMember("TotalScanAngle"))
	{
		config.totalScanAngle = doc["TotalScanAngle"].GetFloat();
	}
	else
	{
		config.totalScanAngle = 360.0f;
	}

	config.detEltCount = doc["DetectorElementCount"].GetInt();
	config.views = doc["Views"].GetInt();

	config.detEltSize = doc["DetectorElementSize"].GetFloat();
	config.detOffCenter = doc["DetectorOffcenter"].GetFloat();

	if (doc.HasMember("OversampleSize"))
	{
		config.oversampleSize = doc["OversampleSize"].GetInt();
	}
#pragma endregion
}

void mango::FpjClass::InitParam()
{
	InitializeU_Agent(u, config.detEltCount*config.oversampleSize, config.detEltSize/config.oversampleSize, config.detOffCenter);
	InitializeBeta_Agent(beta, config.views, config.startAngle,config.totalScanAngle);

	cudaDeviceSynchronize();

	MallocManaged_Agent(image, config.imgDim*config.imgDim*config.sliceCount * sizeof(float));
	MallocManaged_Agent(sinogram, config.detEltCount*config.views*config.sliceCount * sizeof(float));
	MallocManaged_Agent(sinogram_large, config.detEltCount * config.oversampleSize * config.views * config.sliceCount * sizeof(float));

}

void mango::FpjClass::ReadImageFile(const char * filename)
{
#pragma warning (disable : 4996)

	FILE* fp = fopen(filename, "rb");
	if (fp == NULL)
	{
		fprintf(stderr, "Cannot open file %s!\n", filename);
		exit(3);
	}

	fread(image, sizeof(float), config.imgDim*config.imgDim*config.sliceCount, fp);

	fclose(fp);
}

void mango::FpjClass::SaveSinogram(const char * filename)
{
#pragma warning (disable : 4996)

	FILE* fp = fopen(filename, "wb");
	if (fp == NULL)
	{
		fprintf(stderr, "Cannot save to file %s!\n", filename);
		exit(4);
	}
	fwrite(sinogram, sizeof(float), config.detEltCount * config.views * config.sliceCount, fp);

	fclose(fp);
}

void mango::FpjClass::ForwardProjectionBilinear()
{
	ForwardProjectionBilinear_Agent(image, sinogram_large, u, beta, config);

	BinSinogram(sinogram_large, sinogram, config);
}
