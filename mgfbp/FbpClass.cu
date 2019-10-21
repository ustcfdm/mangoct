#include "stdafx.h"
#include "FbpClass.cuh"

#include "FbpClass_Agent.cuh"


mango::Config mango::FbpClass::config;
float* mango::FbpClass::u = nullptr;
float* mango::FbpClass::beta = nullptr;
float* mango::FbpClass::reconKernel = nullptr;

mango::FbpClass::FbpClass()
{
}


mango::FbpClass::~FbpClass()
{
	FreeMemory_Agent(u);
	FreeMemory_Agent(beta);
	FreeMemory_Agent(reconKernel);

	FreeMemory_Agent(sinogram);
	FreeMemory_Agent(sinogram_filter);
	FreeMemory_Agent(image);

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
void mango::FbpClass::ReadConfigFile(const char * filename)
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

	// ========================================================================================
	// input and output directory and files
	// ========================================================================================

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

	// get value of saveFilterSinogram
	config.saveFilterSinogram = doc["SaveFilteredSinogram"].GetBool();


	// ========================================================================================
	// sinogram and slice parameters
	// ========================================================================================

	config.sgmWidth = doc["SinogramWidth"].GetUint();
	config.sgmHeight = doc["SinogramHeight"].GetUint();
	config.views = doc["Views"].GetUint();
	config.sliceCount = doc["SliceCount"].GetUint();

	config.detEltSize = doc["DetectorElementSize"].GetFloat();
	config.detOffCenter = doc["DetectorOffcenter"].GetFloat();

	config.sid = doc["SourceIsocenterDistance"].GetFloat();
	config.sdd = doc["SourceDetectorDistance"].GetFloat();

	// ========================================================================================
	// reconstruction parameters
	// ========================================================================================

	//BeamHardeningCorrection
	config.doBeamHardeningCorr = false;
	if (doc.HasMember("BeamHardeningCorrection"))
	{
		const js::Value& bh = doc["BeamHardeningCorrection"];
		if (bh.Size() > 10)
		{
			fprintf(stderr, "You have more than 10 beam hardening correction parameters.\n");
			exit(1);
		}

		for (js::SizeType i = 0; i < bh.Size(); i++)
		{
			config.beamHardening[i] = bh[i].GetFloat();
		}
		config.doBeamHardeningCorr = true;
	}


	config.imgDim = doc["ImageDimension"].GetUint();

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

	config.imgRot = doc["ImageRotation"].GetFloat();

	config.xCenter = doc["ImageCenter"][0].GetFloat();
	config.yCenter = doc["ImageCenter"][1].GetFloat();

	// reconstruction kernel
	if (doc.HasMember("HammingFilter"))
	{
		config.kernelName = "HammingFilter";

		config.kernelParam.push_back(doc["HammingFilter"].GetFloat());

	}
	else if (doc.HasMember("QuadraticFilter"))
	{
		config.kernelName = "QuadraticFilter";

		for (js::SizeType i = 0; i < doc["QuadraticFilter"].Size(); i++)
		{
			config.kernelParam.push_back(doc["QuadraticFilter"][i].GetFloat());
		}
	}
	else if (doc.HasMember("Polynomial"))
	{
		config.kernelName = "Polynomial";

		if (doc["Polynomial"].Size() >7)
		{
			fprintf(stderr, "Do not support more than 7 paramters of Polynomial kernel!\n");
			exit(1);
		}

		for (js::SizeType i = 0; i < doc["Polynomial"].Size(); i++)
		{
			config.kernelParam.push_back(doc["Polynomial"][i].GetFloat());
		}
	}
	else if (doc.HasMember("Delta"))
	{
		config.kernelName = "Delta";
		config.kernelParam.push_back(doc["Delta"].GetFloat());
	}
	else if (doc.HasMember("Hilbert"))
	{
		config.kernelName = "Hilbert";
		config.kernelParam.push_back(doc["Hilbert"].GetFloat());
	}
	else if (doc.HasMember("Hilbert_angle"))
	{
		config.kernelName = "Hilbert_angle";
		config.kernelParam.push_back(doc["Hilbert_angle"].GetFloat());
	}
	else
	{
		fprintf(stderr, "Did not find reconstruction kernel! Please check your config file: %s.\n", filename);
		exit(1);
	}

}

void mango::FbpClass::InitParam()
{
	InitializeU_Agent(u, config.sgmWidth, config.detEltSize, config.detOffCenter);

	InitializeBeta_Agent(beta, config.views, config.imgRot);

	InitializeReconKernel_Agent(reconKernel, config.sgmWidth, config.detEltSize, config.kernelName, config.kernelParam);

	cudaDeviceSynchronize();

	MallocManaged_Agent(sinogram, config.sgmWidth*config.sgmHeight*config.sliceCount * sizeof(float));
	MallocManaged_Agent(sinogram_filter, config.sgmWidth*config.views*config.sliceCount * sizeof(float));
	MallocManaged_Agent(image, config.imgDim*config.imgDim*config.sliceCount * sizeof(float));
}

void mango::FbpClass::ReadSinogramFile(const char * filename)
{
#pragma warning (disable : 4996)

	FILE* fp = fopen(filename, "rb");
	if (fp==NULL)
	{
		fprintf(stderr, "Cannot open file %s!\n", filename);
		exit(3);
	}

	fread(sinogram, sizeof(float), config.sgmWidth*config.sgmHeight*config.sliceCount, fp);

	fclose(fp);
}

void mango::FbpClass::SaveFilteredSinogram(const char * filename)
{
#pragma warning (disable : 4996)

	FILE* fp = fopen(filename, "wb");
	if (fp==NULL)
	{
		fprintf(stderr, "Cannot save to file %s!\n", filename);
		exit(4);
	}
	fwrite(sinogram_filter, sizeof(float), config.sgmWidth * config.views * config.sliceCount, fp);
	fclose(fp);
}

void mango::FbpClass::CorrectBeamHardening()
{
	CorrectBeamHardening_Agent(sinogram, config);
}

void mango::FbpClass::SaveImage(const char * filename)
{
#pragma warning (disable : 4996)

	FILE* fp = fopen(filename, "wb");
	if (fp == NULL)
	{
		fprintf(stderr, "Cannot save to file %s!\n", filename);
		exit(4);
	}
	fwrite(image, sizeof(float), config.imgDim * config.imgDim * config.sliceCount, fp);
	fclose(fp);
}

void mango::FbpClass::FilterSinogram()
{
	FilterSinogram_Agent(sinogram, sinogram_filter, reconKernel, u, config);
}

void mango::FbpClass::BackprojectPixelDriven()
{
	BackprojectPixelDriven_Agent(sinogram_filter, image, u, beta, config);
}
