#include "FpjClass.h"
#include "FpjClass_Agent.cuh"
#include "stdafx.h"

mango::Config mango::FpjClass::config;
float* mango::FpjClass::sdd_array = nullptr;
float* mango::FpjClass::sid_array = nullptr;
float* mango::FpjClass::offcenter_array = nullptr;
float* mango::FpjClass::v = nullptr;
float* mango::FpjClass::u = nullptr;
float* mango::FpjClass::beta = nullptr;
float* mango::FpjClass::swing_angle_array = nullptr;


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
	namespace fs = std::filesystem;

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
	namespace fs = std::filesystem;
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

	if (doc.HasMember("ConeBeam"))
	{
		config.coneBeam = doc["ConeBeam"].GetBool();
	}

#pragma endregion


#pragma region  geometry and detector parameters
	config.sid = doc["SourceIsocenterDistance"].GetFloat();

	if (doc.HasMember("SIDFile"))
	{
		printf("--nonuniform SID--\n");
		config.nonuniformSID = true;
		config.sidFile = doc["SIDFile"].GetString();
	}
	else
	{
		config.nonuniformSID = false;
	}

	config.sdd = doc["SourceDetectorDistance"].GetFloat();

	if (doc.HasMember("SDDFile"))
	{
		printf("--nonuniform SDD--\n");
		config.nonuniformSDD = true;
		config.sddFile = doc["SDDFile"].GetString();
	}
	else
	{
		config.nonuniformSDD = false;
	}

	if (doc.HasMember("SwingAngleFile"))
	{
		printf("--nonzero swing angle--\n");
		config.nonZeroSwingAngle = true;
		config.swingAngleFile = doc["SwingAngleFile"].GetString();
	}
	else
	{
		config.nonZeroSwingAngle = false;
	}

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

	if (doc.HasMember("ScanAngleFile"))
	{
		printf("--nonuniform scan angle--\n");
		config.nonuniformScanAngle = true;
		config.scanAngleFile = doc["ScanAngleFile"].GetString();
	}
	else
	{
		config.nonuniformScanAngle = false;
	}

	if (abs(config.totalScanAngle - 360.0f) < 0.001f)
	{
		printf("--FULL scan--\n");
	}
	else
	{
		printf("--SHORT scan (%.1f degrees)--\n", config.totalScanAngle);
	}

	config.detEltCount = doc["DetectorElementCount"].GetInt();
	config.views = doc["Views"].GetInt();

	config.detEltSize = doc["DetectorElementSize"].GetFloat();
	config.detOffCenter = doc["DetectorOffcenter"].GetFloat();

	if (doc.HasMember("DetectorOffCenterFile"))
	{
		printf("--nonuniform offcenter--\n");
		config.nonuniformOffCenter = true;
		config.offCenterFile = doc["DetectorOffCenterFile"].GetString();
	}
	else
	{
		config.nonuniformOffCenter = false;
	}

	if (doc.HasMember("OversampleSize"))
	{
		config.oversampleSize = doc["OversampleSize"].GetInt();
	}
#pragma endregion

#pragma region cone beam parameters
	if (config.coneBeam == true)
	{
		printf("--CONE beam--\n");
		if (doc.HasMember("ImageSliceThickness"))
		{
			config.sliceThickness = doc["ImageSliceThickness"].GetFloat();
		}
		if (doc.HasMember("DetectorZElementCount"))
		{
			config.detZEltCount = doc["DetectorZElementCount"].GetInt();
		}
		if (doc.HasMember("DetectorElementHeight"))
		{
			config.detEltHeight = doc["DetectorElementHeight"].GetFloat();
		}
		if (doc.HasMember("DetectorZOffcenter"))
		{
			config.detZoffCenter = doc["DetectorZOffcenter"].GetFloat();
		}
	}
	else
	{
		config.detZEltCount = config.sliceCount;
	}
#pragma endregion

#pragma region water mu parameters
	if (doc.HasMember("WaterMu"))
	{
		printf("--Images are in HU values--\n");
		config.converToHU = true;
		config.waterMu = doc["WaterMu"].GetFloat();
	}
	else
	{
		config.converToHU = false;
	}
#pragma endregion
}

void mango::FpjClass::InitParam()
{
	if (config.nonuniformSDD == true)
	{
		InitializeNonuniformPara_Agent(sdd_array, config.views, config.sddFile);
	}
	else
	{
		InitializeDistance_Agent(sdd_array, config.sdd, config.views);
	}

	if (config.nonuniformSID == true)
	{
		InitializeNonuniformPara_Agent(sid_array, config.views, config.sidFile);
	}
	else
	{
		InitializeDistance_Agent(sid_array, config.sid, config.views);
	}

	if (config.nonuniformOffCenter == true)
	{
		InitializeNonuniformPara_Agent(offcenter_array, config.views, config.offCenterFile);
		float*offcenter_array_cpu = new float[config.views];
		cudaMemcpy(offcenter_array_cpu, offcenter_array, sizeof(float)*config.views, cudaMemcpyDeviceToHost);
		config.detOffCenter = offcenter_array_cpu[0];
	}
	else
	{
		InitializeDistance_Agent(offcenter_array, config.detOffCenter, config.views);
	}

	if (config.nonZeroSwingAngle == true)
	{
		InitializeNonuniformPara_Agent(swing_angle_array, config.views, config.swingAngleFile);
		// the same function of offcenter can be used to import swing angles
	}
	else
	{
		InitializeDistance_Agent(swing_angle_array, 0.0f, config.views);
	}


	if (config.coneBeam == true)
	{
		InitializeU_Agent(v, config.detZEltCount, config.detEltHeight, config.detZoffCenter);
	}

	InitializeU_Agent(u, config.detEltCount*config.oversampleSize, config.detEltSize/config.oversampleSize, config.detOffCenter);
	if (config.nonuniformScanAngle == true)
	{
		InitializeNonuniformBeta_Agent(beta, config.views, config.startAngle, config.scanAngleFile);
		float*beta_cpu = new float[config.views];
		cudaMemcpy(beta_cpu, beta, sizeof(float)*config.views, cudaMemcpyDeviceToHost);
		config.totalScanAngle = (beta_cpu[config.views - 1] - beta_cpu[0] + beta_cpu[1] - beta_cpu[0]) / 3.1415926f * 180;
	}
	else
	{
		InitializeBeta_Agent(beta, config.views, config.startAngle, config.totalScanAngle);
	}

	cudaDeviceSynchronize();

	MallocManaged_Agent(image, config.imgDim*config.imgDim*config.sliceCount * sizeof(float));
	MallocManaged_Agent(sinogram, config.detEltCount*config.views * sizeof(float));// the size of the sinogram is limited to one slice and will be saved slice by slice
	MallocManaged_Agent(sinogram_large, config.detEltCount * config.oversampleSize *config.views  * sizeof(float));// the size of the sinogram is limited to one slice and will be saved slice by slice

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

	// if the image has been converted to HU values
	// we need to convert it back to mu values
	if (config.converToHU)
	{
		for (int idx = 0; idx < config.imgDim*config.imgDim*config.sliceCount; idx++)
		{
			image[idx] = (image[idx] + 1000.0f) / 1000.0f*config.waterMu;
		}
	}

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
	fwrite(sinogram, sizeof(float), config.detEltCount * config.views * config.detZEltCount, fp);

	fclose(fp);
}

void mango::FpjClass::ForwardProjectionBilinear()
{
	ForwardProjectionBilinear_Agent(image, sinogram_large,sid_array,sdd_array, offcenter_array, u, v, beta, swing_angle_array, config, 0);

	BinSinogram(sinogram_large, sinogram, config);
}


void mango::FpjClass::ForwardProjectionBilinearAndSave(const char* filename)
{
	printf("\nProcessing slice# ");
	for (int z_idx = 0; z_idx < config.detZEltCount; z_idx++)
	{
		if (z_idx != 0)
		{
			printf("\b\b\b\b\b\b\b");
		}
		printf("%3d/%3d", z_idx + 1, config.detZEltCount);
		ForwardProjectionBilinear_Agent(image, sinogram_large, sid_array, sdd_array, offcenter_array, u, v, beta, swing_angle_array, config, z_idx);

		BinSinogram(sinogram_large, sinogram, config);

		cudaDeviceSynchronize();

		SaveSinogramSlice(filename, sinogram, z_idx, config);
	}	
}