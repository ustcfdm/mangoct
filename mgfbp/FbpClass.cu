#include "stdafx.h"
#include "FbpClass.cuh"

#include "FbpClass_Agent.cuh"


mango::Config mango::FbpClass::config;
float* mango::FbpClass::sid_array = nullptr;
float* mango::FbpClass::sdd_array = nullptr;
float* mango::FbpClass::offcenter_array = nullptr;
float* mango::FbpClass::pmatrix_array = nullptr;
float* mango::FbpClass::u = nullptr;
float* mango::FbpClass::v = nullptr;
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
				fprintf(stderr, "Did not find substring \"%s\" to be replaced!\n (Every sinogram file in the input folder should have this substring)\n", replace[2 * i].c_str());
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

	if (doc.HasMember("TotalScanAngle"))
		config.totalScanAngle = doc["TotalScanAngle"].GetFloat();
	else
		config.totalScanAngle = 360.0f;

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

	// for recon with pmatrix
	if (doc.HasMember("PMatrixFile"))
	{
		printf("--pmatrix applied--\n");
		config.pmatrixFlag = true;
		config.pmatrixFile = doc["PMatrixFile"].GetString();
	}
	else
	{
		config.pmatrixFlag = false;
	}

	// for cone beam reconstruction
	if (doc.HasMember("SliceThickness"))
		config.sliceThickness = doc["SliceThickness"].GetFloat();
	else
		config.sliceThickness = 0;

	if (doc.HasMember("SliceOffCenter"))
		config.sliceOffcenter = doc["SliceOffCenter"].GetFloat();
	else
		config.sliceOffcenter = 0;

	// scan angle file info for nonuniform scan angles
	if (doc.HasMember("ScanAngleFile"))
	{
		config.nonuniformScanAngle = true;
		config.scanAngleFile = doc["ScanAngleFile"].GetString();
	}
	else
	{
		config.nonuniformScanAngle = false;
	}

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

	if (doc.HasMember("ConeBeam"))
		config.coneBeam = doc["ConeBeam"].GetBool();
	else
		config.coneBeam = false; 

	if (config.coneBeam)
		printf("--CONE beam--\n");
	else
		printf("--FAN beam--\n");

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

	if (doc.HasMember("ImageRotation"))
		config.imgRot = doc["ImageRotation"].GetFloat();
	else
		config.imgRot = 0;

	config.xCenter = doc["ImageCenter"][0].GetFloat();
	config.yCenter = doc["ImageCenter"][1].GetFloat();

	// for cone beam reconstruction
	if (doc.HasMember("ImageSliceThickness") && config.coneBeam)
		config.imgSliceThickness = doc["ImageSliceThickness"].GetFloat();
	else
		config.imgSliceThickness = config.sliceThickness; 

	if (doc.HasMember("ImageSliceCount") && config.coneBeam)
	{
		config.imgSliceCount = doc["ImageSliceCount"].GetUint();
	}	
	else
	{
		config.imgSliceCount = config.sliceCount;
	}
		
	
	if (doc.HasMember("ImageCenterZ") && config.coneBeam)
		config.zCenter = doc["ImageCenterZ"].GetFloat();
	else
		config.zCenter = config.sliceOffcenter;


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
	else if (doc.HasMember("GaussianApodizedRamp"))
	{
		config.kernelName = "GaussianApodizedRamp";
		config.kernelParam.push_back(doc["GaussianApodizedRamp"].GetFloat());
	}
	else
	{
		fprintf(stderr, "Did not find reconstruction kernel! Please check your config file: %s.\n", filename);
		exit(1);
	}

	if (doc.HasMember("WaterMu"))
	{
		printf("--Convert To HU--\n");
		config.converToHU = true;
		config.waterMu = doc["WaterMu"].GetFloat();
	}
	else
	{
		config.converToHU = false;
	}

}

void mango::FbpClass::InitParam()
{

	if (config.nonuniformSDD == true)
	{
		InitializeNonuniformSDD_Agent(sdd_array, config.views, config.sddFile);
	}
	else
	{
		InitializeDistance_Agent(sdd_array, config.sdd, config.views);
	}

	if (config.nonuniformSID == true)
	{
		InitializeNonuniformSID_Agent(sid_array, config.views, config.sidFile);
	}
	else
	{
		InitializeDistance_Agent(sid_array, config.sid, config.views);
	}

	if (config.nonuniformOffCenter == true)
	{
		InitializeNonuniformOffCenter_Agent(offcenter_array, config.views, config.offCenterFile);
	}
	else
	{
		InitializeDistance_Agent(offcenter_array, config.detOffCenter, config.views);
	}

	if (config.pmatrixFlag == true)
	{
		InitializePMatrix_Agent(pmatrix_array,config.views, config.pmatrixFile);
	}
	{
		;
	}

	float* offcenter_array_cpu = new float[config.views];
	cudaMemcpy(offcenter_array_cpu, offcenter_array, config.views * sizeof(float), cudaMemcpyDeviceToHost);

	InitializeU_Agent(u, config.sgmWidth, config.detEltSize, offcenter_array_cpu[0]);

	InitializeU_Agent(v, config.sliceCount, config.sliceThickness, config.sliceOffcenter);

	if (config.nonuniformScanAngle == true)
	{
		printf("--nonuniform scan angle--\n");
		InitializeNonuniformBeta_Agent(beta, config.views, config.imgRot, config.scanAngleFile);

		float *beta_cpu = new float[config.views];
		cudaMemcpy(beta_cpu,beta , sizeof(float)*config.views, cudaMemcpyDeviceToHost);

		config.totalScanAngle = (beta_cpu[config.views - 1] - beta_cpu[0])/float(config.views)*float(config.views + 1)/ 3.1415926f * 180;
		//It is not easy to define the total scan angle for a non uniform scan. 
		//This equation is just one method. 
	}
	else
	{
		InitializeBeta_Agent(beta, config.views, config.imgRot, config.totalScanAngle);
	}

	if (360.0f - abs(config.totalScanAngle) < 0.01f)
	{
		config.shortScan = false;
		printf("--FULL scan--\n");
	}
	else
	{
		config.shortScan = true;
		printf("--SHORT scan (scan angle = %.2f degrees)--\n", abs(config.totalScanAngle));
	}


	InitializeReconKernel_Agent(reconKernel, config.sgmWidth, config.detEltSize, config.kernelName, config.kernelParam);

	cudaDeviceSynchronize();

	MallocManaged_Agent(sinogram, config.sgmWidth*config.sgmHeight*config.sliceCount * sizeof(float));
	MallocManaged_Agent(sinogram_filter, config.sgmWidth*config.views*config.sliceCount * sizeof(float));
	//MallocManaged_Agent(image, config.imgDim*config.imgDim*config.imgSliceCount * sizeof(float));
	MallocManaged_Agent(image, config.imgDim*config.imgDim* sizeof(float));// recon slice by slice
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
	fwrite(image, sizeof(float), config.imgDim * config.imgDim * config.imgSliceCount, fp);
	fclose(fp);
}

void mango::FbpClass::FilterSinogram()
{
	FilterSinogram_Agent(sinogram, sinogram_filter, reconKernel, u, config,beta,sdd_array, offcenter_array);
}

void mango::FbpClass::BackprojectPixelDriven()
{
	BackprojectPixelDriven_Agent(sinogram_filter, image, sdd_array,sid_array, offcenter_array, pmatrix_array, u, v, beta, config, 0);
}

//this function is to reconstruct and save the image slice by slice
void mango::FbpClass::BackprojectPixelDrivenAndSave(const char* filename)

{
	printf(" slice# ");
	for (int z_idx = 0; z_idx < config.imgSliceCount; z_idx++)
	{
		if (z_idx != 0)
		{
			printf("\b\b\b\b\b\b\b");
		}
		printf("%3d/%3d", z_idx + 1, config.imgSliceCount);
		BackprojectPixelDriven_Agent(sinogram_filter, image, sdd_array, sid_array, offcenter_array, pmatrix_array, u, v, beta, config, z_idx);

		cudaDeviceSynchronize();
		//printf("%f", image[0]);
		//convert to HU
		
		if (config.converToHU)
		{
			for (int row_idx = 0; row_idx < config.imgDim; row_idx++)
			{
				for (int col_idx = 0; col_idx < config.imgDim; col_idx++)
				{
					image[row_idx*config.imgDim + col_idx] = (image[row_idx*config.imgDim + col_idx] - config.waterMu) / config.waterMu  * 1000.0f;
				}
			}
		}

		SaveReconImageSlice(filename, image, z_idx, config);

		
	}
}
