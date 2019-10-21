// mglog.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "mglog.h"

namespace fs = std::experimental::filesystem;
namespace mg = mango;
namespace js = rapidjson;

void LoadConfigFile(js::Document& d, Config& config);
std::vector<mg::Matrix> GetBkgData(const Config& config, const int fileId = 1);
std::vector<mg::Matrix> GetPrelogSinogram(const mg::Matrix& obj, const Config& config);
void SmoothoutMargianlData(std::vector<mg::Matrix>& sgm, const Config& config);

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
		if (!ifs)
		{
			printf("Cannot open config file '%s'!\n", argv[configIdx]);
			exit(-2);
		}
		js::IStreamWrapper isw(ifs);
		js::Document doc;
		doc.ParseStream<js::kParseCommentsFlag | js::kParseTrailingCommasFlag>(isw);

		LoadConfigFile(doc, config);

		//////////////////////////////////////////////////////////////////////////////
		// Step 2: acquire and process bakcground data
		//////////////////////////////////////////////////////////////////////////////

		printf("Processing background file: %s\n", config.bkgFile.c_str());
		std::vector<mg::Matrix> bkg_TE = GetBkgData(config);

		std::vector<mg::Matrix> bkg_HE;
		if (config.dualEnergySubtraction)
		{
			printf("Processing background file: %s\n", config.bkgFile2.c_str());
			bkg_HE = GetBkgData(config, 2);
		}


		//////////////////////////////////////////////////////////////////////////////
		// Step 2: acquire and process object data
		//////////////////////////////////////////////////////////////////////////////

		// If perform dual energy subtraction
		if (config.dualEnergySubtraction)
		{
			mg::Matrix obj_TE(config.detectorHeight, config.detectorWidth, config.objViews);
			mg::Matrix obj_LE(config.detectorHeight, config.detectorWidth, config.objViews);
			mg::Matrix obj_HE(config.detectorHeight, config.detectorWidth, config.objViews);


			// get background file
			std::vector<mg::Matrix> bkg_LE;
			for (size_t i = 0; i < bkg_TE.size(); i++)
			{
				bkg_LE.push_back(bkg_TE[i] - bkg_HE[i]);
			}

			// repeat for each output file
			for (size_t i = 0; i < config.outputFiles3.size(); i++)
			{
				printf("    Processing %s and %s ...", config.inputFiles[i].c_str(), config.inputFiles2[i].c_str());

				// read the evi file
				obj_TE.ReadEviFile((fs::path(config.inputDir) / config.inputFiles[i]).string().c_str(), config.offsetToFirstImage, config.gap);
				obj_HE.ReadEviFile((fs::path(config.inputDir) / config.inputFiles2[i]).string().c_str(), config.offsetToFirstImage, config.gap);

				// get the pre-log sinogram
				std::vector<mg::Matrix> sgm_TE = GetPrelogSinogram(obj_TE, config);
				std::vector<mg::Matrix> sgm_HE = GetPrelogSinogram(obj_HE, config);

				std::vector<mg::Matrix> sgm_LE;
				for (size_t i = 0; i < sgm_TE.size(); i++)
				{
					sgm_LE.push_back(sgm_TE[i] - sgm_HE[i]);
				}

				//-----------------------------------------------------
				// for "LE"
				//-----------------------------------------------------
				if (std::find(config.outputTypes.begin(),config.outputTypes.end(),"LE") != config.outputTypes.end())
				{
					// get the post-log sinogram
					for (size_t i = 0; i < sgm_LE.size(); i++)
					{
						sgm_LE[i] = (bkg_LE[i] / sgm_LE[i]).Log();
						sgm_LE[i].SetNanOrInf(0.0f);
					}


					// interpolate white lines

					if (config.interpolateWhiteLines)
					{
						for (size_t k = 0; k < sgm_LE.size(); k++)
						{
							for (unsigned row = 0; row < config.objViews; row++)
							{
								for (unsigned col = 255; col < 5119; col += 256)
								{
									sgm_LE[k](row, col) = sgm_LE[k](row, col - 1) * 2 / 3 + sgm_LE[k](row, col + 2) / 3;
									sgm_LE[k](row, col + 1) = sgm_LE[k](row, col - 1) / 3 + sgm_LE[k](row, col + 2) * 2 / 3;
								}
							}
						}
					}

					// rebin sinogram data
					for (size_t k = 0; k < sgm_LE.size(); k++)
					{
						sgm_LE[k].Rebin(config.rebinSize, mg::Axis::Col, true);
					}

					// smoothout marginal data
					if (config.smoothoutMarginalData)
					{
						SmoothoutMargianlData(sgm_LE, config);
					}

					// save to file
					std::string saveFullName = (fs::path(config.outputDir) / config.outputFiles3[i]).string();
					sgm_LE[0].SaveRawFile(saveFullName.c_str());
					for (size_t k = 1; k < sgm_LE.size(); k++)
					{
						sgm_LE[k].AppendRawFile(saveFullName.c_str());
					}
					printf("\t->\tSave to %s\n", config.outputFiles3[i].c_str());
				}

				//-----------------------------------------------------
				// for "TE"
				//-----------------------------------------------------
				if (std::find(config.outputTypes.begin(), config.outputTypes.end(), "TE") != config.outputTypes.end())
				{
					// Get post-log sinogram
				// get the post-log sinogram
					for (size_t i = 0; i < sgm_TE.size(); i++)
					{
						sgm_TE[i] = (bkg_TE[i] / sgm_TE[i]).Log();
						sgm_TE[i].SetNanOrInf(0.0f);
					}


					// interpolate white lines

					if (config.interpolateWhiteLines)
					{
						for (size_t k = 0; k < sgm_TE.size(); k++)
						{
							for (unsigned row = 0; row < config.objViews; row++)
							{
								for (unsigned col = 255; col < 5119; col += 256)
								{
									sgm_TE[k](row, col) = sgm_TE[k](row, col - 1) * 2 / 3 + sgm_TE[k](row, col + 2) / 3;
									sgm_TE[k](row, col + 1) = sgm_TE[k](row, col - 1) / 3 + sgm_TE[k](row, col + 2) * 2 / 3;
								}
							}
						}
					}

					// rebin sinogram data
					for (size_t k = 0; k < sgm_TE.size(); k++)
					{
						sgm_TE[k].Rebin(config.rebinSize, mg::Axis::Col, true);
					}

					// smoothout marginal data
					if (config.smoothoutMarginalData)
					{
						SmoothoutMargianlData(sgm_TE, config);
					}

					// save to file
					std::string saveFullName = (fs::path(config.outputDir) / config.outputFiles[i]).string();
					sgm_TE[0].SaveRawFile(saveFullName.c_str());
					for (size_t k = 1; k < sgm_TE.size(); k++)
					{
						sgm_TE[k].AppendRawFile(saveFullName.c_str());
					}
					printf("\t\t\t->\tSave to %s\n", config.outputFiles[i].c_str());
				}
				
				//-----------------------------------------------------
				// for "HE"
				//-----------------------------------------------------
				if (std::find(config.outputTypes.begin(), config.outputTypes.end(), "HE") != config.outputTypes.end())
				{
					// Get post-log sinogram
				// get the post-log sinogram
					for (size_t i = 0; i < sgm_HE.size(); i++)
					{
						sgm_HE[i] = (bkg_HE[i] / sgm_HE[i]).Log();
						sgm_HE[i].SetNanOrInf(0.0f);
					}


					// interpolate white lines

					if (config.interpolateWhiteLines)
					{
						for (size_t k = 0; k < sgm_HE.size(); k++)
						{
							for (unsigned row = 0; row < config.objViews; row++)
							{
								for (unsigned col = 255; col < 5119; col += 256)
								{
									sgm_HE[k](row, col) = sgm_HE[k](row, col - 1) * 2 / 3 + sgm_HE[k](row, col + 2) / 3;
									sgm_HE[k](row, col + 1) = sgm_HE[k](row, col - 1) / 3 + sgm_HE[k](row, col + 2) * 2 / 3;
								}
							}
						}
					}

					// rebin sinogram data
					for (size_t k = 0; k < sgm_HE.size(); k++)
					{
						sgm_HE[k].Rebin(config.rebinSize, mg::Axis::Col, true);
					}

					// smoothout marginal data
					if (config.smoothoutMarginalData)
					{
						SmoothoutMargianlData(sgm_HE, config);
					}

					// save to file
					std::string saveFullName = (fs::path(config.outputDir) / config.outputFiles2[i]).string();
					sgm_HE[0].SaveRawFile(saveFullName.c_str());
					for (size_t k = 1; k < sgm_HE.size(); k++)
					{
						sgm_HE[k].AppendRawFile(saveFullName.c_str());
					}
					printf("\t\t\t->\tSave to %s\n", config.outputFiles2[i].c_str());
				}
			}
		}
		// do not perform dual energy subtraction
		else
		{
			mg::Matrix obj(config.detectorHeight, config.detectorWidth, config.objViews);

			// repeat for each output file
			for (size_t i = 0; i < config.outputFiles.size(); i++)
			{
				printf("    Processing %s...", config.inputFiles[i].c_str());

				// read the evi file
				obj.ReadEviFile((fs::path(config.inputDir) / config.inputFiles[i]).string().c_str(), config.offsetToFirstImage, config.gap);

				// get the pre-log sinogram
				std::vector<mg::Matrix> sgm = GetPrelogSinogram(obj, config);

				// pre-log sinogram to post-log sinogram
				for (size_t k = 0; k < sgm.size(); k++)
				{
					sgm[k] = (bkg_TE[k] / sgm[k]).Log();
					sgm[k].SetNanOrInf(0.0f);
				}

				// interpolate white lines

				if (config.interpolateWhiteLines)
				{
					for (size_t k = 0; k < sgm.size(); k++)
					{
						for (unsigned row = 0; row < config.objViews; row++)
						{
							for (unsigned col = 255; col < 5119; col += 256)
							{
								sgm[k](row, col) = sgm[k](row, col - 1) * 2 / 3 + sgm[k](row, col + 2) / 3;
								sgm[k](row, col + 1) = sgm[k](row, col - 1) / 3 + sgm[k](row, col + 2) * 2 / 3;
							}
						}
					}
				}

				// rebin sinogram data
				for (size_t k = 0; k < sgm.size(); k++)
				{
					sgm[k].Rebin(config.rebinSize, mg::Axis::Col, true);
				}

				// smoothout marginal data
				if (config.smoothoutMarginalData)
				{
					SmoothoutMargianlData(sgm, config);
				}

				// save to file
				std::string saveFullName = (fs::path(config.outputDir) / config.outputFiles[i]).string();
				sgm[0].SaveRawFile(saveFullName.c_str());
				for (size_t k = 1; k < sgm.size(); k++)
				{
					sgm[k].AppendRawFile(saveFullName.c_str());
				}
				printf("\t->\tSave to %s\n", config.outputFiles[i].c_str());

			}
		}


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
		fprintf(stderr, "Output directory %s does not exist!\n", config.outputDir.c_str());
		exit(1);
	}

	// Get background files
	if (d["BackgroundFile"].IsArray())
	{
		const js::Value& bkgFiles = d["BackgroundFile"];
		config.bkgFile = bkgFiles[0].GetString();
		if (bkgFiles.Size() >= 2)
		{
			config.bkgFile2 = bkgFiles[1].GetString();
		}
	}
	else
	{
		config.bkgFile = d["BackgroundFile"].GetString();
	}

	// check background file
	if (!fs::exists(inDir / config.bkgFile))
	{
		fprintf(stderr, "Cannot find backgorund file %s!\n", config.bkgFile.c_str());
		exit(1);
	}
	if (config.bkgFile2 != "" && !fs::exists(inDir / config.bkgFile2))
	{
		fprintf(stderr, "Cannot find backgorund file %s!\n", config.bkgFile2.c_str());
		exit(1);
	}

	if (d["InputFiles"].IsArray())
	{
		config.inputFiles = GetInputFileNames(config.inputDir, d["InputFiles"][0].GetString());
		if (d["InputFiles"].Size() >= 2)
		{
			config.inputFiles2 = GetInputFileNames(config.inputDir, d["InputFiles"][1].GetString());
		}
	}
	else
	{
		config.inputFiles = GetInputFileNames(config.inputDir, d["InputFiles"].GetString());
	}

	// check input file names
	if (d["InputFiles"].IsArray())
	{
		if (config.inputFiles.size() == 0)
		{
			fprintf(stderr, "No file name like %s!\n", d["InputFiles"][0].GetString());
			exit(1);
		}
		if (d["InputFiles"].Size() >= 2 && config.inputFiles2.size() == 0)
		{
			fprintf(stderr, "No file name like %s!\n", d["InputFiles"][1].GetString());
			exit(1);
		}
	}
	else
	{
		if (config.inputFiles.size() == 0)
		{
			fprintf(stderr, "No file name like %s!\n", d["InputFiles"].GetString());
			exit(1);
		}
	}


	// Dual energy subtraction or not?
	if (d.HasMember("DualEnergySubtraction"))
		config.dualEnergySubtraction = d["DualEnergySubtraction"].GetBool();
	else
		config.dualEnergySubtraction = false;

	// Output types (i.e. "TE", "LE", or "HE")
	if (d.HasMember("OutputTypes"))
	{
		for (auto& v : d["OutputTypes"].GetArray())
		{
			config.outputTypes.push_back(v.GetString());
		}
	}

	// Subtraction file name replace
	std::vector<std::string> subtractionReplace;
	if (d.HasMember("SubtractionFileNameReplace"))
	{
		for (auto& v : d["SubtractionFileNameReplace"].GetArray())
		{
			subtractionReplace.push_back(v.GetString());
		}
	}

	// get output file names
	const js::Value& replaceJs = d["OutputFileReplace"];
	std::vector<std::string> replace;
	for (js::SizeType i = 0; i < replaceJs.Size(); i++)
	{
		replace.push_back(replaceJs[i].GetString());
	}

	config.outputFiles = GetOutputFileNames(config.inputFiles, replace, d["OutputFilePrefix"].GetString());
	if (config.dualEnergySubtraction)
	{
		// For "HE" files
		config.outputFiles2 = GetOutputFileNames(config.inputFiles2, replace, d["OutputFilePrefix"].GetString());
		// For "LE" files
		config.outputFiles3 = GetOutputFileNames(config.outputFiles, subtractionReplace, "");
	}


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

// acquire background file data, and process the data according to config infomation
std::vector<mg::Matrix> GetBkgData(const Config& config, const int fileId)
{
	// read background data file
	std::string tmpBkgFile = fileId == 1 ? config.bkgFile : config.bkgFile2;
	mg::Matrix bkg = mg::Matrix::ReadEviFile((fs::path(config.inputDir) / tmpBkgFile).string().c_str(), config.detectorHeight, config.detectorWidth, config.bkgViews, config.offsetToFirstImage, config.gap);

	// take the average along views(pages) direction
	bkg.Average(mg::Axis::Page);

	// the bkg data to be returned
	std::vector<mg::Matrix> b;

	unsigned idx = config.sliceStartIdx;
	for (unsigned i = 0; i < config.sliceCount; i++)
	{
		b.push_back(mg::Matrix::Average(bkg, mg::Axis::Row, idx, idx + config.sliceThickness));
		idx += config.sliceThickness;
	}

	return b;
}

std::vector<mg::Matrix> GetPrelogSinogram(const mg::Matrix& obj, const Config& config)
{
	std::vector<mg::Matrix> sgm;

	unsigned idx = config.sliceStartIdx;
	for (unsigned i = 0; i < config.sliceCount; i++)
	{
		sgm.push_back(mg::Matrix::Average(obj, mg::Axis::Row, idx, idx + config.sliceThickness).Reshape(config.objViews, config.detectorWidth));
		idx += config.sliceThickness;
	}

	return sgm;
}

void SmoothoutMargianlData(std::vector<mg::Matrix>& sgm, const Config& config)
{
	for (size_t slice = 0; slice < sgm.size(); slice++)
	{
		for (unsigned row = 0; row < config.objViews; row++)
		{
			// left part
			for (unsigned col = 0; col < config.smoothoutLeftIdx; col++)
				sgm[slice](row, col) = 0;

			// right part
			for (unsigned col = config.smoothoutRightIdx; col < sgm[slice].Cols(); col++)
				sgm[slice](row, col) = 0;
		}
	}
}