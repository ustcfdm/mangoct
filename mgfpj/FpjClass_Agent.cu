#include "FpjClass_Agent.cuh"
#include <stdio.h>
#include "stdafx.h"

#define PI 3.1415926536f
#define STEPSIZE 0.2f

__global__ void InitDistance(float *distance_array, const float distance, const int V)
{
	int tid = threadIdx.x + blockDim.x * blockIdx.x;
	if (tid < V)
	{
		distance_array[tid] = distance;
	}
}

__global__ void InitU(float* u, const int N, const float du, const float offcenter)
{
	int tid = threadIdx.x + blockDim.x * blockIdx.x;
	if (tid < N)
	{
		u[tid] = (tid - (N - 1) / 2.0f) * du + offcenter;
	}
}

__global__ void InitBeta(float* beta, const int V, const float startAngle, const float totalScanAngle)
{
	int tid = threadIdx.x + blockDim.x * blockIdx.x;
	if (tid<V)
	{
		beta[tid] = (totalScanAngle / V * tid + startAngle) * PI / 180.0f;
	}
}

// img: image data
// sgm: sinogram data
// u: array of each detector element position
// beta: array of each view angle [radian]
// M: image dimension
// S: number of image slices
// N_z: number of detector elements in Z direction
// N: number of detector elements (sinogram width)
// V: number of views (sinogram height)
// dx: image pixel size [mm]
// dz: image slice thickness [mm]
// sid: source to isocenter distance
// sdd: source to detector distance
__global__ void ForwardProjectionBilinear_device(float* img, float* sgm, const float* u, const float *v, const float* offcenter_array, const float* beta, int M, int S,\
	int N, int N_z, int V, float dx, float dz, const float* sid_array, const float* sdd_array, bool conebeam, \
	int z_element_begin_idx, int z_element_end_idx)
{
	int col = threadIdx.x + blockDim.x * blockIdx.x;
	int row = threadIdx.y + blockDim.y * blockIdx.y;


	if (col < N && row < V && z_element_end_idx <= N_z)
	{
		// half of image side length
		float D = M * dx / 2.0f;
		// half of image thickness
		float D_z = 0.0f;
		if (conebeam)
		{
			D_z = float(S) * dz / 2.0f;
			
		}
		else
		{
			dz = 0;
		}
		//get the sid and sdd for a given view
		float sid = sid_array[row];
		float sdd = sdd_array[row];

		// current source position
		float xs = sid * cosf(beta[row]);
		float ys = sid * sinf(beta[row]);
		float zs = 0;

		// calculate offcenter bias
		float offcenter_bias = offcenter_array[row] - offcenter_array[0];

		// current detector element position
		float xd = -(sdd - sid) * cosf(beta[row]) + (u[col]+ offcenter_bias) * sinf(beta[row]);
		float yd = -(sdd - sid) * sinf(beta[row]) - (u[col]+ offcenter_bias) * cosf(beta[row]);
		float zd = 0;

		// step point region
		float L_min = sid - sqrt(2 * D * D + D_z * D_z);
		float L_max = sid + sqrt(2 * D * D + D_z * D_z);

		// source to detector element distance
		float sed = sqrtf((xs - xd)*(xs - xd) + (ys - yd)*(ys - yd));// for fan beam case

		// the point position
		float x, y, z;
		// the point index
		int kx, ky, kz;
		// weighting factor for linear interpolation
		float wx, wy, wz;

		// the most upper left image pixel position
		float x0 = -D + dx / 2.0f;
		float y0 = D - dx / 2.0f;
		float z0 = 0;
		if (conebeam)
		{
			z0 = -D_z + dz / 2.0f;// first slice is at the bottom
		}

		// repeat for each slice
		for (int slice = z_element_begin_idx; slice < z_element_end_idx; slice++)
		{
			// initialization
			//sgm[row*N + col + N * V * slice] = 0;
			sgm[row*N + col] = 0;
			if (conebeam)
			{
				
				zd = v[slice];
				
				sed = sqrtf((xs - xd)*(xs - xd) + (ys - yd)*(ys - yd) + (zs - zd)*(zs - zd));
			}

			// calculate line integration
			for (float L = L_min; L <= L_max; L+= STEPSIZE*sqrt(dx*dx+dz*dz/2.0f))
			{
				// get the current point position 
				x = xs + (xd - xs) * L / sed;
				y = ys + (yd - ys) * L / sed;
				if (conebeam)
				{
					z = zs + (zd - zs) * L / sed;
				}

				// get the current point index
				kx = floorf((x - x0) / dx);
				ky = floorf((y0 - y) / dx);

				if (conebeam)
					kz = floorf((z - z0) / dz);

				// get the image pixel value at the current point
				if(kx>=0 && kx+1<M && ky>=0 && ky+1<M)
				{
					// get the weighting factor
					wx = (x - kx * dx - x0) / dx;
					wy = (y0 - y - ky * dx) / dx;

					// perform bilinear interpolation
					if (conebeam == false)
					{
						sgm[row*N + col] += (1 - wx)*(1 - wy)*img[ky*M + kx + M * M*slice] // upper left
							+ wx * (1 - wy) * img[ky*M + kx + 1 + M * M*slice] // upper right
							+ (1 - wx) * wy * img[(ky + 1)*M + kx + M * M*slice] // bottom left
							+ wx * wy * img[(ky + 1)*M + kx + 1 + M * M*slice];	// bottom right
					}
					else if (conebeam == true && kz >= 0 && kz + 1 < S)
					{
						wz = (z - kz * dz - z0) / dz;
						float sgm_val_lowerslice = (1 - wx)*(1 - wy)*img[ky*M + kx + M * M*kz] // upper left
							+ wx * (1 - wy) * img[ky*M + kx + 1 + M * M*kz] // upper right
							+ (1 - wx) * wy * img[(ky + 1)*M + kx + M * M*kz] // bottom left
							+ wx * wy * img[(ky + 1)*M + kx + 1 + M * M*kz];	// bottom right
						float sgm_val_upperslice = (1 - wx)*(1 - wy)*img[ky*M + kx + M * M*(kz+1)] // upper left
							+ wx * (1 - wy) * img[ky*M + kx + 1 + M * M*(kz + 1)] // upper right
							+ (1 - wx) * wy * img[(ky + 1)*M + kx + M * M*(kz + 1)] // bottom left
							+ wx * wy * img[(ky + 1)*M + kx + 1 + M * M*(kz + 1)];	// bottom right

						sgm[row*N + col] += (1 - wz)*sgm_val_lowerslice + wz * sgm_val_upperslice;
					}
					
				}
			}

			sgm[row*N + col] *= STEPSIZE * sqrt(dx*dx + dz * dz);

		}
	}
}

// sgm_large: sinogram data before binning
// sgm: sinogram data after binning
// N: number of detector elements (after binning)
// V: number of views
// S: number of slices
// binSize: bin size
__global__ void BinSinogram(float* sgm_large, float* sgm, int N, int V, int S, int binSize)
{
	int col = threadIdx.x + blockDim.x * blockIdx.x;
	int row = threadIdx.y + blockDim.y * blockIdx.y;

	if (col < N && row < V)
	{
		// repeat for each slice
		for (int slice = 0; slice < S; slice++)
		{
			// initialization
			sgm[row * N + col + N * V * slice] = 0;

			// sum over each bin
			for (int i = 0; i < binSize; i++)
			{
				sgm[row * N + col + N * V * slice] += sgm_large[row * N * binSize + col*binSize + i + slice * N * binSize * V];
			}
			// take average
			sgm[row * N + col + N * V * slice] /= binSize;
		}
	}
}

void InitializeDistance_Agent(float* &distance_array, const float distance, const int V)
{
	if (distance_array != nullptr)
		cudaFree(distance_array);

	cudaMalloc((void**)&distance_array, V * sizeof(float));
	InitDistance << <(V + 511) / 512, 512 >> > (distance_array, distance, V);
}

void InitializeNonuniformSDD_Agent(float* &distance_array, const int V, const std::string& distanceFile)
{
	namespace fs = std::experimental::filesystem;
	namespace js = rapidjson;

	if (distance_array != nullptr)
		cudaFree(distance_array);

	cudaMalloc((void**)&distance_array, V * sizeof(float));

	float* distance_array_cpu = new float[V];
	std::ifstream ifs(distanceFile);
	if (!ifs)
	{
		printf("Cannot find SDD information file '%s'!\n", distanceFile.c_str());
		exit(-2);
	}
	rapidjson::IStreamWrapper isw(ifs);
	rapidjson::Document doc;
	doc.ParseStream<js::kParseCommentsFlag | js::kParseTrailingCommasFlag>(isw);
	js::Value distance_jsonc_value;
	if (doc.HasMember("SourceDetectorDistance"))
	{

		distance_jsonc_value = doc["SourceDetectorDistance"];

		if (distance_jsonc_value.Size() != V)
		{
			printf("Number of sdd values is %d while the number of Views is %d!\n", distance_jsonc_value.Size(), V);
			exit(-2);
		}

		for (unsigned i = 0; i < distance_jsonc_value.Size(); i++)
		{
			distance_array_cpu[i] = distance_jsonc_value[i].GetFloat();
		}

	}
	else
	{
		printf("Did not find SourceDetectorDistance member in jsonc file!\n");
		exit(-2);
	}

	cudaMemcpy(distance_array, distance_array_cpu, sizeof(float)*V, cudaMemcpyHostToDevice);

	cudaDeviceSynchronize();
}

void InitializeNonuniformSID_Agent(float* &distance_array, const int V, const std::string& distanceFile)
{
	namespace fs = std::experimental::filesystem;
	namespace js = rapidjson;

	if (distance_array != nullptr)
		cudaFree(distance_array);

	cudaMallocManaged((void**)&distance_array, V * sizeof(float));
	std::ifstream ifs(distanceFile);
	if (!ifs)
	{
		printf("Cannot find SID information file '%s'!\n", distanceFile.c_str());
		exit(-2);
	}
	rapidjson::IStreamWrapper isw(ifs);
	rapidjson::Document doc;
	doc.ParseStream<js::kParseCommentsFlag | js::kParseTrailingCommasFlag>(isw);
	js::Value distance_jsonc_value;
	if (doc.HasMember("SourceIsocenterDistance"))
	{

		distance_jsonc_value = doc["SourceIsocenterDistance"];

		if (distance_jsonc_value.Size() != V)
		{
			printf("Number of sid values is %d while the number of Views is %d!\n", distance_jsonc_value.Size(), V);
			exit(-2);
		}

		for (unsigned i = 0; i < distance_jsonc_value.Size(); i++)
		{
			distance_array[i] = distance_jsonc_value[i].GetFloat();
		}

	}
	else
	{
		printf("Did not find SourceIsocenterDistance member in jsonc file!\n");
		exit(-2);
	}

	cudaDeviceSynchronize();
}

void InitializeNonuniformOffCenter_Agent(float* &offcenter_array, const int V, const std::string& offCenterFile)
{
	namespace fs = std::experimental::filesystem;
	namespace js = rapidjson;

	if (offcenter_array != nullptr)
		cudaFree(offcenter_array);

	cudaMallocManaged((void**)&offcenter_array, V * sizeof(float));
	std::ifstream ifs(offCenterFile);
	if (!ifs)
	{
		printf("Cannot find Offcenter information file '%s'!\n", offCenterFile.c_str());
		exit(-2);
	}
	rapidjson::IStreamWrapper isw(ifs);
	rapidjson::Document doc;
	doc.ParseStream<js::kParseCommentsFlag | js::kParseTrailingCommasFlag>(isw);
	js::Value distance_jsonc_value;
	if (doc.HasMember("OffcenterArray"))
	{

		distance_jsonc_value = doc["OffcenterArray"];

		if (distance_jsonc_value.Size() != V)
		{
			printf("Number of offcenter values is %d while the number of Views is %d!\n", distance_jsonc_value.Size(), V);
			exit(-2);
		}

		for (unsigned i = 0; i < distance_jsonc_value.Size(); i++)
		{
			offcenter_array[i] = distance_jsonc_value[i].GetFloat();
		}

	}
	else
	{
		printf("Did not find OffcenterArray member in jsonc file!\n");
		exit(-2);
	}

	cudaDeviceSynchronize();
}

void InitializeU_Agent(float* &u, const int N, const float du, const float offcenter)
{
	if (u != nullptr)
		cudaFree(u);

	cudaMalloc((void**)&u, N * sizeof(float));
	InitU <<<(N + 511) / 512, 512 >>> (u, N, du, offcenter);
}

void InitializeBeta_Agent(float *& beta, const int V, const float startAngle, const float totalScanAngle)
{
	if (beta != nullptr)
		cudaFree(beta);

	cudaMalloc((void**)&beta, V * sizeof(float));
	InitBeta <<< (V + 511) / 512, 512 >>> (beta, V, startAngle, totalScanAngle);
}

void InitializeNonuniformBeta_Agent(float* &beta, const int V, const float rotation, const std::string& scanAngleFile)
{
	namespace fs = std::experimental::filesystem;
	namespace js = rapidjson;

	if (beta != nullptr)
		cudaFree(beta);

	cudaMallocManaged((void**)&beta, V * sizeof(float));
	std::ifstream ifs(scanAngleFile);
	if (!ifs)
	{
		printf("Cannot find angle information file '%s'!\n", scanAngleFile.c_str());
		exit(-2);
	}
	rapidjson::IStreamWrapper isw(ifs);
	rapidjson::Document doc;
	doc.ParseStream<js::kParseCommentsFlag | js::kParseTrailingCommasFlag>(isw);
	js::Value scan_angle_jsonc_value;
	if (doc.HasMember("ScanAngle"))
	{

		scan_angle_jsonc_value = doc["ScanAngle"];

		if (scan_angle_jsonc_value.Size() != V)
		{
			printf("Number of scan angles is %d while the number of Views is %d!\n",scan_angle_jsonc_value.Size(),V);
			exit(-2);
		}

		for (unsigned i = 0; i < scan_angle_jsonc_value.Size(); i++)
		{
			beta[i] = rotation / 180.0f*PI + scan_angle_jsonc_value[i].GetFloat() / 180.0*PI;
		}

	}
	else
	{
		printf("Did not find ScanAngle member in jsonc file!\n");
		exit(-2);
	}

	cudaDeviceSynchronize();
}

void ForwardProjectionBilinear_Agent(float *& image, float * &sinogram, const float* sid_array, const float* sdd_array, const float* offcenter_array,\
	const float* u, const float *v, const float* beta, const mango::Config & config, int z_element_idx)
{
	dim3 grid((config.detEltCount*config.oversampleSize + 7) / 8, (config.views + 7) / 8);
	dim3 block(8, 8);

	ForwardProjectionBilinear_device<<<grid, block>>>(image, sinogram, u, v, offcenter_array, beta, config.imgDim, config.sliceCount,\
		config.detEltCount*config.oversampleSize, config.detZEltCount, config.views, config.pixelSize, config.sliceThickness, sid_array, sdd_array, config.coneBeam, z_element_idx, z_element_idx+1);

	cudaDeviceSynchronize();
}

void BinSinogram(float* &sinogram_large, float* &sinogram, const mango::Config& config)
{
	dim3 grid((config.detEltCount + 7) / 8, (config.views + 7) / 8);
	dim3 block(8, 8);
	
	BinSinogram <<<grid, block >>> (sinogram_large, sinogram, config.detEltCount, config.views, 1, config.oversampleSize);
	// since the sinogram has only one slice, the z_element count is 1

	cudaDeviceSynchronize();
}

void SaveSinogramSlice(const char * filename, float*&sinogram_slice, int z_element_idx, const mango::Config& config)
{
	FILE* fp = NULL;
	if (z_element_idx == 0)
		fp = fopen(filename, "wb");
	else
		fp = fopen(filename, "ab");

	if (fp == NULL)
	{
		fprintf(stderr, "Cannot save to file %s!\n", filename);
		exit(4);
	}
	fwrite(sinogram_slice, sizeof(float), config.detEltCount * config.views, fp);
	fclose(fp);
}

void MallocManaged_Agent(float * &p, const int size)
{
	cudaMallocManaged((void**)&p, size);
}

void FreeMemory_Agent(float* &p)
{
	cudaFree(p);
	p = nullptr;
}
