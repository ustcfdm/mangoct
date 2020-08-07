#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "FpjClass.h"

// Initialize sdd or sid, the array of sdd or sid across views
// V: number of views
void InitializeDistance_Agent(float* &distance_array, const float distance, const int V);

// Initialize nonuniform sdd, the array of sdd across views
// V: number of views
void InitializeNonuniformSDD_Agent(float* &distance_array, const int V, const std::string& distanceFile);

// Initialize nonuniform sid, the array of sid across views
// V: number of views
void InitializeNonuniformSID_Agent(float* &distance_array, const int V, const std::string& distanceFile);

void InitializeNonuniformOffCenter_Agent(float* &offcenter_array, const int V, const std::string& offCenterFile);

void InitializeNonuniformPara_Agent(float* &para_array, const int V, const std::string& paraFile);

// Initialize u, the array of each detector element coordiante
// u: array of detector elements
// N: number of detector elements
// du: detector element size [mm]
// offcenter: detector off-center [mm]
void InitializeU_Agent(float* &u, const int N, const float du, const float offcenter);

// Initialize beta, the array of each view angle
// beta: array of view angles [radius]
// V: number of views
void InitializeBeta_Agent(float* &beta, const int V, const float startAngle, const float totalScanAngle);

// Initialize beta from an external jsonc file
// The purpose is for non uniform beta distribution
// V: number of views
// rotation: rotate the reconstructed image [degree]
// scanAngleFile: name of the jsonc file to save the angles 
void InitializeNonuniformBeta_Agent(float* &beta, const int V, const float rotation, const std::string& scanAngleFile);

// Forward projection, using bilinear interpolation
void ForwardProjectionBilinear_Agent(float* &image, float* &sinogram, const float* sid_array, const float* sdd_array, const float* offcenter_array,\
	const float* u,const float* v, const float* beta, const float* swing_angle_array, const mango::Config& config, int z_element_idx);

// Bin the sinogram data along detector direction
void BinSinogram(float* &sinogram_large, float* &sinogram, const mango::Config& config);

// Save one slice of the sinogram data
void SaveSinogramSlice(const char * filename, float*&sinogram_slice, int z_element_idx, const mango::Config& config);

// Malloc the memory as a given size
void MallocManaged_Agent(float* &p, const int size);

// free memory
void FreeMemory_Agent(float* &p);