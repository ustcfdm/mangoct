#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <string>
#include <vector>

#include "FbpClass.cuh"

// free memory
void FreeMemory_Agent(float* &p);

// Initialize sdd or sid, the array of sdd or sid across views
// V: number of views
void InitializeDistance_Agent(float* &distance_array, const float distance, const int V);

// Initialize nonuniform sdd, the array of sdd across views
// V: number of views
void InitializeNonuniformSDD_Agent(float* &distance_array, const int V, const std::string& distanceFile);

// Initialize nonuniform sid, the array of sid across views
// V: number of views
void InitializeNonuniformSID_Agent(float* &distance_array, const int V, const std::string& distanceFile);

void InitializeNonuniformOffCenter_Agent(float* &distance_array, const int V, const std::string& distanceFile);

void InitializePMatrix_Agent(float* &pmatrix_array, const int V, const std::string& pmatrixFile);

// Initialize u, the array of each detector element coordiante
// u: array of detector elements
// N: number of detector elements
// du: detector element size [mm]
// offcenter: detector off-center [mm]
void InitializeU_Agent(float* &u, const int N, const float du, const float offcenter);

// Initialize beta, the array of each view angle
// beta: array of view angles [radius]
// V: number of views
// rotation: rotate the reconstructed image [degree]
// totalScanAngle: total scan angle for short scan [degree]
void InitializeBeta_Agent(float* &beta, const int V, const float rotation, const float totalScanAngle);

// Initialize beta from an external jsonc file
// The purpose is for non uniform beta distribution
// V: number of views
// rotation: rotate the reconstructed image [degree]
// scanAngleFile: name of the jsonc file to save the angles 
void InitializeNonuniformBeta_Agent(float* &beta, const int V, const float rotation, const std::string& scanAngleFile);

// Initialize reconstruction kernel
// reconKernel: array of reconstruction kernel
// N: number of detector elements
// du: detector element size
void InitializeReconKernel_Agent(float* &reconKernel, const int N, const float du, const std::string& kernelName, const std::vector<float>& kernelParam);

// Malloc the memory as a given size
void MallocManaged_Agent(float* &p, const int size);

// Perform beam hardening correction
void CorrectBeamHardening_Agent(float* sgm, mango::Config & config);

// Filter the sinogram data
void FilterSinogram_Agent(float* sgm, float* sgm_flt, float* reconKernel, float* u, mango::Config& config, float*beta, float*sdd_array, float * offcenter_array);

// Backproject the image using pixel-driven method
void BackprojectPixelDriven_Agent(float* sgm_flt, float* img, float* sdd_array, float* sid_array, float* offcenter_array, float* pmatrix_array, float* u, float *v, float* beta, mango::Config& config, int z_idx);

// Save reconstructed images slice by slice
void SaveReconImageSlice(const char* filename, float* rec_image, int z_idx, const mango::Config& config);
