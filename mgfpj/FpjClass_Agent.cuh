#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "FpjClass.h"

// Initialize u, the array of each detector element coordiante
// u: array of detector elements
// N: number of detector elements
// du: detector element size [mm]
// offcenter: detector off-center [mm]
void InitializeU_Agent(float* &u, const int N, const float du, const float offcenter);

// Initialize beta, the array of each view angle
// beta: array of view angles [radius]
// V: number of views
void InitializeBeta_Agent(float* &beta, const int V, const float startAngle);

// Forward projection, using bilinear interpolation
void ForwardProjectionBilinear_Agent(float* &image, float* &sinogram, const float* u, const float* beta, const mango::Config& config);

// Bin the sinogram data along detector direction
void BinSinogram(float* &sinogram_large, float* &sinogram, const mango::Config& config);

// Malloc the memory as a given size
void MallocManaged_Agent(float* &p, const int size);

// free memory
void FreeMemory_Agent(float* &p);