#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <string>
#include <vector>

// free memory
void FreeMemory_Agent(float* p);

// Initialize u, the array of each detector element coordiante
// u: array of detector elements
// N: number of detector elements
// du: detector element size [mm]
// offcenter: detector off-center [mm]
void InitializeU_Agent(float* u, const int N, const float du, const float offcenter);

// Initialize beta, the array of each view angle
// beta: array of view angles [radius]
// V: number of views
// rotation: rotate the reconstructed image [degree]
void InitializeBeta_Agent(float* beta, const int V, const float rotation);

// Initialize reconstruction kernel
// reconKernel: array of reconstruction kernel
// N: number of detector elements
// du: detector element size
void InitializeReconKernl_Agent(float* reconKernel, const int N, const float du, const std::string& kernelName, const std::vector<float>& kernelParam);