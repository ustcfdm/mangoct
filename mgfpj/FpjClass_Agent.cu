#include "FpjClass_Agent.cuh"

#define PI 3.1415926536f

__global__ void InitU(float* u, const int N, const float du, const float offcenter)
{
	int tid = threadIdx.x + blockDim.x * blockIdx.x;
	if (tid < N)
	{
		u[tid] = (tid - (N - 1) / 2.0f) * du + offcenter;
	}
}

__global__ void InitBeta(float* beta, const int V)
{
	int tid = threadIdx.x + blockDim.x * blockIdx.x;
	if (tid<V)
	{
		beta[tid] = (360.0f / V * tid) * PI / 180;
	}
}

void InitializeU_Agent(float* &u, const int N, const float du, const float offcenter)
{
	if (u != nullptr)
		cudaFree(u);

	cudaMalloc((void**)&u, N * sizeof(float));
	InitU << <(N + 511) / 512, 512 >> > (u, N, du, offcenter);
}

void InitializeBeta_Agent(float *& beta, const int V)
{
	if (beta != nullptr)
		cudaFree(beta);

	cudaMalloc((void**)&beta, V * sizeof(float));
	InitBeta << < (V + 511) / 512, 512 >> > (beta, V);
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
