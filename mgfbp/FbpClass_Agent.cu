#include "FbpClass_Agent.cuh"
#include <stdio.h>

#define PI 3.1415926536f

__global__ void InitU(float* u, const int N, const float du, const float offcenter)
{
	int tid = threadIdx.x + blockDim.x * blockIdx.x;
	if (tid < N)
	{
		u[tid] = (tid - (N - 1) / 2.0f) * du + offcenter;
	}
}

__global__ void InitBeta(float* beta, const int V, const float rotation)
{
	int tid = threadIdx.x + blockDim.x * blockIdx.x;
	if (tid<V)
	{
		beta[tid] = (360.0f / V * tid + rotation) * PI / 180;
	}
}

__global__ void InitReconKernel_Hamming(float* reconKernel, const int N, const float du, const float t)
{
	int tid = threadIdx.x + blockDim.x * blockIdx.x;
	if (tid < 2* N - 1)
	{
		// the center element index is N-1
		int n = tid - (N - 1);

		// ramp part
		if (n==0)
			reconKernel[tid] = t / (4 * du*du);
		else if (n%2 ==0)
			reconKernel[tid] = 0;
		else
			reconKernel[tid] = -t / (n*n * PI*PI * du*du);

		// cosine part
		int sgn = n % 2 == 0 ? 1 : -1;

		reconKernel[tid] += (1 - t)* (sgn / (2 * PI*du*du) * (1.0f / (1 + 2 * n) + 1.0f / (1 - 2 * n))
			- 1 / (PI*PI*du*du) * (1.0f / (1 + 2 * n) / (1 + 2 * n) + 1.0f / (1 - 2 * n) / (1 - 2 * n)));
	}
}

__global__ void InitReconKernel_Quadratic(float* reconKernel, const int N, const float du, const int paramNum, const float p1, const float p2, const float p3)
{
	int idx = threadIdx.x + blockDim.x * blockIdx.x;

	if (idx < 2*N )
	{
		float a, b, c;

		float kn = 1 / (2 * du);

		if (paramNum==2)
		{
			// p1 = t, p2 = h, p3 is ignored
			a = (p2 - 1) / (kn*kn * (1 - 2 * p1));
			b = -2 * a*p1*kn;
			c = 1.0f;
		}
		else
		{
			a = p1;
			b = p2;
			c = p3;
		}

		reconKernel[idx] = 0.0f;

		float du2 = du * du;
		float du3 = du2 * du;
		float du4 = du3 * du;

		int n = idx - (N - 1);
		if (n==0)
		{
			// H3(x)
			reconKernel[idx] += a / 32 / du4;
			// H2(x)
			reconKernel[idx] += b / 12 / du3;
			// H1(x)
			reconKernel[idx] += c / 4 / du2;
		}
		else if (n%2==0)
		{
			// H3(x)
			reconKernel[idx] += a * 3 / (8 * n*n * PI*PI * du4);
			// H2(x)
			reconKernel[idx] += b / (2 * n*n * PI*PI * du3);
			// H1(x)
			// do nothing, H1(even) is zero
		}
		else
		{
			// H3(x)
			reconKernel[idx] += a * 3 / (8 * n*n * PI*PI * du4) *  (4 /(n*n*PI*PI) - 1);
			// H2(x)
			reconKernel[idx] += -b / (2 * n*n * PI*PI * du3);
			// H1(x)
			reconKernel[idx] += -c / (n*n * PI*PI * du2);
		}
	}
}



void InitializeU_Agent(float* u, const int N, const float du, const float offcenter)
{
	if (u!=nullptr)
		cudaFree(u);

	cudaMalloc((void**)&u, N * sizeof(float));
	InitU <<<(N + 511) / 512, 512 >>> (u, N, du, offcenter);
}

void InitializeBeta_Agent(float* beta, const int V, const float rotation)
{
	if (beta!=nullptr)
		cudaFree(beta);
	
	cudaMalloc((void**)&beta, V * sizeof(float));
	InitBeta<<< (V+511)/512, 512>>> (beta, V, rotation);
}

void InitializeReconKernl_Agent(float* reconKernel, const int N, const float du, const std::string& kernelName, const std::vector<float>& kernelParam)
{
	if (reconKernel!=nullptr)
		cudaFree(reconKernel);

	cudaMalloc((void**)&reconKernel, (2 * N - 1) * sizeof(float));

	if (kernelName=="HammingFilter")
	{
		InitReconKernel_Hamming << <(2 * N - 1 + 511) / 512, 512 >> > (reconKernel, N, du, kernelParam[0]);
	}
	else if (kernelName=="QuadraticFilter")
	{
		float lastParam = 0.0f;
		if (kernelParam.size() == 3)
			lastParam = kernelParam[2];

		InitReconKernel_Quadratic << <(2 * N - 1 + 511) / 512, 512 >> > (reconKernel, N, du, int(kernelParam.size()), kernelParam[0], kernelParam[1], lastParam);
	}
}

void FreeMemory_Agent(float* p)
{
	cudaFree(p);
	p = nullptr;
}