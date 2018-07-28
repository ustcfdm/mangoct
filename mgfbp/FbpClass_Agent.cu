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

	if (idx < 2*N -1 )
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

__global__ void InitReconKernel_Polynomial(float* reconKernel, const int N, const float du, const float p6, const float p5, const float p4, const float p3, const float p2, const float p1, const float p0)
{
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	if (idx < 2 * N -1)
	{
		int n = idx - (N - 1);
		reconKernel[idx] = 0.0f;
		float kn = 1 / (2 * du);

		float du2 = du * du;
		float du3 = du2 * du;
		float du4 = du3 * du;

		if (n==0)
		{
			// H7(x)
			reconKernel[idx] += p6 * powf(kn, 8) / 4;
			// H6(x)
			reconKernel[idx] += p5 * powf(kn, 7) * 2 / 7;
			// H5(x)
			reconKernel[idx] += p4 * powf(kn, 6) / 3;
			// H4(x)
			reconKernel[idx] += p3 * powf(kn, 5) * 2 /5;
			// H3(x)
			reconKernel[idx] += p2 * powf(kn,4) / 2;
			// H2(x)
			reconKernel[idx] += p1 * 2 * kn*kn*kn / 3;
			// H1(x)
			reconKernel[idx] += p0 * kn*kn;
		}
		else if (n%2==0)
		{
			// H7(x)
			reconKernel[idx] += p6 * 7 * (360 - 30 * n*n*PI*PI + powf(n*PI, 4)) / (128 * du2* powf(du*n*PI, 6));
			// H6(x)
			reconKernel[idx] += p5 * 3 * (120 - 20 * n*n*PI*PI + powf(n*PI, 4)) / (32 * du*powf(du*n*PI, 6));
			// H5(x)
			reconKernel[idx] += p4 * 5 * (n*n*PI*PI - 12) / (32 * du2 *powf(du*n*PI, 4));
			// H4(x)
			reconKernel[idx] += p3 * (n*n*PI*PI - 6) / (4 * du * powf(du*n*PI, 4));
			// H3(x)
			reconKernel[idx] += p2 * 3 / (8 * du4 * n*n * PI*PI);
			// H2(x)
			reconKernel[idx] += p1 / (2 * n*n *PI*PI * du3);
			// H1(x)
			// do nothing, H1(even) is zero
		}
		else
		{
			// H7(x)
			reconKernel[idx] += p6 * 7 * (1440 - 360 * n*n*PI*PI + 30 * powf(n*PI, 4) - powf(n*PI, 6)) / (128 * powf(du*n*PI, 8));
			// H6(x)
			reconKernel[idx] += -p5 * 3 * (120 - 20 * n*n*PI*PI + powf(n*PI, 4)) / (32 * du*powf(du*n*PI, 6));
			// H5(x)
			reconKernel[idx] += -p4 * 5 * (48 - 12 * n*n*PI*PI + powf(n*PI, 4)) / (32 * powf(du*n*PI, 6));
			// H4(x)
			reconKernel[idx] += p3 * (6 - n*n*PI*PI) / (4 * du * powf(du*n*PI, 4));
			// H3(x)
			reconKernel[idx] += p2 * (4 - n * n*PI*PI) * 3 / (8 * powf(du*n*PI, 4));
			// H2(x)
			reconKernel[idx] += -p1 / (2 * n*n *PI*PI * du3);
			// H1(x)
			reconKernel[idx] += -p0 / (n*n *PI*PI * du2);
		}
	}
}


// weight the sinogram data
// sgm: sinogram (width x height x slice)
// N: width
// V: height (views)
// S: slice
// sdd: source to detector distance
__global__ void WeightSinogram_device(float* sgm, const float* u, const int N, const int V, const int S, float sdd)
{
	int col = threadIdx.x + blockDim.x * blockIdx.x;
	int row = threadIdx.y + blockDim.y * blockIdx.y;

	if (col<N && row <V)
	{
		for (int i = 0; i < S; i++)
		{
			sgm[row*N + col + i * N*V] *= sdd * sdd / sqrtf(u[col] * u[col] + sdd * sdd);
		}
	}
}

// convolve the sinogram data
// sgm_flt: sinogram data after convolving
// sgm: initial sinogram data
// reconKernel: reconstruction kernel
// N: sinogram width
// H: sinogram height
// V: number of views
// S: number of slices
// u: the position (coordinate) of each detector element
// du: detector element size [mm]
__global__ void ConvolveSinogram_device(float* sgm_flt, const float* sgm, float* reconKernel, const int N, const int H, const int V, const int S, const float* u, const float du)
{
	int col = threadIdx.x + blockDim.x * blockIdx.x;
	int row = threadIdx.y + blockDim.y * blockIdx.y;

	if (col < N && row<V)
	{
		for (int slice = 0; slice < S; slice++)
		{
			sgm_flt[row*N + col + slice * N*V] = 0;

			for (int i = 0; i < N; i++)
			{
				sgm_flt[row*N + col + slice * N*V] += sgm[row*N + i + slice * N*H] * reconKernel[N - 1 - col + i];
			}
			sgm_flt[row*N + col + slice * N*V] *= du;
		}

	}
}

// backproject the image using pixel-driven method
// sgm: sinogram data
// img: image data
// U: each detector element position [mm]
// beta: view angle [radius]
// N: number of detector elements
// V: number of views
// S: number of slices
// M: image dimension
// sdd: source to detector distance [mm]
// sid: source to isocenter distance [mm]
// du: detector element size [mm]
// dx: image pixel size [mm]
// (xc, yc): image center position [mm, mm]
__global__ void BackprojectPixelDriven_device(float* sgm, float* img, float* u, float* beta, const int N, const int V, const int S, const int M, const float sdd, const float sid, const float du, const float dx, const float xc, const float yc)
{
	int col = threadIdx.x + blockDim.x * blockIdx.x;
	int row = threadIdx.y + blockDim.y * blockIdx.y;

	if (col<M && row<M)
	{
		float x = (col - (M - 1) / 2.0f)*dx + xc;
		float y = ((M - 1) / 2.0f - row)*dx + yc;

		float U, u0;
		float w;
		int k;

		for (int slice = 0; slice < S; slice++)
		{
			img[row*M + col + slice * M*M] = 0;

			for (int view = 0; view < V; view++)
			{
				U = sid - x * cosf(beta[view]) - y * sinf(beta[view]);
				u0 = sdd * (x*sinf(beta[view]) - y * cosf(beta[view])) / U;

				k = floorf((u0 - u[0]) / du);
				if (k<0 || k+1>N-1)
				{
					img[row*M + col + slice * M*M] = 0;
					break;
				}

				w = u0 - u[k];

				img[row*M + col + slice * M*M] += sid / U / U * (w*sgm[view*N + k + 1 + slice * N*V] + (1 - w)*sgm[view*N + k + slice * N*V]);

			}
			img[row*M + col + slice * M*M] *= PI / V;
		}
	}
}


void InitializeU_Agent(float* &u, const int N, const float du, const float offcenter)
{
	if (u!=nullptr)
		cudaFree(u);

	cudaMalloc((void**)&u, N * sizeof(float));
	InitU <<<(N + 511) / 512, 512 >>> (u, N, du, offcenter);
}

void InitializeBeta_Agent(float* &beta, const int V, const float rotation)
{
	if (beta!=nullptr)
		cudaFree(beta);
	
	cudaMalloc((void**)&beta, V * sizeof(float));
	InitBeta<<< (V+511)/512, 512>>> (beta, V, rotation);
}

void InitializeReconKernel_Agent(float* &reconKernel, const int N, const float du, const std::string& kernelName, const std::vector<float>& kernelParam)
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
	else if (kernelName=="Polynomial")
	{
		// TODO: 
		// InitReconKernel_Polynomial <<<...>>> (...);
		float p[7] = { 0 };

		for (size_t i = 0; i < kernelParam.size(); i++)
		{
			p[i] = kernelParam[kernelParam.size() -1 - i];
		}

		//InitReconKernel_Polynomial << <(2 * N - 1 + 511) / 512, 512 >> > (reconKernel, N, du, p[0], p[1], p[2], p[3], p[4], p[5], p[6]);
		InitReconKernel_Polynomial << <(2 * N - 1 + 511) / 512, 512 >> > (reconKernel, N, du, p[6], p[5], p[4], p[3], p[2], p[1], p[0]);
	}
}

void MallocManaged_Agent(float * &p, const int size)
{
	cudaMallocManaged((void**)&p, size);
}

void FilterSinogram_Agent(float * sgm, float* sgm_flt, float* reconKernel, float* u, mango::Config & config)
{
	// Step 1: weight the sinogram
	dim3 grid((config.sgmWidth + 15) / 16, (config.sgmHeight + 15) / 16);
	dim3 block(16, 16);

	WeightSinogram_device << <grid, block >> >(sgm, u, config.sgmWidth, config.sgmHeight, config.sliceCount, config.sdd);

	// Step 2: convolve the sinogram
	ConvolveSinogram_device << <grid, block >> > (sgm_flt, sgm, reconKernel, config.sgmWidth, config.sgmHeight, config.views, config.sliceCount, u, config.detEltSize);

	cudaDeviceSynchronize();
}

void BackprojectPixelDriven_Agent(float * sgm_flt, float * img, float * u, float* beta, mango::Config & config)
{
	dim3 grid((config.imgDim + 15) / 16, (config.imgDim + 15) / 16);
	dim3 block(16, 16);

	BackprojectPixelDriven_device<<<grid,block>>>(sgm_flt, img, u, beta, config.sgmWidth, config.views, config.sliceCount, config.imgDim, config.sdd, config.sid, config.detEltSize, config.pixelSize, config.xCenter, config.yCenter);

	cudaDeviceSynchronize();
}


void FreeMemory_Agent(float* &p)
{
	cudaFree(p);
	p = nullptr;
}