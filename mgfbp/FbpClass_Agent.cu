#include "FbpClass_Agent.cuh"
#include <stdio.h>
#include "stdafx.h"
#define PI 3.1415926536f


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

__global__ void InitBeta(float* beta, const int V, const float rotation, const float totalScanAngle)
{
	int tid = threadIdx.x + blockDim.x * blockIdx.x;
	if (tid < V)
	{
		beta[tid] = (totalScanAngle / V * tid + rotation) * PI / 180;
	}
}

__global__ void InitReconKernel_Hamming(float* reconKernel, const int N, const float du, const float t)
{
	int tid = threadIdx.x + blockDim.x * blockIdx.x;
	if (tid < 2 * N - 1)
	{
		// the center element index is N-1
		int n = tid - (N - 1);

		// ramp part
		if (n == 0)
			reconKernel[tid] = t / (4 * du*du);
		else if (n % 2 == 0)
			reconKernel[tid] = 0;
		else
			reconKernel[tid] = -t / (n*n * PI*PI * du*du);

		// cosine part
		int sgn = n % 2 == 0 ? 1 : -1;

		reconKernel[tid] += (1 - t)* (sgn / (2 * PI*du*du) * (1.0f / (1 + 2 * n) + 1.0f / (1 - 2 * n))
			- 1 / (PI*PI*du*du) * (1.0f / (1 + 2 * n) / (1 + 2 * n) + 1.0f / (1 - 2 * n) / (1 - 2 * n)));
	}
}

__global__ void InitReconKernel_Delta(float* reconKernel, const int N, const float du, const float t)
{
	int tid = threadIdx.x + blockDim.x * blockIdx.x;
	if (tid < 2 * N - 1)
	{
		// the center element index is N-1
		int n = tid - (N - 1);

		if (n == 0)
			reconKernel[tid] = t;
		else
			reconKernel[tid] = 0;
	}
}

//initialize a Gaussian kernel
//This kernel will be used along with the ramp kernel
//delta is in number of pixels, which is the standard deviation of the gaussian
//This kernel is normalized
__global__ void InitReconKernel_GaussianApodized(float* reconKernel, const int N, const float du, const float delta)
{
	int tid = threadIdx.x + blockDim.x * blockIdx.x;
	if (tid < 1)
	{
		// the center element index is N-1
		float temp_sum = 0;
		for (int i = 0; i < 2 * N - 1; i++)
		{
			int n = i - (N - 1);
			reconKernel[i] = exp(-float(n) * float(n) / 2.0 / delta / delta);
			temp_sum = temp_sum + reconKernel[i];
		}

		for (int i = 0; i < 2 * N - 1; i++)
		{
			reconKernel[i] = reconKernel[i] / temp_sum / du;
		}
	}
}


__global__ void InitReconKernel_Quadratic(float* reconKernel, const int N, const float du, const int paramNum, const float p1, const float p2, const float p3)
{
	int idx = threadIdx.x + blockDim.x * blockIdx.x;

	if (idx < 2 * N - 1)
	{
		float a, b, c;

		float kn = 1 / (2 * du);

		if (paramNum == 2)
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
		if (n == 0)
		{
			// H3(x)
			reconKernel[idx] += a / 32 / du4;
			// H2(x)
			reconKernel[idx] += b / 12 / du3;
			// H1(x)
			reconKernel[idx] += c / 4 / du2;
		}
		else if (n % 2 == 0)
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
			reconKernel[idx] += a * 3 / (8 * n*n * PI*PI * du4) *  (4 / (n*n*PI*PI) - 1);
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
	if (idx < 2 * N - 1)
	{
		int n = idx - (N - 1);
		reconKernel[idx] = 0.0f;
		float kn = 1 / (2 * du);

		float du2 = du * du;
		float du3 = du2 * du;
		float du4 = du3 * du;

		if (n == 0)
		{
			// H7(x)
			reconKernel[idx] += p6 * powf(kn, 8) / 4;
			// H6(x)
			reconKernel[idx] += p5 * powf(kn, 7) * 2 / 7;
			// H5(x)
			reconKernel[idx] += p4 * powf(kn, 6) / 3;
			// H4(x)
			reconKernel[idx] += p3 * powf(kn, 5) * 2 / 5;
			// H3(x)
			reconKernel[idx] += p2 * powf(kn, 4) / 2;
			// H2(x)
			reconKernel[idx] += p1 * 2 * kn*kn*kn / 3;
			// H1(x)
			reconKernel[idx] += p0 * kn*kn;
		}
		else if (n % 2 == 0)
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
			reconKernel[idx] += p3 * (6 - n * n*PI*PI) / (4 * du * powf(du*n*PI, 4));
			// H3(x)
			reconKernel[idx] += p2 * (4 - n * n*PI*PI) * 3 / (8 * powf(du*n*PI, 4));
			// H2(x)
			reconKernel[idx] += -p1 / (2 * n*n *PI*PI * du3);
			// H1(x)
			reconKernel[idx] += -p0 / (n*n *PI*PI * du2);
		}
	}
}

__global__ void InitReconKernel_Hilbert(float* reconKernel, const int N, const float du, const float t)
{
	int tid = threadIdx.x + blockDim.x * blockIdx.x;
	if (tid < 2 * N - 1)
	{
		int n = tid - (N - 1);

		if (n % 2 == 0)
			reconKernel[tid] = 0;
		else
		{
			reconKernel[tid] = 1 / (PI * PI * n * du);
			if (t < 0)
				reconKernel[tid] = -reconKernel[tid];
		}
	}
}

// weight the sinogram data
// sgm: sinogram (width x height x slice)
// N: width
// H: height
// V: views
// S: slice
// sliceThickness: mm
// sliceOffcenter: mm
// sdd: source to detector distance
// totalScanAngle
__global__ void WeightSinogram_device(float* sgm, const float* u, const int N, const int H, const int V, \
	const int S, const float sliceThickness, const float sliceOffcenter, float* sdd_array, float totalScanAngle, bool shortScan, float *beta_array, float* offcenter_array)
{
	int col = threadIdx.x + blockDim.x * blockIdx.x;
	int row = threadIdx.y + blockDim.y * blockIdx.y;

	if (col < N && row < V)
	{
		float offcenter_bias = offcenter_array[row] - offcenter_array[0];
		float u_actual = u[col] + offcenter_bias;//actual u value due to non uniform offcenter

		float sdd = sdd_array[row];
		for (int i = 0; i < S; i++)
		{
			float v = sliceThickness * (i - (float(S) / 2.0f + 0.5)) + sliceOffcenter;
			sgm[row*N + col + i * N*H] *= sdd * sdd / sqrtf((u_actual)*(u_actual)+sdd * sdd + v * v);
			//the loop is to include all the slices since there may be more than one slice
		}

		if (shortScan)
		{
			//this beta is different from the beta_array
			//To calculate the parker weighting, beta should begin with zero degree
			//while the betaArray includes the start rotation angle

			//adding abs function to deal with the case when totalScanAngle is negative
			float beta = abs(beta_array[row] - beta_array[0]);
			float rotation_direction = abs(totalScanAngle) / (totalScanAngle);
			float gamma = atan(u_actual / sdd) * rotation_direction;

			//float beta = abs(totalScanAngle/ 180.0f * PI ) / float(V) * float(row) ;
			//float gamma = abs(totalScanAngle) / totalScanAngle * atan(u[col] / sdd);
			float gamma_max = abs(totalScanAngle) * PI / 180.0f - PI;

			//calculation of the parker weighting
			float weighting = 0;
			if (beta >= 0 && beta < gamma_max - 2 * gamma)
			{
				weighting = sin(PI / 2 * beta / (gamma_max - 2 * gamma));
				weighting = weighting * weighting;
			}
			else if (beta >= gamma_max - 2 * gamma && beta < PI - 2 * gamma)
			{
				weighting = 1;
			}
			else if (beta >= PI - 2 * gamma && beta <= PI + gamma_max)
			{
				weighting = sin(PI / 2 * (PI + gamma_max - beta) / (gamma_max + 2 * gamma));
				weighting = weighting * weighting;
			}
			else
			{
				//printf("ERROR!");
			}
			for (int i = 0; i < S; i++)
			{
				sgm[row*N + col + i * N*H] *= weighting;
			}
		}
		else
		{
			;
		}
	}

}


// weight the sinogram data of Hilbert kernel (for phase contrast imaging)
// sgm: sinogram (width x height x slice)
// N: width
// V: height (views)
// S: slice
// sdd: source to detector distance
__global__ void WeightSinogramHilbert_device(float* sgm, const float* u, const int N, const int V, const int S, float sdd)
{
	int col = threadIdx.x + blockDim.x * blockIdx.x;
	int row = threadIdx.y + blockDim.y * blockIdx.y;

	if (col < N && row < V)
	{
		for (int i = 0; i < S; i++)
		{
			sgm[row*N + col + i * N*V] *= sqrtf(u[col] * u[col] + sdd * sdd);
		}
	}
}


// weight the sinogram data of Hilbert kernel (for phase contrast imaging) along angle direction (temporary test)
// sgm: sinogram (width x height x slice)
// N: width
// V: height (views)
// S: slice
// sdd: source to detector distance
__global__ void WeightSinogramHilbert_angle_device(float* sgm, const float* u, const int N, const int V, const int S, float sdd)
{
	int col = threadIdx.x + blockDim.x * blockIdx.x;
	int row = threadIdx.y + blockDim.y * blockIdx.y;

	if (col < N && row < V)
	{
		for (int i = 0; i < S; i++)
		{
			sgm[row*N + col + i * N*V] *= sdd / sqrtf(u[col] * u[col] + sdd * sdd);
		}
	}
}

// perform beam hardening correction of sinogram
// sgm: sinogram (width x height x slice)
// N: width
// V: height (views)
// S: slice
// p0-p9: correction parameters
__global__ void CorrectBeamHardening_device(float* sgm, const int N, const int V, const int S, float p0, float p1, float p2, float p3, float p4, float p5, float p6, float p7, float p8, float p9)
{
	int col = threadIdx.x + blockDim.x * blockIdx.x;
	int row = threadIdx.y + blockDim.y * blockIdx.y;

	if (col < N && row < V)
	{
		for (int i = 0; i < S; i++)
		{
			float oldSgm = sgm[row*N + col + i * N*V];
			sgm[row*N + col + i * N*V] = p0 + p1 * powf(oldSgm, 1) + p2 * powf(oldSgm, 2) + p3 * powf(oldSgm, 3) + p4 * powf(oldSgm, 4) + p5 * powf(oldSgm, 5) + p6 * powf(oldSgm, 6) + p7 * powf(oldSgm, 7) + p8 * powf(oldSgm, 8) + p9 * powf(oldSgm, 9);
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

	if (col < N && row < V)
	{
		for (int slice = 0; slice < S; slice++)
		{
			// temporary variable to speed up
			float sgm_flt_local = 0;

			for (int i = 0; i < N; i++)
			{
				sgm_flt_local += sgm[row*N + i + slice * N*H] * reconKernel[N - 1 - col + i];
			}
			sgm_flt[row*N + col + slice * N*V] = sgm_flt_local * du;
		}

	}
}

// backproject the image using pixel-driven method
// sgm: sinogram data
// img: image data
// U: each detector element position [mm]
// u: detector pixel array
// v: detector pixel array in z axis
// beta: view angle [radius]
// N: number of detector elements
// V: number of views
// totalScanAngle [rads]
// S: number of slices of the sinogram
// coneBeam: whether the recon is conbeam recon or not
// M: image dimension
// imgS: image slice count
// sdd: source to detector distance [mm]
// sid: source to isocenter distance [mm]
// du: detector element size [mm]
// dv: detector element height [mm]
// dx: image pixel size [mm]
// dz: image slice thickness [mm]
// (xc, yc, zc): image center position [mm, mm, mm]
__global__ void BackprojectPixelDriven_device(float* sgm, float* img, float* u, float* v, float* beta, bool shortScan, const int N, const int V, \
	const int S, bool coneBeam, const int M, const int imgS, float* sdd_array, float* sid_array, float* offcenter_array, const float dx, const float dz, \
	const float xc, const float yc, const float zc, int imgS_idx)
{

	int col = threadIdx.x + blockDim.x * blockIdx.x;
	int row = threadIdx.y + blockDim.y * blockIdx.y;

	float du = u[1] - u[0];
	float dv = v[1] - v[0];


	if (col < M && row < M && imgS_idx <= imgS)
	{

		float x = (col - (M - 1) / 2.0f)*dx + xc;
		float y = ((M - 1) / 2.0f - row)*dx + yc;

		float z;

		float U, u0, v0;
		float mag_factor;
		float w;
		int k;

		float w_z;//weight for cbct
		int k_z;//index for cbct
		float delta_beta;// delta_beta for the integral calculation (nonuniform scan angle)

		float lower_row_val, upper_row_val;

		for (int slice = imgS_idx; slice < imgS_idx + 1; slice++)
		{

			z = (slice - (float(imgS) - 1.0f) / 2.0f) * dz + zc;

			// temporary local variable to speed up
			float img_local = 0;

			for (int view = 0; view < V; view++)
			{
				float offcenter_bias = offcenter_array[view] - offcenter_array[0];
				float sid = sid_array[view];
				float sdd = sdd_array[view];
				//calculation of delta_beta for the integral calculation
				if (view == 0)
					delta_beta = abs(beta[1] - beta[0]);
				else if (view == V - 1)
					delta_beta = abs(beta[view] - beta[view - 1]);
				else
					delta_beta = abs(beta[view + 1] - beta[view - 1]) / 2.0f;

				U = sid - x * cosf(beta[view]) - y * sinf(beta[view]);

				//calculate the magnification
				mag_factor = sdd / U;

				// find u0 
				u0 = mag_factor * (x*sinf(beta[view]) - y * cosf(beta[view]));


				k = floorf((u0 - (u[0] + offcenter_bias)) / du);
				if (k<0 || k + 1>N - 1)
				{
					img_local = 0;
					break;
				}

				w = (u0 - (u[k] + offcenter_bias)) / du;

				// for cone beam ct, we also need to find v0
				if (coneBeam && abs(dv) > 0.00001f)
				{
					v0 = mag_factor * z;
					// weight for cbct recon
					k_z = floorf((v0 - v[0]) / dv);
					if (k_z<0 || k_z + 1>S - 1)
					{
						img_local = 0;
						break;
					}

					w_z = (v0 - v[k_z]) / dv;

					lower_row_val = (w*sgm[view*N + k + 1 + k_z * N*V] + (1 - w)*sgm[view*N + k + k_z * N*V]);
					upper_row_val = (w*sgm[view*N + k + 1 + (k_z + 1) * N*V] + (1 - w)*sgm[view*N + k + (k_z + 1) * N*V]);

					img_local += sid / U / U * (w_z*upper_row_val + (1 - w_z)*lower_row_val) * delta_beta;
				}

				else
				{
					img_local += sid / U / U * (w*sgm[view*N + k + 1 + slice * N*V] + (1 - w)*sgm[view*N + k + slice * N*V]) * delta_beta;
				}
			}

			//judge whether the scan is a full scan or a short scan
			if (shortScan)
			{
				//printf("this is a full scan");
				img[row*M + col] = img_local;

			}
			else
				img[row*M + col] = img_local / 2.0f;


		}
	}
}

__global__ void BackprojectPixelDriven_pmatrix_device(float* sgm, float* img, float* u, float* v, float* beta, float* pmatrix, \
	bool shortScan, const int N, const int V, const int S, bool coneBeam, const int M, const int imgS, float* sdd_array, float* sid_array, \
	const float dx, const float dz, const float xc, const float yc, const float zc, int imgS_idx)
{
	int col = threadIdx.x + blockDim.x * blockIdx.x;
	int row = threadIdx.y + blockDim.y * blockIdx.y;

	if (col < M && row < M && imgS_idx < imgS)
	{

		float x_after_rotation = (col - (M - 1) / 2.0f)*dx + xc;
		float y_after_rotation = ((M - 1) / 2.0f - row)*dx + yc;

		float x = x_after_rotation * cos(beta[0]) + y_after_rotation * sin(beta[0]);//(col - (M - 1) / 2.0f)*dx + xc;
		float y = y_after_rotation * cos(beta[0]) - x_after_rotation * sin(beta[0]);//((M - 1) / 2.0f - row)*dx + yc;
		float z;
		float U;
		float w;
		int k;
		float w_z;//weight for cbct
		int k_z;//index for cbct
		float delta_beta;// delta_beta for the integral calculation (nonuniform scan angle)

		float lower_row_val, upper_row_val;

		for (int slice = imgS_idx; slice < imgS_idx + 1; slice++)
		{

			z = (slice - (float(imgS) - 1.0f) / 2.0f) * dz + zc;

			// temporary local variable to speed up
			float img_local = 0;

			for (int view = 0; view < V; view++)
			{
				float sid = sid_array[view];
				float sdd = sdd_array[view];
				//calculation of delta_beta for the integral calculation
				if (view == 0)
					delta_beta = abs(beta[1] - beta[0]);
				else if (view == V - 1)
					delta_beta = abs(beta[view] - beta[view - 1]);
				else
					delta_beta = abs(beta[view + 1] - beta[view - 1]) / 2.0f;



				//use pmatrix to calculate the corresponding index on the detector
				int pos_in_matrix = 12 * view;
				float k_u_divide_mag = pmatrix[pos_in_matrix] * x + pmatrix[pos_in_matrix + 1] * y + pmatrix[pos_in_matrix + 2] * z + pmatrix[pos_in_matrix + 3] * 1;
				float one_divide_mag = pmatrix[pos_in_matrix + 8] * x + pmatrix[pos_in_matrix + 9] * y + pmatrix[pos_in_matrix + 10] * z + pmatrix[pos_in_matrix + 11] * 1;

				float k_f = k_u_divide_mag / one_divide_mag;//float number of k
				k = floorf(k_f);

				//the pmatrix is acquired with beta[0]=0
				//however, in a real recon, the image need to be rotated
				//we need to retrieve the beta value for the pmatrix recon
				//for calculation of U
				float beta_pmatrix = beta[view] - beta[0];
				U = sid - x * cosf(beta_pmatrix) - y * sinf(beta_pmatrix);

				if (k<0 || k + 1>N - 1)
				{
					img_local = 0;
					break;
				}

				w = k_f - k;

				// for cone beam ct, we also need to find v0
				if (coneBeam)
				{
					float k_z_divide_mag = pmatrix[pos_in_matrix + 4] * x + pmatrix[pos_in_matrix + 5] * y + pmatrix[pos_in_matrix + 6] * z + pmatrix[pos_in_matrix + 7] * 1;
					float k_z_f = k_z_divide_mag / one_divide_mag;//float number of k_z
					k_z = floorf(k_z_f);
					if (k_z<0 || k_z + 1>S - 1)
					{
						img_local = 0;
						break;
					}

					w_z = k_z_f - k_z;

					lower_row_val = (w*sgm[view*N + k + 1 + k_z * N*V] + (1 - w)*sgm[view*N + k + k_z * N*V]);
					upper_row_val = (w*sgm[view*N + k + 1 + (k_z + 1) * N*V] + (1 - w)*sgm[view*N + k + (k_z + 1) * N*V]);

					img_local += sid / U / U * (w_z*upper_row_val + (1 - w_z)*lower_row_val) * delta_beta;
				}

				else
				{
					img_local += sid / U / U * (w*sgm[view*N + k + 1 + slice * N*V] + (1 - w)*sgm[view*N + k + slice * N*V]) * delta_beta;
				}
			}

			//judge whether the scan is a full scan or a short scan
			if (shortScan)
			{
				//printf("this is a full scan");
				img[row*M + col] = img_local;

			}
			else
			{
				img[row*M + col] = img_local / 2.0f;
			}


		}
	}
}


// backproject the image using pixel-driven method for Hilbert kernel (for phase contrast imaging)
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
__global__ void BackprojectPixelDrivenHilbert_device(float* sgm, float* img, float* u, float* beta, const int N, const int V, \
	const int S, const int M, const float sdd, const float sid, const float du, const float dx, const float xc, const float yc, int imgS_idx)
{
	int col = threadIdx.x + blockDim.x * blockIdx.x;
	int row = threadIdx.y + blockDim.y * blockIdx.y;

	if (col < M && row < M && imgS_idx < S)
	{
		float x = (col - (M - 1) / 2.0f)*dx + xc;
		float y = ((M - 1) / 2.0f - row)*dx + yc;

		float U, u0;
		float w;
		int k;

		for (int slice = imgS_idx; slice < imgS_idx + 1; slice++)
		{
			img[row*M + col] = 0;

			for (int view = 0; view < V; view++)
			{
				U = sid - x * cosf(beta[view]) - y * sinf(beta[view]);
				u0 = sdd * (x*sinf(beta[view]) - y * cosf(beta[view])) / U;

				k = floorf((u0 - u[0]) / du);
				if (k<0 || k + 1>N - 1)
				{
					img[row*M + col] = 0;
					break;
				}

				w = (u0 - u[k]) / du;

				img[row*M + col] += 1 / U * (w*sgm[view*N + k + 1 + slice * N*V] + (1 - w)*sgm[view*N + k + slice * N*V]);

			}
			img[row*M + col] *= PI / V;
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
		printf("\nCannot find SDD information file '%s'!\n", distanceFile.c_str());
		exit(-2);
	}
	rapidjson::IStreamWrapper isw(ifs);
	rapidjson::Document doc;
	doc.ParseStream<js::kParseCommentsFlag | js::kParseTrailingCommasFlag>(isw);
	js::Value distance_jsonc_value;
	if (doc.HasMember("SourceDetectorDistance"))
	{

		distance_jsonc_value = doc["SourceDetectorDistance"];

	} 
	else if(doc.HasMember("Value"))//a new version of the program uses value as member to avoid complex member names
	{
		distance_jsonc_value = doc["Value"];
	}
	else
	{
		printf("\nDid not find Value member in jsonc file!\n");
		exit(-2);
	}

	if (distance_jsonc_value.Size() != V)
	{
		printf("\nNumber of sdd values is %d while the number of Views is %d!\n", distance_jsonc_value.Size(), V);
		exit(-2);
	}

	for (unsigned i = 0; i < distance_jsonc_value.Size(); i++)
	{
		distance_array_cpu[i] = distance_jsonc_value[i].GetFloat();
	}

	cudaMemcpy(distance_array, distance_array_cpu, V * sizeof(float), cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
}

void InitializeNonuniformSID_Agent(float* &distance_array, const int V, const std::string& distanceFile)
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
		printf("\nCannot find SID information file '%s'!\n", distanceFile.c_str());
		exit(-2);
	}
	rapidjson::IStreamWrapper isw(ifs);
	rapidjson::Document doc;
	doc.ParseStream<js::kParseCommentsFlag | js::kParseTrailingCommasFlag>(isw);
	js::Value distance_jsonc_value;
	if (doc.HasMember("SourceIsocenterDistance"))
	{
		distance_jsonc_value = doc["SourceIsocenterDistance"];
	}
	else if (doc.HasMember("Value"))//a new version of the program uses value as member to avoid complex member names
	{
		distance_jsonc_value = doc["Value"];
	}
	else
	{
		printf("\nDid not find SourceIsocenterDistance or Value member in jsonc file!\n");
		exit(-2);
	}

	if (distance_jsonc_value.Size() != V)
	{
		printf("\nNumber of sid values is %d while the number of Views is %d!\n", distance_jsonc_value.Size(), V);
		exit(-2);
	}

	for (unsigned i = 0; i < distance_jsonc_value.Size(); i++)
	{
		distance_array_cpu[i] = distance_jsonc_value[i].GetFloat();
	}

	cudaMemcpy(distance_array, distance_array_cpu, V * sizeof(float), cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
}

void InitializeNonuniformOffCenter_Agent(float* &distance_array, const int V, const std::string& distanceFile)
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
		printf("\nCannot find Offcenter information file '%s'!\n", distanceFile.c_str());
		exit(-2);
	}
	rapidjson::IStreamWrapper isw(ifs);
	rapidjson::Document doc;
	doc.ParseStream<js::kParseCommentsFlag | js::kParseTrailingCommasFlag>(isw);
	js::Value distance_jsonc_value;
	if (doc.HasMember("OffcenterArray"))
	{
		distance_jsonc_value = doc["OffcenterArray"];
	}
	else if (doc.HasMember("Value"))//a new version of the program uses value as member to avoid complex member names
	{
		distance_jsonc_value = doc["Value"];
	}
	else
	{
		printf("\nDid not find OffcenterArray or Value member in jsonc file!\n");
		exit(-2);
	}

	if (distance_jsonc_value.Size() != V)
	{
		printf("\nNumber of offcenter values is %d while the number of Views is %d!\n", distance_jsonc_value.Size(), V);
		exit(-2);
	}

	for (unsigned i = 0; i < distance_jsonc_value.Size(); i++)
	{
		distance_array_cpu[i] = distance_jsonc_value[i].GetFloat();
	}

	cudaMemcpy(distance_array, distance_array_cpu, V * sizeof(float), cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
}

void InitializePMatrix_Agent(float* &pmatrix_array, const int V, const std::string& pmatrixFile)
{
	namespace fs = std::experimental::filesystem;
	namespace js = rapidjson;

	if (pmatrix_array != nullptr)
		cudaFree(pmatrix_array);

	//cudaMallocManaged((void**)&pmatrix_array, 12 * V * sizeof(float));
	cudaMalloc((void**)&pmatrix_array, 12 * V * sizeof(float));
	//cudaMallocManaged somehow does not work for this function
	//so cudaMalloc and cudaMemcpy is used


	float* pmatrix_array_cpu = new float[12 * V];


	std::ifstream ifs(pmatrixFile);
	if (!ifs)
	{
		printf("\nCannot find pmatrix information file '%s'!\n", pmatrixFile.c_str());
		exit(-2);
	}
	rapidjson::IStreamWrapper isw(ifs);
	rapidjson::Document doc;
	doc.ParseStream<js::kParseCommentsFlag | js::kParseTrailingCommasFlag>(isw);
	js::Value pmatrix_jsonc_value;
	if (doc.HasMember("PMatrix"))
	{
		pmatrix_jsonc_value = doc["PMatrix"];
	}
	else if(doc.HasMember("Value"))
	{
		pmatrix_jsonc_value = doc["Value"];
	}
	else
	{
		printf("\nDid not find PMatrix or Value member in jsonc file!\n");
		exit(-2);
	}
	if (pmatrix_jsonc_value.Size() != 12 * V)
	{
		printf("\nNumber of pmatrix elements is %d while the 12 times number of Views is %d!\n", pmatrix_jsonc_value.Size(), 12 * V);
		exit(-2);
	}

	for (unsigned i = 0; i < 12 * V; i++)
	{
		//printf("\n%d: %f",i, pmatrix_jsonc_value[i].GetFloat());
		pmatrix_array_cpu[i] = pmatrix_jsonc_value[i].GetFloat();
	}
	cudaMemcpy(pmatrix_array, pmatrix_array_cpu, 12 * V * sizeof(float), cudaMemcpyHostToDevice);

	cudaDeviceSynchronize();
}


void InitializeU_Agent(float* &u, const int N, const float du, const float offcenter)
{
	if (u != nullptr)
		cudaFree(u);

	cudaMalloc((void**)&u, N * sizeof(float));
	InitU << <(N + 511) / 512, 512 >> > (u, N, du, offcenter);
}

void InitializeBeta_Agent(float* &beta, const int V, const float rotation, const float totalScanAngle)
{
	if (beta != nullptr)
		cudaFree(beta);

	cudaMalloc((void**)&beta, V * sizeof(float));
	InitBeta << < (V + 511) / 512, 512 >> > (beta, V, rotation, totalScanAngle);
}

void InitializeNonuniformBeta_Agent(float* &beta, const int V, const float rotation, const std::string& scanAngleFile)
//unit of beta is RADs
{
	namespace fs = std::experimental::filesystem;
	namespace js = rapidjson;

	if (beta != nullptr)
		cudaFree(beta);

	cudaMalloc((void**)&beta, V * sizeof(float));
	float* beta_cpu = new float[V];
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
	}
	else if (doc.HasMember("Value"))
	{
		scan_angle_jsonc_value = doc["Value"];
	}
	else
	{
		printf("Did not find ScanAngle or Value member in jsonc file!\n");
		exit(-2);
	}

	if (scan_angle_jsonc_value.Size() != V)
	{
		printf("Number of scan angles is %d while the number of Views is %d!\n", scan_angle_jsonc_value.Size(), V);
		exit(-2);
	}

	for (unsigned i = 0; i < scan_angle_jsonc_value.Size(); i++)
	{
		beta_cpu[i] = rotation / 180.0f*PI + scan_angle_jsonc_value[i].GetFloat() / 180.0*PI;
	}
	cudaMemcpy(beta, beta_cpu, sizeof(float)*V, cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
}

void InitializeReconKernel_Agent(float* &reconKernel, const int N, const float du, const std::string& kernelName, const std::vector<float>& kernelParam)
{
	if (reconKernel != nullptr)
		cudaFree(reconKernel);

	cudaMalloc((void**)&reconKernel, (2 * N - 1) * sizeof(float));

	if (kernelName == "HammingFilter")
	{
		InitReconKernel_Hamming << <(2 * N - 1 + 511) / 512, 512 >> > (reconKernel, N, du, kernelParam[0]);
	}
	if (kernelName == "Delta")
	{
		InitReconKernel_Delta << <(2 * N - 1 + 511) / 512, 512 >> > (reconKernel, N, du, kernelParam[0]);
	}
	else if (kernelName == "QuadraticFilter")
	{
		float lastParam = 0.0f;
		if (kernelParam.size() == 3)
			lastParam = kernelParam[2];

		InitReconKernel_Quadratic << <(2 * N - 1 + 511) / 512, 512 >> > (reconKernel, N, du, int(kernelParam.size()), kernelParam[0], kernelParam[1], lastParam);
	}
	else if (kernelName == "Polynomial")
	{
		// TODO: 
		// InitReconKernel_Polynomial <<<...>>> (...);
		float p[7] = { 0 };

		for (size_t i = 0; i < kernelParam.size(); i++)
		{
			p[i] = kernelParam[kernelParam.size() - 1 - i];
		}

		//InitReconKernel_Polynomial <<<(2 * N - 1 + 511) / 512, 512>>> (reconKernel, N, du, p[0], p[1], p[2], p[3], p[4], p[5], p[6]);
		InitReconKernel_Polynomial << <(2 * N - 1 + 511) / 512, 512 >> > (reconKernel, N, du, p[6], p[5], p[4], p[3], p[2], p[1], p[0]);
	}
	else if (kernelName == "Hilbert" || kernelName == "Hilbert_angle")
	{
		InitReconKernel_Hilbert << <(2 * N - 1 + 511) / 512, 512 >> > (reconKernel, N, du, kernelParam[0]);
	}
	else if (kernelName == "GaussianApodizedRamp")
	{
		InitReconKernel_GaussianApodized << <(2 * N - 1 + 511) / 512, 512 >> > (reconKernel, N, du, kernelParam[0]);
	}
}

void MallocManaged_Agent(float * &p, const int size)
{
	cudaMallocManaged((void**)&p, size);
}


void CorrectBeamHardening_Agent(float* sgm, mango::Config & config)
{
	dim3 grid((config.sgmWidth + 15) / 16, (config.sgmHeight + 15) / 16);
	dim3 block(16, 16);

	CorrectBeamHardening_device << <grid, block >> > (sgm, config.sgmWidth, config.sgmHeight, config.sliceCount, config.beamHardening[0], config.beamHardening[1], config.beamHardening[2], config.beamHardening[3], config.beamHardening[4], config.beamHardening[5], config.beamHardening[6], config.beamHardening[7], config.beamHardening[8], config.beamHardening[9]);

	cudaDeviceSynchronize();

}

void FilterSinogram_Agent(float * sgm, float* sgm_flt, float* reconKernel, float* u, mango::Config & config, float* beta, float * sdd_array, float * offcenter_array)
{
	// Step 1: weight the sinogram
	dim3 grid((config.sgmWidth + 15) / 16, (config.sgmHeight + 15) / 16);
	dim3 block(16, 16);

	// Hilbert kernel for phase contrast imaging
	if (config.kernelName == "Hilbert")
		WeightSinogramHilbert_device << <grid, block >> > (sgm, u, config.sgmWidth, config.sgmHeight, config.sliceCount, config.sdd);
	else if (config.kernelName == "Hilbert_angle")
	{
		printf("Kernel name: %s\n", config.kernelName);
		WeightSinogramHilbert_angle_device << <grid, block >> > (sgm, u, config.sgmWidth, config.sgmHeight, config.sliceCount, config.sdd);
	}
	// Common attenuation imaging
	else
		WeightSinogram_device << <grid, block >> > (sgm, u, config.sgmWidth, config.sgmHeight, config.views, config.sliceCount, \
			config.sliceThickness, config.sliceOffcenter, sdd_array, config.totalScanAngle, config.shortScan, beta, offcenter_array);

	cudaDeviceSynchronize();

	// Step 2: convolve the sinogram
	if (config.kernelName == "GaussianApodizedRamp")
	{
		// if Guassian aposied kernel is used, the sinogram need to be filtered twice
		// first by the ramp filter, then by the gaussian filter
		float du = config.detEltSize;
		float * reconKernel_ramp;
		cudaMalloc((void**)&reconKernel_ramp, (2 * config.sgmWidth - 1) * sizeof(float));
		InitReconKernel_Hamming << <(2 * config.sgmWidth - 1 + 511) / 512, 512 >> > (reconKernel_ramp, config.sgmWidth, du, 1);

		cudaDeviceSynchronize();

		//intermidiate filtration result is saved in sgm_flt_ramp
		float *sgm_flt_ramp;
		//MallocManaged_Agent(sgm_flt_ramp, config.sgmWidth*config.views*config.sliceCount * sizeof(float));
		cudaMalloc((void**)& sgm_flt_ramp, config.sgmWidth * config.views * config.sliceCount * sizeof(float));

		ConvolveSinogram_device << <grid, block >> > (sgm_flt_ramp, sgm, reconKernel_ramp, config.sgmWidth, config.sgmHeight, config.views, config.sliceCount, u, config.detEltSize);
		cudaDeviceSynchronize();
		//the height of the filtered sinogram shrinks to number of views, so the convolution parameters need to be adjusted accordingly
		ConvolveSinogram_device << <grid, block >> > (sgm_flt, sgm_flt_ramp, reconKernel, config.sgmWidth, config.views, config.views, config.sliceCount, u, config.detEltSize);
		cudaDeviceSynchronize();

		// free temporary memory
		cudaFree(reconKernel_ramp);
		cudaFree(sgm_flt_ramp);
	}
	else
	{
		ConvolveSinogram_device << <grid, block >> > (sgm_flt, sgm, reconKernel, config.sgmWidth, config.sgmHeight, config.views, config.sliceCount, u, config.detEltSize);
		cudaDeviceSynchronize();
	}
}

void BackprojectPixelDriven_Agent(float * sgm_flt, float * img, float* sdd_array, float* sid_array, float* offcenter_array, float* pmatrix_array, float * u, float *v, float* beta, mango::Config & config, int z_idx)
{
	dim3 grid((config.imgDim + 15) / 16, (config.imgDim + 15) / 16);
	dim3 block(16, 16);

	// Hilbert kernel for phase contrast imaging
	if (config.kernelName == "Hilbert" || config.kernelName == "Hilbert_angle")
	{
		BackprojectPixelDrivenHilbert_device << <grid, block >> > (sgm_flt, img, u, beta, config.sgmWidth, config.views, \
			config.sliceCount, config.imgDim, config.sdd, config.sid, config.detEltSize, config.pixelSize, config.xCenter, config.yCenter, z_idx);
	}
	// Common attenuation imaging
	else if (config.pmatrixFlag == false)// if pmatrix is not applied
	{
		BackprojectPixelDriven_device << <grid, block >> > (sgm_flt, img, u, v, beta, config.shortScan, config.sgmWidth, config.views, \
			config.sliceCount, config.coneBeam, config.imgDim, config.imgSliceCount, sdd_array, sid_array, offcenter_array, config.pixelSize, config.imgSliceThickness, \
			config.xCenter, config.yCenter, config.zCenter, z_idx);
	}
	else if (config.pmatrixFlag == true)// if pmatrix is applied
	{
		BackprojectPixelDriven_pmatrix_device << <grid, block >> > (sgm_flt, img, u, v, beta, pmatrix_array, config.shortScan, config.sgmWidth, config.views, \
			config.sliceCount, config.coneBeam, config.imgDim, config.imgSliceCount, sdd_array, sid_array, config.pixelSize, config.imgSliceThickness, \
			config.xCenter, config.yCenter, config.zCenter, z_idx);
	}

	cudaDeviceSynchronize();
}

void SaveReconImageSlice(const char* filename, float* rec_image, int z_idx, const mango::Config& config)
{
	FILE* fp = NULL;
	if (z_idx == 0)
		fp = fopen(filename, "wb");
	else
		fp = fopen(filename, "ab");

	if (fp == NULL)
	{
		fprintf(stderr, "Cannot save to file %s!\n", filename);
		exit(4);
	}
	fwrite(rec_image, sizeof(float), config.imgDim*config.imgDim, fp);
	fclose(fp);
}


void FreeMemory_Agent(float* &p)
{
	cudaFree(p);
	p = nullptr;
}