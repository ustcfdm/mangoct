#include "../Matrix/Matrix.h"

using namespace mango;

int main()
{

	Matrix m = Matrix::ReadEviFile("C:\\Users\\mfeng\\Local Files\\For Evan\\Beamhardening water correction\\20170423watercorrection\\EVIdata\\1_TE.EVI", 64u, 5120u, 2720u, 2368u, 64u);

	printf("rows = %d, cols = %d, pages = %d\n", m.Rows(), m.Cols(), m.Pages());


	m.SaveRawFile("C:\\Users\\mfeng\\Local Files\\For Evan\\Beamhardening water correction\\20170423watercorrection\\EVIdata\\1_TE.raw");

	//float one = 1.0f;
	//float zero = 0.0f;
	//float inf = one/zero;

	//float nan_ = zero / zero;

	//if (one<inf)
	//{
	//	printf("one < inf\n");
	//}
	//else
	//{
	//	printf("one >= inf\n");
	//}

	//if (one>nan_)
	//{
	//	printf("one > nan\n");
	//}
	//else if (one<nan_)
	//{
	//	printf("one < nan\n");
	//}
	//else if (one== nan_)
	//{
	//	printf("one == nan\n");
	//}

	//if (isnan(nan_))
	//{
	//	printf("nan is NaN\n");
	//}

	//float log_Inf = log(inf);
	//float log_Nan = log(nan_);
	//printf("log_inf = %f\nlog_nan = %f\n", log_Inf, log_Nan);
	//auto a = nan(ind);


	//mango::Matrix m1(5, 6, 7);
	//float val = 0;
	//for (unsigned page = 0; page < m1.Pages(); page++)
	//{
	//	for (unsigned row = 0; row < m1.Rows(); row++)
	//	{
	//		for (unsigned col = 0; col < m1.Cols(); col++)
	//		{
	//			m1(row, col, page) = val++;
	//		}
	//	}
	//}
	//mango::Matrix m2(m1);
	//mango::Matrix m3(m1);

	//mango::Matrix m11 = mango::Matrix::Average(m1, mango::Axis::Row, 0u, 2u);
	//mango::Matrix m22 = mango::Matrix::Average(m1, mango::Axis::Col, 1u, 4u);
	//mango::Matrix m33 = mango::Matrix::Sum(m1, mango::Axis::Col, 1u, 5u);


	//m1.Average(mango::Axis::Row, 0u, 2u);
	//m2.Average(mango::Axis::Col, 1u, 4u);
	//m3.Average(mango::Axis::Page, 1u, 5u);

	//for (unsigned page = 0; page < m1.Pages(); page++)
	//{
	//	for (unsigned row = 0; row < m1.Rows(); row++)
	//	{
	//		for (unsigned col = 0; col < m1.Cols(); col++)
	//		{
	//			printf("%f\t", m1(row, col, page));
	//		}
	//		printf("\n");
	//	}
	//}

	//m1.SaveRawFile("testfile.raw");

	//Matrix m1(3, 4);
	//for (unsigned row = 0; row < m1.Rows(); row++)
	//{
	//	for (unsigned col = 0; col < m1.Cols(); col++)
	//	{
	//		m1(row, col) = 1;
	//	}
	//}

	//m1.Sum(Axis::Col);

	//for (unsigned row = 0; row < m1.Rows(); row++)
	//{
	//	for (unsigned col = 0; col < m1.Cols(); col++)
	//	{
	//		printf("%f\t", m1(row, col));
	//	}
	//	printf("\n");
	//}



	printf("Press Enter to continue.");
	getchar();
	return 0;
}