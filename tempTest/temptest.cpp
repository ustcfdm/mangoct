#include "../Matrix/Matrix.h"

using namespace mango;

int main()
{
	Matrix m1(3, 4);
	for (unsigned row = 0; row < m1.Rows(); row++)
	{
		for (unsigned col = 0; col < m1.Cols(); col++)
		{
			m1(row, col) = 1;
		}
	}

	m1.Sum(Axis::Col);

	for (unsigned row = 0; row < m1.Rows(); row++)
	{
		for (unsigned col = 0; col < m1.Cols(); col++)
		{
			printf("%f\t", m1(row, col));
		}
		printf("\n");
	}



	printf("Press Enter to continue.");
	getchar();
	return 0;
}