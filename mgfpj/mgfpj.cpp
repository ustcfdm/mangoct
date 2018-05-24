#include "stdafx.h"
#include "FpjClass.cuh"

int main(int argc, char* argv[])
{
	namespace mg = mango;
	namespace fs = std::experimental::filesystem;

	if (argc == 1)
	{
		printf( "This is mangoct forward projection tool, please input config file name as arguments.\n");
		return 0;
	}

	for (int idx = 1; idx < argc; idx++)
	{
		mg::FpjClass fpj;

		printf("Loading config %s...\n", argv[idx]);

		fpj.ReadConfigFile(argv[idx]);
	}


	getchar();
	return 0;
}