#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
	if  (argc != 3)
	{
		fprintf(stderr, "\nusage opf2dist <P1> <P2>\n");
		fprintf(stderr, "P1: output filne name in the OPF binary format\n");
		fprintf(stderr, "P2: input file name in ASCII format\n");
		exit(-1);
	}

	fprintf(stderr, "Program to convert a distance file from OPF to txt format\n");

	FILE *fpIn = NULL, *fpOut = NULL;
	int n_samples;
	float distance;

	fpIn  = fopen(argv[1], "rb");
	fpOut = fopen(argv[2], "w");

	/*reading / writing amount of samples on the file*/
	if (fread(&n_samples, sizeof(int), 1, fpIn) != 1)
	{
		fprintf(stderr, "Could not read number of samples\n");
		exit(-1);
	}

	fprintf(stderr, "Amount of samples: %d\n", n_samples);
	fprintf(fpOut, "%d\n", n_samples);

	/*transfering distances from txt to opf file*/
	for (int i = 0; i < n_samples; i++) {

		for (int j = 0; j < n_samples; j++) {
			if (fread(&distance, sizeof(float), 1, fpIn) != 1) {
				fprintf(stderr, "Could not read distance between sample %d and %d\n", i, j);
				exit(-1);
			}

			fprintf(fpOut, "%f ", distance);
		}
		fprintf(fpOut, "\n");
	}

	fclose(fpIn);
	fclose(fpOut);

	return 0;
}