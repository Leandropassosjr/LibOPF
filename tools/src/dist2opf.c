#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
	if  (argc != 3)
	{
		fprintf(stderr, "\nusage dist2opf <P1> <P2>\n");
		fprintf(stderr, "P1: input file name in ASCII format\n");
		fprintf(stderr, "P2: output file name in the OPF binary format\n");
		exit(-1);
	}

	fprintf(stderr, "Program to convert a distance file into OPF format\n");

	FILE *fpIn = NULL, *fpOut = NULL;
	int n_samples;
	float distance;

	fpIn  = fopen(argv[1], "r");
	fpOut = fopen(argv[2], "wb");

	/*reading / writing amount of samples on the file*/
	if (fscanf(fpIn, "%d", &n_samples) != 1)
	{
		fprintf(stderr, "Could not read number of samples\n");
		exit(-1);
	}
	fprintf(stderr, "Amount of samples: %d\n", n_samples);
	fwrite(&n_samples, sizeof(int), 1, fpOut);

	/*transfering distances from txt to opf file*/
	for (int i = 0; i < n_samples; i++) {

		for (int j = 0; j < n_samples; j++) {
			if (fscanf(fpIn, "%f ", &distance) != 1) {
				fprintf(stderr, "Could not read distance between sample %d and %d\n", i, j);
				exit(-1);
			}

			fwrite(&distance, sizeof(float), 1, fpOut);
		}
	}

	fclose(fpIn);
	fclose(fpOut);

	return 0;
}