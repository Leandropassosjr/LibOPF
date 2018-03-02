#include "OPF.h"
#include <stdio.h>

// Checks if the train/eval/test split uses all samples from the dataset
void CheckInputData(int total, int trainAmount, int evalAmount, int testAmount)
{
	fprintf(stderr, "\nChecking if the amount of samples match the total.");
	fprintf(stderr, "\nTotal amount of samples: %d\n", total);
	if (trainAmount + evalAmount + testAmount != total)
		Error("Amount of samples for train, eval and test do not match the "
				"total amount of samples", "CheckInputData");

	if (trainAmount == 0 || testAmount == 0)
		Error("Amount of train and test samples cannot be 0",
				"CheckInputData");

	fprintf(stderr, "\nOK");
}

int main(int argc, char **argv)
{
	fflush(stdout);
	fprintf(stdout, "\nProgram that generates training, evaluation and test "
			"sets for the OPF classifier\n");
	fprintf(stdout, "\nIf you have any problem, please contact: ");
	fprintf(stdout, "\n- alexandre.falcao@gmail.com");
	fprintf(stdout, "\n- papa.joaopaulo@gmail.com\n");
	fprintf(stdout, "\nLibOPF version 2.0 (2018)\n");
	fprintf(stdout, "\n");
	fflush(stdout);

	if (argc != 6)
	{
		fprintf(stderr, "\nUsage opf_split_hard <P1> <P2> <P3> <P4> <P5>");
		fprintf(stderr, "\nP1: input dataset in the OPF file format");
		fprintf(stderr, "\nP2: raw count of samples for the training set");
		fprintf(stderr, "\nP3: raw count of samples for the evaluation set. "
		    "(leave 0 if the case of no learning)");
		fprintf(stderr, "\nP4: raw count of samples for the test set");
		fprintf(stderr, "\nP5: normalize features? 1 - Yes  0 - No\n\n");
		exit(-1);
	}

	Subgraph *gOriginal = NULL;
	Subgraph *gAux = NULL;
	Subgraph *gTraining = NULL;
	Subgraph *gEvaluating = NULL;
	Subgraph *gTesting = NULL;
	float training_p = atof(argv[2]);
	float evaluating_p = atof(argv[3]);
	float testing_p = atof(argv[4]);
	int normalize = atoi(argv[5]);

	fprintf(stdout, "\nReading data set ...");
	fflush(stdout);
	gOriginal = ReadSubgraph(argv[1]);
	fprintf(stdout, " OK");
	fflush(stdout);

	CheckInputData(gOriginal->nnodes, training_p, evaluating_p, testing_p);

	if (normalize)
		opf_NormalizeFeatures(gOriginal);

	// NOTICE: This part of the code works slightly different from opf_split.
	fprintf(stdout, "\nSplitting data set ...");
	fflush(stdout);
	opf_SplitSubgraphHard(gOriginal, &gTraining, &gAux, training_p);

	fprintf(stderr, "\nAmount of train samples: %d", gTraining->nnodes);

	if (evaluating_p > 0) {
		opf_SplitSubgraphHard(gAux, &gEvaluating, &gTesting, evaluating_p);
	  fprintf(stderr, "\nAmount of evaluation samples: %d", gEvaluating->nnodes);
	}
	else
		gTesting = CopySubgraph(gAux);

	fprintf(stderr, "\nAmount of test samples: %d", gTesting->nnodes);
	fprintf(stdout, " OK");
	fflush(stdout);

	fprintf(stdout, "\nWriting data sets to disk ...");
	fflush(stdout);

	WriteSubgraph(gTraining, "training.dat");
	WriteSubgraph(gTesting, "testing.dat");
	if (evaluating_p > 0)
		WriteSubgraph(gEvaluating, "evaluating.dat");

	fprintf(stdout, " OK");
	fflush(stdout);

	fprintf(stdout, "\nDeallocating memory ...");
	DestroySubgraph(&gOriginal);
	DestroySubgraph(&gAux);
	DestroySubgraph(&gTraining);
	DestroySubgraph(&gEvaluating);
	DestroySubgraph(&gTesting);
	fprintf(stdout, " OK\n");

	return 0;
}
