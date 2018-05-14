#include "OPF.h"
#define _MAJORITY_VOTING

int main(int argc, char **argv)
{
	int i, n, j, op;
    int kmin, kmax;
    int amount_supervised_labels = 0;
    int amount_unsupervised_labels = 0;
	float value;
	char fileName[256];
	FILE *f = NULL;

	fprintf(stdout, "\nProgram that computes clusters by OPF\n");
	fprintf(stdout, "\nIf you have any problem, please contact: ");
	fprintf(stdout, "\n- alexandre.falcao@gmail.com");
	fprintf(stdout, "\n- papa.joaopaulo@gmail.com\n");
	fprintf(stdout, "\nLibOPF version 2.0 (2009)\n");
	fprintf(stdout, "\n");

	if ((argc != 6) && (argc != 7))
	{
		fprintf(stderr, "\nusage opf_cluster <P1> <P2> <P3> <P4> <P5> <P6>");
		fprintf(stderr, "\nP1: unlabeled data set in the OPF file format");
		fprintf(stderr, "\nP2: kmin(minimum degree for the knn graph)");
		fprintf(stderr, "\nP3: kmax(maximum degree for the knn graph)");
		fprintf(stderr, "\nP4: P3 0 (height), 1(area) and 2(volume)");
		fprintf(stderr, "\nP5: value of parameter P3 in (0-1)");
		fprintf(stderr, "\nP6: precomputed distance file (leave it in blank if you are not using this resource");
		exit(-1);
	}

	if (argc == 7)
		opf_PrecomputedDistance = 1;

	fprintf(stdout, "\nReading data file ...");
	Subgraph *g = ReadSubgraph(argv[1]);

	// this field is going to be overwritten during the computation of the optimal clusters
	amount_supervised_labels = g->nlabels;

	if (opf_PrecomputedDistance)
	{
		opf_DistanceValue = opf_ReadDistances(argv[6], &n);
	}

	op = atoi(argv[4]);

	kmin = atoi(argv[2]);
	kmax = atoi(argv[3]);

	if (kmin > kmax || kmin < 1)
	{
		fprintf(stderr, "Kmin must be smaller or equal than Kmax and Kmin must be larger than 0.\n");
		exit(-1);	
	}

	opf_BestkMinCut(g, kmin, kmax); //default kmin = 1

	value = atof(argv[5]);
	if ((value < 1) && (value > 0))
	{
		fprintf(stdout, "\n\n Filtering clusters ... ");
		switch (op)
		{
		case 0:
			fprintf(stdout, "\n by dome height ... ");
			float Hmax = 0;
			for (i = 0; i < g->nnodes; i++)
				if (g->node[i].dens > Hmax)
					Hmax = g->node[i].dens;
			opf_ElimMaxBelowH(g, value * Hmax);
			break;
		case 1:
			fprintf(stdout, "\n by area ... ");
			opf_ElimMaxBelowArea(g, (int)(value * g->nnodes));
			break;
		case 2:
			fprintf(stdout, "\n by volume ... ");
			double Vmax = 0;
			for (i = 0; i < g->nnodes; i++)
				Vmax += g->node[i].dens;
			opf_ElimMaxBelowVolume(g, (int)(value * Vmax / g->nnodes));
			break;
		default:
			fprintf(stderr, "\nInvalid option for parameter P3 ... ");
			exit(-1);
			break;
		}
	}

	fprintf(stdout, "\n\nClustering by OPF ");
	opf_OPFClustering(g);
	printf("num of clusters %d\n", g->nlabels);

    // storing data to-be overwritten
    amount_unsupervised_labels = g->nlabels;

	/* If the training set has true labels, then create a
	   classifier by propagating the true label of each root to
	   the nodes of its tree (cluster). This classifier can be
	   evaluated by running opf_knn_classify on the training set
	   or on unseen testing set. Otherwise, copy the cluster
	   labels to the true label of the training set and write a
	   classifier, which essentially can propagate the cluster
	   labels to new nodes in a testing set. */

    if (amount_supervised_labels == 0) {  // unlabeled training set
	  for (i = 0; i < g->nnodes; i++)
		g->node[i].truelabel = g->node[i].label + 1;
    }
    else {  // labeled training set

#ifdef _MAJORITY_VOTING
		fprintf(stdout, "\nLabeling dataset cluster majority voting.");
        for (i = 0; i < g->nnodes; i++) {
            if (g->node[i].children || g->node[i].root == i) {
				// fprintf(stdout, "\nParent: %d (label: %d)\n\t", i, g->node[i].truelabel);

				// adding 1 because class labels are 1-indexed, need to make room for the last element
				int *frequency = (int*)calloc((amount_supervised_labels + 1), sizeof(int));
				Set *walker = g->node[i].children;

                // taking into account the prototype label when computing the histogram
                frequency[g->node[i].truelabel] += 1;

				// computing label frequency on the neighborhood of the critical point
				while (walker) {
					j = walker->elem;

					frequency[g->node[j].truelabel] += 1;
					walker = walker->next;
				}

				// determining the label mode
				int most_frequent = 0;
				int top_index = 0;

				for (j = 0; j < amount_supervised_labels + 1; j++) {
					if (frequency[j] > most_frequent) {
						most_frequent = frequency[j];
						top_index = j;
					}
				}

				// propagating the most common label to all samples connected to the prototype
				walker = g->node[i].children;
				while (walker) {
					j = walker->elem;
					g->node[j].label = top_index;
					walker = walker->next;
				}

				// the prototype also receives the most common label
				g->node[i].label = top_index;

				// debugging
                /*
				for (j = 0; j < amount_supervised_labels + 1; j++)
					fprintf(stdout, "%d ", frequency[j]);

				fprintf(stdout, "\tMost frequent: %d (%d)", top_index, most_frequent);
				fprintf(stdout, "\n");
                */
            }
        }
#else
		fprintf(stdout, "\nLabeling dataset using the prototype true label.");
	    g->nlabels = 0;
	    for (i = 0; i < g->nnodes; i++){//propagating root labels
	      if (g->node[i].root==i)
	        g->node[i].label = g->node[i].truelabel;
	      else
	        g->node[i].label = g->node[g->node[i].root].truelabel;
	    }

	    for (i = 0; i < g->nnodes; i++){
		  // retrieve the original number of true labels
		  if (g->node[i].label > g->nlabels)
			  g->nlabels = g->node[i].label;
	    }
#endif
    }

	fprintf(stdout, "\nWriting classifier's model file ...");
	fflush(stdout);
	opf_WriteModelFile(g, "classifier.opf");
	fprintf(stdout, " OK");
	fflush(stdout);

	fprintf(stdout, "\nWriting output file ...");
	fflush(stdout);
	sprintf(fileName, "%s.out", argv[1]);
	f = fopen(fileName, "w");
	for (i = 0; i < g->nnodes; i++)
		fprintf(f, "%d\n", g->node[i].label);
	fclose(f);
	fprintf(stdout, " OK");
	fflush(stdout);

    fprintf(stderr,  "\nWriting report.");
    f = fopen("report.txt", "w");
    fprintf(f, "%d %d\n", amount_unsupervised_labels, g->bestk);
    fclose(f);
    fprintf(stderr, "\nDone.");

	fprintf(stdout, "\n\nDeallocating memory ...\n");
	DestroySubgraph(&g);
	if (opf_PrecomputedDistance)
	{
		for (i = 0; i < n; i++)
			free(opf_DistanceValue[i]);
		free(opf_DistanceValue);
	}

	return 0;
}
