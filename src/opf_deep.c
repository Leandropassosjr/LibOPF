#include "OPF.h"

#define ARG_SAMPLES_FILE     1
#define ARG_KMIN_BOTTOM      2
#define ARG_KMAX_BOTTOM      3
#define ARG_KMIN_TOP         4
#define ARG_KMAX_TOP         5
#define ARG_AMOUNT_LEVELS    6
#define ARG_MIN_PROTOTYPES   7
#define ARG_FILTERING_OPTION 8
#define ARG_FILTERING_VALUE  9
#define ARG_DISTANCES        10

void validate(int argc, char **argv) {
    if (argc != (ARG_DISTANCES - 1) && argc != ARG_DISTANCES) {
        fprintf(stderr, "\nusage opf_cluster <P1> <P2> <P3> <P4> <P5> <P6>");
		fprintf(stderr, "\nP1: unlabeled data set in the OPF file format");
        fprintf(stderr, "\nP2: kmin_bottom (minimum degree for the knn graph)");
		fprintf(stderr, "\nP3: kmax_bottom (maximum degree for the knn graph)");
        fprintf(stderr, "\nP4: kmin_top (minimum degree for the knn graph)");
		fprintf(stderr, "\nP5: kmax_top (maximum degree for the knn graph)");
        fprintf(stderr, "\nP6: How many levels to use. Defaults to 1.");
		fprintf(stderr, "\nP7: The least amount of prototypes to generate a new level. Default to 2");
		fprintf(stderr, "\nP8: P3 0 (height), 1(area) and 2 (volume)");
		fprintf(stderr, "\nP9: value of parameter P3 in (0-1)");
		fprintf(stderr, "\nP10: precomputed distance file (leave it in blank if you are not using this resource\n\n");
		exit(-1);
    }
}

void compute_bestk(Subgraph* g, int kmin, int kmax) {
    fprintf(stdout, "Computing k* on [%d %d]", kmin, kmax);
    if (kmin > kmax || kmin < 1) {
		fprintf(stderr, "Kmin must be smaller or equal than Kmax and Kmin must be larger than 0.\n");
		exit(-1);	
    }

    opf_BestkMinCut(g, kmin, kmax);
}

void label_samples(Subgraph *g, int amount_labels) {
    int most_frequent;
    int top_index;
    int i, j;

    if (amount_labels == 0) {
        for (i = 0; i < g->nnodes; i++) {
            g->node[i].label = g->node[i].label + 1;
            g->node[i].truelabel = g->node[i].label;
        }
        return;
    }

    fprintf(stdout, "\nLabeling dataset cluster majority voting.");
    for (i = 0; i < g->nnodes; i++) {
        if (g->node[i].children || g->node[i].root == i) {
            // adding 1 because class labels are 1-indexed, need to make room for the last element
            int *frequency = (int*)calloc((amount_labels + 1), sizeof(int));
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
            most_frequent = 0;
            top_index = 0;
            for (j = 0; j < amount_labels + 1; j++) {
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
        }
    }
}

int count_prototypes(Subgraph *g) {
    int counter = 0;
    for (int i = 0; i < g->nnodes; i++) {
        if (g->node[i].pred == NIL) {
            counter += 1;
        }
    }
    return counter;
}

void write_prototypes(Subgraph *g, char* filepath) {
    FILE* file = fopen(filepath, "w");
    if (!file) {
        fprintf(stderr, "Failed to open file to write prototypes\n");
        exit(-1);
    }
    fprintf(file, "-\n");

    for (int i = 0; i < g->nnodes; i++) {
        fprintf(file, "%d %d ", g->node[i].position, g->node[i].label);
        for (int j = 0; j < g->nfeats; j++) {
            fprintf(file, "%f ", g->node[i].feat[j]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

int main(int argc, char **argv) {
    int amount_samples_dist;
	int init_amount_samples;
    int amount_labels;
    int amount_levels;
    int amount_prototypes = 2;
    int min_amount_prototypes;
    char classifier_name[32];
    char output_filename[32];
	char datafile[32];

    validate(argc, argv);

    amount_levels = atoi(argv[ARG_AMOUNT_LEVELS]);
    if (amount_levels < 1) {
        amount_levels = 1;
    }

    min_amount_prototypes = atoi(argv[ARG_MIN_PROTOTYPES]);
    if (min_amount_prototypes < 2) {
        min_amount_prototypes = 2;
    }
    amount_prototypes = min_amount_prototypes + 1;

    if (argc == 7) {
        // loads the precomputed distances data to global variables
        fprintf(stdout, "Loading precomputed distances.\n");
        opf_PrecomputedDistance = 1;
        opf_DistanceValue = opf_ReadDistances(argv[ARG_DISTANCES], &amount_samples_dist);
    }

    if (atoi(argv[ARG_FILTERING_OPTION]) > 0) {
        fprintf(stderr, "Cluster filtering (parameters P4 and P5) not yet supported.\n");
        exit(-1);
    }

    fprintf(stdout, "Reading the subgraph.\n");
    Subgraph* graph = ReadSubgraph(argv[ARG_SAMPLES_FILE]);
    amount_labels = graph->nlabels;
	init_amount_samples = graph->nnodes;

    for (int level = 0; level < amount_levels; level++) {
        fprintf(stdout, "=============\nLevel %d of %d\n=============\n", level + 1, amount_levels);

        if (amount_prototypes < min_amount_prototypes) {
            amount_levels = level;
            fprintf(stderr, "\nThe dataset contains less than %d prototypes. OPF will stop now.\n", min_amount_prototypes);
            break;
        }

        // level is 0-indexed, the user inputs the amount of levels.
        if (level == amount_levels - 1) {
            compute_bestk(graph, atoi(argv[ARG_KMIN_TOP]), atoi(argv[ARG_KMAX_TOP]));
        }
        else {
            compute_bestk(graph, atoi(argv[ARG_KMIN_BOTTOM]), atoi(argv[ARG_KMAX_BOTTOM]));
        }

        fprintf(stdout, "\nClustering.\n");
        opf_OPFClustering(graph);
        fprintf(stdout, "Amount of clusters %d.\n", graph->nlabels);

        fprintf(stdout, "Labelling the clusters.\n");
        label_samples(graph, amount_labels);

        // writing this level dataset with its corresponding labels         
		fprintf(stdout, "Persisting data partition file\n");
        sprintf(datafile, "data_%d.dat", level);

        WriteSubgraph(graph, datafile);
        fprintf(stdout, "Generating classifier file.\n");

        // writing the classifier of this level to restore later on
        sprintf(classifier_name, "classifier_%d.opf", level);
        opf_WriteModelFile(graph, classifier_name);

        fprintf(stdout, "Generating output files.\n");
        amount_prototypes = count_prototypes(graph);

        fprintf(stdout, "Removing non-prototypes.\n");
        Subgraph* new_graph = CreateGraphFromPrototypes(graph);

        fprintf(stdout, "Deallocating previous graph.\n");
        DestroySubgraph(&graph);

        // making the cycle begin again.
        graph = new_graph;
    }

    write_prototypes(graph, "top_protos.txt");

    DestroySubgraph(&graph);

    // propagating labels from *top* to *bottom*
    Subgraph *classifier, *samples;
    int *label_map = AllocIntArray(init_amount_samples);
    int node_index;

    for (int level = amount_levels - 1; level > 0; level--) {
        fprintf(stdout, "================================================================\n");
        fprintf(stdout, "Loading classifier for level %d - Propagating labels to level %d\n", level + 1, level);
        fprintf(stdout, "================================================================\n");

        sprintf(classifier_name, "classifier_%d.opf", level);
        fprintf(stdout, "%s\n", classifier_name);
        sprintf(output_filename, "data_%d.dat", level - 1);
        fprintf(stdout, "%s\n", output_filename);

        fprintf(stdout, "Loading dataset %s and classifier %s\n", output_filename, classifier_name);
        classifier = opf_ReadModelFile(classifier_name);
        fprintf(stdout, "...\n");
        samples = ReadSubgraph(output_filename);

        // if this is the top level, there's no previous level labels to propagate down
        if (level != amount_levels - 1) {
            fprintf(stdout, "Propagating labels from %d to %d\n", level, level - 1);

            // the current classifier nodes were the sample nodes on the previous iteration. This snippet
            // propagates the top-level labels downwards the pyramid.
            for (int i = 0; i < classifier->nnodes; i++) {
                node_index = classifier->node[i].position;
                classifier->node[i].label = label_map[node_index];

                if (classifier->node[i].label == 0) {
                    fprintf(stderr, "Trying to access node that was not labelled on the previous level!\n");
                    exit(-1);
                }
            }
        }

        fprintf(stdout, "KNN classifying\n");
        opf_OPFknnClassify(classifier, samples);

        // correction factor. Pushing prototype nodes down the tree.
        for (int i = 0; i < classifier->nnodes; i++) {
            int helper = classifier->node[i].position;
            samples->node[helper].label = classifier->node[i].label;
        }

        fprintf(stdout, "Storing labels\n");
        for (int i = 0; i < samples->nnodes; i++) {
            node_index = samples->node[i].position;
            label_map[node_index] = samples->node[i].label;
        }

        fprintf(stdout, "\nDeallocating memory\n");

		if (level > 1) {
			DestroySubgraph(&samples);
			DestroySubgraph(&classifier);
		}
    }

    fprintf(stdout, "Done.\n");

    // deallocating memory
    if (opf_PrecomputedDistance) {
		for (int i = 0; i < amount_samples_dist; i++) {
            free(opf_DistanceValue[i]);
        }
		free(opf_DistanceValue);
	}

    // writing the final level on disk
    sprintf(output_filename, "%s.out", argv[ARG_SAMPLES_FILE]);
	FILE *f = fopen(output_filename, "w");
	for (int i = 0; i < samples->nnodes; i++)
		fprintf(f, "%d\n", samples->node[i].label);
	fclose(f);

    opf_WriteModelFile(classifier, "classifier.opf");

	DestroySubgraph(&samples);
	DestroySubgraph(&classifier);

    return 0;
}
