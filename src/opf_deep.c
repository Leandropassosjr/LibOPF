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
    if (argc != (ARG_DISTANCES + 1) && argc != ARG_DISTANCES) {
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

void label_clusters(Subgraph *g) {
	for (int i = 0; i < g->nnodes; i++) {
		g->node[i].label = g->node[i].label + 1;
		// g->node[i].truelabel = g->node[i].label;
	}
}

// this function is not used for now.
void add_supervised_labels_to_clusters(Subgraph *g, int amount_labels) {
    int most_frequent;
    int top_index;
    int i, j;

    fprintf(stdout, "\nLabeling dataset cluster majority voting.");
    for (i = 0; i < g->nnodes; i++) {
        if (g->node[i].children || g->node[i].pred == NIL) {
            // adding 1 because class labels are 1-indexed, need to make room for the last element
            int *frequency = (int*)calloc((amount_labels + 1), sizeof(int));
            Set *walker = g->node[i].children;

            // taking into account the prototype label when computing the histogram
            frequency[g->node[i].truelabel] += 1;

            // computing label frequency on the neighborhood of the prototype
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
    // this only works because createGraphFromPrototypes is called before, hence
    // *g contains _only_ prototypes at this moment.

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
    int init_amount_labels;
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

	// adding 1 because this checks *length*, not *index*
    if (argc == ARG_DISTANCES + 1) {
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
    init_amount_labels = graph->nlabels;
	init_amount_samples = graph->nnodes;

    for (int level = 1; level <= amount_levels; level++) {
        fprintf(stdout, "=============\nLevel %d of %d\n=============\n", level, amount_levels);

        if (amount_prototypes < min_amount_prototypes) {
            amount_levels = level - 1;
            fprintf(stderr, "\nThe dataset contains less than %d prototypes. OPF will stop now.\n", min_amount_prototypes);
            break;
        }

        // level is 0-indexed, the user inputs the amount of levels.
        if (level == amount_levels) {
            compute_bestk(graph, atoi(argv[ARG_KMIN_TOP]), atoi(argv[ARG_KMAX_TOP]));
        }
        else {
            compute_bestk(graph, atoi(argv[ARG_KMIN_BOTTOM]), atoi(argv[ARG_KMAX_BOTTOM]));
        }

        fprintf(stdout, "\nClustering.\n");
        opf_OPFClustering(graph);
        fprintf(stdout, "Amount of clusters %d.\n", graph->nlabels);

        fprintf(stdout, "Labelling the clusters.\n");
        label_clusters(graph);

        // writing this level dataset with its corresponding labels         
		fprintf(stdout, "Persisting data partition file\n");
        sprintf(datafile, "data_%d.dat", level);

        WriteSubgraph(graph, datafile);
        fprintf(stdout, "\nGenerating classifier file.\n");

        // writing the classifier of this level to restore later on
        sprintf(classifier_name, "classifier_%d.opf", level);
        opf_WriteModelFile_with_children(graph, classifier_name);

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

	fprintf(stdout, "\n================================================\n");
	fprintf(stdout, "Starting to propagate labels from top to bottom.\n");
	fprintf(stdout, "================================================\n");

    Set *walker = NULL;
    for (int level = amount_levels; level > 1; level--) {
        fprintf(stdout, "Propagating labels from classifier %d to %d\n", level, level - 1);

        sprintf(classifier_name, "classifier_%d.opf", level);
        Subgraph* clf_top = opf_ReadModelFile_with_children(classifier_name);

        sprintf(classifier_name, "classifier_%d.opf", level - 1);
        Subgraph* clf_bottom = opf_ReadModelFile_with_children(classifier_name);

        // propagating labels from top to bottom prototypes
        fprintf(stdout, "propagating...\n");
        for (int i = 0; i < clf_top->nnodes; i++) {
            node_index = clf_top->node[i].position;

            for (int j = 0; j < clf_bottom->nnodes; j++) {
                if (clf_bottom->node[j].position == node_index) {
                    clf_bottom->node[j].label = clf_top->node[i].label;
                    break;
                }
            }
        }

        // propagating new labels on the bottom from prototypes to leaves
        // if the node has children, it is a prototype and must have its leaves updated
        fprintf(stdout, "Propagating labels from prototypes to nodes\n");
        for (int i = 0; i < clf_bottom->nnodes; i++) {
            walker = clf_bottom->node[i].children;
            while (walker) {
                node_index = walker->elem;
                clf_bottom->node[node_index].label = clf_bottom->node[i].label;
				walker = walker->next;
            }
        }
        clf_bottom->nlabels = clf_top->nlabels;

        opf_WriteModelFile_with_children(clf_bottom, classifier_name);

        DestroySubgraph(&clf_top);
        DestroySubgraph(&clf_bottom);
    }

    // deallocating memory
    if (opf_PrecomputedDistance) {
		for (int i = 0; i < amount_samples_dist; i++) {
            free(opf_DistanceValue[i]);
        }
		free(opf_DistanceValue);
	}

    // generating the deep classifier
    fprintf(stdout, "Reading Level 1 classifier.\n");
    classifier = opf_ReadModelFile_with_children("classifier_1.opf");

    // writing the final level on disk
    sprintf(output_filename, "%s.out", argv[ARG_SAMPLES_FILE]);
	FILE *f = fopen(output_filename, "w");
	for (int i = 0; i < classifier->nnodes; i++)
		fprintf(f, "%d\n", classifier->node[i].label);
	fclose(f);

    opf_WriteModelFile(classifier, "classifier.opf");
	DestroySubgraph(&classifier);
    fprintf(stdout, "\nDone.\n");

    return 0;
}
