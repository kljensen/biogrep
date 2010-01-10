/*
   (c) Massachusetts of Technology, 2003 
 */

// needed for AIX before pthread 
#ifdef _AIX
struct timespec;
#undef _XOPEN_SOURCE
#define _XOPEN_SOURCE  499
#endif


#include <pthread.h>

#include "fastaSeqIO.h"
#include "patternFunctions.h"
#include "redRegEx.h"
#include <unistd.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#define MYBUFFER 10000


void
usage(char **argv)
{
	fprintf(stdout, "Usage: %s [OPTIONS] -p <pattern file> -i <text file>\n", argv[0]);
	fprintf(stdout, "Try '%s -help' for more information.\n", argv[0]);
}

void
help(char **argv)
{
	fprintf(stdout,
			"\nUsage: %s [OPTIONS] -p <file with patterns> -i <file with FastA>\n", argv[0]);
	fprintf(stdout, "\nSearches for regular expressions in files.");
	fprintf(stdout, "\nFeatures multithreading and FastA-formatted");
	fprintf(stdout, "\ninput options --- see options below.");
	fprintf(stdout, "\n\n");
	fprintf(stdout, "Regex selection and interpretation:");
	fprintf(stdout, "\n  -p FILE1\tget patterns from FILE1, one pattern per line");
	fprintf(stdout, "\n  -P FILE1\tget patterns from Teiresias-formatted FILE1");
	fprintf(stdout, "\n  -c\t\tignore case when matching patterns");
	fprintf(stdout, "\n\n");
	fprintf(stdout, "Query file selection:");
	fprintf(stdout, "\n  -i FILE2\tsearch for patterns in FastA-formatted sequence FILE2");
	fprintf(stdout, "\n  -t\t\ttreat FILE2 as plain text, not FastA-formatted");
	fprintf(stdout, "\n\n");
	fprintf(stdout, "Output format:");
	fprintf(stdout, "\n  -o FILE3\twrite output to FILE3");
	fprintf(stdout, "\n  -T \t\twrite output in Teiresias format");
	fprintf(stdout, "\n  -L \t\talways print FastA sequence labels for matching lines");
	fprintf(stdout, "\n\t\t\t(ignored if either the -t or -T flag is used)");
	fprintf(stdout, "\n  -A N \t\tprint N lines (or sequences if the -t flag is used)");
	fprintf(stdout, "\n\t\t\tafter the matching line (ignored if -T flag is given)");
	fprintf(stdout, "\n  -B N \t\tprint N lines (or sequences if the -t flag is used)");
	fprintf(stdout, "\n\t\t\tbefore the matching line (ignored if -T flag is given)");
	fprintf(stdout, "\n\n");
	fprintf(stdout, "Other options:");
	fprintf(stdout, "\n  -x N\t\tuse N threads (processors)");
	fprintf(stdout, "\n\n");
}

typedef struct {
	int outputFormat;
	int linesBefore;
	int linesAfter;
	int printLabels;
} printFormat_t;

int
processOutput(FILE ** OUTPUT, FILE * OUTPUT_FILE, fSeq_t * mySequences, int numberOfSequences,
			  int numberOfThreads, printFormat_t * myFormat)
{
	int *seqsToPrint = NULL;
	char *myLine1;
	char *myLine2;
	long lineSize;
	int i, j, k;
	int outputFormat = myFormat->outputFormat;
	int linesBefore = myFormat->linesBefore;
	int linesAfter = myFormat->linesAfter;
	int printLabels = myFormat->printLabels;

	// If we're using plain text output, then we need to record
	// which sequences we're going to print out.
	if (outputFormat == 0) {
		seqsToPrint = (int *) malloc(numberOfSequences * sizeof(int));
		if(seqsToPrint == NULL){
			fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
			fflush(stderr);
			exit(EXIT_FAILURE);
		}
		for (i = 0; i < numberOfSequences; i++) {
			seqsToPrint[i] = 0;
		}
	}
	// Look at the output of each thread.  Each file has output like
	// LINE1: 0 1 2 1 3 ....offsets
	// LINE2: support1 support2

	// Allocate some space to hold the lines
	myLine1 = (char *) malloc(MYBUFFER * sizeof(char));
	if(myLine1 == NULL){
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(EXIT_FAILURE);
	}
	myLine2 = (char *) malloc(MYBUFFER * sizeof(char));
	if(myLine2 == NULL){
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(EXIT_FAILURE);
	}
	lineSize = MYBUFFER;


	for (i = 0; i < numberOfThreads; i++) {
		rewind(OUTPUT[i]);

		// Check if the line is longer than the memory we've allocated
		if (lineSize < measureLine(OUTPUT[i])) {
			myLine1 = (char *) realloc(myLine1, measureLine(OUTPUT[i]) * sizeof(char));
			myLine2 = (char *) realloc(myLine2, measureLine(OUTPUT[i]) * sizeof(char));
			lineSize = measureLine(OUTPUT[i]);
		}
		while (fgets(myLine1, (int) lineSize + 1, OUTPUT[i]) != NULL) {
			if (outputFormat == 1) {
				// Check if the line is longer than the memory we've allocated
				if (lineSize < measureLine(OUTPUT[i])) {
					myLine1 = (char *) realloc(myLine1, measureLine(OUTPUT[i]) * sizeof(char));
					myLine2 = (char *) realloc(myLine2, measureLine(OUTPUT[i]) * sizeof(char));
					lineSize = measureLine(OUTPUT[i]);
				}
				fgets(myLine2, (int) lineSize + 1, OUTPUT[i]);
				myLine2[strlen(myLine2) - 1] = '\0';
				/*
				   fprintf(OUTPUT_FILE, "%s\t%s", myLine2, myLine1); 
				 */
				fprintf(OUTPUT_FILE, "%s\t", myLine2);
				fprintf(OUTPUT_FILE, "%s", myLine1);

			} else {
				sscanf(myLine1, "%d", &j);
				for (k = j - linesBefore; k <= j + linesAfter; k++) {
					seqsToPrint[k] = 1;
				}
			}
			// Check if the line is longer than the memory we've allocated
			if (lineSize < measureLine(OUTPUT[i])) {
				myLine1 = (char *) realloc(myLine1, measureLine(OUTPUT[i]) * sizeof(char));
				myLine2 = (char *) realloc(myLine2, measureLine(OUTPUT[i]) * sizeof(char));
				lineSize = measureLine(OUTPUT[i]);
			}
		}
	}
	if (outputFormat == 0) {
		for (i = 0; i < numberOfSequences; i++) {
			if (seqsToPrint[i] == 1) {
				if (printLabels == 1) {
					fprintf(OUTPUT_FILE, "%s\n", mySequences[i].label);
				}
				fprintf(OUTPUT_FILE, "%s\n", mySequences[i].seq);
			}
		}
		free(seqsToPrint);
	}
	free(myLine1);
	free(myLine2);
	return EXIT_SUCCESS;
}

int
main(int argc, char **argv)
{
	int inputOption = 0;
	char *sequenceFile = NULL;
	char *patternFile = NULL;
	char *outputFile = NULL;
	FILE *SEQUENCE_FILE = NULL;
	FILE *PATTERN_FILE = NULL;
	FILE *OUTPUT_FILE = NULL;
	FILE **OUTPUT = NULL;
	int numberOfSequences = 0;
	int numberOfPatterns = 0;
	fSeq_t *mySequences = NULL;

	tPat_t *arrayOfTeiresiasPatterns;
	redRegEx_t *myArrayOfRedRegExs;

	pthread_t *myThreads;
	pthread_attr_t tattr;
	int ret;
	mtMatchPass_t *myMtMatchPass;
	int numberOfThreads = 1;
	int i;
	int regExsPerThread;
	int completedRegExs;
	int (*parseFunct) () = &ParseTxtLine;
	fSeq_t *(*seqReadFunct) () = &ReadTxtSeqs;
	printFormat_t myFormat;
	int ignoreCase = 0;

	myFormat.linesBefore = 0;
	myFormat.linesAfter = 0;
	myFormat.outputFormat = 0;
	myFormat.printLabels = 0;

	/*
	   Get use command-line options 
	 */
	while ((inputOption = getopt(argc, argv, "i:tp:o:P:x:hTcA:B:L")) != EOF) {
		switch (inputOption) {
		case 'i':
			sequenceFile = optarg;
			seqReadFunct = &ReadFSeqs;
			break;
		case 'c':
			ignoreCase = REG_ICASE;
			break;
		case 't':
			seqReadFunct = &ReadTxtSeqs;
			myFormat.printLabels+=2;
			break;
		case 'p':
			patternFile = optarg;
			parseFunct = &ParseTxtLine;
			break;
		case 'P':
			patternFile = optarg;
			parseFunct = &ParseTPatLine;
			break;
		case 'x':
			numberOfThreads = atoi(optarg);
			break;
		case 'o':
			outputFile = optarg;
			break;
		case 'h':
			help(argv);
			return EXIT_SUCCESS;
		case 'B':
			myFormat.linesBefore = atoi(optarg);
			break;
		case 'A':
			myFormat.linesAfter = atoi(optarg);
			break;
		case 'T':
			myFormat.outputFormat = 1;
			break;
		case 'L':
			myFormat.printLabels++;
			break;
		case '?':
			fprintf(stderr, "Unknown option `-%c'.\n", optopt);
			usage(argv);
			return EXIT_SUCCESS;
		default:
			usage(argv);
			return EXIT_SUCCESS;
		}
	}
	if (sequenceFile == NULL || patternFile == NULL) {
		usage(argv);
		return EXIT_SUCCESS;
	}
	// Open the sequence file
	if ((SEQUENCE_FILE = fopen(sequenceFile, "r")) == NULL) {
		fprintf(stderr, "Couldn't open file %s; %s\n", sequenceFile, strerror(errno));
		exit(EXIT_FAILURE);
	}
	// Open the pattern file
	if ((PATTERN_FILE = fopen(patternFile, "r")) == NULL) {
		fprintf(stderr, "Couldn't open file %s; %s\n", patternFile, strerror(errno));
		exit(EXIT_FAILURE);
	}
	// Open the output file
	if (outputFile != NULL) {
		if ((OUTPUT_FILE = fopen(outputFile, "w")) == NULL) {
			fprintf(stderr, "Couldn't open file %s; %s\n", outputFile, strerror(errno));
			exit(EXIT_FAILURE);
		}
	} else {
		OUTPUT_FILE = stdout;
	}

	// Open a series of temporary files
	OUTPUT = (FILE **) malloc(numberOfThreads * sizeof(FILE *));
	if(OUTPUT == NULL){
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(EXIT_FAILURE);
	}
	for (i = 0; i < numberOfThreads; i++) {
		if ((OUTPUT[i] = tmpfile()) == NULL) {
			fprintf(stderr, "Couldn't open file; %s\n", strerror(errno));
			exit(EXIT_FAILURE);
		}
	}


	// Allocate some sequences
	mySequences = seqReadFunct(SEQUENCE_FILE, &numberOfSequences);
	if (mySequences == NULL) {
		fprintf(stderr, "\nError reading your sequences/text.");
		fprintf(stderr, "\nCheck the format/size of the file.");
		fprintf(stderr, "\nERROR:  %s\n", strerror(errno));
		return EXIT_FAILURE;
	}
	// Read in the patterns
	arrayOfTeiresiasPatterns = ReadTPats(PATTERN_FILE, 1, &numberOfPatterns, parseFunct);
	if (arrayOfTeiresiasPatterns == NULL) {
		fprintf(stderr, "\nError reading your patterns.");
		fprintf(stderr, "\nCheck the format/size of the file.");
		fprintf(stderr, "\nERROR:  %s\n", strerror(errno));
		FreeFSeqs(mySequences, numberOfSequences);
		return EXIT_FAILURE;
	}
	// Make array of regexs
	myArrayOfRedRegExs =
		makeRedRegExAFromTPatA(arrayOfTeiresiasPatterns, numberOfPatterns, ignoreCase);
	if (myArrayOfRedRegExs == NULL) {
		fprintf(stderr, "\nError compiling your regular expressions.");
		fprintf(stderr, "\nCheck the format/size of the file.");
		fprintf(stderr, "\nERROR:  %s\n", strerror(errno));
		FreeFSeqs(mySequences, numberOfSequences);
		FreeTPatA(arrayOfTeiresiasPatterns, numberOfPatterns);
		return EXIT_FAILURE;
	}
	// Search sequences for regexs
	myThreads = (pthread_t *) malloc(numberOfThreads * sizeof(pthread_t));
	if(myThreads == NULL){
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(EXIT_FAILURE);
	}
	myMtMatchPass = (mtMatchPass_t *) malloc(numberOfThreads * sizeof(mtMatchPass_t));
	if(myMtMatchPass == NULL){
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(EXIT_FAILURE);
	}
	regExsPerThread = (int) numberOfPatterns / numberOfThreads;
	completedRegExs = 0;
	for (i = 0; i < numberOfThreads; i++) {
		myMtMatchPass[i].myArrayOfRedRegExs = &myArrayOfRedRegExs[completedRegExs];
		myMtMatchPass[i].myArrayOfFastaSequences = mySequences;
		if (i == numberOfThreads - 1 || completedRegExs + regExsPerThread > numberOfPatterns) {
			myMtMatchPass[i].numberOfPatterns = numberOfPatterns - completedRegExs;
		} else {
			myMtMatchPass[i].numberOfPatterns = regExsPerThread;
		}
		myMtMatchPass[i].numberOfSequences = numberOfSequences;
		myMtMatchPass[i].OUTPUT = OUTPUT[i];
		/*
		   myMtMatchPass[i].OUTPUT = stdout;
		 */
		myMtMatchPass[i].threadNumber = i;
		myMtMatchPass[i].outputFormat = myFormat.outputFormat;


		ret = pthread_attr_init(&tattr);
		pthread_create(&myThreads[i], NULL, (void *(*)(void *)) MTmatchRedRedRegExAOnFSeqA,
					   (void *) &myMtMatchPass[i]);
		completedRegExs += myMtMatchPass[i].numberOfPatterns;
	}

	// Wait for the threads to complete
	for (i = 0; i < numberOfThreads; i++) {
		ret = pthread_join(myThreads[i], NULL);
	}
	pthread_attr_destroy(&tattr);
	free(myMtMatchPass);

	// Free up threads
	free(myThreads);

	/*WriteFSeqA(stdout, mySequences, 0, numberOfSequences-1);*/

	processOutput(OUTPUT, OUTPUT_FILE, mySequences, numberOfSequences, numberOfThreads, &myFormat);

	for (i = 0; i < numberOfThreads; i++) {
		fclose(OUTPUT[i]);
	}

	// Free output file handles
	free(OUTPUT);

	// Free up redRegExs
	freeRedRegExA(myArrayOfRedRegExs, numberOfPatterns);

	// Free up fastaSequences
	FreeFSeqs(mySequences, numberOfSequences);

	// Free up patterns
	FreeTPatA(arrayOfTeiresiasPatterns, numberOfPatterns);

	// Close the input files
	fclose(SEQUENCE_FILE);
	fclose(PATTERN_FILE);
	fclose(OUTPUT_FILE);
	return 0;
}
