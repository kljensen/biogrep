/*
   (c) Massachusetts of Technology, 2003
 */
# include <stdlib.h>
# include <string.h>
# include <errno.h>
# include <regex.h>

/*
   needed for AIX before pthread 
 */
struct timespec;
#undef _XOPEN_SOURCE
#define _XOPEN_SOURCE  499
# include <pthread.h>

# include "patternFunctions.h"
# include "fastaSeqIO.h"
# include "redRegEx.h"


void
MTmatchRedRedRegExAOnFSeqA(void *arg)
{
	mtMatchPass_t *input;
	redRegEx_t *myArrayOfRedRegExs;
	fSeq_t *myArrayOfFastaSequences;
	int numberOfPatterns;
	int numberOfSequences;
	FILE *OUTPUT;
	int i;
	int threadNumber;
	int outputFormat;

	input = (mtMatchPass_t *) arg;
	myArrayOfRedRegExs = input->myArrayOfRedRegExs;
	myArrayOfFastaSequences = input->myArrayOfFastaSequences;
	numberOfPatterns = input->numberOfPatterns;
	numberOfSequences = input->numberOfSequences;
	OUTPUT = input->OUTPUT;
	threadNumber = input->threadNumber;
	outputFormat = input->outputFormat;



	for (i = 0; i < numberOfPatterns; i++) {
		if (outputFormat == 1) {
			fprintf(OUTPUT, "%s", myArrayOfRedRegExs[i].string);
		}
		MTMatchRedRedRegExOnFSeqA(&(myArrayOfRedRegExs[i]), myArrayOfFastaSequences,
								  numberOfSequences, OUTPUT, outputFormat);
	}
}

int
MTMatchRedRedRegExOnFSeqA(redRegEx_t * myRedRegEx, fSeq_t * myArrayOfFastaSequences,
						  int numberOfSequences, FILE * OUTPUT, int outputFormat)
{
	int j;
	char outputPrefix[1000];
	int offset;
	regmatch_t *myMatch = NULL;
	size_t no_sub = myRedRegEx->regex->re_nsub + 1;
	int support1, support2;
	int status;


	if ((myMatch = (regmatch_t *) malloc(sizeof(regmatch_t) * no_sub)) == NULL) {
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(EXIT_FAILURE);
	}


	no_sub = myRedRegEx->regex->re_nsub + 1;
	support1 = 0;
	support2 = 0;
	for (j = 0; j < numberOfSequences; j++) {
		offset = 0;
		sprintf(outputPrefix, " %d", j);
		status = 0;
		while (regexec
			   (myRedRegEx->regex, myArrayOfFastaSequences[j].seq + offset, no_sub, myMatch,
				0) == 0) {
			if (outputFormat == 1) {
				fprintf(OUTPUT, "%s %d", outputPrefix, (int) myMatch->rm_so + offset);
			} else {
				fprintf(OUTPUT, "%d\n", j);
				break;
			}
			offset += myMatch->rm_so + 1;
			if (status == 0) {
				support1++;
				support2++;
				status = 1;
			} else {
				support1++;
			}
		}
	}
	if (outputFormat == 1) {
		fprintf(OUTPUT, "\n%d\t%d\n", support1, support2);
	}
	free(myMatch);
	return EXIT_SUCCESS;
}

int
matchRedRedRegExAOnFSeqA(redRegEx_t * myArrayOfRedRegExs, fSeq_t * myArrayOfFastaSequences,
						 int numberOfPatterns, int numberOfSequences, FILE * OUTPUT)
{
	int i, j;
	char outputPrefix[1000];
	for (i = 0; i < numberOfPatterns; i++) {
		fprintf(OUTPUT, "%s", myArrayOfRedRegExs[i].string);
		for (j = 0; j < numberOfSequences; j++) {
			sprintf(outputPrefix, " %d", j);
			matchRedRegExOnFSeq(&(myArrayOfRedRegExs[i]), &(myArrayOfFastaSequences[j]),
								outputPrefix, OUTPUT);
		}
		fprintf(OUTPUT, "\n");
	}
	return EXIT_SUCCESS;
}


int
matchRedRegExOnFSeq(redRegEx_t * myRedRegEx, fSeq_t * myFastaSequence,
					char *prefixString, FILE * OUTPUT)
{
	int offset;
	regmatch_t *myMatch;
	size_t no_sub = myRedRegEx->regex->re_nsub + 1;

	/*
	   Allocate space for the regexmatch 
	 */
	if ((myMatch = (regmatch_t *) malloc(sizeof(regmatch_t) * no_sub)) == NULL) {
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(EXIT_FAILURE);
	}

	offset = 0;
	while (regexec(myRedRegEx->regex, myFastaSequence->seq + offset, no_sub, myMatch, 0) == 0) {
		fprintf(OUTPUT, "%s %d", prefixString, (int) myMatch->rm_so + offset);
		offset += myMatch->rm_so + 1;
	}
	free(myMatch);
	return EXIT_SUCCESS;
}

redRegEx_t *
makeRedRegExAFromTPatA(tPat_t * myArrayOfTeiresiasPatterns, int numberOfPatterns, int ignoreCase)
{
	int i;
	redRegEx_t *myRedRegExArray;
	myRedRegExArray = (redRegEx_t *) malloc(numberOfPatterns * sizeof(redRegEx_t));
	if (myRedRegExArray == NULL) {
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(EXIT_FAILURE);
	}
	for (i = 0; i < numberOfPatterns; i++) {
		myRedRegExArray[i] = makeRedRegExFromTPat(&(myArrayOfTeiresiasPatterns[i]), ignoreCase);
		if (myRedRegExArray[i].regex == NULL) {
			fprintf(stderr, "Check format of pattern number %d!\n", i);
			free(myRedRegExArray);
			myRedRegExArray = NULL;
			return myRedRegExArray;
		}
	}
	return myRedRegExArray;
}

int
freeRedRegExA(redRegEx_t * myArrayOfRedRegExs, int numberOfRegExs)
{
	int i;
	for (i = 0; i < numberOfRegExs; i++) {
		freeRedRegEx(&(myArrayOfRedRegExs[i]));
	}
	free(myArrayOfRedRegExs);
	myArrayOfRedRegExs = NULL;
	return EXIT_SUCCESS;
}

int
freeRedRegEx(redRegEx_t * myRedRegEx)
{
	regfree(myRedRegEx->regex);
	free(myRedRegEx->regex); /* changed to fix memory leak */
	myRedRegEx->regex = NULL;
	return EXIT_SUCCESS;
}

redRegEx_t
makeRedRegExFromTPat(tPat_t * myTeiresiasPattern, int ignoreCase)
{
	int errorNumber = 0;
	size_t errorLength;
	char *myBuffer;
	redRegEx_t myRedRegEx;

	strcpy(myRedRegEx.string, myTeiresiasPattern->string);
	myRedRegEx.regex = (regex_t *) malloc(sizeof(regex_t));
	if (myRedRegEx.regex == NULL) {
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(EXIT_FAILURE);
	}
	memset(myRedRegEx.regex, 0, sizeof(regex_t));

	errorNumber = regcomp(myRedRegEx.regex, myRedRegEx.string, REG_EXTENDED | ignoreCase);
	if (errorNumber != 0) {
		errorLength = regerror(errorNumber, myRedRegEx.regex, (char *) NULL, (size_t) 0);

		myBuffer = (char *) malloc(errorLength);
		if (myBuffer != NULL) {
			regerror(errorNumber, myRedRegEx.regex, myBuffer, errorLength);
			fprintf(stderr, "%s in pattern = %s\n", myBuffer, myRedRegEx.string);
			exit(EXIT_FAILURE);
		}
		regfree(myRedRegEx.regex);
		myRedRegEx.regex = NULL;
		return myRedRegEx;
	}
	return myRedRegEx;
}
