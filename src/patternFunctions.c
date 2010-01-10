/*
   (c) Massachusetts of Technology, 2003
 */
# include <stdlib.h>
# include <string.h>
# include <errno.h>
# include <math.h>
# include "patternFunctions.h"
# include "fastaSeqIO.h"

#define BUFFER 100000
#define BIG_BUFFER 1000000

// This will print a Teiresias pattern object
int
printTPat(FILE * OUTPUT, tPat_t * myTeiresiasPattern, int hasOffsets)
{
	int i;

	fprintf(OUTPUT, "%d\t", myTeiresiasPattern->support);
	fprintf(OUTPUT, "%d\t", myTeiresiasPattern->support2);
	if (myTeiresiasPattern->logOdds != 0) {
		fprintf(OUTPUT, "%f\t", myTeiresiasPattern->logOdds);
	}
	fprintf(OUTPUT, "%s", myTeiresiasPattern->string);
	fflush(OUTPUT);
	if (hasOffsets == 1) {
		for (i = 0; i < myTeiresiasPattern->support - 1; i++) {
			fprintf(OUTPUT, " %d %d", myTeiresiasPattern->offset[i].seq,
					myTeiresiasPattern->offset[i].pos);
		}
	}
	return EXIT_SUCCESS;
}


// This function takes a file handle and a int boolean (0 or 1) for
// if you want offsets or not.  It will read all the patterns in an
// output file from Teiresias, which may or may not have logOdds values
// 
tPat_t *
ReadTPats(FILE * INPUT, int getOffsets, int *numberOfPatterns, int (*parseFunct) ())
{
	int i;
	int countedPatterns;
	char *buf;
	long tls = BUFFER;
	long ls;
	tPat_t *arrayOfTeiresiasPatterns;


	// --- Count the patterns
	*numberOfPatterns = countLines(INPUT);
	rewind(INPUT);

	// --- Allocate an array of teiresiasPattern objects
	arrayOfTeiresiasPatterns = (tPat_t *) malloc((*numberOfPatterns) * sizeof(tPat_t));
	if (arrayOfTeiresiasPatterns == NULL) {
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(EXIT_FAILURE);
	}
	// Allocate the initial buffer space
	buf = (char *) malloc(tls * sizeof(char));
	if (buf == NULL) {
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(EXIT_FAILURE);
	}
	// Check to make sure line is shorter than 
	// the space we've allocated to hold it
	ls = measureLine(INPUT);
	if (ls > tls) {
		tls = ls + 1;
		buf = (char *) realloc(buf, tls * sizeof(char));
		if (buf == NULL) {
			fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
			fflush(stderr);
			free(arrayOfTeiresiasPatterns);
			arrayOfTeiresiasPatterns = NULL;
			return arrayOfTeiresiasPatterns;
		}
	}
	// Read through each line 
	countedPatterns = 0;
	while (fgets(buf, tls, INPUT) != NULL) {

		// Get rid of the newline
		if (buf[strlen(buf) - 1] == '\n') {
			buf[strlen(buf) - 1] = '\0';
		}

		if (buf[0] == '#' || strlen(buf) == 0) {
			continue;
		}
		i = parseFunct(buf, getOffsets, &(arrayOfTeiresiasPatterns[countedPatterns]));

		countedPatterns++;

		// Check to make sure line is shorter than 
		// the space we've allocated to hold it
		ls = measureLine(INPUT);
		if (ls > tls) {
			tls = ls + 1;
			buf = (char *) realloc(buf, tls * sizeof(char));
			if (buf == NULL) {
				fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
				fflush(stderr);
				free(arrayOfTeiresiasPatterns);
				arrayOfTeiresiasPatterns = NULL;
				return arrayOfTeiresiasPatterns;
			}
		}
	}

	// resize array of teiresiasPattern objects
	*numberOfPatterns = countedPatterns;
	arrayOfTeiresiasPatterns =
		(tPat_t *) realloc(arrayOfTeiresiasPatterns, (*numberOfPatterns) * sizeof(tPat_t));
	if (arrayOfTeiresiasPatterns == NULL) {
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		return arrayOfTeiresiasPatterns;
	}

	rewind(INPUT);
	free(buf);
	return arrayOfTeiresiasPatterns;
}


// Counts all the Teiresias patterns in a file given a file handle
int
CountTPats2(FILE * INPUT)
{
	char buffer[BUFFER];
	int numberOfPatterns = 0;
	while (fgets(buffer, BUFFER, INPUT) != NULL) {
		if (buffer[0] != '#') {
			numberOfPatterns++;
		}
	}
	return numberOfPatterns;
}

int
MeasurePattern(char *myPattern)
{
	int j = 0;
	int size = 0;
	while ((unsigned int) j < strlen(myPattern)) {
		if (myPattern[j] == '[') {
			j++;
			while (myPattern[j] != ']') {
				j++;
			}
			size++;
		} else {
			size++;
		}
		j++;
	}
	return size;
}

// Parses a line of a Teireisas file and returns a teiresiasPattern
// object with or without an offset list.
int
ParseTPatLine(char *buffer, int getOffsets, tPat_t * myTeiresiasPattern)
{
	int a, b;
	double c;
	char myPattern[PATTERN_SIZE];
	char *placeHolder;
	int e, f;
	int countedOffsets;


	if (sscanf(buffer, "%*d%*d%lf", &c) == 1 && finite(c) != 0) {
		// there ARE logOdds valuedd
		sscanf(buffer, "%d%d%lf%s", &a, &b, &c, myPattern);

	} else if (sscanf(buffer, "%d%d%*s", &a, &b) == 2) {
		// there are no logOdds values
		sscanf(buffer, "%d%d%s", &a, &b, myPattern);
		c = 0;
	} else {
		sscanf(buffer, "%s", myPattern);
		a = 0;
		b = 0;
		c = 0;
	}

	// Make the teiresiasPattern object
	myTeiresiasPattern->size = MeasurePattern(myPattern);
	myTeiresiasPattern->support = a;
	myTeiresiasPattern->support2 = b;
	myTeiresiasPattern->logOdds = c;
	strcpy(myTeiresiasPattern->string, myPattern);

	if (getOffsets == 1) {

		myTeiresiasPattern->offset = (offsetList_t *) malloc(a * sizeof(offsetList_t));
		if (myTeiresiasPattern == NULL) {
			fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
			fflush(stderr);
			exit(EXIT_FAILURE);
		}

		countedOffsets = 0;
		placeHolder = buffer;

		// while there are spaces after the pattern
		while ((placeHolder = strstr(placeHolder, " "))) {

			// scan in a (seq, pos) pair
			sscanf(placeHolder, "%d%d", &e, &f);

			// add these to the offset list for this pattern
			myTeiresiasPattern->offset[countedOffsets].seq = e;
			myTeiresiasPattern->offset[countedOffsets].pos = f;

			// go to the next pair
			placeHolder++;
			placeHolder = strstr(placeHolder, " ");
			placeHolder++;
			countedOffsets++;
		}
	} else {
		myTeiresiasPattern->offset = NULL;
	}


	return EXIT_SUCCESS;
}


// Parses a line of a Teireisas file and returns a teiresiasPattern
// object with or without an offset list.
int
ParseTxtLine(char *buffer, int getOffsets, tPat_t * myTeiresiasPattern)
{
	char myPattern[PATTERN_SIZE];
	strcpy(myPattern, buffer);
	if (myPattern[strlen(myPattern) - 1] == '\n') {
		myPattern[strlen(myPattern) - 1] = '\0';
	}
	// Make the teiresiasPattern object
	myTeiresiasPattern->size = MeasurePattern(myPattern);
	myTeiresiasPattern->support = 0;
	myTeiresiasPattern->support2 = 0;
	myTeiresiasPattern->logOdds = 0;
	strcpy(myTeiresiasPattern->string, myPattern);
	myTeiresiasPattern->offset = NULL;

	return EXIT_SUCCESS;
}


int
FreeTPatA(tPat_t * arrayOfTeiresiasPatterns, int numberOfPatterns)
{
	int i;
	for (i = 0; i < numberOfPatterns; i++) {
		if (arrayOfTeiresiasPatterns[i].offset != NULL) {
			free(arrayOfTeiresiasPatterns[i].offset);
			arrayOfTeiresiasPatterns[i].offset = NULL;
		}
	}
	free(arrayOfTeiresiasPatterns);
	arrayOfTeiresiasPatterns = NULL;
	return EXIT_SUCCESS;
}
