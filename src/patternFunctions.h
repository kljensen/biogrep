/*(c) Massachusetts of Technology, 2003*/
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <errno.h>

#ifndef __PATTERN_FUNCTIONS_H_
#define __PATTERN_FUNCTIONS_H_

#define PATTERN_SIZE 10000


/*
 *      Structure for a list of offsets, or integer pairs
 */
typedef struct {
    int seq;
    int pos;
} offsetList_t;

typedef struct {
    int support;
    int support2;
    int size;
    double logOdds;
    char string[PATTERN_SIZE];
    offsetList_t *offset;
} tPat_t;



int CountTPats(FILE * INPUT);
int MeasurePattern(char *pattern);
int ParseTPatLine(char *buffer, int getOffsets, tPat_t * myTeiresiasPattern);
int ParseTxtLine(char *buffer, int getOffsets, tPat_t * myTeiresiasPattern);
tPat_t *ReadTPats(FILE * INPUT, int getOffsets, int *numberOfPatterns, int (*parseFunct)() );
int printTPat(FILE * OUTPUT, tPat_t * myTeiresiasPattern, int hasOffsets);
int FreeTPatA(tPat_t * arrayOfTeiresiasPatterns, int numberOfPatterns);

#endif				/* __PATTERN_FUNCTIONS_H_ */
