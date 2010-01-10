/*(c) Massachusetts of Technology, 2003*/
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <errno.h>
# include <regex.h>

#ifdef _AIX
    /*
       needed for AIX before pthread 
     */
    struct timespec;
    #undef _XOPEN_SOURCE
    #define _XOPEN_SOURCE  499
#endif

# include <pthread.h>

# include "patternFunctions.h"
# include "fastaSeqIO.h"


#ifndef __RED_REG_EX_H_
#define __RED_REG_EX_H_


typedef struct {
    char string[PATTERN_SIZE];
    regex_t *regex;
} redRegEx_t;

/*
   for multithreading 
 */
typedef struct {
    redRegEx_t *myArrayOfRedRegExs;
    fSeq_t *myArrayOfFastaSequences;
    int numberOfPatterns;
    int numberOfSequences;
    int threadNumber;
    int outputFormat;
    FILE *OUTPUT;
} mtMatchPass_t;

int MTMatchRedRedRegExOnFSeqA(redRegEx_t * myRedRegEx,
			      fSeq_t * myArrayOfFastaSequences,
			      int numberOfSequences, FILE * OUTPUT, int outputFormat);
int matchRedRedRegExAOnFSeqA(redRegEx_t * myArrayOfRedRegExs,
			     fSeq_t * myArrayOfFastaSequences,
			     int numberOfPatterns, int numberOfSequences, FILE * OUTPUT);
int matchRedRegExOnFSeq(redRegEx_t * myRedRegEx, fSeq_t * myFastaSequence,
			char *prefixString, FILE * OUTPUT);
redRegEx_t makeRedRegExFromTPat(tPat_t * myTeiresiasPattern, int ignoreCase);
int freeRedRegEx(redRegEx_t * myRedRegEx);
redRegEx_t *makeRedRegExAFromTPatA(tPat_t * myArrayOfTeiresiasPatterns, int numberOfPatterns, int ignoreCase);
int freeRedRegExA(redRegEx_t * myArrayOfRedRegExs, int numberOfRegExs);
/*
   multi threaded 
 */
void MTmatchRedRedRegExAOnFSeqA(void *arg);


#endif				/* ___RED_REG_EX_H_ */
