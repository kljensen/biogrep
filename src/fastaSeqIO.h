/*
   (c) Massachusetts of Technology, 2003
 */
#ifndef __FASTA_SEQ_IO_H_
#define __FASTA_SEQ_IO_H_

#include <stdio.h>




typedef struct {
    char *seq;
    char *label;
} fSeq_t;

long measureLine(FILE * INPUT);
long countLines(FILE * INPUT);
long CountFSeqs(FILE * INPUT);
int initAofFSeqs(fSeq_t *aos, int numSeq);
fSeq_t *ReadFSeqs(FILE * INPUT, int *numberOfSequences);
int FreeFSeqs(fSeq_t * arrayOfSequences, int numberOfSequences);
int WriteFSeqA(FILE * MY_FILE, fSeq_t * arrayOfSequences, int start, int stop);
fSeq_t *ReadTxtSeqs(FILE * INPUT, int *numberOfSequences);

#endif				/* __FASTA_SEQ_IO_H_ */
