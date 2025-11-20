#ifndef PARSER_H
#define PARSER_H

#include "basic.h"

char* load_genome_from_fasta(const char* filename, long* genome_length);
FILE* open_fastq_file(const char* filename);
int read_fastq_to_string(FILE* fp, char*header, char*sequence, char*plus, char* quality);

#endif