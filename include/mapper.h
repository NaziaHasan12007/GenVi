#ifndef MAPPER_H
#define MAPPER_H

#include "basic.h"
#include "parser.h"
#include "hashtable.h"
#include "trimmer.h"
#include "trie.h"
#include "aligner.h"
#include "statistics.h"

HashTable* build_genome_index(const char* genome_string);
GenomicStats* process_and_map_reads(const char* fastq_filename, HashTable* index, const char* genome_string, long genome_length);
#endif