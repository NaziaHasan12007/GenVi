#ifndef MAPPER_H
#define MAPPER_H

#include "basic.h"
#include "parser.h"
#include "hashtable.h"
#include "trimmer.h"
#include "trie.h"
#include "aligner.h"
#include "statistics.h"
#include "colors.h"
#include <pthread.h>

#define MAX_SEED_HITS 50      
#define REPETITIVE_SEED_THRESHOLD 20 
#define EXTENSION_MARGIN 1.1   
#define MAX_WINDOW_SIZE 2000
#define NUM_THREADS 6
#define BATCH_SIZE 5000 
#define RESET "\e[0m"

typedef struct {
    char sequence[MAX_LINE_LENGTH];
    char quality[MAX_LINE_LENGTH];
} Read;

typedef struct {
    const char* genome_string;
    long start;
    long end;
    HashTable* index;
} IndexWorkerData;

typedef struct {
    Read* read_batch;
    int count;
    HashTable* index;
    const char* genome_string;
    long genome_length;
    unsigned long* local_coverage; 
    long reads_discarded;
    long reads_failed_map;
} ThreadData;

HashTable* build_genome_index(const char* genome_string);

GenomicStats* process_and_map_reads(const char* fastq_filename, HashTable* index, const char* genome_string, long genome_length);

#endif