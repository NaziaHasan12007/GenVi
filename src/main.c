#include "basic.h"
#include "parser.h"
#include "trimmer.h"
#include "trie.h"
#include "hashtable.h"
#include "aligner.h"
#include "mapper.h"
#include "statistics.h"
#include "eco_efficiency.h"

void print_warnings(){
    fprintf(stderr, "ERROR: Invalid Arguments.\n\n");
    fprintf(stderr, "USAGE: \n  ./bin/genvi <FASTA_Filename> <FASTQ File_name> <Starting time> <Ending time>");
    exit(1);
    
}

int main(int argc, char *argv[]){
    if(argc!=3){
        print_warnings();
    }

    char* fasta_filename= argv[1];
    char* fastq_filename= argv[2];
   // char* start_time= argv[3];
    //char* end_time= argv[4];

    printf("Wecome to GenVi: Genomic Eco-Efficiency Analyzer\n");
    printf("=====================================================\n\n");
    printf("Starting GenVi Analysis\n");
    printf("------------------------");
    printf("[Phase 1] Starting: Building Reference Index from the FASTA file %s\n\n", fasta_filename);

    long genome_length=0;
    char* genome_string=NULL;
    genome_string=load_genome_from_fasta(fasta_filename, &genome_length);
    if(genome_string==NULL){
        fprintf(stderr, "[Phase 1] failed to load genome from FASTA.\n");
        return 1;
    }
 
    HashTable* genome_index=build_genome_index(genome_string);
    if(genome_index==NULL){
        fprintf(stderr, "[Phase 1] failed to build genome index from FASTA.\n");
        return 1;
    }
    printf("[Phase 1] Completed: Genome Index Built Successful\n\n");
    printf("[Phase 2&3] Starting: Processing FASTQ from the file and mapping the reads%s\n\n", fastq_filename);
    GenomicStats* statistics=process_and_map_reads(fastq_filename, genome_index, genome_string, genome_length);
    if(statistics==NULL){
        fprintf(stderr, "[Phase 2] failed to read, processing and mapping");
        ht_free(genome_index);
        free(genome_string);
        return 1;
    }
    long clean_reads=statistics->total_reads_aligned;
    float clean_read_percentage=(float)(((double)clean_reads/statistics->total_reads_processed)*100);
    //Have to call the Fastq Parser and generate k-mer and trie-builder
    printf("[Phase 2&3] Completed: FASTQ Processing Completed. Found %f%% Valid Reads\n\n", clean_read_percentage);
    
    
    //count Average Coverage (Mean)
    //  - Median Coverage
    //  - Coverage Standard Deviation (how "even" is the coverage?)
    //  - Coverage Breadth (% of genome covered at least 1x)
    //  - A Coverage Histogram (how much of the genome is at 0x, 1-5x, etc.)
    /*printf("[Phase 3] Completed: Genome Mapping and Statistics Completed\n\n");
    printf("[Phase 4] Starting: Calculating Carbon Foootprint\n\n");
    //have to eco-efficiency call
    printf("[Phase 4] Completed: carbon Analysis Cimpleted\n\n");
    printf("[Phase 5] Starting: Generating The Final report\n\n");
    //have to call the final report
    printf("--------------------------------\n\n");
    printf("GenVi Analysis Completed\n");
    printf("Cleaning The Memory.....\n");
    printf("Done! Thank You");
    return 0;*/
}