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
    if(argc!=5){
        print_warnings();
    }

    char* fasta_filename= argv[1];
    char* fastq_filename= argv[2];
    char* start_time= argv[3];
    char* end_time= argv[4];

    printf("Wecome to GenVi: Genomic Eco-Efficiency Analyzer\n");
    printf("=====================================================\n\n");
    printf("Starting GenVi Analysis\n");
    printf("------------------------");
    printf("[Phase 1] Starting: Building Reference Index from the FASTA file %s\n\n", fasta_filename);

    //Have to call the Fasta Parser, Hash table builder
    printf("[Phase 1] Completed: Genome Index Built Successful\n\n");
    printf("[Phase 2] Starting: Processing FASTQ from the file %s\n\n", fastq_filename);
    //Have to call the Fastq Parser and generate k-mer and trie-builder
    printf("[Phase 2] Completed: FASTQ Processing Completed. Found %d%% Valid Reads\n\n" /*clean_read counts*/);
    printf("[Phase 3] Starting: Map Clean Reading to Genome Index\n\n");
    //have to call related functions
    printf("[Phase 3]Reading Maps Completed. Now Calculating Statistics.....\n\n");
    //count:Average Coverage (Mean)
    //  - Median Coverage
    //  - Coverage Standard Deviation (how "even" is the coverage?)
    //  - Coverage Breadth (% of genome covered at least 1x)
    //  - A Coverage Histogram (how much of the genome is at 0x, 1-5x, etc.)
    printf("[Phase 3] Completed: Genome Mapping and Statistics Completed\n\n");
    printf("[Phase 4] Starting: Calculating Carbon Foootprint\n\n");
    //eco-efficiency call
    printf("[Phase 4] Completed: carbon Analysis Cimpleted\n\n");
    printf("[Phase 5] Starting: Generating The Final report\n\n");
    //call the final report
    printf("--------------------------------\n\n");
    printf("GenVi Analysis Completed\n");
    printf("Cleaning The Memory.....\n");
    printf("Done! Thank You");
    return 0;
}