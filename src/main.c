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

    printf("Wecome to GenVi: Genomic Environmental Vitality Index\n");
    printf("=====================================================\n\n");
    printf("Starting GenVi Analysis\n");
    printf("------------------------\n");
    printf("[Phase 1] Starting: Building Reference Index from the FASTA file %s\n\n", fasta_filename);

    long genome_length=0;
    char* genome_string=NULL;
    genome_string=load_genome_from_fasta(fasta_filename, &genome_length);
    if(genome_string==NULL){
        fprintf(stderr, "[Phase 1] Failed to load genome from FASTA.\n");
        return 1;
    }
 
    HashTable* genome_index=build_genome_index(genome_string);
    if(genome_index==NULL){
        fprintf(stderr, "[Phase 1] Failed to build genome index from FASTA.\n");
        return 1;
    }
    printf("[Phase 1] Completed: Genome Index Built Successful\n\n");
    printf("[Phase 2&3] Starting: Processing FASTQ from the file and mapping the reads %s\n\n", fastq_filename);
    GenomicStats* statistics=process_and_map_reads(fastq_filename, genome_index, genome_string, genome_length);
    if(statistics==NULL){
        fprintf(stderr, "[Phase 2] failed to read, processing and mapping");
        ht_free(genome_index);
        free(genome_string);
        return 1;
    }
    long clean_reads=statistics->total_reads_aligned;
    float clean_read_percentage=(float)(((double)clean_reads/statistics->total_reads_processed)*100);
    printf("[Phase 2&3] Completed: FASTQ Processing Completed. Found %f%% Valid Reads\n\n", clean_read_percentage);
    printf("Genomic Statistics\n");
    printf("========================\n");
    printf("Genome Length: %ld\n", statistics->genome_length);
    printf("Total Reads Processed: %ld\n", statistics->total_reads_processed);
    printf("Total Reads Aligned: %ld\n", statistics->total_reads_aligned);
    printf("Reads Discarded or Failed to Map: %ld\n", statistics->total_reads_trimmed + statistics->total_reads_failed_to_map);
    printf("Percent Reads Discarded: %.2f%%\n", statistics->percent_reads_discarded);
    printf("Mean Coverage: %.2f\n", statistics->mean_coverage);
    printf("Median Coverage: %.2f\n", statistics->median_coverage);
    printf("Coverage Standard Deviation: %.2f\n", statistics->std_dev_coverage);
    printf("Coverage Breadth: %.2f%%\n", statistics->genome_breadth_percent);
    printf("GC Content: %.2f%%\n", statistics->gc_content_percent);
    printf("N50 metric: %ld\n", statistics->n50_metric);
    printf("Coverage Quartile Ranges\n-------------------------\n");
    printf("Q1: %.2f\n", statistics->coverage_q1);
    printf("Q3: %.2f\n", statistics->coverage_q3);
    printf("IQR(Inter Quartile Range): %.2f\n\n\n", statistics->coverage_iqr);
    printf("Done for now!! Thank you for using GENVI!");
    
}