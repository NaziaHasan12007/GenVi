#include "basic.h"
#include "parser.h"
#include "trimmer.h"
#include "trie.h"
#include "hashtable.h"
#include "aligner.h"
#include "mapper.h"
#include "statistics.h"
#include "eco_efficiency.h"
#include "recommendation.h"
#include "colors.h"

#define RESET "\e[0m"

void print_warnings(){
    fprintf(stderr, BRED "ERROR: Invalid Arguments. " RESET "\n\n");
    fprintf(stderr, YEL"USAGE: \n  ./bin/genvi <FASTA_Filename> <FASTQ File_name> <Sensor_file> <Emmission_file> <Starting time> <Ending time>" RESET);
    exit(1);
    
}

int main(int argc, char *argv[]){
    if(argc!=7){
        print_warnings();
    }
    char* fasta_filename=argv[1];
    char* fastq_filename=argv[2];
    char* sensor_file=argv[3];
    char* emmission_file=argv[4];
    char* start_time=argv[5];
    char* end_time=argv[6];

    printf(BBLU "Wecome to GenVi: Genomic Environmental Vitality Index\n" RESET);
    printf(BBLU"=====================================================\n\n"RESET);
    printf(BWHT"Starting GenVi Analysis\n"RESET);
    printf("------------------------\n");
    printf(CYN"[Phase 1] Starting: Building Reference Index from the FASTA file %s\n\n"RESET, fasta_filename);

    long genome_length=0;
    char* genome_string=NULL;
    genome_string=load_genome_from_fasta(fasta_filename, &genome_length);
    if(genome_string==NULL){
        fprintf(stderr, RED"[Phase 1] Failed to load genome from FASTA.\n"RESET);
        return 1;
    }
 
    HashTable* genome_index=build_genome_index(genome_string);
    if(genome_index==NULL){
        fprintf(stderr, RED"[Phase 1] Failed to build genome index from FASTA.\n"RESET);
        return 1;
    }
    printf(GRN"[Phase 1] Completed: Genome Index Built Successful\n\n"RESET);
    printf(CYN"[Phase 2&3] Starting: Processing FASTQ from the file and mapping the reads %s\n\n"RESET, fastq_filename);
    GenomicStats* stats=process_and_map_reads(fastq_filename, genome_index, genome_string, genome_length);
    CarbonReport report = calculate_carbon_footprint(sensor_file, emmission_file,start_time, end_time);
    if(stats==NULL){
        fprintf(stderr, RED"[Phase 2] failed to read, processing and mapping"RESET);
        ht_free(genome_index);
        free(genome_string);
        return 1;
    }
    printf(GRN "[Phase 2–3] Completed: Reads mapped successfully\n\n" RESET);

    ParameterNode node;
    node.timestamp=time(NULL);
    node.forecast_carbon_intensity=report.carbon_intensity;
    node.forecast_temparature_c=report.average_temp_celsius;
    node.urgency=0;  
    Weights prefs={.carbon_weight  = 0.6, .cooling_weight = 0.3, .urgency_weight = 0.1};
    Recommendation rec;
    get_parameter_recommendation(&node, 1, prefs, &rec);

    print_final_report(stats, &report, &rec);
    ht_free(genome_index);
    free(genome_string);
    free(stats);
}