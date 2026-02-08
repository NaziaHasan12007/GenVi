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

void print_warnings(){
    fprintf(stderr, "ERROR: Invalid Arguments.\n\n");
    fprintf(stderr, "USAGE: \n  ./bin/genvi <FASTA_Filename> <FASTQ File_name> <Sensor_file> <Emmission_file> <Starting time> <Ending time>");
    exit(1);
    
}

int main(int argc, char *argv[]){
    if(argc!=7){
        print_warnings();
    }

    char* fasta_filename= argv[1];
    char* fastq_filename= argv[2];
    char* sensor_file=argv[3];
    char* emmission_file=argv[4];
    char* start_time= argv[5];
    char* end_time= argv[6];

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
    GenomicStats* stats=process_and_map_reads(fastq_filename, genome_index, genome_string, genome_length);
    CarbonReport report = calculate_carbon_footprint(sensor_file, emmission_file,start_time, end_time);
    if(stats==NULL){
        fprintf(stderr, "[Phase 2] failed to read, processing and mapping");
        ht_free(genome_index);
        free(genome_string);
        return 1;
    }
   
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