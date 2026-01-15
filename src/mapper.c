#include "mapper.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define MAX_SEED_HITS 50      
#define REPETITIVE_SEED_THRESHOLD 20 
#define EXTENSION_MARGIN 1.1   
#define MAX_WINDOW_SIZE 2000

HashTable* build_genome_index(const char* genome_string){
    printf("[Phase 1]Creating genome index...\n");
    HashTable* index=ht_create();  
    if(!index){
        fprintf(stderr, "[Phase 1]Failed to create hash table.\n");
        return NULL;
    }
    long genome_length=strlen(genome_string);
    const int k=SEED_KMER_SIZE;
    char* kmer_buffer=malloc(k+1);
    if(!kmer_buffer){
        fprintf(stderr, "[Phase 1] Failed to allocate kmer buffer.\n");
        ht_free(index);
        return NULL;
    }
    unsigned long repetitive_kmers = 0;
    for(long i=0; i<=genome_length-k; i++){
        strncpy(kmer_buffer, genome_string+i, k);
        kmer_buffer[k]='\0';
        Ht_Node* existing=ht_search(index, kmer_buffer);
        if(existing && existing->location_count>=REPETITIVE_SEED_THRESHOLD){
            repetitive_kmers++;  
            continue;
        }
        ht_insert(index, kmer_buffer, (unsigned long)i);
        if(i>0 && i%100000==0){
            printf("[Phase 1]Indexed %ld/%ld k-mers (%.2f%%)\r", i, genome_length - k, (double)i / (genome_length - k) * 100.0);
            fflush(stdout);
        }
    }
    printf("\n[Phase 1] Genome indexing complete. %lu repetitive seeds skipped.\n", repetitive_kmers);
    free(kmer_buffer);
    return index;
}

GenomicStats* process_and_map_reads(const char* fastq_filename, HashTable* index, const char* genome_string, long genome_length){
    FILE* fp=open_fastq_file(fastq_filename);
    if(!fp){
      return NULL;
    } 
    unsigned long* coverage_map=(unsigned long*)calloc(genome_length, sizeof(unsigned long));
    if(!coverage_map){
        fprintf(stderr, "[Mapper ERROR]Failed to allocate coverage_map\n");
        fclose(fp);
        return NULL;
    }

    char header[MAX_LINE_LENGTH];
    char sequence[MAX_LINE_LENGTH];
    char plus[MAX_LINE_LENGTH];
    char quality[MAX_LINE_LENGTH];
    const int k=SEED_KMER_SIZE;
    char* seed_buffer=malloc(k+1);

    if(!seed_buffer){
        fprintf(stderr, "[Mapper ERROR] Failed to allocate seed_buffer\n");
        free(coverage_map);
        fclose(fp);
        return NULL;
    }

    int* buf_pre=malloc((MAX_WINDOW_SIZE+1)*sizeof(int));
    int* buf_cur=malloc((MAX_WINDOW_SIZE+1)*sizeof(int));
    if(!buf_pre || !buf_cur){
        fprintf(stderr, "[Mapper ERROR] Failed to allocate DP buffers\n");
        free(seed_buffer);
        free(coverage_map);
        fclose(fp);
        return NULL;
    }

    long total_reads=0;
    long reads_discarded=0;
    long reads_failed_map=0;
    unsigned long total_quality_sum=0;
    unsigned long total_bases_sequenced=0;

    while(read_fastq_to_string(fp, header, sequence, plus, quality)){
        total_reads++;
        int read_len=strlen(sequence);

        for(int i=0; i<read_len; i++) {
            total_quality_sum+=((int)quality[i]-33);
        }
        total_bases_sequenced+=read_len;

        int trimmed_len=trim_read_by_quality(sequence, quality);
        if(trimmed_len<MIN_READ_LENGTH){
            reads_discarded++;
            continue;
        }

        int best_score=0;
        long best_loc=-1;
        int max_window_len=(int)(trimmed_len*EXTENSION_MARGIN);

        if(max_window_len>MAX_WINDOW_SIZE){
           max_window_len=MAX_WINDOW_SIZE;
        } 

        for(int i=0; i<=trimmed_len-k; i+=k){
            strncpy(seed_buffer, sequence+i, k);
            seed_buffer[k]='\0';
            Ht_Node* hits=ht_search(index, seed_buffer);
            if(!hits){
             continue;
            } 
            if(hits->location_count > REPETITIVE_SEED_THRESHOLD || hits->location_count > MAX_SEED_HITS){
              continue;
            } 

            for(size_t h= 0; h< hits->location_count; h++) {
                long loc=hits->locations[h];
                long genome_start=loc - i;

                if(genome_start<0||genome_start+max_window_len>genome_length){
                   continue;
                } 

                int score=run_dp_sw(sequence, genome_string + genome_start, trimmed_len, max_window_len, buf_pre, buf_cur);
                if(score>best_score){
                    best_score=score;
                    best_loc=genome_start;
                }
            }
        }
        const int ALIGNMENT_SCORE_THRESHOLD=(int)(trimmed_len * MATCH_SCORE * 0.7);
        if(best_score>=ALIGNMENT_SCORE_THRESHOLD && best_loc>=0){
            for(int j=0; j<trimmed_len; j++){
                if(best_loc + j < genome_length)
                    coverage_map[best_loc+j]++;
            }
        } 
        else{
            reads_failed_map++;
        }
        if (total_reads%10000==0) {
            printf("[Mapper] Processed %ld reads... \r", total_reads);
            fflush(stdout);
        }
    }
    fclose(fp);
    free(buf_pre);
    free(buf_cur);
    free(seed_buffer);
    printf("\n[Mapper] Completed processing %ld reads.\n", total_reads);
    printf("[Mapper] Reads discarded: %ld, failed to map: %ld\n", reads_discarded, reads_failed_map);
    GenomicStats* stats=calculate_genomic_stats(coverage_map, genome_length, total_reads, reads_discarded, reads_failed_map, total_quality_sum, total_bases_sequenced, genome_string);
    free(coverage_map);
    return stats;
}