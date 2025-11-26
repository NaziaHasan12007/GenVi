#include "mapper.h"

HashTable* build_genome_index(const char* genome_string){
    printf("[Phase 1]Creating new hash table index......\n");
    HashTable* index=ht_create();
    if(index==NULL){
        fprintf(stderr, "[Phase 1]Hashtable build failed\n");
        return NULL;
    }
    long genome_length=strlen(genome_string); 
    const int k=SEED_KMER_SIZE;
    char kmer_buffer[k + 1];
    printf("[Phase 1]Genome Length: %ld bases\n", genome_length);
    printf("[Phase 1]Indexing %ld k-mers (k=%d). This may take a moment...\n", genome_length-(k-1), k);
    for(long i=0; i<=genome_length-k; i++){
        strncpy(kmer_buffer, genome_string+i, k);
        kmer_buffer[k]='\0';
        ht_insert(index, kmer_buffer, (unsigned long)i);
        if(i>0 && i%100000==0){
            printf("[Phase 1]Progress:%.2f%% \r", (double)i/genome_length*100.0);
            fflush(stdout);
        }
    }
    printf("\n[Phase 1]Indexing complete.\n");
    return index;

}

GenomicStats* process_and_map_reads(const char* fastq_filename, HashTable* index, const char* genome_string, long genome_length){
    printf("  [Phase 2/3] Starting read processing and mapping...\n");
    FILE* fp = open_fastq_file(fastq_filename);
    if(fp==NULL){ return NULL; }
    
    unsigned long* coverage_map=(unsigned long*)calloc(genome_length, sizeof(unsigned long));
    if(coverage_map==NULL){
        fprintf(stderr,"[Mapper ERROR] Failed to allocate coverage map.\n");
        fclose(fp);
        return NULL;
    }
    TrieNode* trie=trie_create();
    if(trie==NULL){
        fprintf(stderr,"[Mapper ERROR] Failed to create Trie.\n");
        fclose(fp);
        free(coverage_map);
        return NULL;
    }
    char header[MAX_LINE_LENGTH];
    char sequence[MAX_LINE_LENGTH];
    char plus[MAX_LINE_LENGTH];
    char quality[MAX_LINE_LENGTH];
    
    long total_reads=0;
    long reads_discarded=0;
    long reads_failed_map=0;
    unsigned long total_quality_sum=0;
    unsigned long total_bases_sequenced=0;
    
    const int k=SEED_KMER_SIZE;
    char kmer_buffer[k+1];
    char seed_buffer[k+1];
    const int ALIGNMENT_SCORE_THRESHOLD=(k*MATCH_SCORE)/2;
 
    while(read_fastq_to_string(fp, header, sequence, plus, quality)){
        total_reads++;
        int raw_read_length=strlen(quality);
        for(int i=0; i<raw_read_length; i++) {
            total_quality_sum+=((int)quality[i]-33);
        }
        total_bases_sequenced+=raw_read_length;

        int trimmed_length=trim_read_by_quality(sequence, quality);
        if(trimmed_length<MIN_READ_LENGTH){
            reads_discarded++;
            continue;
        }
        if(trimmed_length>=k){
            for(int i=0; i<=trimmed_length-k; i++){
                strncpy(kmer_buffer, sequence+i, k);
                kmer_buffer[k]='\0';
                trie_insert(trie, kmer_buffer);
            }
            analyze_and_correct_read(trie, sequence);
        }
        char* clean_read_sequence=sequence;
        if(trimmed_length<k){
            reads_failed_map++;
            continue;
        }
        strncpy(seed_buffer, clean_read_sequence, k);
        seed_buffer[k]='\0';
        Ht_Node* hits=ht_search(index, seed_buffer);
        if(hits==NULL){
            reads_failed_map++;
            continue;
        }
        int alignment_found=0;
        for(size_t i=0; i<hits->location_count; i++){
            unsigned long location=hits->locations[i];
            if(location+trimmed_length > genome_length){
                continue; 
            }
            const char* genome_segment=genome_string+location;
            int score=run_smith_waterman(clean_read_sequence, genome_segment);
            if(score>ALIGNMENT_SCORE_THRESHOLD){
                alignment_found=1;
                for(int read_pos=0; read_pos<trimmed_length; read_pos++){
                    coverage_map[location+read_pos]++;
                }
                break;
            }
        }
        if(!alignment_found){ 
            reads_failed_map++; 
        }
        if(total_reads % 10000==0){
            printf("[Phase 2/3] Processed %ld reads... \r", total_reads);
            fflush(stdout);
        }
    }
    
    fclose(fp);
    printf("\n[Phase 2/3]All %ld reads processed. Calculating final statistics...\n", total_reads);
    trie_free(trie);
    GenomicStats* stats=calculate_genomic_stats(coverage_map, genome_length, total_reads, reads_discarded, reads_failed_map, total_quality_sum, total_bases_sequenced, genome_string);
    free(coverage_map);
    printf("[Phase 2/3]Mapping and statistics complete.\n");
    return stats;
}