#include "mapper.h"

void* indexing_worker(void* arg) {
    IndexWorkerData* data=(IndexWorkerData*)arg;
    int k=SEED_KMER_SIZE;
    char kmer[SEED_KMER_SIZE+1];

    for (long i=data->start; i<=data->end; i++){
       
        strncpy(kmer, data->genome_string + i, k);
        kmer[k]='\0';
        
        ht_insert(data->index, kmer, (unsigned long)i);
        
        if (i>0&&i%200000==0) {
            printf(CYN "[Phase 1] Indexing progress: %ld bases\r" RESET, i);
            fflush(stdout);
        }
    }
    return NULL;
}

HashTable* build_genome_index(const char* genome_string) {
    printf(CYN "[Phase 1] Initializing multithreaded indexing...\n" RESET);
    HashTable* index=ht_create();
    if (!index) return NULL;

    long genome_length=strlen(genome_string);
    pthread_t threads[NUM_THREADS];
    IndexWorkerData tdata[NUM_THREADS];

    long total_kmers=genome_length-SEED_KMER_SIZE;
    long chunk=total_kmers/NUM_THREADS;

    for(int t=0; t<NUM_THREADS; t++){
        tdata[t].genome_string=genome_string;
        tdata[t].index=index;
        tdata[t].start=t*chunk;
        tdata[t].end=(t==NUM_THREADS-1)?total_kmers:(t+1)*chunk-1;

        pthread_create(&threads[t], NULL, indexing_worker, &tdata[t]);
    }

    for(int t=0; t<NUM_THREADS; t++){
        pthread_join(threads[t], NULL);
    }
    printf(GRN "\n[Phase 1] Indexing complete.\n" RESET);
    return index;
}

void* mapper_worker(void* arg) {
    ThreadData* data=(ThreadData*)arg;
    const int k=SEED_KMER_SIZE;

    char* seed_buffer=malloc(k+1);
    int* buf_pre=malloc((MAX_WINDOW_SIZE+1)*sizeof(int));
    int* buf_cur=malloc((MAX_WINDOW_SIZE+1)*sizeof(int));

    for(int r=0; r<data->count; r++){
        const char* sequence=data->read_batch[r].sequence;
        const char* quality=data->read_batch[r].quality;

        int trimmed_len=trim_read_by_quality((char*)sequence, (char*)quality);
        if(trimmed_len<MIN_READ_LENGTH){
            data->reads_discarded++;
            continue;
        }

        int best_score=0;
        long best_loc=-1;
        int max_window_len=(int)(trimmed_len*EXTENSION_MARGIN);
        if(max_window_len>MAX_WINDOW_SIZE)max_window_len=MAX_WINDOW_SIZE;

        for(int i=0; i<=trimmed_len-k; i+=k) {
            strncpy(seed_buffer, sequence+i, k);
            seed_buffer[k]='\0';

            Ht_Node* hits=ht_search(data->index, seed_buffer);
            if(!hits||hits->location_count>REPETITIVE_SEED_THRESHOLD) continue;

            for(size_t h = 0; h < hits->location_count && h < MAX_SEED_HITS; h++) {
                long loc = hits->locations[h];
                long genome_start = loc - i;

                if(genome_start < 0 || genome_start + max_window_len > data->genome_length)
                    continue;

                int score = run_dp_sw(sequence, data->genome_string + genome_start, trimmed_len, max_window_len, buf_pre, buf_cur);
                if(score > best_score) {
                    best_score = score;
                    best_loc = genome_start;
                }
                
                if(best_score > (trimmed_len * MATCH_SCORE * 0.95)) break;
            }
        }

        const int threshold = (int)(trimmed_len * MATCH_SCORE * 0.7);
        if(best_score >= threshold && best_loc >= 0) {
            for(int j = 0; j < trimmed_len; j++) {
                if(best_loc + j < data->genome_length) {
                    data->local_coverage[best_loc + j]++;
                }
            }
        } else {
            data->reads_failed_map++;
        }
    }

    free(seed_buffer);
    free(buf_pre);
    free(buf_cur);
    return NULL;
}

GenomicStats* process_and_map_reads(const char* fastq_filename, HashTable* index, const char* genome_string, long genome_length) {

    FILE* fp=open_fastq_file(fastq_filename);
    if(!fp){ 
        return NULL;
    }
    unsigned long* global_coverage = calloc(genome_length, sizeof(unsigned long));
    Read* batch=malloc(sizeof(Read) * BATCH_SIZE);
    
    long total_processed=0, total_discarded=0, total_failed=0;

    printf(CYN "[Phase 2&3] Starting multithreaded mapping...\n" RESET);

    while(1) {
        int loaded=0;
        char h[MAX_LINE_LENGTH], p[MAX_LINE_LENGTH], s[MAX_LINE_LENGTH], q[MAX_LINE_LENGTH];
        
        while(loaded<BATCH_SIZE&&read_fastq_to_string(fp, h, s, p, q)){
            strncpy(batch[loaded].sequence, s, MAX_LINE_LENGTH);
            strncpy(batch[loaded].quality, q, MAX_LINE_LENGTH);
            loaded++;
        }

        if(loaded==0) break;

        pthread_t threads[NUM_THREADS];
        ThreadData tdata[NUM_THREADS];

        int base_chunk=loaded/NUM_THREADS;
        int remainder=loaded%NUM_THREADS;
        int current_start=0;

        for(int t=0; t<NUM_THREADS; t++){
            int current_chunk=base_chunk+(t<remainder?1:0);
            
            tdata[t].read_batch=&batch[current_start];
            tdata[t].count=current_chunk;
            tdata[t].index=index;
            tdata[t].genome_string=genome_string;
            tdata[t].genome_length=genome_length;
            tdata[t].local_coverage=calloc(genome_length, sizeof(unsigned long));
            tdata[t].reads_discarded=0;
            tdata[t].reads_failed_map=0;

            pthread_create(&threads[t], NULL, mapper_worker, &tdata[t]);
            current_start += current_chunk;
        }

        for(int t=0; t<NUM_THREADS; t++){
            pthread_join(threads[t], NULL);
            total_discarded+=tdata[t].reads_discarded;
            total_failed+=tdata[t].reads_failed_map;

            for(long i=0; i<genome_length; i++){
                global_coverage[i]+=tdata[t].local_coverage[i];
            }
            free(tdata[t].local_coverage);
        }

        total_processed+=loaded;
        printf(CYN "[Mapper] Processed %ld reads...\r" RESET, total_processed);
        fflush(stdout);
    }

    GenomicStats* stats=calculate_genomic_stats(global_coverage, genome_length, total_processed, total_discarded, total_failed, 0, 0, genome_string);

    fclose(fp);
    free(batch);
    free(global_coverage);
    return stats;
}