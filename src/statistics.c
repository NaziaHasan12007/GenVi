#include "statistics.h"
#include <math.h>

static double calculate_mean(unsigned long*coverage_map,long genome_length,unsigned long*total_depth_sum){
    if(genome_length==0){
        return 0.0;
    }
    unsigned long sum=0;
    for(long i=0;i<genome_length;i++){
        sum+=coverage_map[i];
    }
    *total_depth_sum=sum;
    return (double)sum/genome_length;
}

static double calculate_std_dev(unsigned long*coverage_map,long genome_length,double mean){
    if(genome_length==0){
        return 0.0;
    }
    double variance_sum=0;
    for(long i=0;i<genome_length;i++){
        double diff=(double)coverage_map[i]-mean;
        variance_sum+=diff*diff;
    }
    return sqrt(variance_sum/genome_length);
}

static double calculate_breadth(unsigned long*coverage_map,long genome_length){
    if(genome_length==0){
        return 0.0;
    }
    long covered_bases=0;
    for(long i=0;i<genome_length;i++){
        if(coverage_map[i]>0){
            covered_bases++;
        }
    }
    return (double)covered_bases/genome_length*100.0;
}

static double calculate_gc_content(const char*genome_string,long genome_length){
    if(genome_length==0){
        return 0.0;
    }
    long gc_count=0;
    for(long i=0;i<genome_length;i++){
        char c=genome_string[i];
        if(c=='G'||c=='C'||c=='g'||c=='c'){
            gc_count++;
        }
    }
    return (double)gc_count/genome_length*100.0;
}

static void calculate_histogram(unsigned long *coverage_map,long genome_length,double *histogram_out){
    if(genome_length==0){
        return;
    }
    long histogram_counts[HISTOGRAM_BINS]={0};
    long max_bin_index=HISTOGRAM_BINS-1;
    for(long i=0;i<genome_length;i++){
        unsigned long coverage=coverage_map[i];
        int bin_index;

        if(coverage>=HISTOGRAM_BINS){
            bin_index=max_bin_index;
        }
        else{
            bin_index=(int)coverage;
        }
        histogram_counts[bin_index]++;
    }
    printf("\nCoverage Histogram:\n");
    printf("-------------------\n");
    for(int i=0;i<HISTOGRAM_BINS;i++){
        double percentage=(double)histogram_counts[i]/genome_length*100.0;
        histogram_out[i]=percentage;
        if(i==max_bin_index){
            printf("%2d+  | ",i);
        }
        else{
            printf("%2d   | ",i);
        }
        int bar_len=histogram_counts[i]/(genome_length/50+1);
        for(int j=0;j<bar_len;j++){
            putchar('#');
        }
        printf(" (%ld, %.2f%%)\n",histogram_counts[i],percentage);
    }
}

void merge(double*arr,long l,long m,long r){
    long n1=m-l+1,n2=r-m;
    double*L=malloc(n1*sizeof(double));
    double*R=malloc(n2*sizeof(double));
    if(L==NULL||R==NULL){
        free(L);
        free(R);
        return;
    }
    for(long i=0;i<n1;i++){
        L[i]=arr[l+i];
    }
    for(long j=0;j<n2;j++){
        R[j]=arr[m+1+j];
    }
    long i=0,j=0,k=l;
    while(i<n1&&j<n2){
        if(L[i]<=R[j]){
            arr[k++]=L[i++];
        }else{
            arr[k++]=R[j++];
        }
    }
    while(i<n1){
        arr[k++]=L[i++];
    }
    while(j<n2){
        arr[k++]=R[j++];
    }
    free(L);
    free(R);
}

void merge_sort(double*arr,long l,long r){
    if(l<r){
        long m=l+(r-l)/2;
        merge_sort(arr,l,m);
        merge_sort(arr,m+1,r);
        merge(arr,l,m,r);
    }
}

void merge_long(long*arr,long l,long m,long r){
    long n1=m-l+1,n2=r-m;
    long*L=malloc(n1*sizeof(long));
    long*R=malloc(n2*sizeof(long));
    if(L==NULL||R==NULL){
        free(L);
        free(R);
        return;
    }
    for(long i=0;i<n1;i++){
        L[i]=arr[l+i];
    }
    for(long j=0;j<n2;j++){
        R[j]=arr[m+1+j];
    }
    long i=0,j=0,k=l;
    while(i<n1&&j<n2){
        if(L[i]<=R[j]){
            arr[k++]=L[i++];
        }else{
            arr[k++]=R[j++];
        }
    }
    while(i<n1){
        arr[k++]=L[i++];
    }
    while(j<n2){
        arr[k++]=R[j++];
    }
    free(L);
    free(R);
}

void merge_sort_long(long*arr,long l,long r){
    if(l<r){
        long m=l+(r-l)/2;
        merge_sort_long(arr,l,m);
        merge_sort_long(arr,m+1,r);
        merge_long(arr,l,m,r);
    }
}

static void compute_basic_metrics(GenomicStats*stats,unsigned long*coverage_map,const char*genome_string){
    unsigned long total_depth_sum=0;
    stats->gc_content_percent=calculate_gc_content(genome_string,stats->genome_length);
    stats->mean_coverage=calculate_mean(coverage_map,stats->genome_length,&total_depth_sum);
    stats->std_dev_coverage=calculate_std_dev(coverage_map,stats->genome_length,stats->mean_coverage);
    stats->genome_breadth_percent=calculate_breadth(coverage_map,stats->genome_length);
    calculate_histogram(coverage_map,stats->genome_length,stats->histogram);
}

static void compute_quartiles_and_median(GenomicStats*stats,unsigned long*coverage_map){
    long len=stats->genome_length;
    if(len==0){
        return;
    }
    double*temp_array=malloc(len*sizeof(double));
    if(temp_array==NULL){
        return;
    }
    for(long i=0;i<len;i++){
        temp_array[i]=(double)coverage_map[i];
    }
    merge_sort(temp_array,0,len-1);
    long idx_q1=len/4,idx_med=len/2,idx_q3=len*3/4;
    stats->coverage_q1=temp_array[idx_q1];
    stats->coverage_q3=temp_array[idx_q3];
    stats->coverage_iqr=stats->coverage_q3-stats->coverage_q1;
    if(len%2==0){
        stats->median_coverage=(temp_array[idx_med-1]+temp_array[idx_med])/2.0;
    }else{
        stats->median_coverage=temp_array[idx_med];
    }
    free(temp_array);
}

static void compute_n50_metric(GenomicStats*stats,unsigned long*coverage_map){
    long len=stats->genome_length;
    long*contig_lengths=malloc((len/2+1)*sizeof(long));
    if(contig_lengths==NULL){
        return;
    }
    long contig_count=0,current_length=0,total_covered_bases=0;
    for(long i=0;i<len;i++){
        if(coverage_map[i]>0){
            current_length++;
        }else{
            if(current_length>0){
                contig_lengths[contig_count++]=current_length;
                total_covered_bases+=current_length;
                current_length=0;
            }
        }
    }
    if(current_length>0){
        contig_lengths[contig_count++]=current_length;
        total_covered_bases+=current_length;
    }
    if(total_covered_bases==0){
        free(contig_lengths);
        stats->n50_metric=0;
        return;
    }
    merge_sort_long(contig_lengths,0,contig_count-1);
    long threshold=total_covered_bases/2,running_sum=0;
    for(long i=contig_count-1;i>=0;i--){
        running_sum+=contig_lengths[i];
        if(running_sum>=threshold){
            stats->n50_metric=contig_lengths[i];
            break;
        }
    }
    free(contig_lengths);
}

GenomicStats* calculate_genomic_stats(unsigned long*coverage_map,long genome_length,long total_reads,long reads_discarded,long reads_failed_map,unsigned long total_quality_sum,unsigned long total_bases_sequenced,const char*genome_string){
    GenomicStats*stats=calloc(1,sizeof(GenomicStats));
    if(stats==NULL){
        fprintf(stderr,"[Statistics ERROR] Failed to allocate memory for stats report.\n");
        return NULL;
    }
    stats->genome_length=genome_length;
    stats->total_reads_processed=total_reads;
    stats->total_reads_trimmed=reads_discarded;
    stats->total_reads_failed_to_map=reads_failed_map;
    stats->total_reads_aligned=total_reads-reads_discarded-reads_failed_map;
    if(total_reads>0){
        stats->percent_reads_discarded=(double)(reads_discarded+reads_failed_map)/total_reads*100.0;
    }
    if(total_bases_sequenced>0){
        stats->avg_read_quality=(double)total_quality_sum/total_bases_sequenced;
    }else{
        stats->avg_read_quality=0.0;
    }
    compute_basic_metrics(stats,coverage_map,genome_string);
    compute_quartiles_and_median(stats,coverage_map);
    compute_n50_metric(stats,coverage_map);
    return stats;
}