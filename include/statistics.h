#ifndef STATISTICS_H
#define STATISTICS_H
#include "basic.h"

typedef struct{
    long total_reads_processed;
    long total_reads_trimmed;
    long total_reads_failed_to_map;
    long total_reads_aligned;
    double percent_reads_discarded;
    double avg_read_quality;
    long genome_length;
    double gc_content_percent;
    double mean_coverage;
    double median_coverage;
    double std_dev_coverage;
    double coverage_q1;      
    double coverage_q3;      
    double coverage_iqr;
    long n50_metric;
    double genome_breadth_percent;
    double histogram[HISTOGRAM_BINS];
}GenomicStats;

GenomicStats* calculate_genomic_stats(unsigned long* coverage_map, long genome_length, long total_reads, long reads_discarded, long reads_failed_to_map, unsigned long total_quality_sum, unsigned long total_bases_sequenced, const char* genome_string);
void merge_sort(double* arr, long l, long r);
void merge(double* arr, long l, long m, long r);
void merge_sort_long(long* arr, long l, long r);
void merge_long(long* arr, long l, long m, long r);

#endif