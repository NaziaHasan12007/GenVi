#ifndef ECO_EFFICIENCY_H
#define ECO_EFFICIENCY_H
#include <time.h>
#ifndef MAX_LINE_LENGTH
#define MAX_LINE_LENGTH 1024
#endif
typedef struct GenomicStats GenomicStats;

typedef struct {
    double total_kwh;                 
    double total_carbon_kg;            
    double average_temp_celsius;       
    int    sensor_readings_count;     
    int    skipped_rows;              
} CarbonReport;
CarbonReport calculate_carbon_footprint(
    const char *sensor_file,
    const char *emission_file,
    const char *start_time,
    const char *end_time
);

void print_final_report(
    const GenomicStats *stats,
    const CarbonReport *report
);

#endif