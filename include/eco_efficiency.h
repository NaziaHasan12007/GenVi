#ifndef ECO_EFFICIENCY_H
#define ECO_EFFICIENCY_H
#include "statistics.h" 
#include "recommendation.h"
#include <time.h>
#ifndef MAX_LINE_LENGTH
#define MAX_LINE_LENGTH 2048
#endif

typedef struct {
    double total_kwh;                 
    double total_carbon_kg;            
    double average_temp_celsius;       
    int    sensor_readings_count;     
    int    skipped_rows;  
    double carbon_intensity;            
} CarbonReport;
CarbonReport calculate_carbon_footprint( const char *sensor_file, const char *emmission_file, const char *start_time, const char *end_time);

void print_final_report(const GenomicStats *stats, const CarbonReport *report, const Recommendation *rec);

#endif