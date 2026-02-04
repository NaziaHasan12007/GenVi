#include "eco_efficiency.h"
#include "statistics.h"
#include <assert.h>

#define MAX_INTERVALS 128

typedef struct {
    time_t start;
    time_t end;
    double kg_co2_per_kwh;
}EmissionInterval;

static time_t parse_timestamp_utc(const char *s);
static int load_emission_profile(const char *file, EmissionInterval *arr, int max);
static double lookup_emission_factor(time_t t, const EmissionInterval *arr, int count);

static time_t parse_timestamp_utc(const char *s){
    struct tm tm={0};
    if (!strptime(s, "%Y-%m-%dT%H:%M:%S", &tm)){
       return (time_t)-1;
    }
   tm.tm_isdst=-1;

   #if defined(__unix__)||defined(__APPLE__)
       return timegm(&tm);   
   #else
        return mktime(&tm);  
   #endif
}

static int load_emission_profile(const char *file, EmissionInterval *arr,int max){
    FILE *fp=fopen(file, "r");
    if(!fp){
      return -1;
    } 
    char line[MAX_LINE_LENGTH];
    int count=0;
    fgets(line, sizeof(line), fp); 
    while(fgets(line, sizeof(line), fp) && count < max){
        line[strcspn(line, "\r\n")]='\0';
        char *saveptr;
        char *tok=strtok_r(line, ",", &saveptr);
        if(!tok){
          continue;
        }
        time_t start=parse_timestamp_utc(tok);
        tok=strtok_r(NULL, ",", &saveptr);
        if(!tok){
          continue;
        } 
        time_t end=parse_timestamp_utc(tok);
        tok=strtok_r(NULL, ",", &saveptr);
        if(!tok){
          continue;
        } 
        char *endptr;
        double factor=strtod(tok, &endptr);
        if(*endptr){
          continue;
        } 
        arr[count++]=(EmissionInterval){ start, end, factor };
    }
    fclose(fp);
    return count;
}

static double lookup_emission_factor(time_t t, const EmissionInterval *arr, int count){
    for(int i=0; i<count; i++){
        if(t>=arr[i].start&&t<arr[i].end){
          return arr[i].kg_co2_per_kwh;
        }      
    }
    return -1.0;
}

CarbonReport calculate_carbon_footprint(
    const char *sensor_file,
    const char *emission_file,
    const char *start_time,
    const char *end_time)
{
    CarbonReport report={0};

    time_t start_t=parse_timestamp_utc(start_time);
    time_t end_t=parse_timestamp_utc(end_time);

    if(start_t==-1||end_t==-1||end_t<start_t)
        return report;

    EmissionInterval intervals[MAX_INTERVALS];
    int interval_count=load_emission_profile(emission_file, intervals, MAX_INTERVALS);
    if (interval_count<=0){
       return report;
    }
    FILE *fp=fopen(sensor_file, "r");
    if(!fp){
        return report;
    }
    char line[MAX_LINE_LENGTH];
    double temp_sum=0.0;
    int skipped=0;
    fgets(line, sizeof(line), fp); 
    while (fgets(line, sizeof(line), fp)) {
        line[strcspn(line, "\r\n")] = '\0';
        char *saveptr;
        char *tok=strtok_r(line, ",", &saveptr);
        if(!tok){ 
            skipped++; 
            continue; 
        }
        time_t t=parse_timestamp_utc(tok);
        if(t<start_t||t>end_t){
           continue;
        } 
        tok=strtok_r(NULL, ",", &saveptr);
        if(!tok){ 
            skipped++; 
            continue; 
        }
        char *endptr;
        double kwh=strtod(tok, &endptr);
        if(*endptr){ 
            skipped++; 
            continue; 
        }

        tok=strtok_r(NULL, ",", &saveptr);
        if(!tok){ 
            skipped++; 
            continue; 
        }
        double temp=strtod(tok, &endptr);
        if(*endptr){ 
            skipped++; 
            continue; 
        }
        double factor=lookup_emission_factor(t, intervals, interval_count);
        if(factor<0){ 
            skipped++; 
            continue; 
        }
        report.total_kwh+=kwh;
        report.total_carbon_kg+=kwh*factor;
        temp_sum+=temp;
        report.sensor_readings_count++;
    }
    fclose(fp);
    if (report.sensor_readings_count>0){
       report.average_temp_celsius=temp_sum/report.sensor_readings_count;
    }
    report.skipped_rows=skipped;
    return report;
}

void print_final_report(const GenomicStats *stats, const CarbonReport *report){
    assert(stats && report);

    printf("\n========== ECO-EFFICIENCY REPORT ==========\n\n");

    printf("Energy Consumed:      %.2f kWh\n", report->total_kwh);
    printf("Carbon Emitted:       %.3f kg CO₂\n", report->total_carbon_kg);
    printf("Avg Lab Temperature:  %.1f °C\n", report->average_temp_celsius);
    printf("Skipped Sensor Rows:  %d\n\n", report->skipped_rows);
    if(report->total_carbon_kg>0.0){
        double eco_score=stats->mean_coverage/report->total_carbon_kg;
        printf("------------------------------------------\n");
        printf("ECO-EFFICIENCY SCORE: %.2f x/kg\n", eco_score);
        printf("------------------------------------------\n");
    } 
    else {
        printf("ECO-EFFICIENCY SCORE: N/A (zero carbon)\n");
    }
    printf("\n==========================================\n");
}