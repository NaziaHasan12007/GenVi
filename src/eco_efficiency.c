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
    static int last=0;
    if(count<=0){
        return -1.0;
    }
    if(t>=arr[last].start&&t<arr[last].end){
        return arr[last].kg_co2_per_kwh;
    }
    for(int i=0;i<count;i++){
        if(t>=arr[i].start&&t<arr[i].end){
            last=i;
            return arr[i].kg_co2_per_kwh;
        }
    }
    return -1.0;
}

CarbonReport calculate_carbon_footprint(const char *sensor_file, const char* emmission_file, const char *start_time, const char *end_time)
{
    CarbonReport report={0};

    time_t start_t=parse_timestamp_utc(start_time);
    time_t end_t=parse_timestamp_utc(end_time);

    if(start_t==-1||end_t==-1||end_t<start_t)
        return report;

    EmissionInterval intervals[MAX_INTERVALS];
    int interval_count=load_emission_profile(emmission_file, intervals, MAX_INTERVALS);
    if (interval_count<=0){
       return report;
    }
    FILE *fp=fopen(sensor_file, "r");
    if(!fp){
        return report;
    }
    char line[MAX_LINE_LENGTH];
    double temp_weighted_sum=0.0;
    double temp_time_sum=0.0;
    int skipped=0;
    time_t prev_t=0;

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

        if(prev_t==0){
            prev_t=t;
            continue;
        }

        double delta_hours=difftime(t, prev_t)/3600.0;
        if(delta_hours<=0){
            prev_t=t;
            continue;
        }

        tok=strtok_r(NULL, ",", &saveptr);
        if(!tok){ 
            skipped++; 
            prev_t=t;
            continue; 
        }
        char *endptr;
        double energy_kwh=strtod(tok, &endptr);
        if(*endptr){ 
            skipped++; 
            prev_t=t;
            continue; 
        }

        tok=strtok_r(NULL, ",", &saveptr);
        if(!tok){ 
            skipped++; 
            prev_t=t;
            continue; 
        }
        double temp=strtod(tok, &endptr);
        if(*endptr){ 
            skipped++; 
            prev_t=t;
            continue; 
        }
        double factor=lookup_emission_factor(t, intervals, interval_count);
        if(factor<0){ 
            skipped++; 
            prev_t=t;
            continue; 
        }
        report.total_kwh+=energy_kwh;
        report.total_carbon_kg+=energy_kwh*factor;

        temp_weighted_sum+=temp*delta_hours;
        temp_time_sum+=delta_hours;

        report.sensor_readings_count++;
        prev_t=t;
    }
    fclose(fp);
    if (temp_time_sum>0.0){
       report.average_temp_celsius=temp_weighted_sum/temp_time_sum;
    }
    report.skipped_rows=skipped;

    if(report.total_kwh>0.0){
        report.carbon_intensity=report.total_carbon_kg/report.total_kwh;
    }
    return report;
}

void print_final_report(const GenomicStats *stats, const CarbonReport *report, const Recommendation *rec){
   if(!stats || !report){
      return;
    } 
    printf("[Phase 4&5] Starting: Genome statistics, Carbon report, Eco-efficiciency score and Recomendations...\n\n");
    printf("\n========== GENOMIC STATISTICS ==========\n\n");
    printf("Genome Length: %ld\n", stats->genome_length);
    printf("Total Reads Processed: %ld\n", stats->total_reads_processed);
    printf("Total Reads Aligned: %ld\n", stats->total_reads_aligned);
    printf("Reads Discarded or Failed to Map: %ld\n", stats->total_reads_trimmed + stats->total_reads_failed_to_map);
    printf("Percent Reads Discarded: %.2f%%\n", stats->percent_reads_discarded);
    printf("Mean Coverage: %.2f\n", stats->mean_coverage);
    printf("Median Coverage: %.2f\n", stats->median_coverage);
    printf("Coverage Std Dev: %.2f\n", stats->std_dev_coverage);
    printf("Coverage Breadth: %.2f%%\n", stats->genome_breadth_percent);
    printf("GC Content: %.2f%%\n", stats->gc_content_percent);
    printf("N50 Metric: %ld\n", stats->n50_metric);
    printf("Coverage Q1: %.2f\nQ3: %.2f\nIQR: %.2f\n", stats->coverage_q1, stats->coverage_q3, stats->coverage_iqr);

    printf("\n========== ECO-EFFICIENCY REPORT ==========\n\n");
    printf("Energy Consumed:      %.2f kWh\n", report->total_kwh);
    printf("Carbon Emitted:       %.3f kg CO₂\n", report->total_carbon_kg);
    printf("Carbon Intensity:     %.3f kg CO₂/kWh\n", report->carbon_intensity);
    printf("Avg Lab Temperature:  %.1f °C\n", report->average_temp_celsius);
    printf("Skipped Sensor Rows:  %d\n\n", report->skipped_rows);

    if(report->total_carbon_kg > 0.0){
        double eco_score = stats->mean_coverage / report->total_carbon_kg;
        printf("------------------------------------------\n");
        printf("ECO-EFFICIENCY SCORE: %.2f x/kg\n", eco_score);
        printf("------------------------------------------\n");
    } else {
        printf("ECO-EFFICIENCY SCORE: N/A (zero carbon)\n");
    }
    if(rec!=NULL) {
        printf("----- Environmental Recommendation -----\n");
        printf("Predicted Score        : %.4f (lower is better)\n", rec->predicted_score);
        printf("Suggested Carbon Level : %.2f gCO2/kWh\n", rec->carbon_reccomedation);
        printf("Suggested Cooling Adj. : %.2f C\n", rec->cooling_reccomedation);
        printf("Allowed Delay (hrs)    : %d\n\n", rec->allowed_delay);
    }
    printf("[Phase 4&5] Completed : Genome statistics, Carbon report, Eco-efficiciency score and Recomendations");
    printf("Done Computing!! See you again!\n");
    printf("\n==========================================\n");
}