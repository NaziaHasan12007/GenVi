#include "eco_efficiency.h"
#include<string.h>
#include<stdlib.h>

static double estimated_duration(int data_point){
    if(data_point>0){
        return (double)data_point/3600.0;
    }
    return 0.0;
}

EcoEffieciency carbon_footprint(const char* filename, const char* start_time, const char* end_time){
    printf("[Eco Effeciency] Reading sensor log: %s\n", filename);
    printf("[Eco Efficiency] Time Window: %s to %s\n", start_time, end_time);

    EcoEffieciency report={0};
    FILE* fp= fopen(filename, "r");

    if(fp==NULL){
       fprintf(stderr, "[Eco Efficiency-Warning] Could not open sensor log.");
       report.avg_power = E_CPU_POWER;
       report.data_point_count= 0;
       report.duration = 0.25;
    }
    else{
        char line[256];
        double total_watts = 0;
        int count = 0;

        fgets(line, sizeof(line), fp);

        while (fgets(line, sizeof(line), fp)) {
            char timestamp[64];
            double power = 0.0;

            char* token = strtok(line, ",");
            if (token != NULL) {
                strncpy(timestamp, token, 63);
                token = strtok(NULL, ",");
                if (token != NULL) {
                    power = atof(token);
                }
            }

            if (strcmp(timestamp, start_time) >= 0 && strcmp(timestamp, end_time) <= 0) {
                total_watts += power;
                count++;
            }
        }
        fclose(fp);

        report.data_point_count = count;
        if (count > 0) {
            report.avg_power = total_watts / count;
            report.duration = (double)count / 3600.0; 
        } else {
            report.avg_power = 0;
            report.duration = 0;
        }
    }
    report.total_energy = (report.avg_power / 1000.0) * report.duration * PUE_FACTOR;
    report.total_CO2 = report.total_energy * CARBON_INTENSITY_FACTOR;

    return report;
}
