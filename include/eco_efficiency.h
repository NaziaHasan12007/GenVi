#ifndef ECO_EFFICIENCY_H
#define ECO_EFFICIENCY_H

#include "basic.h"
#include"statistics.h"

#define E_CPU_POWER 45.0
#define CARBON_INTENSITY_FACTOR 475.0
#define PUE_FACTOR 1.67 

typedef struct {
   double total_energy;
   double total_CO2;
   double avg_power;
   double duration;
   int data_point_count;
}EcoEffieciency;

EcoEffieciency carbon_footprint(const char* filename, const char* start_time, const char* end_time);

void print_final_report(GenomicStats* stats, EcoEffieciency* carbon);

#endif