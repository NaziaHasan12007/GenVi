#ifndef RECOMMENDATION_H
#define RECOMMENDATION_H

#include<time.h>
#include "basic.h"

typedef struct{
    double carbon_weight;
    double cooling_weight;
    double urgency_weight;
}Weights;

typedef struct{
    time_t timestamp;
    double forecast_carbon_intensity;
    double forecast_temparature_c;
    int urgency;
}ParameterNode;

typedef struct{
   time_t best_time;
   double predicted_score;
   double carbon_reccomedation;
   double cooling_reccomedation;
   int allowed_delay;
}Recommendation;

void get_parameter_recommendation(
    const ParameterNode *nodes,
    int count,
    Weights prefs,
    Recommendation *out_rec
);

#endif
