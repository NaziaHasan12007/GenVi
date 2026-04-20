#include "recommendation.h"
#include<float.h>

static double normalize(double val, double ideal, double max_dev) {
    double diff=val-ideal;
    if(diff<0){
      diff=0.0;
    }  
    double n=diff/max_dev;
    if(n>1.0){
       n=1.0;
    } 
    return n;
}

static double calculate_cooling_penalty(double temp_c) {
    if(temp_c>21.0){
      return (temp_c-21.0)*0.05;
    } 
    if(temp_c<18.0){
      return(18.0-temp_c)*0.02;
    } 
    return 0.0;
}

void get_parameter_recommendation(const ParameterNode *nodes, int count, Weights prefs, Recommendation *out_rec){
    if(count<=0||!out_rec){
       return;
    } 

    double best_score=DBL_MAX;
    int best_index=0;

    for (int i=0; i<count; i++) {
        double carbon_norm=normalize(nodes[i].forecast_carbon_intensity, 200.0, 600.0);
        double cooling_penalty=calculate_cooling_penalty(nodes[i].forecast_temparature_c);
        double cooling_norm=cooling_penalty/1.0; 
        double urgency_norm=nodes[i].urgency?1.0:0.0;
        double score=carbon_norm*prefs.carbon_weight+cooling_norm*prefs.cooling_weight+urgency_norm*prefs.urgency_weight;
        if(score<best_score){
            best_score=score;
            best_index=i;
        }
    }
    const ParameterNode *best_node=&nodes[best_index];
    out_rec->best_time=best_node->timestamp;
    out_rec->predicted_score=best_score;
    out_rec->carbon_reccomedation=best_node->forecast_carbon_intensity>300.0?300.0:best_node->forecast_carbon_intensity;
    out_rec->cooling_reccomedation=best_node->forecast_temparature_c>21.0?21.0:best_node->forecast_temparature_c;
    out_rec->allowed_delay=best_node->urgency?0:1;
}