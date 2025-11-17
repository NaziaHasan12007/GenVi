#include<trimmer.h>

int trim_read_by_quality(char* sequence, char* quality){
    int og_length=strlen(sequence);
    if(og_length<TRIM_WINDOW_SIZE){
        return og_length;
    }
    for(int i=0; i<og_length-TRIM_WINDOW_SIZE; i++){
        int window_quality_sum=0;
        for(int j=0; j<TRIM_WINDOW_SIZE; j++){
            char char_quality= quality[i+j];
            int score=(int)char_quality-33;
            window_quality_sum+=score;
        }
        float window_quality_avg=(float)window_quality_sum/TRIM_WINDOW_SIZE;
        if(window_quality_avg<TRIM_QUALITY_THRESHOLD){
            sequence[i]='\0';
            quality[i]='\0';     
        }
    }
    return og_length;
}