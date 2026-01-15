#include "aligner.h"
#include <string.h>
#include <stdlib.h>

#define MAX(a,b) (((a)>(b))?(a):(b))

int run_dp_sw(const char* read, const char* genome_segment, int read_len, int max_window_len, int* buf_pre, int* buf_cur){
    int window_len=(int)(read_len*1.1); 
    if(window_len>max_window_len){
       window_len=max_window_len;
    } 
    memset(buf_pre, 0, (window_len+1)*sizeof(int));
    memset(buf_cur, 0, (window_len+1)*sizeof(int));
    int max_score=0;
    for(int i=1; i<=read_len; i++){
        for(int j=1; j<=window_len; j++){
            int match=(read[i-1]==genome_segment[j-1])?MATCH_SCORE:MISMATCH_PENALTY;
            int diag=buf_pre[j-1]+match;
            int up=buf_pre[j]+GAP_OPEN_PENALTY;
            int left=buf_cur[j-1]+GAP_OPEN_PENALTY;
            
            int score=MAX(0, MAX(diag, MAX(up, left)));
            
            buf_cur[j]=score;
            if (score > max_score) max_score = score;
        }
        int* temp=buf_pre; 
        buf_pre=buf_cur; 
        buf_cur=temp;
    }
    return max_score;
}