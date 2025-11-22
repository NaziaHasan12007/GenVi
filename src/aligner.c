#include "aligner.h"

static int max(int a, int b, int c, int d){
    int max=a;
    if(b>max){
        max=b;
    }
    if(c>max){
        max=c;
    }
    if(d>max){
        max=d;
    }
    return max;
}

int run_dp_sw(const char* sample_seq, const char* genome_seq){
    int sample_len=strlen(sample_seq);
    int genome_len=strlen(genome_seq);
    int* pre_row=(int*)calloc(genome_len+1, sizeof(int));
    int* cur_row=(int*)calloc(genome_len+1, sizeof(int));
    
    if(pre_row==NULL || cur_row==NULL){
        fprintf(stderr, "[Aligner error] calloc failed for DP matrix allocation");
        free(pre_row);
        free(cur_row);
        return 0;
    }
    int max_score=0;
    int match_mismatch_score;
    
    for(int i=1; i<=sample_len; i++){
        for(int j=1; j<=genome_len; j++){
            if(sample_seq[i-1]==genome_seq[i-1]){
               match_mismatch_score=MATCH_SCORE;
            }
            else{
               match_mismatch_score=MISMATCH_PENALTY;
            }
            int diagonal_score=pre_row[j-1]+match_mismatch_score;
            int up_score=pre_row[j]+GAP_OPEN_PENALTY;
            int left_score=cur_row[j-1]+GAP_OPEN_PENALTY;
            int score=max(0, diagonal_score, up_score, left_score);
            cur_row[j]=score;
            if(score>max_score){
                max_score=score;
            }

        }
        int* temp=pre_row;
        pre_row=cur_row;
        cur_row=temp;
    }
    free(pre_row);
    free(cur_row);
    return max_score;

}