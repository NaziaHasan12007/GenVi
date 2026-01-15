#ifndef ALIGNER_H
#define ALIGNER_H

#include "basic.h" 

int run_dp_sw(const char* read, const char* genome_segment, int read_len, int max_window_len, int* buf_pre, int* buf_cur);

#endif