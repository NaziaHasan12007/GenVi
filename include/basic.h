#ifndef BASIC_H
#define BASIC_H

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<limits.h>
#include<math.h>

#define TRIM_WINDOW_SIZE 5
#define MIN_READ_LENGTH 35
#define TRIM_QUALITY_THRESHOLD 20
#define SEED_KMER_SIZE 21
#define MAX_LINE_LENGTH 2048
#define HASH_TABLE_SIZE 4999999UL
#define WEAK_THRESHOLD 3
#define STRONG_THRESHOLD 5
#define STRONG_MARGIN 1
#define MATCH_SCORE 5
#define MISMATCH_PENALTY -4
#define GAP_OPEN_PENALTY -10

#endif