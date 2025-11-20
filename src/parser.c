#include "parser.h"

char* load_genome_from_fasta(const char* filename, long* genome_length){
    FILE* fp= fopen(filename, "r");
    if(fp==NULL){
        fprintf(stderr, "[Parser Error] Can't open FASTA file: %s\n", filename);
        return NULL;
    }
    char buffer_line[MAX_LINE_LENGTH];
    long buffer_size=1048576;
    char* genome_string=(char*)calloc(buffer_size, sizeof(char));
    if(genome_string==NULL){
        fprintf(stderr, "[Parser ERROR] Initial calloc failed for genome string (out of memory).\n");
        fclose(fp);
        return NULL;
    }
    long current_position=0;
    if(fgets(buffer_line, MAX_LINE_LENGTH, fp)==NULL){
        fprintf(stderr, "[Parser Error] FASTA file is not in the correct format or is empty.\n");
        fclose(fp);
        free(genome_string);
        return NULL;
    }
    while(fgets(buffer_line, MAX_LINE_LENGTH, fp)!=NULL){
        if(buffer_line[0]=='>'){
            continue;
        }
        buffer_line[strcspn(buffer_line, "\n")]='\0';
        buffer_line[strcspn(buffer_line, "\r")]='\0';
        int line_length=strlen(buffer_line);
        if(line_length==0){
            continue;
        }
        if(current_position+line_length+1>buffer_size){
            long new_buffer_size=buffer_size*2;
            while(current_position+line_length+1>new_buffer_size){
                new_buffer_size*=2;
            }
            char* new_buffer=(char*)realloc(genome_string, new_buffer_size);
            if(new_buffer==NULL){
                fprintf(stderr, "[Parser Error] Realloc failed(out of memory)");
                fclose(fp);
                free(genome_string);
                return NULL;
            }
            genome_string=new_buffer;
            buffer_size=new_buffer_size;
        }
        memcpy(genome_string+current_position, buffer_line, line_length);
        current_position+=line_length;
    }
    fclose(fp);
    genome_string[current_position]='\0';
    char*final_genome_string=(char*)realloc(genome_string, current_position+1);
    if(final_genome_string==NULL){
        fprintf(stderr, "[Parser Warning] Could not shrink-to-fit.Using oversized buffer.\n");
        genome_length=current_position;
        return genome_string; 
    }
    
    genome_length=current_position;
    return final_genome_string;
}

FILE* open_fastq_file(const char* filename){
    FILE* fp=fopen(filename, "r");
    if(fp==NULL){
        fprintf(stderr,"[Parser Error] Could not open FASTQ file: %s\n", filename);
        return NULL;
    }
    return fp;
}

int read_fastq_to_string(FILE* fp, char*header, char*sequence, char*plus, char* quality){
    if(fgets(header, MAX_LINE_LENGTH, fp)==NULL){
        return 0;
    }
    if(fgets(sequence, MAX_LINE_LENGTH, fp)==NULL){
        return 0;
    }
    if(fgets(plus, MAX_LINE_LENGTH, fp)==NULL){
        return 0;
    }
    if(fgets(quality, MAX_LINE_LENGTH, fp)==NULL){
        return 0;
    }
    header[strcspn(header, "\n")]='\0';
    header[strcspn(header, "\r")]='\0';
    sequence[strcspn(sequence, "\n")]='\0';
    sequence[strcspn(sequence, "\r")]='\0';
    plus[strcspn(plus, "\n")]='\0';
    plus[strcspn(plus, "\r")]='\0';
    quality[strcspn(quality, "\n")]='\0';
    quality[strcspn(quality, "\r")]='\0';
    return 1;
}