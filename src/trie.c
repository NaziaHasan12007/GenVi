#include "trie.h"
#include <ctype.h>
#include <limits.h>

static int dna_to_index(char c){
    switch(toupper((unsigned char)c)){
        case 'A':return 0;
        case 'C':return 1;
        case 'T':return 2;
        case 'G':return 3;
        default :return 4;
    }
}

TrieNode* trie_create(){
    TrieNode* n=(TrieNode*)calloc(1, sizeof(TrieNode));
    if(n==NULL){
       fprintf(stderr,"[Trie] Calloc failed");
       return NULL;
    }
    return n;
}

void trie_free(TrieNode* root){
    if(!root){
        return;
    }
    for(int i=0; i<5; i++){
        if(root->children[i]){
            trie_free(root->children[i]);
        }
    }
    free(root);
}

void trie_insert(TrieNode* root, const char* kmer){
    if(!root || !kmer){
        return;
    }
    TrieNode* current=root;
    size_t len=strlen(kmer);

    for(size_t i=0; i<len; i++){
       int index=dna_to_index(kmer[i]);
       if((current->children[index])==NULL){
          current->children[index]=trie_create();
          if(current->children[index] == NULL){
                return;
            }
       }
       current=current->children[index];
       if(current->count != UINT_MAX){
          current->count++;
       }
    }
}

static TrieNode* trie_search_node(TrieNode* root, const char* kmer){
    if(!root || !kmer){
        return NULL;
    }
    size_t len=strlen(kmer);
    TrieNode* current=root;
    for(size_t i=0; i<len; i++){
       int index=dna_to_index(kmer[i]);
       if((current->children[index])==NULL){
          return NULL;
          current=current->children[index];
       }
    }
    if(current && current->kmer_end){
        return current;
    }  
   return NULL;
}
unsigned int trie_search(TrieNode* root, const char* kmer){
    TrieNode* node=trie_search_node(root, kmer);
    if(node!=NULL){
      return node->count;
    }
    return 0u;
}

static char find_strong_neighbor(TrieNode* root, char* kmer_buffer, int k){
    if(!root || !kmer_buffer || k<=0){
      return '\0';
    } 
    const char bases[4]={'A','C','G','T'};
    unsigned int best_count=0;
    unsigned int second_best_count=0;
    char best_base='\0';
    for(int i=0; i<k; i++){
        char original=kmer_buffer[i];
        for(int j=0; j<4; j++){
            char nj=bases[j];
            if(toupper((unsigned char)original)==nj){
              continue;
            } 
            kmer_buffer[i]=nj; 
            TrieNode* node=trie_search_node(root, kmer_buffer);
            unsigned int count;
            if(node!=NULL){
              count=node->count;
            }
            else{
              count=0u;
            }
            if(count>(unsigned int)STRONG_THRESHOLD){
                if(count>best_count){
                    second_best_count=best_count;
                    best_count=count;
                    best_base=nj;
                } 
                else if(count>second_best_count) {
                    second_best_count=count;
                }
            }
        }
        kmer_buffer[i]=original; 
    }
    if(best_count>0 && (best_count>=(second_best_count+(unsigned int)STRONG_MARGIN))){
        return best_base;
    }
    return '\0';
}

void analyze_and_correct_kmer(TrieNode* root, char* sequence){
    if(!root || !sequence){
      return;
    } 
    int k=SEED_KMER_SIZE;
    size_t sequence_len=strlen(sequence);
    if(sequence_len<(size_t)k){
      return;
    } 
    char* kmer_buffer=(char*)malloc((size_t)k+1);
    if(!kmer_buffer){
        fprintf(stderr, "[Trie] malloc failed for kmer_buffer\n");
        return;
    }
    kmer_buffer[k]='\0';
    char* original_copy=(char*)malloc(sequence_len+1);
    if(!original_copy){
        fprintf(stderr, "[Trie] malloc failed for original_copy\n");
        free(kmer_buffer);
        return;
    }
    memcpy(original_copy, sequence, sequence_len+1);
    for(size_t i=0; i<=sequence_len-(size_t)k; i++){
        memcpy(kmer_buffer, original_copy+i, k);
        kmer_buffer[k]='\0';
        TrieNode* node=trie_search_node(root, kmer_buffer);
        unsigned int count;
        if(node != NULL){
          count = node->count;
        }
        else{
           count=0u;
        }
        if(count<(unsigned int)WEAK_THRESHOLD){
            char correction=find_strong_neighbor(root, kmer_buffer, k);
            if(correction!='\0'){
                size_t correction_pos=i+(size_t)((k-1)/2);
                if(correction_pos<sequence_len){
                    sequence[correction_pos]=correction;
                }
                
            }
        }
    }
    free(kmer_buffer);
    free(original_copy);
}
