#ifndef TRIE_H
#define TRIE_H

#include "basic.h"

typedef struct TrieNode{
    int kmer_end;
    unsigned int count;
    struct TrieNode* children[5];
}TrieNode;

TrieNode* trie_create();
void trie_insert(TrieNode* root, const char* kmer);
unsigned int trie_search(TrieNode* root, const char* kmer);
void trie_free(TrieNode* root);
void analyze_and_correct_kmer(TrieNode* root, char* sequence);

#endif