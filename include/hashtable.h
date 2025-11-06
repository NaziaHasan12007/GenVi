#ifndef HASHTABLE_H
#define HASHTABLE_H

#include "basic.h"

typedef struct Ht_Node{
    char* key;
    unsigned long* locations;
    size_t location_count;
    size_t location_capacity;
    struct Ht_Node* next;
}Ht_Node;

typedef struct HashTable{
    size_t size;
    Ht_Node** items;

}HashTable;

unsigned long hash_function(const char* key);
HashTable* ht_create(void);
void ht_free(HashTable* table);
void ht_insert(HashTable* table, const char* key, unsigned long location);
Ht_Node* ht_search(HashTable* table, const char* key);

#endif