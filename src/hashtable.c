#include "hashtable.h"

static Ht_Node* create_new_node(const char* key, unsigned long location);
static void add_location_to_node(Ht_Node* node, unsigned long location);

unsigned long hash_function(const char* str){
    unsigned long hash=5381;
    int c;
    while((c=*str++)){
        hash=((hash<<5)+hash)+c; 
    }
    return hash;
}

HashTable* ht_create(void){
    HashTable* table=(HashTable*)malloc(sizeof(HashTable));
    if(table==NULL){
        fprintf(stderr, "[Hash Table ERROR] malloc failed for HashTable struct.\n");
        return NULL;
    }
    table->size=(size_t)HASH_TABLE_SIZE;
    table->items=(Ht_Node**)calloc(table->size, sizeof(Ht_Node*));
    if(table->items==NULL){
        fprintf(stderr, "[Hash Table ERROR] calloc failed for HashTable items array.\n");
        free(table);
        return NULL;
    }
    return table;
}

void ht_free(HashTable* table){
    if(table==NULL) return;
    for(size_t i=0; i<table->size; i++){
        Ht_Node* node=table->items[i];
        
        while(node!=NULL){
            Ht_Node* next=node->next; 
            free(node->key);         
            free(node->locations);   
            free(node);              
            node = next;
        }
    }
    free(table->items);
    free(table);
}

void ht_insert(HashTable* table, const char* key, unsigned long location){
    unsigned long hash=hash_function(key);
    size_t slot=hash%table->size; 
    Ht_Node* node=table->items[slot];
    Ht_Node* prev=NULL;
    while(node!=NULL){
        if(strcmp(node->key, key)==0){
            add_location_to_node(node, location);
            return; 
        }
        prev=node;
        node=node->next;
    }
    Ht_Node* new_node=create_new_node(key, location);
    if(new_node==NULL){
        return; 
    }
    if(prev==NULL){
        table->items[slot]=new_node;
    } 
    else {
        prev->next=new_node;
    }
}

Ht_Node* ht_search(HashTable* table, const char* key){
    unsigned long hash=hash_function(key);
    size_t slot=hash%table->size;
    Ht_Node* node=table->items[slot];

    while(node!=NULL){
        if(strcmp(node->key, key)==0){
            return node;
        }
        node=node->next;
    }
    return NULL;
}

static Ht_Node* create_new_node(const char* key, unsigned long location){
    Ht_Node* node=(Ht_Node*)malloc(sizeof(Ht_Node));
    if(node==NULL){
        fprintf(stderr, "[Hash Table ERROR] malloc failed for Ht_Node.\n");
        return NULL;
    }
    node->key=(char*)malloc(strlen(key)+1);
    if(node->key==NULL){
        fprintf(stderr, "[Hash Table ERROR] malloc failed for node key.\n");
        free(node);
        return NULL;
    }
    strcpy(node->key, key);
    node->location_capacity=8; 
    node->locations=(unsigned long*)malloc(node->location_capacity*sizeof(unsigned long));
    if(node->locations==NULL){
        fprintf(stderr, "[Hash Table ERROR] malloc failed for locations array.\n");
        free(node->key);
        free(node);
        return NULL;
    }
    node->locations[0]=location;
    node->location_count=1;
    node->next=NULL;
    return node;
}

static void add_location_to_node(Ht_Node* node, unsigned long location){
    if(node->location_count==node->location_capacity){ 
        size_t new_capacity=node->location_capacity*2;
        unsigned long* new_locations=(unsigned long*)realloc(node->locations, new_capacity*sizeof(unsigned long));
        if(new_locations==NULL){
            fprintf(stderr, "[Hash Table ERROR] realloc failed for locations array.\n");
            return; 
        }
        node->locations=new_locations;
        node->location_capacity=new_capacity;
    }
    node->locations[node->location_count]=location;
    node->location_count++;
}
