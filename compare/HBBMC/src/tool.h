#ifndef MCE_TOOL_H
#define MCE_TOOL_H


#include <cstdio>
#include <cstdlib>
#include <sys/resource.h>
#include <functional>
#include <vector>
#include "hash.hpp"

#define N_EDGES 60000000 // 6M
#define N_NODES 200000000 // 5M

using namespace std;

class Edge_t {
public:
    int u = 0;
    int v = 0;

    Edge_t();
    Edge_t(int u, int v);

    bool operator<(const Edge_t& e) const;
    bool operator==(const Edge_t& e) const;

    struct Hash_Edge_t {
        size_t operator()(const Edge_t& e) const;
    };
};

class Graph_t {
public:
    int v_size = 0;
    int e_size = 0;

    Edge_t* edges = nullptr;
    int* new2old = nullptr;

    std::vector<CuckooHash> cuhash;

    Graph_t();
    ~Graph_t();

    bool connect(int u, int v);
    void read_graph(const char* r_file);
    void clean_edges(const char* r_file, const char* w_file);
    void er_model(int n, int rho);
};


class HashMap_t {

    struct HashItem_t {
        Edge_t key;
        int val = -1;
    };

private:
    const size_t TAB_SIZE = 0x7fffff;
    const size_t MAX_COLL = 500;
    int* table_size;
    HashItem_t** table;

public:
    HashMap_t();
    ~HashMap_t();

    int exist(const Edge_t& key);
    void insert(const Edge_t& key, int val);
    void remove(const Edge_t& key);
};

class KeyVal_t {
public:
    int key = 0;
    int val = 0;

    KeyVal_t();
    KeyVal_t(int key, int val);
    bool operator<(const KeyVal_t& kv) const;
};

class Heap_t {
private:
    unsigned n{};
    KeyVal_t *kv_list = nullptr;
    int *pt{};

    void swap(unsigned i, unsigned j);
    void bubble_up(unsigned i);
    void bubble_down(unsigned i);

public:
    Heap_t();
    ~Heap_t();

    bool empty();
    void insert(KeyVal_t kv);
    KeyVal_t pop();
    void update(unsigned key);
    KeyVal_t min_element();

    void make_heap(const int *v_list, unsigned v_size);
};

class Stack_t {
public:
    unsigned n;
    int* v_list = nullptr;
    int top;

    Stack_t();
    ~Stack_t();

    bool empty();
    void pop();
    void push(int v);
    void print();

    void make_stack(int size);
};

class Queue_t {
public:
    unsigned n;
    int* v_list = nullptr;
    int begin, end;

    Queue_t();
    ~Queue_t();

    bool empty();
    bool full();
    int pop();
    void push(int v);

    void make_queue(int size);
};

class Stat_t {
private:
    struct rusage s;
    struct rusage t;

public:
    double runtime;

    Stat_t();

    void setup();
    void collapse();
};



#endif //MCE_TOOL_H
