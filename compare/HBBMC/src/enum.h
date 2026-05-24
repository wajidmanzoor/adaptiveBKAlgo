#ifndef MCE_TOMITA_H
#define MCE_TOMITA_H

#include "tool.h"

class Tomita_t {
private:

    Graph_t g;

    Stack_t R;

    int d_max = 0;
    int core_num = 0;

    int* d = nullptr;  
    int** adj = nullptr; 
    int* rank = nullptr;

    int** P = nullptr;
    int* P_size = nullptr;
    int** X = nullptr;
    int* X_size = nullptr;
    int* lab = nullptr;
    int* stat_v = nullptr;

    void print_maximal();

    bool is_plex = false;

public:

    unsigned long long maximal_clique = 0;

    Tomita_t(const char* r_file);
    ~Tomita_t();

    void Tomita_pivot();
    void Tomita_pivot_rec(int l);
    void Tomita_opt(); 
    void Tomita_opt_rec(int l);
    void Tomita_pivot_mat();
    void Tomita_pivot_mat_rec(int l);
    void Tomita_opt_mat();
    void Tomita_opt_mat_rec(int l);
    void Tomit_rev_mat();
    void Tomita_rev_mat_rec(int l);
    void list_in_plex(int l);
};


#endif //MCE_TOMITA_H
