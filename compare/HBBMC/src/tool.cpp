#include "tool.h"
#include <cassert>
#include <set>
#include <random>

size_t h(const Edge_t& e) {
    int s_ = e.u < e.v ? e.u : e.v;
    int t_ = e.u < e.v ? e.v : e.u;
    size_t hash = 1;
    std::hash<int> H;

    hash ^= H(s_) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
    hash ^= H(t_) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
    return hash;
}

Edge_t::Edge_t() = default;

Edge_t::Edge_t(int u, int v) {
    this->u = u;
    this->v = v;
}

bool Edge_t::operator<(const Edge_t &e) const {
    return (this->u < e.u) || (this->u == e.u && this->v < e.v);
}

bool Edge_t::operator==(const Edge_t &e) const {
    return this->u == e.u && this->v == e.v;
}

size_t Edge_t::Hash_Edge_t::operator()(const Edge_t &e) const {
    return h(e);
}

Graph_t::Graph_t() {
    v_size = 0;
    e_size = 0;

    new2old = nullptr;
    edges = nullptr;
};

Graph_t::~Graph_t() {
    if (new2old) {
        free(new2old);
        new2old = nullptr;
    }

    if (edges) {
        free(edges);
        edges = nullptr;
    }
}

void Graph_t::read_graph(const char *r_file) {
    FILE *fp;
    static bool is_read = false;

    if ((fp = fopen(r_file, "r")) == nullptr) {
        printf("Cannot open file %s.\n", r_file);
        exit(0);
    }

    Edge_t e;
    int u, v, i;
    int *old2new = (int*) malloc(N_NODES * sizeof(int));
    for (i = 0; i < N_NODES; i++) old2new[i] = -1;

    new2old = (int*) malloc(N_NODES * sizeof(int));
    edges = new Edge_t [N_EDGES];

    if (!is_read) printf("Reading the graph.\n");

    while (fscanf(fp, "%d %d%*[^\n]%*c", &u, &v) == 2) {

        if (u > N_NODES || v > N_NODES) {
            printf("Enlarge N_NODES to at least %u.\n", (u > v ? u : v));
            exit(0);
        }
        if (old2new[u] == -1) {
            new2old[v_size] = u;
            old2new[u] = v_size++;
        }
        if (old2new[v] == -1) {
            new2old[v_size] = v;
            old2new[v] = v_size++;
        }

        e = Edge_t(old2new[u], old2new[v]);
        edges[e_size++] = e;
    }

    edges = (Edge_t*) realloc(edges, e_size * sizeof(Edge_t));

    if (!is_read) {
        printf("Finish reading the graph. |V| = %d, |E| = %d\n", v_size, e_size);
        is_read = true;
    }

    fclose(fp);
    free(old2new);

}

void Graph_t::clean_edges(const char *r_file, const char *w_file) {
    set<Edge_t> E;
    int u, v;

    FILE *f_r = fopen(r_file, "r");
    while (fscanf(f_r, "%u %u%*[^\n]%*c", &u, &v) == 2) {
        if (u == v) {
            printf("%u Self loop.\n", u);
            continue;
        }
        if (E.find(Edge_t(u, v)) != E.end()) {
            printf("(%u, %u) duplicates.\n", u, v);
            continue;
        }

        E.insert(Edge_t(u, v));

    }
    fclose(f_r);

    FILE *f_w = fopen(w_file, "w");
    for (Edge_t e : E) {
        fprintf(f_w, "%u %u\n", e.u, e.v);
    }
    fclose(f_w);
}

bool Graph_t::connect(int u, int v) {
    return cuhash[u].find(v);
}

void Graph_t::er_model(int n, int rho) {
    int s, t;
    long long k, e = n * rho;
    double p = (double) 2 * rho / (n - 1);
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<> dist(0, 1);

    FILE *fp = fopen((to_string((double)n/1000000) + "x" + to_string(rho) + ".clean").c_str(), "w");
    this->v_size = n;
    this->e_size = 0;

    k = 1 + log(dist(gen)) / log(1 - p);

    while (log(k) <= 2 * log(n)) {
        s = (k - 1) / n;
        t = (k - 1) % n;

        if (s < t) {
            this->e_size++;
            fprintf(fp, "%u %u\n", s, t);
        }

        k += 1 + log(dist(gen)) / log(1 - p);
    }

    printf("|V| = %d, |E| = %d\n", v_size, e_size);

    fclose(fp);
}

HashMap_t::HashMap_t() {
    table = new HashMap_t::HashItem_t* [TAB_SIZE];
    table_size = new int [TAB_SIZE]();
}

HashMap_t::~HashMap_t() {
    delete [] table;
    delete [] table_size;
}

int HashMap_t::exist(const Edge_t &key) {
    size_t pos = h(key) % TAB_SIZE;
    int i, val;
    Edge_t e;

    if (table_size[pos] > 0) {
        for (i = 0; i < table_size[pos]; i++) {
            e = table[pos][i].key;
            val = table[pos][i].val;
            if ((e.u == key.u && e.v == key.v) || (e.u == key.v && e.v == key.u)) return val;
        }
    }

    return -1;
}

void HashMap_t::insert(const Edge_t &key, int val) {
    size_t pos = h(key) % TAB_SIZE;

    if (table_size[pos] == 0)
        table[pos] = new HashMap_t::HashItem_t [MAX_COLL];

    assert(table_size[pos] < MAX_COLL);     // Enlarge max_coll or redefine hash function h(e)

    table[pos][table_size[pos]].key = key;
    table[pos][table_size[pos]++].val = val;
}

void HashMap_t::remove(const Edge_t &key) {
    size_t pos = h(key) % TAB_SIZE;
    int i;
    Edge_t e;

    if (table_size[pos] > 0) {
        for (i = 0; i < table_size[pos]; i++) {
            e = table[pos][i].key;
            if ((e.u == key.u && e.v == key.v) || (e.u == key.v && e.v == key.u)) {
                table[pos][i] = table[pos][--table_size[pos]];
                if (table_size[pos] == 0) delete [] table[pos];
                break;
            }
        }
        return;
    }

    exit(-1);
}

KeyVal_t::KeyVal_t() = default;

KeyVal_t::KeyVal_t(int key, int val) {
    this->key = key;
    this->val = val;
}

bool KeyVal_t::operator<(const KeyVal_t &kv) const {
    return this->val > kv.val || (this->val == kv.val && this->key < kv.key);
}

Heap_t::Heap_t() = default;

Heap_t::~Heap_t() {
    delete [] this->kv_list;
    delete [] this->pt;
}

void Heap_t::swap(unsigned i, unsigned j) {
    KeyVal_t kv_tmp = this->kv_list[i];
    int pt_tmp = this->pt[kv_tmp.key];
    this->pt[this->kv_list[i].key] = this->pt[this->kv_list[j].key];
    this->kv_list[i] = this->kv_list[j];
    this->pt[this->kv_list[j].key] = pt_tmp;
    this->kv_list[j] = kv_tmp;
}

void Heap_t::bubble_up(unsigned int i) {
    unsigned j = (i - 1) >> 1;
    while (i > 0) {
        if (this->kv_list[j].val > this->kv_list[i].val) {
            this->swap(i, j);
            i = j;
            j = (i - 1) >> 1;
        }
        else break;
    }
}

void Heap_t::bubble_down(unsigned int i) {
    unsigned l = (i << 1) + 1, r = l + 1, j;
    while (l < this->n) {
        j = ((r < this->n) && (this->kv_list[r].val < this->kv_list[l].val)) ? r : l;
        if (this->kv_list[i].val > this->kv_list[j].val) {
            this->swap(i, j);
            i = j;
            l = (i << 1) + 1;
            r = l + 1;
        }
        else break;
    }
}

bool Heap_t::empty() {
    return this->n == 0;
}

void Heap_t::insert(KeyVal_t kv) {
    this->pt[kv.key] = this->n;
    this->kv_list[this->n] = kv;
    this->bubble_up(this->n);
    this->n++;
}

KeyVal_t Heap_t::pop() {
    assert(!this->empty());
    KeyVal_t min = this->kv_list[0];
    this->pt[min.key] = -1;
    this->kv_list[0] = this->kv_list[--(this->n)];
    this->pt[this->kv_list[0].key] = 0;
    this->bubble_down(0);
    return min;
}

KeyVal_t Heap_t::min_element() {
    assert(!this->empty());
    return this->kv_list[0];
}

void Heap_t::update(unsigned int key) {
    int i = this->pt[key];
    if (i != -1) {
        ((this->kv_list[i]).val)--;
        this->bubble_up(i);
    }
}

void Heap_t::make_heap(const int *v_list, unsigned int v_size) {
    unsigned i;
    KeyVal_t kv;

    this->n = 0;
    this->pt = new int [v_size];
    for (i = 0; i < v_size; i++) this->pt[i] = -1;
    this->kv_list = new KeyVal_t [v_size];

    for (i = 0; i < v_size; i++) {
        kv.key = i;
        kv.val = v_list[i];
        this->insert(kv);
    }
}

extern int p;

Stack_t::Stack_t() = default;

Stack_t::~Stack_t() {
    if (v_list) {
        free(this->v_list);
        v_list = nullptr;
    }
}

void Stack_t::make_stack(int size) {
    this->n = size;
    this->v_list = (int*) malloc(size * sizeof(int));
    this->top = 0;
    p = 2;
}

bool Stack_t::empty() {
    return this->top == 0;
}

void Stack_t::pop() {
    --this->top;
}

void Stack_t::push(int v) {
    this->v_list[this->top++] = v;
}

void Stack_t::print() {
    int i;
    for (i = 0; i < this->top; i++) printf("%d ", this->v_list[i]);
    printf("\n");
}

Queue_t::Queue_t() = default;

Queue_t::~Queue_t() {
    if (v_list) {
        free(v_list);
        v_list = nullptr;
    }
}

void Queue_t::make_queue(int size) {
    this->n = size;
    this->v_list = (int*) malloc(size * sizeof(int));
    this->begin = 0;
    this->end = 0;
}

bool Queue_t::empty() {
    if (this->begin == this->end) {
        this->begin = 0;
        this->end = 0;
        return true;
    }

    return false;
}

bool Queue_t::full() {
    return this->end == this->n;
}

int Queue_t::pop() {
    assert(!this->empty());
    return this->v_list[this->begin++];
}

void Queue_t::push(int v) {
    assert(!this->full());
    this->v_list[end++] = v;
}

Stat_t::Stat_t() = default;

void Stat_t::setup() {
    if (getrusage(RUSAGE_THREAD, &(this->s)) != 0) {
        fprintf(stderr, "The running time info couldn't be collected successfully.\n");
        exit(0);
    }
}

void Stat_t::collapse() {
    if (getrusage(RUSAGE_THREAD, &(this->t)) != 0) {
        fprintf(stderr, "The running time info couldn't be collected successfully.\n");
        exit(0);
    }
    this->runtime =  ((float)(t.ru_utime.tv_sec - s.ru_utime.tv_sec)) * 1e3 +
            ((float)(t.ru_utime.tv_usec - s.ru_utime.tv_usec)) * 1e-3;      // ms
}
