#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>
#include <algorithm>
#include <queue>
#include <iostream>

#include "tool.h"
#include "io.hpp"
//#include "heap.hpp"
#include "hash.hpp"

class UniGraph_t {

public:
    int32_t n, m, maxDu, maxDv;

    std::vector <Edge_t> edges;
    std::vector <int32_t> d, cd;

    std::vector <CuckooHash> cuhash;
    std::vector <int32_t> old_lables;

    int32_t lrs[2];
    int32_t k = 0;

    UniGraph_t() = default;

    UniGraph_t(const char* filePath, int mode, const std::string &order) {

        FILE *fp;
        static bool is_read = false;

        if ((fp = fopen(filePath, "r")) == nullptr) {
            printf("Cannot open file %s.\n", filePath);
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


        if (mode == 1) read(filePath, order);
        else readWithOutUVM(filePath, order);

        cuhash[0].resize(n[0]);
        cuhash[1].resize(n[1]);
        for (int32_t t = 0; t <= 1; t++) {
            for (int32_t i = 0; i < n[t]; i++) {
                int32_t d = p[t][i + 1] - p[t][i];
                cuhash[t][i].reserve(d + 1);
                for (int32_t j = p[t][i]; j < p[t][i + 1]; ++j)
                    cuhash[t][i].insert(e[t][j]);
            }
        }
    }

    void read(const std::string &filePath, const std::string &order) {
        fastIO in(filePath, "r");

        n = in.getUInt();
        m = in.getUInt();

        edges.resize(m);
        e1.resize(m);
        e2.resize(m);
        pU.resize(n1 + 5);
        pV.resize(n2 + 5);

        for (int32_t i = 0; i < m; i++) {
            edges[i].u = in.getUInt();
            edges[i].v = in.getUInt();
        }

        cores[0].resize(n1);
        cores[1].resize(n2);
        old_lables[0].resize(n1);
        old_lables[1].resize(n2);
//        if (order == "core") changeToCoreOrderVersion2();
//        else if (order == "two") changeToTwoHopCoreOrder();
//        else {
//            printf("error order\n");
//            exit(1);
//        }

        n[0] = n1;
        n[1] = n2;
        p[0] = std::move(pU);
        p[1] = std::move(pV);
        e[0] = std::move(e1);
        e[1] = std::move(e2);
    }

    void readWithOutUVM(const std::string &filePath, const std::string &order) {
        fastIO in(filePath, "r");

        m = 0;
        int32_t minL = 1 << 30, minR = 1 << 30;
        int32_t maxL = 0, maxR = 0;
        while (!in.empty()) {
            int32_t u = in.getUInt();
            int32_t v = in.getUInt();

            edges.push_back(Edge{u, v});
            m++;
            minL = std::min(minL, u);
            minR = std::min(minR, v);
            maxL = std::max(maxL, u);
            maxR = std::max(maxR, v);
        }

        for (int32_t i = 0; i < m; i++) {
            edges[i].u -= minL;
            edges[i].v -= minR;
        }

        n1 = maxL - minL + 1;
        n2 = maxR - minR + 1;

        e1.resize(m);
        e2.resize(m);
        pU.resize(n1 + 5);
        pV.resize(n2 + 5);

        cores[0].resize(n1);
        cores[1].resize(n2);
        old_lables[0].resize(n1);
        old_lables[1].resize(n2);
//        if (order == "core") changeToCoreOrderVersion2();
//        else if (order == "two") changeToTwoHopCoreOrder();
//        else {
//            printf("order vertex id-is\n");
            VertexIDOrder();
//        }

        n[0] = n1;
        n[1] = n2;
        p[0] = std::move(pU);
        p[1] = std::move(pV);
        e[0] = std::move(e1);
        e[1] = std::move(e2);
    }

//    void changeToTwoHopCoreOrder() {
//        std::vector <int32_t> d1, d2;
//
//        d1.resize(n1);
//        d2.resize(n2);
//
//        for (int32_t i = 0; i < m; i++) {
//            ++d1[edges[i].u];
//            ++d2[edges[i].v];
//        }
//
//        maxDu = 0;
//        for (int32_t i = 0; i < n1; i++) {
//            maxDu = std::max(maxDu, d1[i]);
//        }
//        maxDv = 0;
//        for (int32_t i = 0; i < n2; i++) {
//            maxDv = std::max(maxDv, d2[i]);
//        }
//
//        pU[0] = 0;
//        for (int32_t u = 0; u < n1; u++) {
//            pU[u + 1] = d1[u] + pU[u];
//        }
//        for (int32_t i = 0; i < m; i++) {
//            e1[pU[edges[i].u]++] = edges[i].v;
//        }
//        pU[0] = 0;
//        for (int32_t u = 0; u < n1; u++) {
//            pU[u + 1] = d1[u] + pU[u];
//        }
//
//        pV[0] = 0;
//        for (int32_t v = 0; v < n2; v++) {
//            pV[v + 1] = d2[v] + pV[v];
//        }
//        for (int32_t i = 0; i < m; i++) {
//            e2[pV[edges[i].v]++] = edges[i].u;
//        }
//        pV[0] = 0;
//        for (int32_t v = 0; v < n2; v++) {
//            pV[v + 1] = d2[v] + pV[v];
//        }
//
//        int32_t N = n1 + n2;
//        std::vector <int32_t> dequeue;
//        std::vector <int32_t> cdeg;
//        dequeue.resize(N);
//        cdeg.resize(N, 0);
//        int32_t tsize = 0;
//
//        n[0] = n1;
//        n[1] = n2;
//        for (int32_t tt = 0; tt <= 1; tt++) {
//            for (int32_t i = 0; i < n[tt]; i++) {
//                int32_t d = 0;
//                if (tt == 0) {
//                    d = d1[i];
//                    cdeg[i] = d;
//                    if (d + k < lrs[1])
//                        dequeue[tsize++] = i;
//                } else {
//                    d = d2[i];
//                    cdeg[i + n1] = d;
//                    if (d + k < lrs[0])
//                        dequeue[tsize++] = i + n1;
//                }
//            }
//        }
//
//
//        int32_t rsize = 0;
//        while (tsize > rsize) {
//            int32_t s = rsize;
//            rsize = tsize;
//            for (int32_t i = s; i < rsize; ++i) {
//                int v = dequeue[i];
//                if (v >= n1) {
//                    v -= n1;
//                    for (int32_t j = pV[v]; j < pV[v + 1]; ++j) {
//                        int32_t u = e2[j];
//                        if (cdeg[u] + k >= lrs[1]) {
//                            cdeg[u]--;
//                            if (cdeg[u] + k < lrs[1])
//                                dequeue[tsize++] = u;
//                        }
//                    }
//                } else {
//                    for (int32_t j = pU[v]; j < pU[v + 1]; ++j) {
//                        int32_t u = e1[j];
//                        u += n1;
//                        if (cdeg[u] + k >= lrs[0]) {
//                            cdeg[u]--;
//                            if (cdeg[u] + k < lrs[0])
//                                dequeue[tsize++] = u;
//                        }
//                    }
//                }
//            }
//        }
//
//        int32_t maxTwoHopDegreeU = 0, maxTwoHopDegreeV = 0;
//        d1.clear();
//        std::vector <int32_t> stk;
//        int32_t n = std::max(n1, n2);
//        std::vector <int32_t> ids(n);
//        std::vector <int32_t> keys(n);
//        std::vector <int32_t> labelsL(n1);
//        std::vector <int32_t> labelsR(n2);
//        int32_t l = 0;
//
//        for (int32_t u = 0; u < n1; u++) ids[u] = u;
//        for (int32_t u = 0; u < n1; u++) {
//            if (cdeg[u] + k < lrs[1]) {
//                keys[u] = 0;
//                continue;
//            }
//            stk.clear();
//            int32_t tmp = 0;
//
//            for (int32_t i = pU[u]; i < pU[u + 1]; i++) {
//                int32_t v = e1[i];
//                if (cdeg[v + n1] + k < lrs[0]) continue;
//                tmp++;
//                for (int32_t j = pV[v]; j < pV[v + 1]; j++) {
//                    int32_t w = e2[j];
//                    if (d1[w] == 0 && cdeg[w] + k >= lrs[1]) {
//                        d1[w] = 1;
//                        stk.push_back(w);
//                    }
//                }
//            }
//            tmp += stk.size();
//
//            maxTwoHopDegreeU = std::max(maxTwoHopDegreeU, tmp);
//            keys[u] = tmp;
//
//            for (auto w: stk) d1[w] = 0;
//        }
//
//        Heap_t heap(n1, maxTwoHopDegreeU + 1);
//        heap.init(n1, maxTwoHopDegreeU + 1, ids.data(), keys.data());
//        core[0] = 0;
//        core[1] = 0;
//
//        for (int32_t i = 0; i < n1; i++) {
//            int32_t u, degU;
//
//            if (!heap.pop_min(u, degU)) printf("errorLheap\n");
//            old_lables[0][l] = u;
//            labelsL[u] = l++;
//            core[0] = std::max(core[0], degU);
//            cores[0][l - 1] = core[0];
//            if (core[0] == 0) continue;
//            stk.clear();
//            for (int32_t i = pU[u]; i < pU[u + 1]; i++) {
//                int32_t v = e1[i];
//                if (cdeg[v + n1] + k < lrs[0]) continue;
//                for (int32_t j = pV[v]; j < pV[v + 1]; j++) {
//                    int32_t w = e2[j];
//                    if (d1[w] == 0 && cdeg[w] + k >= lrs[1]) {
//                        d1[w] = 1;
//                        stk.push_back(w);
//                        heap.decrement(w, 1);
//                    }
//                }
//            }
//
//            for (auto w: stk) d1[w] = 0;
//        }
//
//        d2.clear();
//        l = 0;
//
//        for (int32_t u = 0; u < n2; u++) ids[u] = u;
//        for (int32_t u = 0; u < n2; u++) {
//            if (cdeg[u + n1] + k < lrs[0]) {
//                keys[u] = 0;
//                continue;
//            }
//            stk.clear();
//            int32_t tmp = 0;
//
//            for (int32_t i = pV[u]; i < pV[u + 1]; i++) {
//                int32_t v = e2[i];
//                if (cdeg[v] + k < lrs[1]) continue;
//                tmp++;
//                for (int32_t j = pU[v]; j < pU[v + 1]; j++) {
//                    int32_t w = e1[j];
//                    if (d2[w] == 0 && cdeg[w + n1] + k >= lrs[0]) {
//                        d2[w] = 1;
//                        stk.push_back(w);
//                    }
//                }
//            }
//            tmp += stk.size();
//
//            maxTwoHopDegreeV = std::max(maxTwoHopDegreeV, tmp);
//            keys[u] = tmp;
//
//            for (auto w: stk) d2[w] = 0;
//        }
//
//        Heap_t rheap(n2, maxTwoHopDegreeV + 1);
//        rheap.init(n2, maxTwoHopDegreeV + 1, ids.data(), keys.data());
//
//        for (int32_t i = 0; i < n2; i++) {
//            int32_t u, degU;
//
//            if (!rheap.pop_min(u, degU)) printf("errorLheap\n");
//            old_lables[1][l] = u;
//            labelsR[u] = l++;
//            core[1] = std::max(core[1], degU);
//            cores[1][l - 1] = core[1];
//            if (core[1] == 0) continue;
//            stk.clear();
//            for (int32_t i = pV[u]; i < pV[u + 1]; i++) {
//                int32_t v = e2[i];
//                if (cdeg[v] + k < lrs[1]) continue;
//                for (int32_t j = pU[v]; j < pU[v + 1]; j++) {
//                    int32_t w = e1[j];
//                    if (d2[w] == 0 && cdeg[w + n1] + k >= lrs[0]) {
//                        d2[w] = 1;
//                        stk.push_back(w);
//                        heap.decrement(w, 1);
//                    }
//                }
//            }
//
//            for (auto w: stk) d2[w] = 0;
//        }
//        for (int32_t i = 0; i < n1; ++i) {
//            int v = labelsL[i];
//            cores[0][v] = cdeg[i];
//        }
//        for (int32_t i = 0; i < n2; ++i) {
//            int v = labelsR[i];
//            cores[1][v] = cdeg[i + n1];
//        }
//
//        for (int32_t i = 0; i < m; i++) {
//            edges[i].u = labelsL[edges[i].u];
//            edges[i].v = labelsR[edges[i].v];
//        }
//
//        std::fill(d1.begin(), d1.begin() + n1, 0);
//        std::fill(d2.begin(), d2.begin() + n2, 0);
//        std::fill(pU.begin(), pU.begin() + n1 + 1, 0);
//        std::fill(pV.begin(), pV.begin() + n2 + 1, 0);
//
//        for (int32_t i = 0; i < m; i++) {
//            ++d1[edges[i].u];
//            ++d2[edges[i].v];
//        }
//
//        for (int32_t i = 0; i < n1; i++) {
//            pU[i + 1] = pU[i] + d1[i];
//        }
//        for (int32_t i = 0; i < n2; i++) {
//            pV[i + 1] = pV[i] + d2[i];
//        }
//
//        for (int32_t i = 0; i < m; i++) {
//            e1[pU[edges[i].u]++] = edges[i].v;
//        }
//        for (int32_t i = 0; i < m; i++) {
//            e2[pV[edges[i].v]++] = edges[i].u;
//        }
//
//        pU[0] = pV[0] = 0;
//        for (int32_t i = 0; i < n1; i++) {
//            pU[i + 1] = pU[i] + d1[i];
//        }
//        for (int32_t i = 0; i < n2; i++) {
//            pV[i + 1] = pV[i] + d2[i];
//        }
//
//        for (int32_t i = 0; i < n1; i++) {
//            std::sort(e1.begin() + pU[i], e1.begin() + pU[i + 1]);
//        }
//        for (int32_t i = 0; i < n2; i++) {
//            std::sort(e2.begin() + pV[i], e2.begin() + pV[i + 1]);
//        }
//    }
//
//    void coreReduction(int p, int q) {
//        std::queue <int32_t> qL, qR;
//        std::vector <int32_t> d1(n1 + 1), d2(n2 + 1);
//        std::vector <int32_t> labelsL(n1 + 1), labelsR(n2 + 1);
//        std::vector<bool> visL(n1 + 1), visR(n2 + 1);
//
//        for (int32_t i = 0; i < n1; i++) {
//            d1[i] = deg1(i);
//            if (deg1(i) < q) {
//                qL.push(i);
//                visL[i] = true;
//            }
//        }
//        for (int32_t i = 0; i < n2; i++) {
//            d2[i] = deg2(i);
//            if (deg2(i) < p) {
//                qR.push(i);
//                visR[i] = true;
//            }
//        }
//
//        while (!qL.empty() || !qR.empty()) {
//            while (!qL.empty()) {
//                int32_t u = qL.front();
//                qL.pop();
//
//                for (int32_t i = 0; i < d1[u]; i++) {
//                    int32_t v = e1[pU[u] + i];
//                    // if(d2[v] < q) continue;
//
//                    for (int32_t j = pV[v]; j < pV[v] + d2[v]; j++) {
//                        if (e2[j] == u) {
//                            --d2[v];
//                            std::swap(e2[j], e2[pV[v] + d2[v]]);
//
//                            if (d2[v] == p - 1 && !visR[v]) {
//                                qR.push(v);
//                                visR[v] = true;
//                            }
//                            break;
//                        }
//                    }
//                }
//            }
//
//            while (!qR.empty()) {
//                int32_t v = qR.front();
//                qR.pop();
//
//                for (int32_t i = 0; i < d2[v]; i++) {
//                    int32_t u = e2[pV[v] + i];
//                    // if(d1[u] < p) continue;
//
//                    for (int32_t j = pU[u]; j < pU[u] + d1[u]; j++) {
//                        if (e1[j] == v) {
//                            --d1[u];
//                            std::swap(e1[j], e1[pU[u] + d1[u]]);
//
//                            if (d1[u] == q - 1 && !visL[u]) {
//                                qL.push(u);
//                                visL[u] = true;
//                            }
//                            break;
//                        }
//                    }
//                }
//            }
//        }
//
//        int32_t pL = 1, pR = 1;
//        for (int32_t u = 0; u < n1; u++) {
//            if (!visL[u]) labelsL[u] = pL++;
//        }
//        for (int32_t v = 0; v < n2; v++) {
//            if (!visR[v]) labelsR[v] = pR++;
//        }
//
//        int32_t pm = 0;
//        for (int32_t u = 0; u < n1; u++) {
//            if (visL[u]) continue;
//            for (int32_t i = pU[u]; i < pU[u + 1]; i++) {
//                int32_t v = e1[i];
//                if (!visR[v]) {
//                    edges[pm].u = labelsL[u] - 1;
//                    edges[pm].v = labelsR[v] - 1;
//                    ++pm;
//                }
//            }
//        }
//        m = pm;
//
//        n1 = pL - 1;
//        n2 = pR - 1;
//
//        std::fill(pU.begin(), pU.begin() + n1 + 1, 0);
//        std::fill(pV.begin(), pV.begin() + n2 + 1, 0);
//
//        std::fill(d1.begin(), d1.begin() + n1 + 1, 0);
//        std::fill(d2.begin(), d2.begin() + n2 + 1, 0);
//        changeToDegreeOrder();
//
//    }

    void VertexIDOrder() {
        printf("here VertexIDOrder\n");
        fflush(stdout);
        std::vector <int32_t> d1, d2;

        d1.resize(n1);
        d2.resize(n2);

        for (int32_t i = 0; i < m; i++) {
            ++d1[edges[i].u];
            ++d2[edges[i].v];
        }

        maxDu = 0;
        for (int32_t i = 0; i < n1; i++) {
            maxDu = std::max(maxDu, d1[i]);
        }
        maxDv = 0;
        for (int32_t i = 0; i < n2; i++) {
            maxDv = std::max(maxDv, d2[i]);
        }

        pU[0] = 0;
        for (int32_t u = 0; u < n1; u++) {
            pU[u + 1] = d1[u] + pU[u];
        }
        for (int32_t i = 0; i < m; i++) {
            e1[pU[edges[i].u]++] = edges[i].v;
        }
        pU[0] = 0;
        for (int32_t u = 0; u < n1; u++) {
            pU[u + 1] = d1[u] + pU[u];
        }

        pV[0] = 0;
        for (int32_t v = 0; v < n2; v++) {
            pV[v + 1] = d2[v] + pV[v];
        }
        for (int32_t i = 0; i < m; i++) {
            e2[pV[edges[i].v]++] = edges[i].u;
        }
        pV[0] = 0;
        for (int32_t v = 0; v < n2; v++) {
            pV[v + 1] = d2[v] + pV[v];
        }

        int32_t maxTwoHopDegreeU = 0, maxTwoHopDegreeV = 0;
        d1.clear();
        std::vector <int32_t> labelsL(n1);
        std::vector <int32_t> labelsR(n2);


        int32_t N = n1 + n2;
        std::vector <int32_t> dequeue;
        std::vector <int32_t> cdeg;
        dequeue.resize(N);
        cdeg.resize(N, 0);
        int32_t tsize = 0;
        n[0] = n1;
        n[1] = n2;
        for (int32_t tt = 0; tt <= 1; tt++) {
            for (int32_t i = 0; i < n[tt]; i++) {
                int32_t d = 0;
                if (tt == 0) {
                    d = d1[i];
                    cdeg[i] = d;
                    if (d + k < lrs[1])
                        dequeue[tsize++] = i;
                } else {
                    d = d2[i];
                    cdeg[i + n1] = d;
                    if (d + k < lrs[0])
                        dequeue[tsize++] = i + n1;
                }
            }
        }

        int32_t rsize = 0;
        while (tsize > rsize) {
            int32_t s = rsize;
            rsize = tsize;
            for (int32_t i = s; i < rsize; ++i) {
                int v = dequeue[i];
                if (v >= n1) {
                    v -= n1;
                    for (int32_t j = pV[v]; j < pV[v + 1]; ++j) {
                        int32_t u = e2[j];
                        if (cdeg[u] + k >= lrs[1]) {
                            cdeg[u]--;
                            if (cdeg[u] + k < lrs[1])
                                dequeue[tsize++] = u;
                        }
                    }
                } else {
                    for (int32_t j = pU[v]; j < pU[v + 1]; ++j) {
                        int32_t u = e1[j];
                        u += n1;
                        if (cdeg[u] + k >= lrs[0]) {
                            cdeg[u]--;
                            if (cdeg[u] + k < lrs[0])
                                dequeue[tsize++] = u;
                        }
                    }
                }
            }
        }

        std::vector <int32_t> sn1, sn2;
        sn1.reserve(n1);
        sn2.reserve(n2);
        for (int32_t i = 0; i < n1; ++i) sn1.push_back(i);
        for (int32_t i = 0; i < n2; ++i) sn2.push_back(i);
        sort(sn1.begin(), sn1.end(), [&](int32_t x, int32_t y) -> bool { return cdeg[x] < cdeg[y]; });
        sort(sn2.begin(), sn2.end(), [&](int32_t x, int32_t y) -> bool { return cdeg[x + n1] < cdeg[y + n1]; });
        for (int32_t i = 0; i < n1; ++i) {
            int32_t v = sn1[i];
            labelsL[v] = i;
            cores[0][i] = cdeg[v];
        }
        for (int32_t i = 0; i < n2; ++i) {
            int32_t v = sn2[i];
            labelsR[v] = i;
            cores[1][i] = cdeg[v + n1];
        }


        for (int32_t i = 0; i < m; i++) {
            edges[i].u = labelsL[edges[i].u];
            edges[i].v = labelsR[edges[i].v];
        }

        std::fill(d1.begin(), d1.begin() + n1, 0);
        std::fill(d2.begin(), d2.begin() + n2, 0);
        std::fill(pU.begin(), pU.begin() + n1 + 1, 0);
        std::fill(pV.begin(), pV.begin() + n2 + 1, 0);

        for (int32_t i = 0; i < m; i++) {
            ++d1[edges[i].u];
            ++d2[edges[i].v];
        }

        for (int32_t i = 0; i < n1; i++) {
            pU[i + 1] = pU[i] + d1[i];
        }
        for (int32_t i = 0; i < n2; i++) {
            pV[i + 1] = pV[i] + d2[i];
        }

        for (int32_t i = 0; i < m; i++) {
            e1[pU[edges[i].u]++] = edges[i].v;
        }
        for (int32_t i = 0; i < m; i++) {
            e2[pV[edges[i].v]++] = edges[i].u;
        }

        pU[0] = pV[0] = 0;
        for (int32_t i = 0; i < n1; i++) {
            pU[i + 1] = pU[i] + d1[i];
        }
        for (int32_t i = 0; i < n2; i++) {
            pV[i + 1] = pV[i] + d2[i];
        }

        for (int32_t i = 0; i < n1; i++) {
            std::sort(e1.begin() + pU[i], e1.begin() + pU[i + 1]);
        }
        for (int32_t i = 0; i < n2; i++) {
            std::sort(e2.begin() + pV[i], e2.begin() + pV[i + 1]);
        }
    }


    void changeToCoreOrderVersion2() {
        std::vector <int32_t> d1, d2;

        d1.resize(n1);
        d2.resize(n2);

        for (int32_t i = 0; i < m; i++) {
            ++d1[edges[i].u];
            ++d2[edges[i].v];
        }

        maxDu = 0;
        for (int32_t i = 0; i < n1; i++) {
            maxDu = std::max(maxDu, d1[i]);
        }
        maxDv = 0;
        for (int32_t i = 0; i < n2; i++) {
            maxDv = std::max(maxDv, d2[i]);
        }

        pU[0] = 0;
        for (int32_t u = 0; u < n1; u++) {
            pU[u + 1] = d1[u] + pU[u];
        }
        for (int32_t i = 0; i < m; i++) {
            e1[pU[edges[i].u]++] = edges[i].v;
        }
        pU[0] = 0;
        for (int32_t u = 0; u < n1; u++) {
            pU[u + 1] = d1[u] + pU[u];
        }

        pV[0] = 0;
        for (int32_t v = 0; v < n2; v++) {
            pV[v + 1] = d2[v] + pV[v];
        }
        for (int32_t i = 0; i < m; i++) {
            e2[pV[edges[i].v]++] = edges[i].u;
        }
        pV[0] = 0;
        for (int32_t v = 0; v < n2; v++) {
            pV[v + 1] = d2[v] + pV[v];
        }

        Heap_t l_heap(n1, maxDu + 1), r_heap(n2, maxDv + 1);
        int32_t n = std::max(n1, n2);
        std::vector <int32_t> ids(n);
        std::vector <int32_t> keys(n);
        std::vector <int32_t> labelsL(n1);
        std::vector <int32_t> labelsR(n2);
        int32_t l1 = 0, l2 = 0;

        for (int32_t i = 0; i < n1; i++) {
            ids[i] = i;
            keys[i] = d1[i] + 1;
        }
        l_heap.init(n1, maxDu + 1, ids.data(), keys.data());
        for (int32_t i = 0; i < n2; i++) {
            ids[i] = i;
            keys[i] = d2[i] + 1;
        }
        r_heap.init(n2, maxDv + 1, ids.data(), keys.data());

        core[0] = 0;
        core[1] = 0;
        for (int32_t i = 0; i < n1 + n2; i++) {
            int32_t u, degU = n2 + 11;
            int32_t v, degV = n1 + 11;

            if (!l_heap.empty() && !l_heap.pop_min(u, degU)) printf("error l_heap\n");

            if (!r_heap.empty() && !r_heap.pop_min(v, degV)) printf("error r_heap\n");

            if (degU <= degV) {
                if (degV != n1 + 11) r_heap.insert(v, degV);

                for (int32_t j = pU[u]; j < pU[u + 1]; j++) r_heap.decrement(e1[j]);

                old_lables[0][l1] = u;
                labelsL[u] = l1++;
                core[0] = std::max(core[0], degU);
                cores[0][l1 - 1] = core[0];
            } else {
                if (degU != n2 + 11) l_heap.insert(u, degU);

                for (int32_t j = pV[v]; j < pV[v + 1]; j++) l_heap.decrement(e2[j]);

                old_lables[1][l2] = v;
                labelsR[v] = l2++;
                core[1] = std::max(core[1], degV);
                cores[1][l2 - 1] = core[1];
            }
        }

        for (int32_t i = 0; i < m; i++) {
            edges[i].u = labelsL[edges[i].u];
            edges[i].v = labelsR[edges[i].v];
        }

        std::fill(d1.begin(), d1.begin() + n1, 0);
        std::fill(d2.begin(), d2.begin() + n2, 0);
        std::fill(pU.begin(), pU.begin() + n1 + 1, 0);
        std::fill(pV.begin(), pV.begin() + n2 + 1, 0);

        for (int32_t i = 0; i < m; i++) {
            ++d1[edges[i].u];
            ++d2[edges[i].v];
        }

        for (int32_t i = 0; i < n1; i++) {
            pU[i + 1] = pU[i] + d1[i];
        }
        for (int32_t i = 0; i < n2; i++) {
            pV[i + 1] = pV[i] + d2[i];
        }

        for (int32_t i = 0; i < m; i++) {
            e1[pU[edges[i].u]++] = edges[i].v;
        }
        for (int32_t i = 0; i < m; i++) {
            e2[pV[edges[i].v]++] = edges[i].u;
        }

        pU[0] = pV[0] = 0;
        for (int32_t i = 0; i < n1; i++) {
            pU[i + 1] = pU[i] + d1[i];
        }
        for (int32_t i = 0; i < n2; i++) {
            pV[i + 1] = pV[i] + d2[i];
        }

        for (int32_t i = 0; i < n1; i++) {
            std::sort(e1.begin() + pU[i], e1.begin() + pU[i + 1]);
        }
        for (int32_t i = 0; i < n2; i++) {
            std::sort(e2.begin() + pV[i], e2.begin() + pV[i + 1]);
        }

    }
//
//    void changeToCoreOrder() {
//        std::vector <int32_t> d1, d2;
//
//        d1.resize(n1);
//        d2.resize(n2);
//
//        for (int32_t i = 0; i < m; i++) {
//            ++d1[edges[i].u];
//            ++d2[edges[i].v];
//        }
//
//        maxDu = 0;
//        for (int32_t i = 0; i < n1; i++) {
//            maxDu = std::max(maxDu, d1[i]);
//        }
//        maxDv = 0;
//        for (int32_t i = 0; i < n2; i++) {
//            maxDv = std::max(maxDv, d2[i]);
//        }
//
//        pU[0] = 0;
//        for (int32_t u = 0; u < n1; u++) {
//            pU[u + 1] = d1[u] + pU[u];
//        }
//        for (int32_t i = 0; i < m; i++) {
//            e1[pU[edges[i].u]++] = edges[i].v;
//        }
//        pU[0] = 0;
//        for (int32_t u = 0; u < n1; u++) {
//            pU[u + 1] = d1[u] + pU[u];
//        }
//
//        pV[0] = 0;
//        for (int32_t v = 0; v < n2; v++) {
//            pV[v + 1] = d2[v] + pV[v];
//        }
//        for (int32_t i = 0; i < m; i++) {
//            e2[pV[edges[i].v]++] = edges[i].u;
//        }
//        pV[0] = 0;
//        for (int32_t v = 0; v < n2; v++) {
//            pV[v + 1] = d2[v] + pV[v];
//        }
//
//        Heap_t l_heap(n1, maxDu + 1), r_heap(n2, maxDv + 1);
//        int32_t n = std::max(n1, n2);
//        std::vector <int32_t> ids(n);
//        std::vector <int32_t> keys(n);
//        std::vector <int32_t> labelsL(n1);
//        std::vector <int32_t> labelsR(n2);
//
//        for (int32_t i = 0; i < n1; i++) {
//            ids[i] = i;
//            keys[i] = d1[i] + 1;
//        }
//        l_heap.init(n1, maxDu + 1, ids.data(), keys.data());
//
//        for (int32_t i = 0; i < n2; i++) {
//            ids[i] = i;
//            keys[i] = d2[i] + 1;
//        }
//        r_heap.init(n2, maxDv + 1, ids.data(), keys.data());
//
//        int32_t minN = std::min(n1, n2);
//
//        for (int32_t i = 0; i < minN; i++) {
//            int32_t u, degU;
//            int32_t v, degV;
//
//            if (!l_heap.pop_min(u, degU)) printf("error l_heap\n");
//            if (!r_heap.pop_min(v, degV)) printf("error r_heap\n");
//
//            labelsL[u] = i;
//            labelsR[v] = i;
//            for (int32_t j = pU[u]; j < pU[u + 1]; j++) r_heap.decrement(e1[j]);
//
//            for (int32_t j = pV[v]; j < pV[v + 1]; j++) l_heap.decrement(e2[j]);
//
//        }
//        if (n1 < n2) {
//            for (int32_t j = n1; j < n2; j++) {
//                int32_t v, degV;
//                if (!r_heap.pop_min(v, degV)) printf("error r_heap\n");
//                labelsR[v] = j;
//            }
//        }
//        else {
//            for (int32_t j = n2; j < n1; j++) {
//                int32_t u, degU;
//                if (!l_heap.pop_min(u, degU)) printf("error l_heap\n");
//                labelsL[u] = j;
//            }
//        }
//
//        for (int32_t i = 0; i < m; i++) {
//            edges[i].u = labelsL[edges[i].u];
//            edges[i].v = labelsR[edges[i].v];
//        }
//
//        std::fill(d1.begin(), d1.begin() + n1, 0);
//        std::fill(d2.begin(), d2.begin() + n2, 0);
//        std::fill(pU.begin(), pU.begin() + n1 + 1, 0);
//        std::fill(pV.begin(), pV.begin() + n2 + 1, 0);
//
//        for (int32_t i = 0; i < m; i++) {
//            ++d1[edges[i].u];
//            ++d2[edges[i].v];
//        }
//
//        for (int32_t i = 0; i < n1; i++) {
//            pU[i + 1] = pU[i] + d1[i];
//        }
//        for (int32_t i = 0; i < n2; i++) {
//            pV[i + 1] = pV[i] + d2[i];
//        }
//
//        for (int32_t i = 0; i < m; i++) {
//            e1[pU[edges[i].u]++] = edges[i].v;
//        }
//        for (int32_t i = 0; i < m; i++) {
//            e2[pV[edges[i].v]++] = edges[i].u;
//        }
//
//        pU[0] = pV[0] = 0;
//        for (int32_t i = 0; i < n1; i++) {
//            pU[i + 1] = pU[i] + d1[i];
//        }
//        for (int32_t i = 0; i < n2; i++) {
//            pV[i + 1] = pV[i] + d2[i];
//        }
//
//        for (int32_t i = 0; i < n1; i++) {
//            std::sort(e1.begin() + pU[i], e1.begin() + pU[i + 1]);
//        }
//        for (int32_t i = 0; i < n2; i++) {
//            std::sort(e2.begin() + pV[i], e2.begin() + pV[i + 1]);
//        }
//    }
//
//    void rawOrder() {
//        std::vector <int32_t> d1, d2;
//
//        d1.resize(n1);
//        d2.resize(n2);
//
//        for (int32_t i = 0; i < m; i++) {
//            ++d1[edges[i].u];
//            ++d2[edges[i].v];
//        }
//
//        maxDu = 0;
//        for (int32_t i = 0; i < n1; i++) {
//            maxDu = std::max(maxDu, d1[i]);
//        }
//        maxDv = 0;
//        for (int32_t i = 0; i < n2; i++) {
//            maxDv = std::max(maxDv, d2[i]);
//        }
//
//        pU[0] = 0;
//        for (int32_t u = 0; u < n1; u++) {
//            pU[u + 1] = d1[u] + pU[u];
//        }
//        for (int32_t i = 0; i < m; i++) {
//            e1[pU[edges[i].u]++] = edges[i].v;
//        }
//        pU[0] = 0;
//        for (int32_t u = 0; u < n1; u++) {
//            pU[u + 1] = d1[u] + pU[u];
//        }
//
//        pV[0] = 0;
//        for (int32_t v = 0; v < n2; v++) {
//            pV[v + 1] = d2[v] + pV[v];
//        }
//        for (int32_t i = 0; i < m; i++) {
//            e2[pV[edges[i].v]++] = edges[i].u;
//        }
//        pV[0] = 0;
//        for (int32_t v = 0; v < n2; v++) {
//            pV[v + 1] = d2[v] + pV[v];
//        }
//
//        for (int32_t i = 0; i < n1; i++) {
//            std::sort(e1.begin() + pU[i], e1.begin() + pU[i + 1]);
//        }
//        for (int32_t i = 0; i < n2; i++) {
//            std::sort(e2.begin() + pV[i], e2.begin() + pV[i + 1]);
//        }
//    }
//
//    void changeToDegreeOrder() {
//        std::vector <int32_t> d1, d2;
//        std::vector <int32_t> label1, label2;
//
//        d1.resize(std::max(n1, n2) + 1);
//        d2.resize(std::max(n1, n2) + 1);
//        label1.resize(n1);
//        label2.resize(n2);
//
//        for (int32_t i = 0; i < m; i++) {
//            ++d1[edges[i].u];
//            ++d2[edges[i].v];
//        }
//
//        maxDu = 0;
//        for (int32_t i = 0; i < n1; i++) {
//            maxDu = std::max(maxDu, d1[i]);
//        }
//
//        maxDv = 0;
//        for (int32_t i = 0; i < n2; i++) {
//            maxDv = std::max(maxDv, d2[i]);
//        }
//
//        for (int32_t i = 0; i < n1; i++) {
//            pV[d1[i] + 1]++;
//        }
//        for (int32_t i = 0; i < maxDu; i++) {
//            pV[i + 1] += pV[i];
//        }
//        for (int32_t i = 0; i < n1; i++) {
//            label1[i] = pV[d1[i]]++;
//        }
//        for (int32_t i = 0; i < n2; i++) {
//            pU[d2[i] + 1]++;
//        }
//        for (int32_t i = 0; i < maxDv; i++) {
//            pU[i + 1] += pU[i];
//        }
//
//        for (int32_t i = 0; i < n2; i++) {
//            label2[i] = pU[d2[i]]++;
//        }
//
//        for (int32_t i = 0; i < m; i++) {
//            edges[i].u = label1[edges[i].u];
//            edges[i].v = label2[edges[i].v];
//        }
//
//        std::fill(d1.begin(), d1.begin() + std::max(n1, n2) + 1, 0);
//        std::fill(d2.begin(), d2.begin() + std::max(n1, n2) + 1, 0);
//        std::fill(pU.begin(), pU.begin() + n1 + 1, 0);
//        std::fill(pV.begin(), pV.begin() + n2 + 1, 0);
//
//        for (int32_t i = 0; i < m; i++) {
//            ++d1[edges[i].u];
//            ++d2[edges[i].v];
//        }
//
//        for (int32_t i = 0; i < n1; i++) {
//            pU[i + 1] = pU[i] + d1[i];
//        }
//        for (int32_t i = 0; i < n2; i++) {
//            pV[i + 1] = pV[i] + d2[i];
//        }
//
//        for (int32_t i = 0; i < m; i++) {
//            e1[pU[edges[i].u]++] = edges[i].v;
//        }
//        for (int32_t i = 0; i < m; i++) {
//            e2[pV[edges[i].v]++] = edges[i].u;
//        }
//        pU[0] = pV[0] = 0;
//        for (int32_t i = 0; i < n1; i++) {
//            pU[i + 1] = pU[i] + d1[i];
//        }
//        for (int32_t i = 0; i < n2; i++) {
//            pV[i + 1] = pV[i] + d2[i];
//        }
//
//        for (int32_t i = 0; i < n1; i++) {
//            std::sort(e1.begin() + pU[i], e1.begin() + pU[i + 1]);
//        }
//        for (int32_t i = 0; i < n2; i++) {
//            std::sort(e2.begin() + pV[i], e2.begin() + pV[i + 1]);
//        }
//    }

    void swapUV() {
        std::swap(n[0], n[1]);
        std::swap(maxDu, maxDv);

        p[0].swap(p[1]);
        e[0].swap(e[1]);
    }

    int32_t deg1(int32_t u) {
        return p[0][u + 1] - p[0][u];
    }

    int32_t deg2(int32_t v) {
        return p[1][v + 1] - p[1][v];
    }

    bool connect(int32_t u, int32_t v, int32_t t) {
        return cuhash[t][u].find(v);
    }

    void print() {
        printf("U:\n");
        for (int32_t u = 0; u < n[0]; u++) {
            printf("%u:", u);
            for (int32_t i = p[0][u]; i < p[0][u + 1]; i++) {
                printf("%u ", e[0][i]);
            }
            printf("\n");
        }
        printf("V:\n");
        for (int32_t v = 0; v < n[1]; v++) {
            printf("%u:", v);
            for (int32_t i = p[1][v]; i < p[1][v + 1]; i++) {
                printf("%u ", e[1][i]);
            }
            printf("\n");
        }
    }

    int32_t deg(int32_t u, int32_t t) {
        return p[t][u + 1] - p[t][u];
    }

};


#endif