#pragma once

#include "../graph/graph.h"
#include "../tools/types.hpp"
#include "../tools/linearSet.hpp"
#include <iosfwd>
#include <string>

class MCE {
private:
    Graph g;
    ui cnt = 0;
    ui threshold = 0;

    LinearSet V;
    void BKTomitaRec(ui deep, std::vector<ui> & C, ui P, ui SX, ui EX);

    std::vector<std::vector<ui>> adj;
    std::vector<std::vector<ui>> edAdj, stAdj;
    std::vector<ui> reId;
    std::vector<bool> inC;
    //void BKTomitaRecCCR(ui deep, std::vector<ui> & C, ui SC, ui EC, ui SP, ui SX, ui EX);

   // void BKTomitaRecCCRV2(ui deep, std::vector<ui> & C, ui SC, ui EC, ui SP, ui SX, ui EX);

    LinearSet Cv3, Rv3, Pv3, Xv3;
    ui numOfNonCNodes = 0;
    ui currentRoot = 0;
    std::vector<std::vector<ui>> allCliques;

    void reportCurrentClique(ui EC, ui ER);
    void BKTomitaRecCCRV3(ui deep, ui ER, ui EC, ui EP, ui SX, ui EX);

    //void BKTomitaRecCCRV4(ui deep, ui ER, ui EC, ui EP, ui SX, ui EX);

   // void BKTomitaRecCCRV5(ui deep, ui ER, ui EC, ui EP, ui SX, ui EX);
#define COUNT_NUM_SEARCH_NODES
#ifdef COUNT_NUM_SEARCH_NODES
    ull numOfDFSSearchNodes = 0;
#endif

public:
    MCE(Graph && g) : g(std::move(g)) {
        V.resize(this->g.n);
    }

    // Every reported clique is retained, even when it is not printed. Vertex
    // labels are restored to the labels used in the input graph.
    const std::vector<std::vector<ui>> & cliques() const { return allCliques; }
    void saveCliques(const std::string & path) const;
    void printStoredCliques(std::ostream & output) const;

    ui run();

    //ui runCoreCliqueRemoval();

   // ui runCoreCliqueRemovalV2();

    ui runCoreCliqueRemovalV3();

    //ui runCoreCliqueRemovalV4();
    
    //ui runCoreCliqueRemovalV5();
};
