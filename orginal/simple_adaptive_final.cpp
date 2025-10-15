#include "inc/common.h"
#include "inc/graph.h"

class SimpleAdaptiveEnumeration {
private:
    ui n;
    vector<vector<ui>> adjList;
    ui cliqueCount;
    ui maxCliqueSize;
    
    // Skip-mask and reordering components
    vector<bool> skip_mask;     // marks vertices already covered
    vector<ui> global_order;    // dynamic vertex processing order
    
public:
    SimpleAdaptiveEnumeration(Graph &g) {
        n = g.n;
        adjList.resize(n);
        cliqueCount = 0;
        maxCliqueSize = 0;
        skip_mask.resize(n, false);
        
        // Fill adjacency lists
        for (ui i = 0; i < n; i++) {
            for (ui j = g.offset[i]; j < g.offset[i + 1]; j++) {
                ui neighbor = g.neighbors[j];
                if (neighbor < n) {
                    adjList[i].push_back(neighbor);
                }
            }
            sort(adjList[i].begin(), adjList[i].end());
        }
        
        // Initialize global order
        global_order.resize(n);
        iota(global_order.begin(), global_order.end(), 0);
    }
    
    bool isConnected(ui u, ui v) const {
        return binary_search(adjList[u].begin(), adjList[u].end(), v);
    }
    
    bool isClique(const vector<ui> &R) const {
        for (size_t i = 0; i < R.size(); i++) {
            for (size_t j = i + 1; j < R.size(); j++) {
                if (!isConnected(R[i], R[j])) return false;
            }
        }
        return true;
    }
    
    void handleClique(const vector<ui> &R) {
        cliqueCount++;
        maxCliqueSize = max(maxCliqueSize, (ui)R.size());
        
        cout << "Found maximal clique: { ";
        for (ui v : R) cout << v << " ";
        cout << "}" << endl;
        
        // Mark covered vertices
        for (ui v : R) skip_mask[v] = true;
        
        // Reorder: uncovered first, covered last
        vector<ui> uncovered, covered;
        for (ui v : global_order) {
            if (skip_mask[v]) covered.push_back(v);
            else uncovered.push_back(v);
        }
        global_order.clear();
        global_order.insert(global_order.end(), uncovered.begin(), uncovered.end());
        global_order.insert(global_order.end(), covered.begin(), covered.end());
    }
    
    void enumerate(vector<ui> &R, ui start_idx) {
        bool expanded = false;
        
        for (ui idx = start_idx; idx < global_order.size(); idx++) {
            ui v = global_order[idx];
            if (skip_mask[v]) continue; // skip redundant vertex
            
            vector<ui> R_next = R;
            R_next.push_back(v);
            
            if (isClique(R_next)) {
                expanded = true;
                enumerate(R_next, idx + 1);
            }
        }
        
        // Leaf node → maximal clique
        if (!expanded && !R.empty()) handleClique(R);
    }
    
    void findAllMaximalCliques() {
        vector<ui> R;
        cliqueCount = 0;
        maxCliqueSize = 0;
        
        cout << "Starting Simple Adaptive Enumeration..." << endl;
        enumerate(R, 0);
        
        cout << "Total Maximal Cliques Found: " << cliqueCount << endl;
        cout << "Maximum Clique Size: " << maxCliqueSize << endl;
    }
    
    ui getCliqueCount() const { return cliqueCount; }
    ui getMaxCliqueSize() const { return maxCliqueSize; }
};