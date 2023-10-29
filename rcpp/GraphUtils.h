#ifndef GRAPHUTILS_H
#define GRAPHUTILS_H

#include <Rcpp.h>
#include <vector>
#include <set>
#include <unordered_map>

class Graph {
public:
    int n;
    std::vector<std::vector<int>> adj;
    std::vector<std::vector<double>> weights;

    Graph(int n) : n(n), adj(n), weights(n) {}

    void addEdge(int u, int v, double w) {
        adj[u].push_back(v);
        weights[u].push_back(w);
    }
};


int countEdges(const Graph& C, const Graph& D) {
    // edge count number;
    int E = 0;
    
    
    std::unordered_map<int, std::set<int>> neighborsC;  // map from node ID to set of neighbors
    std::unordered_map<int, std::vector<double>> weightMapC;  // map from node ID to vector of edge weights
    std::unordered_map<int, std::set<int>> neighborsD;  // map from node ID to set of neighbors

    
    // Populate the maps
    for (int i = 0; i < C.n; i++) { // for each node, i in C
        for (int j = 0; j < C.adj[i].size(); j++) { // for each neighbor, j of node i
            int neighbor = C.adj[i][j];
            double weight = C.weights[i][j];
            neighborsC[i].insert(neighbor);  // add neighbor to set of neighbors
            weightMapC[i].push_back(weight);  // add weight of edge from i to neighbor
        }
    }

    for (int i = 0; i < D.n; i++) {  // for each node, i in D
        for (int j = 0; j < D.adj[i].size(); j++) {  // for each neighbor, j of node i
            int neighbor = D.adj[i][j];  // get neighbor
            neighborsD[i].insert(neighbor);  // add neighbor to set of neighbors
        }
    }

    for (int i = 0; i < C.n; i++) {  // for each node, i in C
        for (int neighbor : neighborsC[i]) {  // for each neighbor of node i
            if (neighborsD[i].count(neighbor) > 0) {  // if neighbor is in D
                auto& weightsC = weightMapC[i];  // get weights of edges in C
                for (int j = 0; j < C.adj[i].size(); j++) {  
                    if (C.adj[i][j] == neighbor) {  // if neighbor is found in adjacency of i
                    E += weightsC[j];  // add weight of edge from i to neighbor
                    }
                }
            }
        }
    }

    return E;

}

Graph createGraph(const Rcpp::NumericMatrix& edgeWeights) {
    // edge weights is an n x n matrix of edge weights connecting genomic cooridanates
    // n is the number of genomic coordinates
    // to save memory, we will consider any 0 weight edge to be non-existent
    matrix edgeWeights = as<mat>(edgeWeights);
    // remove rows and columns with all 0s
    for (int i = 0; i < edgeWeights.n_rows; i++) {
        if (sum(edgeWeights.row(i)) == 0) {
            edgeWeights.shed_row(i);
            edgeWeights.shed_col(i);
            i--;
        }
    }
    // create graph
    Graph G(edgeWeights.n_rows);
    for (int i = 0; i < edgeWeights.n_rows; i++) {
        for (int j = 0; j < edgeWeights.n_cols; j++) {
            if (edgeWeights(i, j) != 0) {
                G.addEdge(i, j, edgeWeights(i, j));
            }
        }
    }
    return G;
}

#endif // GRAPHUTILS_H