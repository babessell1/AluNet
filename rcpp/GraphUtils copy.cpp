#include "C:/Users/bbessell/Miniconda3/envs/r_leiden/lib/R/library/Rcpp/include/Rcpp.h"
using namespace Rcpp;
#include <iostream>
#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>

class Graph {
    // adj is  size n vector of vectors of ints: adj[i] is vector of neighbors of node i 
    // weights is size n vector of vectors of doubles: weights[i] is vector of weights of edges from node i
    // n is number of nodes in graph
public:
    int n; // number of nodes
    std::vector<std::vector<int> > adj; // adjacency list
    std::vector<std::vector<double> > weights; // weights of edges
    
    Graph(int n) : n(n), adj(n), weights(n) {}
    
    void addEdge(int u, int v, double w) {
        adj[u].push_back(v); // add v to u's list of neighbors
        weights[u].push_back(w); // add weight of edge from u to v
    }
};

// [[Rcpp::export]]
int countEdges(const Graph& C, const Graph& D) {
    int E = 0;
    std::unordered_map<int, std::set<int>> neighborsC;  
    std::unordered_map<int, std::vector<double>> weightMapC;  
    std::unordered_map<int, std::set<int>> neighborsD;
    
    for (int i = 0; i < C.n; i++) { 
        for (int j = 0; j < C.adj[i].size(); j++) { 
            int neighbor = C.adj[i][j];
            double weight = C.weights[i][j];
            neighborsC[i].insert(neighbor);  
            weightMapC[i].push_back(weight);  
        }
    }

    for (int i = 0; i < D.n; i++) {  
        for (int j = 0; j < D.adj[i].size(); j++) {  
            int neighbor = D.adj[i][j];  
            neighborsD[i].insert(neighbor); 
        }
    }

    for (int i = 0; i < C.n; i++) {  
        for (int neighbor : neighborsC[i]) {  
            if (neighborsD[i].count(neighbor) > 0) {  
                auto& weightsC = weightMapC[i];  
                for (int j = 0; j < C.adj[i].size(); j++) {  
                    if (C.adj[i][j] == neighbor) {  
                    E += weightsC[j];  
                    }
                }
            }
        }
    }

    return E;
}


// [[Rcpp::export]]
Graph createGraph(const Rcpp::NumericMatrix& edgeWeights) {
    matrix edgeWeights = as<mat>(edgeWeights);

    for (int i = 0; i < edgeWeights.n_rows; i++) {
        if (sum(edgeWeights.row(i)) == 0) {
            edgeWeights.shed_row(i);
            edgeWeights.shed_col(i);
            i--;
        }
    }

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