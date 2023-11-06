#include <Rcpp.h>
using namespace Rcpp;
#include "GraphUtils.h"


Graph::Graph(int n) : n(n), adj(n), weights(n) {}

void Graph::addEdge(int u, int v, double w) {
    adj[u].push_back(v); // add v to u's list of neighbors
    weights[u].push_back(w); // add weight of edge from u to v
}


int countEdges(const Graph& C, const Graph& D) {
    int E = 0;
    std::unordered_map<int, std::set<int>> neighborsC;  
    std::unordered_map<int, std::vector<double>> weightMapC;  
    std::unordered_map<int, std::set<int>> neighborsD;
    
    for (int i = 0; i < C.n; i++) { 
        for (std::vector<int>::size_type j = 0; j < C.adj[i].size(); j++) { 
            int neighbor = C.adj[i][j];
            double weight = C.weights[i][j];
            neighborsC[i].insert(neighbor);  
            weightMapC[i].push_back(weight);  
        }
    }
    for (int i = 0; i < D.n; i++) {  
        for (std::vector<int>::size_type j = 0; j < D.adj[i].size(); j++) {  
            int neighbor = D.adj[i][j];  
            neighborsD[i].insert(neighbor); 
        }
    }
    for (int i = 0; i < C.n; i++) {  
        for (int neighbor : neighborsC[i]) {  
            if (neighborsD[i].count(neighbor) > 0) {  
                auto& weightsC = weightMapC[i];  
                for (std::vector<int>::size_type j = 0; j < C.adj[i].size(); j++) {  
                    if (C.adj[i][j] == neighbor) {  
                    E += weightsC[j];  
                    }
                }
            }
        }
    }

    return E;
}


Rcpp::List graphToRList(const Graph& G) {
    int n = G.n;
    
    Rcpp::NumericVector from;
    Rcpp::NumericVector to;
    Rcpp::NumericVector weight;

    for (int i = 0; i < n; i++) {
        for (std::vector<int>::size_type j = 0; j < G.adj[i].size(); j++) {
            from.push_back(i + 1); // Adding 1 to convert to 1-based indexing
            to.push_back(G.adj[i][j] + 1); // Adding 1 to convert to 1-based indexing
            weight.push_back(G.weights[i][j]);
        }
    }

    // Create a named list with 'from', 'to', and 'weight' components
    Rcpp::List result = Rcpp::List::create(
        Rcpp::Named("from") = from,
        Rcpp::Named("to") = to,
        Rcpp::Named("weight") = weight
    );

    return result;
}

// [[Rcpp::export]]
Rcpp::List createGraph(const Rcpp::NumericMatrix& edgeWeights) {

    Graph G(edgeWeights.nrow());
    for (int i = 0; i < edgeWeights.nrow(); i++) {
        for (int j = 0; j < edgeWeights.ncol(); j++) {
            if (edgeWeights(i, j) != 0) {
                G.addEdge(i, j, edgeWeights(i, j));
            }
        }
    }
    return graphToRList(G);
}
