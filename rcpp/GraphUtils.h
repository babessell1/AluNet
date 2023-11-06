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

    Graph(int n);
    void addEdge(int u, int v, double w);
};

int countEdges(const Graph& C, const Graph& D);
Rcpp::List graphToRList(const Graph& G);
Rcpp::List createGraph(const Rcpp::NumericMatrix& edgeWeights);

#endif