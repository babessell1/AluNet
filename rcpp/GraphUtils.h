#ifndef GRAPHUTILS_H
#define GRAPHUTILS_H

#include <Rcpp.h>
#include <vector>
#include <set>
#include <unordered_map>

class Graph {
public:
    int n;
    std::map<std::string, int> nodeIndexMap;
    std::vector<std::vector<int>> adj;
    std::vector<std::vector<double>> weights;

    Graph(int n);
    void addEdge(const std::string& u, const std::string& v, double w);
    int getNodeIndex(const std::string& node) const;
    std::string getNodeName(int index) const;
    std::vector<int> getNodes() const;
};

Rcpp::List graphToRList(const Graph& G);
// listToGraph
Graph listToGraph(const Rcpp::List& graphList);


#endif