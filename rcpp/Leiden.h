#include <Rcpp.h>
#include "GraphUtils.h"

class Community {
public:
    std::vector<std::string> nodes; // nodes of a graph

    Community(const std::vector<std::string>& nodes); 
    std::vector<int> getIndices();  // returns indices of nodes in community
};

class Partition {
public:
    std::vector<Community> communities;  // communities

    Partition(const std::vector<int>& C);
    std::vector<int> getIndices();  // returns indices of communities in partition
};

class Optimizer {
public:
    Graph G; // graph
    int iter; // number of iterations
    std::vector<std::vector<int>> P; // partition

    Optimizer(const Graph& G, int iterations);
    Partition optimize();  // Leiden iteration, returns new partition
    Partition moveNodesFast();  // Leiden move, updates partition
    Partition refinePartition();  // creates new partition from partition
    Partition mergeNodesSubset(const Community& subset);  // refine community
    Graph aggregateGraph(const Graph& G, const Partition& P);  // aggregate graph
};