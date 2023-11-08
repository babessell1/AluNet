#include <Rcpp.h>
#include "GraphUtils.h"

// inherit partition from graph object?
class Community {
public:
    std::vector<std::string> nodes; // nodes of a graph
    
    Community(const std::vector<std::string>& nodes); 
    std::vector<int> getIndices();  // returns indices of nodes in community

private:
    size_t number_of_nodes;
    size_t number_of_communities;
};

class Partition {
public:
    Partition(const std::vector<int>& C); // should define more features
    std::vector<int> getIndices();  // returns indices of communities in partition
    inline size_t number_of_nodes() { return input_nodes_vector.size(); } // should revise the name
    Graph get_graph(); // should link it to the graph object

private:
    std::vector<Community> communities;  // communities
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
