#ifndef LEIDEN_H
#define LEIDEN_H

#include "GraphUtils.h"

class Optimizer {
public:
    // properties
    Graph& G;
    Partition& P;
    double gamma;
    double temperature;

    // methods
    Optimizer(Graph& G, Partition& P, double gamma, double temperature);
    void optimize();
    bool moveNodesFast();
    //Partition refinePartition() const;
    //Partition mergeNodesSubset(const Community& subset);
    //Graph aggregateGraph(const Graph& G, const Partition& P);
    double calcQuality(double gamma);
    double deltaQuality(int n_idx, int new_c_idx, double gamma, bool recalculate) const;
    void refinePartition(const Partition& P_original);
    std::vector<int> getWellConnectedNodes(const Community& B) const;
    std::vector<Community> getWellConnectedCommunities(const Community& B) const;
    void mergeNodesSubset(Community& S);
    Graph aggregateGraph();
    
};

Partition initializePartition(Graph& G);

#endif