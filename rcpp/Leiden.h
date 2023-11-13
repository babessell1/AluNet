#ifndef LEIDEN_H
#define LEIDEN_H

#include "GraphUtils.h"

class Optimizer {
public:
    // properties
    Graph& G;
    Partition& P;
    double gamma;

    // methods
    Optimizer(Graph& G, Partition& P, double gamma);
    void optimize();
    bool moveNodesFast();
    //Partition refinePartition() const;
    //Partition mergeNodesSubset(const Community& subset);
    //Graph aggregateGraph(const Graph& G, const Partition& P);
    double calcQuality(double gamma);
    double deltaQuality(int n_idx, int new_c_idx, double gamma) const;
};

Partition initializePartition();

#endif