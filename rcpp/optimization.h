#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>

#include <vector>
#include <set>
#include <map>
#include <cfloat>
/*
    These packages should be included in our code,
    but we should modify it if Rcpp not supports these.
*/

// we define a new class that utilizing 

class Optimization{
public:
    double optimization_partition(Rcpp::vector<Partition*> partitions, Rcpp::vector<double> weights, Rcpp::vector<bool> const& is_memeber_fixed, size_t max_common_size);
    double move_nodes(Rcpp::vector<Partition*> partitions, Rcpp::vector<double> weights, Rcpp::vector<bool> const& is_member_fixed, size_t max_community_size); // add more features if possible and necessary

    double move_nodes(Rcpp::vector<Partition*> partitions, Rcpp::vector<double> weights, Rcpp::vector<bool> const& is_member_fixed, size_t max_community_size, ); // use reconstruct method to add a reconstruction of the move_node method
    // to incorporate the move_node method within the community.

    double merge_nodes(Rcpp::vector<Partition*> partitions, Rcpp::vector<double> weights, Rcpp::vector<bool> const& is_member_fixed, size_t max_community_size); // add more features if possible and necessary

    
private:
    // thinking what are needed....
    bool refine_parition = true; // flag that indicates whether we should proceed optimizing the graph;
};

#endif // OPTIMIZATOM_H