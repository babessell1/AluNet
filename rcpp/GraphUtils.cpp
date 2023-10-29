#include <Rcpp.h>
using namespace Rcpp;

// Custom function to create a graph structure
List createGraph(List nodes, List edges, NumericVector weights) {
  int numNodes = nodes.size();
  int numEdges = edges.size();

  List graph;
  graph["nodes"] = nodes;
  graph["edges"] = edges;
  graph["weights"] = weights;

  return graph;
}

// Export the function to R using Rcpp attributes
// [[Rcpp::export]]
List createGraphR(List nodes, List edges, NumericVector weights) {
  return createGraph(nodes, edges, weights);
}
