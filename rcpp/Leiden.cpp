//#include <chrono>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>

//#include <RcppEigen.h>
#include <Rcpp.h>
//#include <progress.hpp>

// [[Rcpp::plugins(cpp11)]]

#include "Leiden.h"

// using namespace Leiden;
//using namespace std::chrono;
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
