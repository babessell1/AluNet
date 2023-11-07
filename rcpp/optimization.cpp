#include "optimization.h"

/*****************************
// the partition object here refers to a particular graph are divided into separate communities or groups, and here the objects hold information about the nodes are grouped in each layer of the multiplex network.

// what we need to do: revise this file so that all the corresponding codes are complementary.

// note that this file is not testable, and we shall add more features to make it possible to run.
*****************************/

/*
    We should define the following operations:

    * Basic Graph object (included in GraphUtils.cpp and GraphUtils.h) [done, but need to revise]
    * Partitions, which contains the method that we classify the communities [haven't work] - like it stores the weight, stores the size and stores the relationships
    * Optimization class that is the core of the algorithm: which contains the major optimization methods and leiden algorithm
    
*/

// this file contains the optimization_partition function
// which is the core of the function.

// get graphs of all layers
// initialize the graphs and layers
// declare which nodes are fixed
// do optimization
//  move_node;
//  merge_subset;
//  

/********************
 * We should add a list of is_member_fixed, which indicates whether there are some nodes that are not able to move,
 * this would be useful for judging whether our function should terminate. 
********************/

// [[Rcpp::export]]
double Optimization::optimization_partition(
    Rcpp::vector<Partition*> partitions, 
    Rcpp::vector<double> weights, 
    Rcpp::vector<bool> const& is_memeber_fixed, 
    size_t max_common_size){


    // define the score result of current partition
    double quality_score = 0.0;

    // get the number of multiplex layers
    size_t number_of_classifications = partitions.size();
    // actually when initializing, we should have $n$ partitions (itself partition)
    if(number_of_classifications == 0){
        throw Exception("No partitions provided.");
    }

    // get the resulting graph of the layers
    Rcpp::vector<Graph*> graphs(number_of_classifications);
    for(size_t i = 0; i < number_of_classifications; i++){
        graphs[i] = partitions[i]->get_graph(); // return the graphs of the specific partition i;
    }
    // now we should return the singleton graph
    // since all of the subgraph have the same number of nodes, 
    // we only take the fist one

    // when initializing, we have number_of_nodes == 1.
    size_t number_of_nodes = graphs[0]->number_of_nodes(); 
    // check and if not, we should throw an error
    
    for(size_t i = 0; i < number_of_classifications; i++){
        if(graphs[i]->number_of_nodes() != number_of_nodes)
            throw Exception("Number of nodes are not equal for all graphs.");
    }

    // we should fix some membership for fixed nodes
    Rcpp::vector<size_t> fixed_nodes;
    Rcpp::vector<size_t> fixed_membership(number_of_nodes);

    // is_member_fixed is a vector that contains whether j is fixed
    for(size_t j = 0; j < number_of_nodes; j++){
        if(is_memeber_fixed[j]){ // if fixed
            fixed_nodes.push_back(j);
            fixed_membership[j] = partitions[0]->membership(j);
        }
    }

    // initialize the vector of collapsed graphs for all layers
    // collapsing graphs is a part of a layer-by-layer optimization process
    // where the partition is optimized at each level.
    Rcpp::vector<Graph* > collapsed_graphs(number_of_classifications);
    Rcpp::vector<Partitions*> collapsed_partitions(number_of_classifications);

    for(size_t layer = 0; layer < number_of_classifications; layer++){
        collapsed_graphs[layer] = graphs[layer];
        collapsed_partitions[layer] = partitions[layer];
    }

    // Declare which nodes in the collapsed graph are fixed, which to start is
    // simply equal to is_membership_fixed
    Rcpp::vector<bool> is_collpased_membership_fixed(is_memeber_fixed);

    // start a vector indicating the aggregating node (initilizing with n)
    // and it maps from the original node to the corresponding aggregated nodes.
    Rcpp::vector<size_t> aggregate_node_per_individual_node = range(number_of_nodes);

    // now we define whether it is possible to further refine
    bool aggregate_further = true; // give a flag

    // as long as there remains improvement iterate
    double improvement = 0.0;

    // this reflects the aggregate node;
    while(aggregate_further)
    {
        // move nodes
        if()
        improve += move_nodes(); // this should be further defined
        else if()
        // merge nodes
        improve += merge_nodes(); // this should be further defined
        // make sure imrpovement 

        // cluster_partition should be specified in the partition cpp file instead of
        // optimization file since this involves the partition method.
        for(size_t i = 0; i < number_of_classifications; i++){
            if(collapsed_partitions[i] != partitions[i]){
                if(this->refine_parition){ // a flag indicating whether to revise the partition
                    partitions[i]->cluster_partition(collapsed_partitions[i], aggregate_node_per_individual_node); //
                }else{
                    partitions[i]->cluster_partition(collapsed_partitions[i]);
                }
            }
        }

        // collapse graph
        // if refine the graph, we separate communities in slightly more
        // fine-grained parts for which we collapse the graph.
        Rcpp::vector<Paritions* > sub_collapsed_partitions(number_of_classifications);

        Rcpp::vector<Graph*> new_collapsed_graph(number_of_classifications);
        Rcpp::vector<Partitions* > new_collapsed_partitions(number_of_classifications);
    
        if(this->refine_partition){
            // first create a new partition, which should be a sub partition
            // of the collapsed partition

            for(size_t i = 0; i < number_of_classifications; i++){
                // create a new partition
                new_collapsed_partitions[i] = collapsed_partitions[i]->create_partition(collapsed_graphs[i]);
            }


        // Then move around nodes but restrict movement to within original communities.
        if(){
            move_nodes(); // reconstruct the move_nodes with 
        }
        else if(){
            merge_nodes_within_community();
        }
        
        // determine new aggregate node per individual nodes

        for(int j = 0; j < number_of_nodes; j++){
            size_t aggregate_node = aggregate_node_per_individual_node[j];
            aggregate_node_per_individual_node[j] = sub_collapsed_partitions[0]->membership(aggregate_node);
        }

        // collapse graph basd on sub collapsed partition
        for(int i = 0; i < number_of_classifications; i++){
            new_collapsed_graph[i] = collapsed_graphs[i]->collapse_graph(sub_collapsed_partitions[i]);
        }

        // determine the membership for the collapsed graph
        Rcpp::vector<size_t> new_collapsed_membership(new_collapsed_graph[0]->number_of_nodes());

        for(size_t i = 0; i < collapsed_graphs[0]->number_of_nodes(); i++){
            size_t new_aggregate_node = sub_collapsed_partitions[0]->membership(i);
            new_collapsed_membership[new_aggregate_node] = collapsed_partitions[0]->membership(i);
        }

        // determine which collapsed nodes are fixed
        is_collapsed_membership_fixed.clear();
        is_collapsed_membership_fixed.resize(new_collapsed_graphs[0]->number_of_nodes(), false);
        for(size_t i = 0; i < number_of_nodes; i++){
            if(is_memeber_fixed[i]){
                is_collapsed_membership_fixed[aggregate_node_per_individual_node[i]] = true;
            }
        }

        // determine whether to aggregate further
        for(size_t j = 0; j < number_of_classifications; j++){
            delete sub_collapsed_partitions[j];
            new_collapsed_partitions[j] = collapsed_partitions[j]->create(new_collapsed_graphs[j], new_collapsed_membership);
        }
    }
    else{
        for(size_t j = 0; j < number_of_classifications; j++){
            new_collapsed_graphs[j] = collapsed_graphs[j]->collapse_graph(collapsed_partitions[j]);
        // Create collapsed partition (i.e. default partition of each node in its own community).
        new_collapsed_partitions[j] = collapsed_partitions[j]->create(new_collapsed_graphs[j]);
        }
    }

    // Determine whether to aggregate further
    // If all memberships are fixed, there is no need to aggregate more.
    aggregate_further = false;
    for(const bool& fixed: is_collapsed_membership_fixed){
        if(!fixed){
            aggregate_further = true;
            break;
        }
    }

    // If not all memberships are fixed, then check if there has been any changes since the last time
    aggregate_further &= (new_collapsed_graphs[0]->number_of_nodes()) && (collapsed_graphs[0]->number_of_nodes() > collapsed_partitions[0]->number_of_communities());
    
    // delete the collapsed partitions and graphs
    for(size_t j = 0; j < number_of_classifications; j++){
        if(collapsed_partitions[j] != partitions[j])
            delete collapsed_partitions[j];
        if(collapsed_graphs[j] != graphs[j]){
            delete collapsed_graphs[j];
        }
    }

     // and set them to the new one.
     // After deletion, set the current collapsed partitions and graphs to the new ones.
    collapsed_partitions = new_collapsed_partitions;
    collapsed_graphs = new_collapsed_graphs;

    }

    // Clean up memory after use.
  for (size_t layer = 0; layer < nb_layers; layer++)
  {
    if (collapsed_partitions[layer] != partitions[layer])
      delete collapsed_partitions[layer];

    if (collapsed_graphs[layer] != graphs[layer])
      delete collapsed_graphs[layer];
  }

  quality_score = 0.0;

  // Renumber communities in the partitions
  partitions[0]->renumber_communities(); // should be included in the partitions methods
  partitions[0]->renumber_communities(fixed_nodes, fixed_membership);

  // Get the membership from partition[0]
  vector<size_t> const& membership = partitions[0]->membership();
 
    // Loop over all other partitions, set their membership to the one from partition[0] and add their quality to the total quality
  for (size_t j = 1; j < number_of_classifications; j++)
  {
    partitions[j]->set_membership(membership);
    quality_score += partitions[j]->quality()*layer_weights[j];
  }
  return improv;
}

double Optimization::move_nodes(Rcpp::vector<Partition*> partitions, Rcpp::vector<double> weights, Rcpp::vector<bool> const& is_member_fixed, size_t max_community_size)
{

    // the first criterion is to find all the unchangeable items
    // i.e., fixed nodes are stable ones

    // skip these nodes when using gready algorithm to find improvements

    // delete empty communities so that we won't have uncessary commons

    // mark all changable neighbors as unstable and rerun the function
    // to find all chengable.

    return 
}

