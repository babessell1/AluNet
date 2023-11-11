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
    std::vector<Partition*> partitions, 
    std::vector<double> weights, 
    std::vector<bool> const& is_memeber_fixed, 
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
    // initialization of the graph structure
    std::vector<Graph*> graphs(number_of_classifications);
    for(size_t i = 0; i < number_of_classifications; i++){
        graphs[i] = partitions[i]->get_graph(); // return the graphs of the specific partition i;
    }
    // now we should return the singleton graph
    // since all of the subgraph have the same number of nodes, 
    // we only take the fist one

    // when initializing, we have number_of_nodes == 1.
    std::vector<int> get_nodes = graphs[0]->getNodes();
    
    size_t number_of_nodes = get_nodes.size(); 
    // check and if not, we should throw an error
    
    for(size_t i = 0; i < number_of_classifications; i++){
        if(graphs[i]->getNodes().size() != number_of_nodes)
            throw Exception("Number of nodes are not equal for all graphs.");
    }

    // we should fix some membership for fixed nodes
    std::vector<size_t> fixed_nodes;
    std::vector<size_t> fixed_membership(number_of_nodes);

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
    std::vector<Graph* > collapsed_graphs(number_of_classifications);
    std::vector<Partitions*> collapsed_partitions(number_of_classifications);

    for(size_t j = 0; j < number_of_classifications; j++){
        collapsed_graphs[j] = graphs[j];
        collapsed_partitions[j] = partitions[j];
    }

    // Declare which nodes in the collapsed graph are fixed, which to start is
    // simply equal to is_membership_fixed
    std::vector<bool> is_collpased_membership_fixed(is_memeber_fixed);

    // start a vector indicating the aggregating node (initilizing with n)
    // and it maps from the original node to the corresponding aggregated nodes.
    std::vector<size_t> aggregate_node_per_individual_node = range(number_of_nodes);

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
        std::vector<Paritions* > sub_collapsed_partitions(number_of_classifications);

        std::vector<Graph*> new_collapsed_graph(number_of_classifications);
        std::vector<Partitions* > new_collapsed_partitions(number_of_classifications);
    
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
        std::vector<size_t> new_collapsed_membership(new_collapsed_graph[0]->getNodes().size());

        for(size_t i = 0; i < collapsed_graphs[0]->getNodes().size(); i++){
            size_t new_aggregate_node = sub_collapsed_partitions[0]->membership(i);
            new_collapsed_membership[new_aggregate_node] = collapsed_partitions[0]->membership(i);
        }

        // determine which collapsed nodes are fixed
        is_collapsed_membership_fixed.clear();
        is_collapsed_membership_fixed.resize(new_collapsed_graphs[0]->getNodes().size(), false);
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
    aggregate_further &= (new_collapsed_graphs[0]->getNodes().size()) && (collapsed_graphs[0]->getNodes().size() > collapsed_partitions[0]->number_of_communities());
    
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
    partitions[j]->updateCommunityMembership(membership);
    quality_score += partitions[j]->quality()*layer_weights[j];
  }
  return improv;
}

double Optimization::move_nodes(std::vector<Partition*> partitions, std::vector<double> weights, std::vector<bool> const& is_member_fixed, size_t max_community_size, bool renumber_fixed_nodes)
{
    // check the number of partitions
    // Define number_of_classifications as the size of the partitions vector
    size_t number_of_classifications = partitions.size();
    if (number_of_classifications == 0) {
        return -1.0;
    }
    // Create a vector of Graph pointers with the size of the number of partitions
    std::vector<Graph*> graphs(number_of_classifications);
    for (size_t i = 0; i < number_of_classifications; i++) {
        graphs[i] = partitions[i]->get_graph();
    }

    // number of nodes in the graph
    // just like the explannation in the optimization function
    // the number of nodes should be equal
    // Define the number_of_nodes as the size of the node set in the first graph
    size_t number_of_nodes = graphs[0]->getNodes().size();

    // get the fixed membership from the iterations
    // we should fix some membership for fixed nodes
    // Initialize a vector for fixed_nodes and fixed_membership using number_of_nodes
    std::vector<size_t> fixed_nodes;
    std::vector<size_t> fixed_membership(number_of_nodes);

    // is_member_fixed is a vector that contains whether j is fixed
    // Populate fixed_nodes and fixed_membership for each node that is fixed
    if(renumber_fixed_nodes){
        for (size_t j = 0; j < number_of_nodes; j++) {
            if (is_memeber_fixed[j]) { // if fixed
                fixed_nodes.push_back(j);
                fixed_membership[j] = partitions[0]->membership(j);
            }
        }
    }

    double total_improvement = 0.0;
    // Assert the node size of each graph is equal to number_of_nodes
    for(Graph* graph: graphs){
        if(graph->getNodes().size() != n){
            throw Exception("Number of nodes are not equal for all graphs.");
        }
    }

    int number_of_moves = 0; // initialize number_of_moves to 0
    // Copy is_member_fixed into is_node_stable
    std::vector<bool> is_node_stable(is_member_fixed);

    // now we first establish the number of the 
    // we should shuffle the order which is greedy
    // and we should use queue to operate the result
    // Create a vector to hold nodes that are not fixed
    std::vector<size_t> nodes;
    for(size_t i = 0; i != is_member_fixed.size(); i++){
        if(!is_member_fixed[i]){
            nodes.push_back(i);
        }
    }

    // then we should resample the nodes
    resample()
    // Initialize a queue using the nodes in the current order
    std::deque<size_t> vertex_order(nodes.begin(), nodes.end());

    // itialize a vector for storing the combined community and initialize a vector for communities
    std::vector<bool> combined_community(partitions[0]->number_of_communities(), false);
    std::vector<size_t> communities;

    // repeat the process until the queue is empty
    // Process nodes in the order from the queue
    while(!vertex_order.empty()){
        size_t first_vertex = vertex_order.front(); vertex_order.pop_front();
        
        // Clear comms
        // Clear previous community assignments
        for (size_t comm : communities)
            combined_community[comm] = false;
        communities.clear();

        // current community
        // Assign current node's community
        size_t first_vertex_community = partitions[0]->membership(first_vertex);
        // Initialize current node's community
        for(size_t community = 0; community < partitions[0]->number_of_communities(); community++){
            for(size_t i = 0; i < number_of_classifications; i++){
                if(partitions[i]->number_of_nodes_in_community(community) > 0 && !combined_community[community]){
                    communities.push_back(community);
                    combined_community[community] = true;
                    break; // Break from for loop in layer
                }
            }
        }
        // Initialize max_community and max_improvement
        size_t max_community = first_vertex_community;
        double max_improvement = (0 < max_community_size && max_community_size < partitions[0]->community_size(first_vertex_community) ? -INFINITY : 10 * DBL_EPSILON);
        double v_size = graphs[0]->size_of_nodes(first_vertex);
        
        // this is the key part
        // Evaluate best community to move current node
        for (size_t comm : communities){
            // If the community is too large, continue
            if(0 < max_community_size && max_community_size < partitions[0]->community_size(comm) + v_size){
                continue;
            }

            double possible_improv = 0.0;

            // Evaluate improvement of moving for each partition
            for (size_t i = 0; i < number_of_classifications; i++){
                // Make sure to multiply it by the weight per layer
                // use the function that returns the difference of moving this node to another community.
                possible_improv += weights[i] * partitions[i]->diff_move(first_vertex, comm);

            }
            // Update maximum improvement and community
            if (possible_improv > max_improv)
            {
                max_comm = comm;
                max_improv = possible_improv;
            }
        }
        // Clear communities for next node
        // Clear comms
        communities.clear();
        is_node_stable[first_vertex] = true;

        // If we actually plan to move the nove
        if (max_comm != first_vertex_community)
        {
        // Keep track of improvement
            total_improvement += max_improv;

            for (size_t i = 0; i < number_of_classifications; i++){
                Partitions* partition = partitions[i];
                // actually move the nodes
                partition->move_node(first_vertex, max_comm);
            }

            // Mark neighbours as unstable (if not in new community and not fixed)
            // Update stability of neighboring nodes
             for (Graph* graph : graphs)
            {
                for (size_t u : graph->get_neighbours(first_vertex)){
                    
                    if (is_node_stable[u] && partitions[0]->membership(u) != max_comm && !is_membership_fixed[u])
                    {
                        vertex_order.push_back(u);
                        is_node_stable[u] = false;
                    }
                }
            }
            // Keep track of number of moves
            number_of_moves += 1;
        }
    }

//Renumber the communities so that they are numbered 0,...,q-1 where q is
// the number of communities. This also removes any empty communities, as they
// will not be given a new number.
    // Renumber communities to remove any gaps
    partitions[0]->renumber_communities();
    if(renumber_fixed_nodes)
        partitions[0]->renumber_communities(fixed_nodes, fixed_membership);
    vector<size_t> const& membership = partitions[0]->membership();
    // Replicate membership across partitions
    for (size_t i = 1; i < number_of_classifications; i++)
    {
        partitions[i]->updateCommunityMembership(membership);
    }
    return total_improvement;
}


double Optimiser::merge_nodes(std::vector<Partition*> partitions, 
    std::vector<double> weights, 
    std::vector<bool> const& is_memeber_fixed, 
    size_t max_common_size,
    int consider_commons,
    bool renumber_fixed_nodes)
{
    // number of classifications
    size_t number_of_classifications = partitions.size();
    if(number_of_classifications == 0){
        return -1.0;
    }

    // get the graph and initialization of the graph using 0
    std::vector<Graph*> graphs(number_of_classifications);
    for(size_t i = 0; i < number_of_classifications; i++){
        graphs[i] = partitions[i]->get_graph();
    }

    // number of nodes in the graph
    size_t number_of_nodes = graphs[0]->getNodes().size(); 
    // check and if not, we should throw an error

    // we should fix some membership for fixed nodes
    std::vector<size_t> fixed_nodes;
    std::vector<size_t> fixed_membership(number_of_nodes);

    // is_member_fixed is a vector that contains whether j is fixed
    if(renumber_fixed_nodes){
        for (size_t j = 0; j < number_of_nodes; j++) {
            if (is_memeber_fixed[j]) { // if fixed
                fixed_nodes.push_back(j);
                fixed_membership[j] = partitions[0]->membership(j);
            }
        }
    }

    double total_improvement = 0.0;
    // Assert the node size of each graph is equal to number_of_nodes
    for(Graph* graph: graphs){
        if(graph->getNodes().size() != n){
            throw Exception("Number of nodes are not equal for all graphs.");
        }
    }

    // Establish vertex order, skipping fixed nodes
    std::vector<size_t> nodes;
    for(size_t j = 0; j != number_of_nodes; j++){
        if(!is_membership_fixed[j])
            nodes.push_back(j);
    }

    // make it random just like the algorithm descriptions
    shuffle();

    // itialize a vector for storing the combined community and initialize a vector for communities
    std::vector<bool> combined_community(partitions[0]->number_of_communities(), false);
    std::vector<size_t> communities;

    // Iterate over all nodes
    for( size_t i : nodes){
        size_t i_within_community = partition[0]->membership(i);
        // Clear comms
        for (size_t comm : communities){
            combined_community[comm] = false;
        }
        communities.clear();

        if(partitions[0]->number_of_nodes_in_community(comm) == 1){
            for(size_t comm = 0; comm < partitions[0]->number_of_communities(); comm++){
                for(size_t i = 0; i < number_of_classifications; i++){
                    if(partitions[i]->number_of_nodes_in_community(comm) > 0 && !combined_community[comm]){
                        communities.push_back(comm);
                        combined_community[comm] = true;
                        break; // Break from for loop in layer
                    }
                }
            }

            size_t max_community = i_within_community;
            double max_improvement = (0 < max_community_size && max_community_size < partitions[0]->community_size(i_within_community) ? -INFINITY : 0);
            double v_size = graphs[0]->size_of_nodes(i);
        
            for(size_t comm : communities){
                // Do not create too-large communities.
                if(0 < max_community_size && max_community_size < partititions[0]->community_size(comm) + v_size){
                    continue;
                }
                double possible_improv = 0;

                for(size_t i = 0; i < number_of_classifications; i++){
                    possible_improv += weights[i] * partitions[i]->diff_move(i_within_community, comm);
                }

                if(possible_improv >= max_improvement){
                    max_community = comm;
                    max_improvement = possible_improv;
                }
            }

            if(max_community != i_within_community){
                // Keep track of improvement
                total_improv += max_improvement;

                for(size_t i = 0; i < number_of_classifications; i++){
                    Partition* partition = partitions[i];

                    // actually move the node
                    partition->move_node(v, max_comm);
                }
            }
        }
    }

    partitions[0]->renumber_communities();
    if (renumber_fixed_nodes)
        partitions[0]->renumber_communities(fixed_nodes, fixed_membership);
    std::vector<size_t> const& membership = partitions[0]->membership();
    for(size_t i = 0; i < number_of_classifications; i++){
        partitions[i]->updateCommunityMembership(membership);
    }

    return total_improvement;

}
