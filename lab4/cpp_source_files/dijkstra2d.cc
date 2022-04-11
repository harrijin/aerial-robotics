#include <queue>

#include "dijkstra2d.h"

namespace game_engine {
  // Anonymous namespace. Put any file-local functions or variables in here
  namespace {
    // Helper struct that functions as a linked list with data. The linked
    // list represents a path. Data members are a node and a cost to reach
    // that node.
    struct NodeWrapper {
      std::shared_ptr<struct NodeWrapper> parent;
      std::shared_ptr<Node2D> node_ptr;
      double cost;

      // Equality operator
      bool operator==(const NodeWrapper& other) const {
        return *(this->node_ptr) == *(other.node_ptr);
      }
    };

    // Helper function. Compares the values of two NodeWrapper pointers.
    // Necessary for the priority queue.
    bool NodeWrapperPtrCompare(
        const std::shared_ptr<NodeWrapper>& lhs, 
        const std::shared_ptr<NodeWrapper>& rhs) {
      return lhs->cost > rhs->cost;
    }

    bool is_explored(std::vector<std::shared_ptr<NodeWrapper>>& explored, std::shared_ptr<Node2D>& node){
      for(auto& exploredNodePtr: explored){
        if(*(exploredNodePtr->node_ptr) == *node){
          return true;
        }
        // if(node->Data()[0] == exploredNodePtr->node_ptr->Data()[0] && node->Data()[1] == exploredNodePtr->node_ptr->Data()[1]){
        //   std::cout << "False negative" << std::endl;
        // }
      }
      return false;
    }
  }

  PathInfo Dijkstra2D::Run(
      const Graph2D& graph, 
      const std::shared_ptr<Node2D> start_ptr, 
      const std::shared_ptr<Node2D> end_ptr) {
    using NodeWrapperPtr = std::shared_ptr<NodeWrapper>;

    ///////////////////////////////////////////////////////////////////
    // SETUP
    // DO NOT MODIFY THIS
    ///////////////////////////////////////////////////////////////////
    Timer timer;
    timer.Start();

    // Use these data structures
    std::priority_queue<
      NodeWrapperPtr,
      std::vector<NodeWrapperPtr>,
      std::function<bool(
          const NodeWrapperPtr&, 
          const NodeWrapperPtr& )>> 
        to_explore(NodeWrapperPtrCompare);

    std::vector<NodeWrapperPtr> explored;

    ///////////////////////////////////////////////////////////////////
    // YOUR WORK GOES HERE
    // SOME EXAMPLE CODE INCLUDED BELOW
    ///////////////////////////////////////////////////////////////////

    // Create a NodeWrapperPtr
    NodeWrapperPtr nw_ptr = std::make_shared<NodeWrapper>();
    nw_ptr->parent = nullptr;
    nw_ptr->node_ptr = start_ptr;
    nw_ptr->cost = 0;
    to_explore.push(nw_ptr);
    NodeWrapperPtr final_ptr;
    while(!to_explore.empty()){
      NodeWrapperPtr cur_node_ptr = to_explore.top();
      to_explore.pop();
      if(is_explored(explored, cur_node_ptr->node_ptr)){
        continue;
      }
      if(*(cur_node_ptr->node_ptr) == *(end_ptr)){
        final_ptr = cur_node_ptr;
        break;
      }
      std::vector<std::shared_ptr<Node2D>> neighbors = graph.Neighbors(cur_node_ptr->node_ptr);
      std::vector<DirectedEdge2D> edges = graph.Edges(cur_node_ptr->node_ptr);
      for(auto &neighbor: neighbors){
        // node has not been explored yet
        NodeWrapperPtr neighbor_ptr = std::make_shared<NodeWrapper>();
        neighbor_ptr->node_ptr = neighbor;
        neighbor_ptr->parent = cur_node_ptr;
        for(auto &edge: edges){
          if(edge.Sink() == neighbor){
            neighbor_ptr->cost = edge.Cost()+cur_node_ptr->cost;
            to_explore.push(neighbor_ptr);
            break;
          }
        }
      }
      explored.push_back(cur_node_ptr);
    }
    // Create a PathInfo
    PathInfo path_info;
    path_info.details.num_nodes_explored = explored.size()+1; // +1 to account for final node
    path_info.details.path_length = 0;
    path_info.details.path_cost = final_ptr->cost;
    path_info.details.run_time = timer.Stop();
    path_info.path = {};

    // Traverse linked list path
    NodeWrapperPtr cur_ptr = final_ptr;
    while(cur_ptr){
      path_info.details.path_length++;
      path_info.path.push_back(cur_ptr->node_ptr);
      cur_ptr = cur_ptr->parent;
    }
    std::reverse(path_info.path.begin(), path_info.path.end());

    // You must return a PathInfo
    return path_info;
  }
  
}
