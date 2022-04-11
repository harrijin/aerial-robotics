#include <queue>
#include <cmath>
#include "a_star2d.h"

namespace game_engine {
  // Anonymous namespace. Put any file-local functions or variables in here
  namespace {
    // Helper struct that functions as a linked list with data. The linked
    // list represents a path. Data members are a node, a cost to reach that
    // node, and a heuristic cost from the current node to the destination.
    struct NodeWrapper {
      std::shared_ptr<struct NodeWrapper> parent;
      std::shared_ptr<Node2D> node_ptr;

      // True cost to this node
      double cost;

      // Heuristic to end node
      double heuristic;

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
      return lhs->cost + lhs->heuristic > rhs->cost + rhs->heuristic;
    }

    bool is_explored(std::vector<std::shared_ptr<NodeWrapper>>& explored, std::shared_ptr<Node2D>& node){
      for(auto& exploredNodePtr: explored){
        if(*(exploredNodePtr->node_ptr) == *node){
          return true;
        }
      }
      return false;
    }

    ///////////////////////////////////////////////////////////////////
    // EXAMPLE HEURISTIC FUNCTION
    // YOU WILL NEED TO MODIFY THIS OR WRITE YOUR OWN FUNCTION
    ///////////////////////////////////////////////////////////////////
    double ZeroHeuristic(
        const std::shared_ptr<Node2D>& current_ptr,
        const std::shared_ptr<Node2D>& end_ptr) {
      return 0;
    }
    double RowHeuristic(
        const std::shared_ptr<Node2D>& current_ptr,
        const std::shared_ptr<Node2D>& end_ptr) {
      return abs(current_ptr->Data()[0]-end_ptr->Data()[0]);
    }
    double EuclideanHeuristic(
        const std::shared_ptr<Node2D>& current_ptr,
        const std::shared_ptr<Node2D>& end_ptr) {
      return sqrt(pow(current_ptr->Data()[0]-end_ptr->Data()[0],2)+
                  pow(current_ptr->Data()[1]-end_ptr->Data()[1],2));
    }
    double OverestimatedHeuristic(
        const std::shared_ptr<Node2D>& current_ptr,
        const std::shared_ptr<Node2D>& end_ptr) {
      return abs(current_ptr->Data()[0]-end_ptr->Data()[0])+
             abs(current_ptr->Data()[1]-end_ptr->Data()[1])+1;
    }

  }

  PathInfo AStar2D::Run(
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
    // Choose Heuristic function
    double (*Heuristic) (const std::shared_ptr<Node2D>& current_ptr, const std::shared_ptr<Node2D>& end_ptr);
    Heuristic = ZeroHeuristic;
    // Create a NodeWrapperPtr
    NodeWrapperPtr nw_ptr = std::make_shared<NodeWrapper>();
    nw_ptr->parent = nullptr;
    nw_ptr->node_ptr = start_ptr;
    nw_ptr->cost = 0;
    nw_ptr->heuristic = Heuristic(start_ptr, end_ptr);
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
      explored.push_back(cur_node_ptr);
      std::vector<std::shared_ptr<Node2D>> neighbors = graph.Neighbors(cur_node_ptr->node_ptr);
      std::vector<DirectedEdge2D> edges = graph.Edges(cur_node_ptr->node_ptr);
      for(auto &neighbor: neighbors){
        NodeWrapperPtr neighbor_ptr = std::make_shared<NodeWrapper>();
        neighbor_ptr->node_ptr = neighbor;
        neighbor_ptr->parent = cur_node_ptr;
        neighbor_ptr->heuristic = Heuristic(neighbor, end_ptr);
        for(auto &edge: edges){
          if(edge.Sink() == neighbor){
            neighbor_ptr->cost = edge.Cost()+cur_node_ptr->cost;
            break;
          }
        }
        to_explore.push(neighbor_ptr);
      }
    }
    // Create a PathInfo
    PathInfo path_info;
    path_info.details.num_nodes_explored = explored.size(); // +1 to account for final node
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
