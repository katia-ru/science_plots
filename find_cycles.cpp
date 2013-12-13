// Copyright 2013 Katia Rubtcova
#include <iostream>
#include <queue>
#include <set>
#include <vector>

using std::cin;
using std::cout;
using std::queue;
using std::set;
using std::vector;

const int MAX_CYCLE_LENGTH = 7;

class CyclesFinder {
 public:
  set<set<int> > FindCycles(const vector<vector<int> >& graph) const {
    set<set<int> > cycles;
    for (size_t vertex = 0; vertex < graph.size(); ++vertex) {
      FindCyclesFromVertex(graph, vertex, &cycles);
    }
    return cycles;
  }

 private:
  void FindCyclesFromVertex(const vector<vector<int> >& graph,
                            int vertex,
                            set<set<int> >* cycles) const {
    vector<int> depthes(graph.size(), 0);
    vector<int> parents(graph.size(), -1);
    queue<int> vertices_queue;
    depthes[vertex] = 1;
    vertices_queue.push(vertex);
    while (!vertices_queue.empty()) {
      int current_vertex = vertices_queue.front();
      vertices_queue.pop();
      for (size_t index = 0; index < graph[current_vertex].size(); ++index) {
        int adjacent_vertex = graph[current_vertex][index];
        if (depthes[current_vertex] > 1 && depthes[adjacent_vertex] > 0 &&
            depthes[adjacent_vertex] + 1 != depthes[current_vertex]) {
          set<int> new_cycle;
          AddVerticesToCycle(parents, adjacent_vertex, &new_cycle);
          AddVerticesToCycle(parents, current_vertex, &new_cycle);
          cycles->insert(new_cycle);
        } else if (depthes[current_vertex] < MAX_CYCLE_LENGTH &&
                   depthes[adjacent_vertex] == 0) {
          depthes[adjacent_vertex] = depthes[current_vertex] + 1;
          parents[adjacent_vertex] = current_vertex;
          vertices_queue.push(adjacent_vertex);
        }
      }
    }
  }

  void AddVerticesToCycle(const vector<int>& parents,
                          int vertex,
                          set<int>* cycle) const {
    while (vertex != -1) {
      cycle->insert(vertex);
      vertex = parents[vertex];
    }
  }
};

vector<vector<int> > Input() {
  vector<vector<int> > graph;
  int vertices_number;
  int edges_number;
  cin >> vertices_number >> edges_number;
  graph.resize(vertices_number, vector<int>());
  for (int number = 0; number < edges_number; ++number) {
    int source_vertex;
    int target_vertex;
    cin >> source_vertex >> target_vertex;
    graph[source_vertex].push_back(target_vertex);
    graph[target_vertex].push_back(source_vertex);
  }
  return graph;
}

void Output(const set<set<int> >& cycles) {
  if (cycles.size() > 0) {
    for (set<set<int> >::const_iterator cycle_iterator = cycles.begin();
         cycle_iterator != cycles.end(); ++cycle_iterator) {
      for (set<int>::const_iterator vertex_iterator =
           cycle_iterator->begin(); vertex_iterator != cycle_iterator->end();
           ++vertex_iterator) {
        cout << *vertex_iterator << " ";
      }
      cout << "\n";
    }
  } else {
    cout << "No cycles found\n";
  }
}

int main() {
  vector<vector<int> > graph = Input();
  CyclesFinder cycles_finder;
  Output(cycles_finder.FindCycles(graph));
  return 0;
}
