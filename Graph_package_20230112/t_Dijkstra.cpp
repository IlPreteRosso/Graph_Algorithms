#include "Graph_Vertex_20230112.hpp"

int main() {
    Graph G; G.build_from_edge_list_with_weight("archive/data/Dijkstra.txt");
    pair<map<const int, int>, map<const int, vector<int>*>> Pair;
    Pair = G.Dijkstra_with_path(0);

    for (auto i : Pair.second) {
        cout << "0 -> " << i.first << ": ";
        for (auto j : *i.second) {
            cout << j << ' ';
        }
        cout << endl;
    }

    return 0;
}