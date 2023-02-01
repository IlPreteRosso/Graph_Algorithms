#include "Graph_Vertex_20230112.hpp"
#include <time.h>


int main() {
    Graph G; G.build_from_edge_list_with_weight("archive/data/EdmondsKarp.txt");
    int t;

    t = G.EdmondsKarp(0, 5);
    cout << t << endl;

    return 0;
}