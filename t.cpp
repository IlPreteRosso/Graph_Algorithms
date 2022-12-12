#include "./Graph_package_20221211/Graph_Vertex_20221211.hpp"
#include <ctime>


using namespace std;

int main() {
    /*
    int a; char b, c;
    string s;
    stringstream ss;
    map<const int, int> SCCs;

    Graph G;

    auto start = time(NULL);

    G.build_from_adj_list_without_weight("/Users/ilpreterosso/GitHub/VSCode/C++/Week4/data/SCC_test.txt");
    
    for (auto i : G.get_v_keys()) {
        cout << i << endl;
    }
    
    
    
    cout << time(NULL) - start << endl;
    
    SCCs = G.SCC();

    for (auto i : SCCs) {
        cout << i.first << ' ' << i.second << endl;
    }
    */

    Graph G; G.build_from_edge_list_with_weight("./data/FordFulkerson_test.txt");
    //Graph another_G = G.FordFulkerson();
    int MAX_FLOW = G.FordFulkerson(0, 5);
    cout << MAX_FLOW << endl;

    /*
    cout << another_G.get_size() << endl;
    
    for (auto i : another_G.get_v_keys()) {
        cout << i << ' ';
        
        auto temp = *(another_G.get_v(i));
        
        for (auto j : temp.get_head_keys()) {
            cout << ' ' << j;
        }
        cout << endl;
        
    }
    */
    
    

    

    return 0;
}