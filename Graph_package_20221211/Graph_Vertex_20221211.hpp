/*
    Created: 20221211
    Definition of advenced Graph algorithms
*/

#include "Graph_Vertex_basics_20221211.hpp"

inline void Graph::DFS_TS(int s, map<const int, bool> &E, deque<int>&TO, bool reversed=false) {
    E[s] = true;

    if (reversed) {
        for (auto v : this->get_v(s)->get_tail_keys()) {
            if (!E[v]) {
                DFS_TS(v, E, TO, reversed);
            }
        }
    }
    else {
        for (auto v : this->get_v(s)->get_head_keys()) {
            if (!E[v]) {
                DFS_TS(v, E, TO, reversed);
            }
        }
    }
    TO.push_front(s);
}

inline void Graph::DFS_SCC(int s, map<const int, bool> &E, deque<int>&TO, map<const int, int> &SCCs, int &numSCC) {
    E[s] = 1;
    SCCs[s] = numSCC;

    for (auto v : this->get_v(s)->get_head_keys()) {
        if (!E[v]) {
            DFS_SCC(v, E, TO, SCCs, numSCC);
        }
    }
}

inline void Graph::update_Res_Aug_Graph(Graph &Res_Aug_Graph, int &min_cut, vector<int> &Path_clean) {
    int path_tail, path_head;
    while(true) {
        path_tail = Path_clean.back(); Path_clean.pop_back();
        // if statement has the same functionality as above while loop
        if (Path_clean.empty()) {
            break;
        }
        path_head = Path_clean.back(); //Path_clean.pop_back(); Commented for the same reason as before

        //========= path_tail -> path_head
        // Remove saturated residual
        if (Res_Aug_Graph.get_e(path_tail, path_head) == min_cut) {
            Res_Aug_Graph.remove_e(path_tail, path_head);
        }
        // decrese unsaturated residual
        else {  
            Res_Aug_Graph.update_e(path_tail, path_head, Res_Aug_Graph.get_e(path_tail, path_head) - min_cut);
        }

        //========= path_head -> path_tail
        // Adding/Update augmented path
        if (Res_Aug_Graph.has_e(path_head, path_tail)) {
            Res_Aug_Graph.update_e(path_head, path_tail, Res_Aug_Graph.get_e(path_head, path_tail) + min_cut);
        }
        else {
            Res_Aug_Graph.add_e(path_head, path_tail, min_cut);
        }
    }
}

inline void Graph::retrieve_FordFulkerson(Graph &Res_Aug_Graph, int &min_cut, vector<int> &Path) {
    // Path without deviation
    vector<int> Path_clean;

    // Update residual and augmented path
    int path_tail, path_head;
    while(true) {

        path_head = Path.back(); Path.pop_back();
        // actuated when the final sink node is popped
        // then path_tail = path_head = s
        if (Path.empty()) {
            break;
        }
        path_tail = Path.back(); //Path.pop_back(); Need to hold this for next iteration (preserve consecutive edges)
        // Backwards searching to tackle with invalid DFSs (Paths not leading to t)
        /*
        Well, after a second thoughts, since I've limited the graph to have only 1 source and 1 sink, 
        There shouldn't be any DFS path that doesn't end with the t node (the unique sink node)
        So probably this backward searching can be neglected.
        Ohhhh nonono, after a third though, loops in the original graph leads to a DFS deviation! 
        So this backward search is still necessary
        */
        while (!Res_Aug_Graph.has_e(path_tail, path_head)) { 
            Path.pop_back();
            path_tail = Path.back();
        }

        // Construct the valid DFS path
        Path_clean.push_back(path_head); //Path_clean.push_back(path_tail); Hold this to prevent extra nodes been appended

        // Update min_cut in the backward DFS
        if (Res_Aug_Graph.get_e(path_tail, path_head) < min_cut) {
            min_cut = Res_Aug_Graph.get_e(path_tail, path_head);
        }
    }
    // Now add the final tail (which is the sink node)
    Path_clean.push_back(path_tail);

    /* For debugging
    cout << "Path_clean: ";
    for (auto i : Path_clean) {
        cout << i << ' ';
    }
    cout << endl;
    */

    update_Res_Aug_Graph(Res_Aug_Graph, min_cut, Path_clean);
}

inline void Graph::DFS_FordFulkerson(Graph &Res_Aug_Graph, int &v, const int s, const int t, int& min_cut, map<const int, bool> &E, vector<int> &Stack, vector<int> &Path, bool &no_s2t_path_remaining) {
            
    Path.push_back(v);
    E[v] = 1;

    for (auto w : Res_Aug_Graph.get_v(v)->get_head_keys()) {

        if (Res_Aug_Graph.get_e(v, w) <= 0 || E[w]) {
            continue;
        }

        // found one path from s to t
        if (w == t) {

            //cout << "c2: " << ++c2 << endl;

            // finalize the path form s to t
            Path.push_back(w); 

            Res_Aug_Graph.retrieve_FordFulkerson(Res_Aug_Graph, min_cut, Path);
            // clear up for next iteration (Path should be empty by now naturally tho, cuz we've iterated the entire path)
            Stack.clear();// Path.clear();  

            // Marking we've found an s-t path, and helps to break the inner while loop
            no_s2t_path_remaining = false;

            // Break this for loop (and the inner while loop), difference paths might share come same nodes, 
            // so if we don't start a new DFS on the updated graph, some paths may be ignored
            break;
        }

        Stack.push_back(w);
    }
}

inline void Graph::build_from_adj_list_without_weight(string dir) {
    //typedef stringstream SS;
    string s;

    if (this->get_size() != 0) {
        cout << "non-empty graph" << endl;
        return;
    }

    fstream fp; fp.open(dir, ios::in);
    if (!fp.is_open()) {
        cout << "File failed to open" << endl;
        return;
    }

    int tail, head;
    while (!fp.eof()) {
        getline(fp, s);
        stringstream ss(s); ss >> tail;     // Write by token

        while (ss >> head) {
            this->add_e(tail, head);
        }
    }

    return;            
}

inline void Graph::build_from_edge_list_with_weight(string dir) {
    string s;

    if (this->get_size() != 0) {
        cout << "non-empty graph" << endl;
        return;
    }

    fstream fp; fp.open(dir, ios::in);
    if (!fp.is_open()) {
        cout << "File failed to open" << endl;
        return;
    }

    int tail, head, weight;
    while(fp >> tail >> head >> weight) {
        this->add_e(tail, head, weight);
    }

    return;
}


inline Graph Graph::undirectedfy() {
    // Since the Verteices are stored as **pointers** in Graph, 
    // So even if G is call by value, the Vertices inside the outside G will be altered
    // therefore we take the copy of G and modify upon that here
    Graph G_temp = this->copy();

    for (int i : this->get_v_keys()) {
        for (int j : this->get_v(i)->get_head_keys()) {
            // assuming no bi-directional path with different weights between every pairs of nodes 
            G_temp.add_e(i, j, this->get_e(i, j), false);
        }
    }

    return G_temp.copy();
}

inline bool Graph::is_connected(bool do_undirectify=true, int fixed_sourse=-1) {
    // Graph G_temp; This raises error of some undifined type conversion things, must do Graph G_temp = this->copy(), i.e. assign the Graph right after it is declared,
    // But not assign value (memory) to the variable in an if-else's scope. 
    /*
    Okay, let's put it this way. When the Graph (G_temp) is declared outside of the if-else's scope (marked by {}), 
    Cpp compiler is supposed to call the variable G_temp automaticly by referece, i.e. Graph &G_temp.
    HOWEVER! The Graph datatype is not native in Cpp so the compiler calls the private variables inside the Graph by reference,
    i.e. map<const int &, Vertex &>, thus causing a conversion error!

    The exact error description is as follows

        candidate function (the implicit copy assignment operator) not viable: 
        no known conversion from 'const pair<const std::__value_type<const int, Vertex *>::key_type, std::__value_type<const int, Vertex *>::mapped_type>' 
        to 'const pair<const int &, Vertex *&>' for 1st argument
        struct _LIBCPP_TEMPLATE_VIS pair

    OMG this is ultterly annoying.

    Maybe using a temporal foo fucntion or bypass the use of of-else loop would solve this
    */

    
    /*
    if (directed) {
        //G_temp = this->undirectedfy();
        Graph G_temp = this->undirectedfy();
    }
    else {
        //G_temp = this->copy();
        Graph G_temp = this->copy();
    }
    */

    // Option 1: traverse all vertecies
    // Option 2: traverse in an undirected graph
    // to prevent missing nodes due to directed edges in DFS. This keeps everything in the same scope (unlike if-else {})
    // I'm a genius
    Graph G_temp = do_undirectify? this->undirectedfy() : this->copy();

    map<const int, bool> E;
    vector<int> vertieces;

    // BFS
    // baically this compares the number of nodes obtained by BFS outsourced
    // from an arbitary node with the number of nodes in the graph
    // so first map every node to false, checking whether there is any remaining false 
    // results in conneted or not

    // initialize Explored map
    if (fixed_sourse != -1) {
        vertieces.push_back(fixed_sourse);
    }
    else {
        vertieces = G_temp.get_v_keys();
    }
    
    for (int i : vertieces) {
        E[i] = false;
    }

    int counter = G_temp.get_size();

    int s = vertieces.front();
    E[s] = true;
    --counter;

    queue<int> Q; Q.push(s);

    while (!Q.empty()) {
        // safe implementation of pop!
        int v = Q.front(); Q.pop();

        for (int e : G_temp.get_v(v)->get_head_keys()) {
            if (!E[e]) {
                E[e] = true;
                --counter;
                Q.push(e);
            }
        }

        // All nodes have been visited, thus the graph is connected
        if (!counter) {
            return true;
        }
    }
    return false;
}

inline map<const int, int> Graph::Dijkstra(const int s) {
    map<const int, int> L;   // distance from s to else
    map<const int, bool> E;     // explored for not
    vector<int> V = this->get_v_keys(); // V is not a pointer to _vertices

    for (int i : V) {
        E[i] = false;
        L[i] = numeric_limits<int>::max();   // numeric_limits<int>::infinity() not working properly
    }

    L[s] = 0;

    int curr_min, u;
    while (V.size()) {
        curr_min = numeric_limits<int>::max();
        
        // find the shortest path (global?)
        for (int v : V) {
            if (L[v] < curr_min) {
                curr_min = L[v];
                u = v;
            }
        }
        
        int index = 0;
        for (int i : V) {
            if (i == u) {
                V.erase(V.begin() + index);
            }
            ++index;
        }
        
        E[u] = true;

        for (int v : this->get_v(u)->get_head_keys()) {
            if (E[v]) {continue;}

            int l = this->get_e(u, v);
            if (L[u] + l < L[v]) {
                L[v] = L[u] + l;
            }
        }
    }
    
    return L;
}

inline deque<int> Graph::TS(map<const int, bool> &E, bool reversed=false) {
    // subfunction -> DFS_TS
    // variables here is passed to DFS_TS
    deque<int> TO;
    vector<int> nodes;
    
    nodes = this->get_v_keys();

    for (auto i : nodes) {
        E[i] = false;
    }

    for (auto i : nodes) {
        if (!E[i]) {
            DFS_TS(i, E, TO, reversed);
        }
    }

    //vector<int> TO_vec = {TO.begin(), TO.end()};

    return TO;
}

inline map<const int, int> Graph::SCC() {
    // subfunction -> DFS_SCC
    // variable here is passed into DFS_SCC
    map<const int, bool> E;
    map<const int, int> SCCs;
    deque<int> TO = this->TS(E, true);
    int numSCC = 0;

    for(auto i : this->get_v_keys()) {
        E[i] = 0;
        SCCs[i] = 0;
    }
    
    for (auto v : TO) {
        if (!E[v]) {
            ++numSCC;
            DFS_SCC(v, E, TO, SCCs, numSCC);
        }
    }

    return SCCs;
}

inline map<const int, bool> Graph::get_init_explored_map() {
    map<const int, bool> E_init;
    for (auto temp_pair : this->_vertices) {
        /*
        if (temp_pair.second->get_tail_keys().empty()) {
            // source node
            ++n_s;
            //memcpy(&s, a_pair.second, sizeof(Vertex));
            s = temp_pair.first;
        }
        if (temp_pair.second->get_head_keys().empty()) {
            // sink node (terminal)
            ++n_t;
            //memcpy(&t, a_pair.second, sizeof(Vertex));
            t = temp_pair.first;
        }
        */

        E_init[temp_pair.first] = 0;
    }

    return E_init;
}

inline int Graph::FordFulkerson(const int s, const int t) {
    int v;    // source, terminate (sink node), regular vertex
    int n_s = 1/*0*/, n_t = 1/*0*/;
    
    map<const int, bool> E_init = get_init_explored_map();
    
    
    // For debugging
    //cout << "s: " << s << " t: " << t << endl; 
    
    // _key = -1 meaning the node is not properly initialized
    //if (s.get_key() == -1 || t.get_key() == -1) {
    //    cout << "sourse/sink node not found" << endl;
    //    return -1;
    //}
    /*
    if ((n_s-1) || (n_t-1)) {
        cout << "sourse/sink node not found or more than one are found" << endl;
        return -1;// *(new Graph);
    }
    */

    Graph Res_Aug_Graph = this->copy();
    vector<int> Stack, Path;  // for adding augmented edge

    // see any reamining path from s -> t
    bool no_s2t_path_remaining = false;

    // Checkpoints for debugging
    //int c0 = 0, c1 = 0, c2 = 0;

    // Exlaination:
    // If the no_s2t_path_remaining stays true, i.e. the if statement for where no_s2t_path_remaining is assigned with false is never entered, 
    // Then the Stack will first be emptyed, breaking the inner loop, and the no_s2t_path_remaining = true breaks the outer while loop.
    while (!no_s2t_path_remaining){//Res_Aug_Graph.FordFulkerson_is_connected(s, t)) {
        //cout << "c0: " << ++c0 << endl;

        // Stack should be empty right now
        // Start a new DFS on the updated Graph
        Stack.push_back(s);
        no_s2t_path_remaining = true;

        map<const int, bool> E (E_init);

        // Note that the min_cut shouldn't be updated in the process of forward DFS
        // (cuz there will be invalid deviation from the final path)
        // Update this in the backward retrieval
        int min_cut = numeric_limits<int>::max();

        // Exlaination:
        // !Stack.empty(): on the last iteration (where there is no s-t path remaining, Stack will be empty (else won't))
        // no_s2t_path_remaining: no_s2t_path_remaining? assume it is true and let's find out in this while loop!
        while (no_s2t_path_remaining && !Stack.empty()) {

            //cout << "c1: " << ++c1 << endl;

            v = Stack.back(); Stack.pop_back(); 
            // Not necessary, because in the for loop below, 
            // an already visited node will never be appened to the Stack
            /*
            if (E[v]) {
                continue;
            }
            */

            Res_Aug_Graph.DFS_FordFulkerson(Res_Aug_Graph, v, s, t, min_cut, E, Stack, Path, no_s2t_path_remaining);
        }
    }

    // And by now the Res_Aug_Graph is completed, without any remaining s-t path.
    // Update the original Graph only on edges that is included in the original Graph (augmented edges excluded)

    // But let's First have a bit of test on the Res_Aug_Graph.

    //cout << "Fini" << endl;

    // if don't return .copy, all the memory of Vertex* will be lost
    //return Res_Aug_Graph.copy();

    //===================
    // NOW! Let's fucking find the actual MAX FLOW!!!
    // Flow = Capacity - Residule
    Graph Flow_Graph = this->copy();

    for (auto v : Flow_Graph.get_v_keys()) {
        for (auto w : Flow_Graph.get_v(v)->get_head_keys()) {
            if (Res_Aug_Graph.has_e(v, w)) {
                Flow_Graph.update_e(v, w, Flow_Graph.get_e(v, w) - Res_Aug_Graph.get_e(v, w));
            }
        }
    }

    int MAX_FLOW = 0;
    for (auto w : Flow_Graph.get_v(s)->get_head_keys()) {
        MAX_FLOW += Flow_Graph.get_e(s, w);
    }

    return MAX_FLOW;
}


