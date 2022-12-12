/*
    Created: 20221211
    This update focues on second version of FordFulkerson's algorithm.
    Mainly on reduce nesting depth and remove redundant functionality/loops
*/

#include <iostream>
#include <map>
#include <vector>
#include <queue>
#include <deque>
#include <limits>
#include <fstream>
#include <sstream>

using namespace std;


class Vertex {
    // __init__ (python)
    private:
        const int _key;
        map<const int, int> _head_of, _tail_of; 


    public:
        // constructor (initialize _key)
        Vertex(int key = -1): _key(key) {}
        // destructor
        ~Vertex() {_head_of.clear(); _tail_of.clear();}

        inline const int get_key() {
            return this->_key;
        }

        inline void add_head(const int head_key, int weight=1) {
            if (weight > -1) {
                _tail_of[head_key] = weight;
            }
            else {
                cout << "Add head error: Invalid weight" << endl;
            }
        }

        inline void add_tail(const int tail_key, int weight=1) {
            if (weight > -1) {
                _head_of[tail_key] = weight;
            }
            else {
                cout << "Add tail error: Invalid weight" << endl;
            }
        }

        void __str__() {
            cout << "key: " << _key << endl;
            cout << "(key)" << "\t" << "(weight)" << endl;
            cout << "tail of: " << endl;
            
            if (_tail_of.begin() == _tail_of.end()) {
                cout << "None" << endl;
            }
            else{
                for (map<const int, int>::iterator i = _tail_of.begin(); i != _tail_of.end(); ++i) {
                cout << i->first << "\t" << i->second << endl;
                }
            }
            
            cout << "head of: " << endl;
            if (_head_of.begin() == _head_of.end()) {
                cout << "None" << endl;
            }
            else {
                for (map<const int, int>::iterator i = _head_of.begin(); i != _head_of.end(); ++i) {
                cout << i->first << "\t" << i->second << endl;
                }
            }
            cout << endl;
        }

        inline bool tail_of(const int head_key) {
            return (_tail_of.find(head_key) != _tail_of.end());
        }

        inline bool head_of(const int tail_key) {
            return (_head_of.find(tail_key) != _head_of.end());
        }

        inline void remove_tail(const int tail_key) {
            if (head_of(tail_key)) {
                _head_of.erase(tail_key);
            }
        }

        inline void remove_head(const int head_key) {
            if (tail_of(head_key)) {
                _tail_of.erase(head_key);
            }
        }

        // return the pointer to the map with tails/heads: weights
        inline vector<int> get_head_keys() {
            vector<int> head_keys;
            for (map<const int, int>::iterator iter_map_int_int = _tail_of.begin(); iter_map_int_int != _tail_of.end(); ++iter_map_int_int) {
                head_keys.push_back(iter_map_int_int->first);
            }

            return head_keys;
        }

        inline vector<int> get_tail_keys() {
            vector<int> tail_keys;
            for (map<const int, int>::iterator iter_map_int_int = _head_of.begin(); iter_map_int_int != _head_of.end(); ++iter_map_int_int) {
                tail_keys.push_back(iter_map_int_int->first);
            }

            return tail_keys;
        }

        inline int get_tail_weight(int tail_key) {
            if (head_of(tail_key)) {
                return _head_of[tail_key];
            }
            throw invalid_argument("Invalid weight look up: No such tail");
        }

        inline int get_head_weight(int head_key) {
            if (tail_of(head_key)) {
                return _tail_of[head_key];
            }
            throw invalid_argument("Invalid weight look up: No such head");
        }

        inline map<const int, int>* get_head_keys_map() {
            return &_tail_of;
        }

        inline map<const int, int>* get_tail_keys_map() {
            return &_head_of;
        }
};



class Graph {
    private:
        map<const int, Vertex*> _vertices;

        // Used in destructure, static member funciton is required
        // argument called by reference in for_each()
        static void delete_ptr_in_map_values(pair<const int, Vertex*> &pair_with_ptr_value)
        {
            delete pair_with_ptr_value.second;
        }
        

        void DFS_TS(int s, map<const int, bool> &E, deque<int>&TO, bool reversed=false) {
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

        void DFS_SCC(int s, map<const int, bool> &E, deque<int>&TO, map<const int, int> &SCCs, int &numSCC) {
            E[s] = 1;
            SCCs[s] = numSCC;

            for (auto v : this->get_v(s)->get_head_keys()) {
                if (!E[v]) {
                    DFS_SCC(v, E, TO, SCCs, numSCC);
                }
            }
        }
        
        // Returns ture if still path form s to t in Res_Aug_Graph
        /* Not needed anymore, integrated in the iteration of FordFulkerson algorithm
        bool FordFulkerson_is_connected(int s, int t) {
            vector<int> Stack {t};
            map<const int, int> E;
            int head;

            for (auto v : this->get_v_keys()) {
                E[v] = 0;
            }
            
            while (!Stack.empty()){
                head = Stack.back(); Stack.pop_back(); E[head] = 1;
                for (auto tail : this->get_v(head)->get_tail_keys()) {
                    if (tail == s) {
                        return true;
                    }
                    if (!E[tail]) {
                        Stack.push_back(tail);
                    }
                }
            }

            // s is never visited
            return false;
        }
        */


        void update_Res_Aug_Graph(Graph &Res_Aug_Graph, int &min_cut, vector<int> &Path_clean) {
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


        void retrieve_FordFulkerson(Graph &Res_Aug_Graph, int &min_cut, vector<int> &Path) {
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


        void DFS_FordFulkerson(Graph &Res_Aug_Graph, int &v, const int s, const int t, int& min_cut, map<const int, bool> &E, vector<int> &Stack, vector<int> &Path, bool &no_s2t_path_remaining) {
            
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


    public:
        ~Graph() {
            for_each(_vertices.begin(), _vertices.end(), delete_ptr_in_map_values);
            _vertices.clear();
        }


        inline int get_size() {
            return _vertices.size();
        }

        inline bool has_v(const int key) {
            return _vertices.find(key) != _vertices.end();
        }

        inline void add_v(Vertex *v) {
            if (this->has_v(v->get_key())) {
                cout << "Error: duplicated vertex" <<endl;
                delete v;   // prevent memory leakage
                return;
            }
            _vertices[v->get_key()] = v;
        }

        inline Vertex* get_v(const int key) {
            if (_vertices.find(key) != _vertices.end()) {
                return _vertices[key];
            }
            throw invalid_argument("Get vertex error: No such Node");
        }

        inline vector<int> get_v_keys() {
            vector<int> v_keys;
            for (map<const int, Vertex*>::iterator iter_int_pVertex = _vertices.begin(); iter_int_pVertex != _vertices.end(); ++iter_int_pVertex) {
                v_keys.push_back(iter_int_pVertex->first);
            }
            return v_keys;
        }

        inline void add_e(const int tail_key, const int head_key, int weight=1, bool directed=true, bool ingore_weight_error=false) {
            if (weight < 0 && !ingore_weight_error) {
                cout << "Add edge error: Invalid weight" << endl;
                return;
            }

            if (!has_v(head_key)) {
                add_v(new Vertex(head_key));
            }
            if (!has_v(tail_key)) {
                add_v(new Vertex(tail_key));
            }

            // add_tail / add_head is assigning key-value pair to dictionary
            // so there's no need to be clear about whether it has the tail / head or not first
            // as long as the tail / head is in the Graph already
            get_v(head_key)->add_tail(tail_key, weight);
            get_v(tail_key)->add_head(head_key, weight);            
            if (!directed) {
                get_v(head_key)->add_head(tail_key, weight);
                get_v(tail_key)->add_tail(head_key, weight);
            }
        }

        inline bool has_e(const int tail_key, const int head_key) {
            if (!has_v(tail_key) || !has_v(head_key)) {
                return false;
            }
            else if (!get_v(head_key)->head_of(tail_key) || !get_v(tail_key)->tail_of(head_key)) {
                return false;
            }

            return true;
        }

        inline int get_e(const int tail_key, const int head_key) {
            if (!has_v(tail_key) || !has_v(head_key)) {
                throw invalid_argument("Get edge error: No such edge: Missing nodes");
            }
            else if (!get_v(head_key)->head_of(tail_key) || !get_v(tail_key)->tail_of(head_key)) {
                cout << "tail_key: " << tail_key << " " << "head_key: " <<head_key << endl;
                throw invalid_argument("Get edge error: No such edge: Non-neighboring nodes");
            }

            return get_v(tail_key)->get_head_weight(head_key);
        }

        inline void update_e(const int tail_key, const int head_key, int weight=1, bool directed=true) {
            if (!has_v(tail_key) || !has_v(head_key)) {
                throw invalid_argument("Update edge error: No such edge: Missing nodes");
            }
            else if (!get_v(head_key)->head_of(tail_key) || !get_v(tail_key)->tail_of(head_key)) {
                throw invalid_argument("Update edge error: No such edge: Non-neighboring nodes");
            }

            get_v(tail_key)->get_head_keys_map()->find(head_key)->second = weight;
            get_v(head_key)->get_tail_keys_map()->find(tail_key)->second = weight;
            if (!directed) {
                get_v(tail_key)->get_tail_keys_map()->find(head_key)->second = weight;
                get_v(head_key)->get_head_keys_map()->find(tail_key)->second = weight;
            }
        }

        inline void remove_e(const int tail_key, const int head_key, bool directed=true) {
            remove_head_of_tail(tail_key, head_key);
            remove_tail_of_head(tail_key, head_key);
            // Change tail to head and head to tail in the input
            if (!directed) {
                remove_head_of_tail(head_key, tail_key);
                remove_tail_of_head(head_key, tail_key);
            }
        }

        inline void remove_v(const int key, bool directed=true) {
            if (has_v(key)) {
                vector <int> head_keys = get_v(key)->get_head_keys();
                vector <int> tail_keys = get_v(key)->get_tail_keys();
                for (int i : head_keys) {
                    // this modifies the data in the head node
                    remove_tail_of_head(key, i);
                    if (!directed) {
                        // this also modifies only the data in the head node
                        remove_head_of_tail(i, key);
                    }
                }
                

                // this step spares the need for deleting heads/tails stores in the node to be deleted
                delete this->get_v(key);
                _vertices.erase(key);
            }
        }

        inline void remove_tail_of_head(const int tail_key, const int head_key) {
            if (has_v(head_key)) {
                get_v(head_key)->remove_tail(tail_key);
            }
            else {
                cout << "Remove tail of head error: No head node" << endl;
            }
        }

        inline void remove_head_of_tail(const int tail_key, const int head_key) {
            if (has_v(tail_key)) {
                get_v(tail_key)->remove_head(head_key);
            }
            else {
                cout << "Remove head of tail error: No tail node" << endl;
            }
        }

        //=================== The main structure of the Graph/Vertex is coded in C, but from now on shifted to C++, so C++ syntax is adopted for future works

        // deep copy (modification-safe)
        // This is very important for member function to return a Graph
        // Because of I stored Vertex* instead of Vertex itself in the Graph._vertices
        Graph copy() {
            Graph G_temp;
            // traverse though all nodes in case the graph is disconnected so that 
            // add_e won't include all the vertices
            // PS: now add_e can automatically add new vertex to the Graph, so the !G_temp.has_v(i) isn't necessary anymore
            
            for (int i : this->get_v_keys()) {
                /*
                if (!G_temp.has_v(i)) {
                    G_temp.add_v(new Vertex(i));  
                }
                */

                for (int j : this->get_v(i)->get_head_keys()) {
                    G_temp.add_e(i, j, this->get_e(i, j));
                }
                for (int k : this->get_v(i)->get_tail_keys()) {
                    G_temp.add_e(k, i, this->get_e(k, i));
                }
            }

            return G_temp;
        }

        // This works really slow, but indeed works lol    (TO_BE_REDONE)
        // Downwards supported for building from edge list
        void build_from_adj_list_without_weight(string dir) {
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

        void build_from_edge_list_with_weight(string dir) {
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

        // This function returns an undirected copy of the original Graph
        Graph undirectedfy() {
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

        
        bool is_connected(bool do_undirectify=true, int fixed_sourse=-1) {
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
        
        
        


        //=================== Dijkstra's algorithm

        // Graph has to be connected else results in infinite loop
        map<const int, int> Dijkstra(const int s) {
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

        //=================== End of Dijkstra


        //=================== SCC algorithm

        // The DFSs are coded with recursion (tho I prefer iteration with Stack lol)
        // On 20221211 I changed my mind, iteration causes too many nesting loops, looks sucks

        //private: void DFS_TS(int s, map<const int, bool> &E, deque<int>&TO, bool reversed=false)
        // for reading convenience
        /*
        void DFS_TS(int s, map<const int, bool> &E, deque<int>&TO, bool reversed=false) {
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
        */


        // Use deque to insert in front
        deque<int> TS(map<const int, bool> &E, bool reversed=false) {
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

        // private: void DFS_SCC(int s, map<const int, bool> &E, deque<int>&TO, map<const int, int> &SCCs, int &numSCC)
        // for reading convenience
        /*
        void DFS_SCC(int s, map<const int, bool> &E, deque<int>&TO, map<const int, int> &SCCs, int &numSCC) {
            E[s] = 1;
            SCCs[s] = numSCC;

            for (auto v : this->get_v(s)->get_head_keys()) {
                if (!E[v]) {
                    DFS_SCC(v, E, TO, SCCs, numSCC);
                }
            }
        }
        */

        map<const int, int> SCC() {
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

        //=================== End of SCC algorithm


        /*
        Graph MST_Prim(const int s) {
            Graph MST = this->copy();

            
        }
        */


        //=================== FordFulkerson
        /*
        Only for directed graphs, no back and forward edges
        The weight is treated as the capacity of the edge
        */

        // Helper functions for FordFulkerson
        map<const int, bool> get_init_explored_map() {
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


        int FordFulkerson(const int s, const int t) {
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

        //=================== End of FordFulkerson
};