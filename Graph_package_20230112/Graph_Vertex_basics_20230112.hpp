/*
    Created: 20221211
    Last modification: 20230105
    This package includes the definitions of basic structure and 
    declerations of algorithms
*/

#include "Graph_Vertex_includes_20230112.hpp"


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

        //=================== 
        // declaration of helper functions
        //===================
        void DFS_TS(int s, map<const int, bool> &E, deque<int>&TO, bool reversed);

        void DFS_SCC(int s, map<const int, bool> &E, deque<int>&TO, map<const int, int> &SCCs, int &numSCC);
        
        // returns the default E dict for each DFS iteration in FordFulkerson
        map<const int, bool> get_init_explored_map();
        // execute the reverse DFS in FordFulkerson
        // obtain clean s-t path, record min_cut and update the Res_Aug_Graph
        void retrieve_FordFulkerson(Graph &Res_Aug_Graph, int &min_cut, vector<int> &Path);
        void update_Res_Aug_Graph(Graph &Res_Aug_Graph, int &min_cut, vector<int> &Path_clean);
        void DFS_FordFulkerson(Graph &Res_Aug_Graph, int &v, const int s, const int t, int& min_cut, map<const int, bool> &E, vector<int> &Stack, vector<int> &Path, bool &no_s2t_path_remaining);
        void update_path(const int u, const int v, map<const int, vector<int>*> &P);

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

        /*=================== helpers
            The main structure of the Graph/Vertex is coded in C, but from now on shifted to C++, so C++ syntax is adopted for future works
            All helper functions should be void and do not call variables throuth this->.. 
        ===================*/

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


        //=================== declarance of advanced funcitons

        // This works really slow, but indeed works lol    (TO_BE_REDONE)
        // Downwards supported for building from edge list
        void build_from_adj_list_without_weight(string dir);

        void build_from_edge_list_with_weight(string dir);

        // This function returns an undirected copy of the original Graph
        Graph undirectedfy();
        
        bool is_connected(bool do_undirectify, int fixed_sourse);
        vector<int> get_connected_vertices_with_fixed_source(int fixed_sourse);

        //=================== Dijkstra's algorithm

        // Graph has to be connected else results in infinite loop
        map<const int, int> Dijkstra(const int s);
        pair<map<const int, int>, map<const int, vector<int>*>> Dijkstra_with_path(const int s);
        //=================== End of Dijkstra


        //=================== SCC algorithm

        // Use deque to insert in front
        deque<int> TS(map<const int, bool> &E, bool reversed);

        map<const int, int> SCC();

        //=================== End of SCC algorithm

        
        //=================== FordFulkerson
        // I coded this up really chunckily, just use EdmondsKarp
        int FordFulkerson(const int s, const int t);
        //=================== End of FordFulkerson


        //=================== Edmonds–Karp
        int EdmondsKarp(const int s, const int t);
        //=================== End of Edmonds–Karp


        //=================== Hungarian Algorithm
        // Operates on adj matrix
        vector<vector<int>> get_adj_matrix(vector<int> column, vector<int> row_v);
        int Hungarian();
        //=================== End of Hungarian Algorithm
};

