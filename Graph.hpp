#ifndef GRAPH_HPP_INCLUDED
#define GRAPH_HPP_INCLUDED

/*dependencies*/
#include <unordered_set>
#include <unordered_map>
#include <list>
#include <vector>
#include <utility> //pair
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

typedef unsigned int Vertex;

class Graph {

    public:

        /*intermediate structures*/

         //edge in an adjacency list (one end-vertex is implicit)
        struct Edge {
            double _weight; //edge weight
            list<Vertex>::iterator _pos; //position in the adjacency list
        };

        //neighborhood
        struct AdjList {
            list<Vertex> _list; //adjacency list
            unordered_map<Vertex, Edge> _map; //index
        };


    private:

        /*flags*/

        bool directed; //directed flag (the graph is directed if set to true)
        bool weighted; //weighted flag (the graph is edge-weighted if set to true)

        /*internal structure*/

        //number of edges (or arcs, if directed)
        int edges_nr;

        //out-neighborhood
        unordered_map<Vertex, AdjList> AdjOut;
        //in-neighborhood (unused if undirected)
        unordered_map<Vertex, AdjList> AdjIn;

        //set of all vertices
        unordered_set<Vertex> vertices;

        //garbage: removed vertices
        list<Vertex> dead_vertices;
        //generation of a default ID
        Vertex gen_id(void);

    public:

        //constructors
		Graph(bool,bool = 0); //construction of a undirected/directed unweighted/weighted null graph
		Graph(); //construction of a null undirected unweighted graph
		Graph(const Graph&); //copy constructor

        //graph information
		inline unsigned int nr_vertices(void) const { return vertices.size(); } //number of nodes
		inline unsigned int nr_edges(void) const { return edges_nr; } //number of edges
		inline bool is_directed(void) const { return directed; }
		inline bool is_weighted(void) const { return weighted; }

        //edge information
		list<pair<Vertex, Vertex>> list_edges(void) const; //enumerate all arcs
		double weight(Vertex, Vertex) const; //arc weight (1 if undirected, undefined if arc does not exist)
		bool is_adjacent(Vertex, Vertex) const;

        //vertex information
		list<Vertex> vertex_list(void) const; //enumerate all vertices
		int degree(Vertex) const; //total degree of a node (-1 if node does not exist)
		int out_degree(Vertex) const;
		int in_degree(Vertex) const;
		list<Vertex> out_neighbours(Vertex) const; //enumerate all (out-)neighbours of a node
		list<Vertex> in_neighbours(Vertex) const;

        //dynamic operations
		bool add_edge(Vertex, Vertex, double = 1); //adds a new edge, if it is not already present
		bool remove_edge(Vertex, Vertex); //removes an edge, if it is already present

		bool add_node(Vertex); //inserts a new isolated vertex, if it does not already exist
		Vertex add_node(void); //inserts a node, with an automatically generated id, and outputs this id
		bool remove_node(Vertex); //removes a node, being given its ID

		void reorder(list<Vertex>); //reorder adjacency list according to some vertex permutation

		bool contract(Vertex,Vertex,Vertex); //contract two adjacent vertices
		int contract(Vertex,Vertex); //contract two adjacent vertices, and outputs the ID of the resulting node

    private:

        void insert(AdjList&, Vertex, double);
		void erase(AdjList&, Vertex);
		void make_first(AdjList&, Vertex);

    friend class ParallelGraphAlgorithms;

    public:

        Graph read_file(string);
        void print_graph(void);
        void print_edges(void);
        void print_weighted_edges(void);

};

/*---------------------------------------- Parallel Algorithms ----------------------------------*/

class ParallelGraphAlgorithms {

	private:
		ParallelGraphAlgorithms();

	private:
		static void make_directed(Graph&); // make Graph directed
		static void isolate(unsigned int, Graph&); // isolate vertex; remove edgesIn from his neighbours
		static void pdfs_rec(unsigned int, Graph&);

	public:

		/*-------------- Edge Isolations -----------------------*/

		static void parallel_dfs(unsigned int, Graph);
		static void parallel_bfs(unsigned int, Graph);

		/*-------------- PBFS --------------------------------
		@ssumption 1: vertices are indexed from 0 to n-1
		@ssumption 2: vectors parent and distance are pre-constructed.
		*/

		struct bfs_t {

			vector<int> _p, _d; //parent, distance

			bfs_t(unsigned n): _p(n), _d(n)	{}
		};

		static void simple_pbfs(unsigned int, Graph, bfs_t&);

		static void D_1_distributed_memory_BFS(unsigned int, Graph, bfs_t&);

		static void a_star_hash_distributed(unsigned int, Graph, bfs_t&);




    /*------------------- Connected Components ----------------------
    @ssumption 1: vertices are indexed from 0 to n-1
    @ssumption 2: vectors p,o are pre-constructed.
    */

    private:
        static bool my_direct_connect(Graph&, vector<unsigned int>&, vector<unsigned int>&); //Lock-free
		static bool my_parent_connect(Graph&, vector<unsigned int>&, vector<unsigned int>&);
		static bool my_direct_root_connect(Graph&, vector<unsigned int>&, vector<unsigned int>&); //Lock-free
		static bool my_parent_root_connect(Graph&, vector<unsigned int>&, vector<unsigned int>&);

    public:
        static void my_algorithm_S(Graph, vector<unsigned int>&, vector<unsigned int>&);
        static void my_algorithm_RA(Graph, vector<unsigned int>&, vector<unsigned int>&);
        static void my_algorithm_A(Graph, vector<unsigned int>&, vector<unsigned int>&);

        static bool shortcut(Graph&, vector<unsigned int>&, vector<unsigned int>&); //Lock-free

		static void alter(Graph&, vector<unsigned int>&);

        static unordered_map<unsigned int, Graph> my_algorithm_R(Graph, vector<unsigned int>&, vector<unsigned int>&);

};

#endif // GRAPH_HPP_INCLUDED
