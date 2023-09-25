#include "Graph.hpp"

//constructors

Graph::Graph(bool d, bool w) : directed(d), weighted(w), edges_nr(0) {}

Graph::Graph() : directed(0), weighted(0), edges_nr(0) {}

Graph::Graph(const Graph& G) : directed(G.is_directed()), weighted(G.is_weighted()), edges_nr(0) {
	for(Vertex v : G.vertex_list()) add_node(v);
	for(auto e : G.list_edges()) add_edge(e.first, e.second, G.weight(e.first, e.second));
}

//edge access

list<pair<Vertex, Vertex>> Graph::list_edges(void) const {
	list<pair<Vertex,Vertex>> _all;
	for(auto Adj: this->AdjOut){
		Vertex v = Adj.first; //non-isolated vertices
		for(Vertex u: Adj.second._list){ //neighbours, in order
			if(is_directed() || v < u)
			 _all.push_back(pair<Vertex,Vertex>(v,u));
		}
	}
	return _all;
}

bool Graph::is_adjacent(Vertex u, Vertex v) const {
	if(AdjOut.find(u) == AdjOut.end()) return 0;
	return (AdjOut.at(u)._map.find(v) != AdjOut.at(u)._map.end());
}

double Graph::weight(Vertex u, Vertex v) const {
	if(!is_adjacent(u, v)) return numeric_limits<double>::infinity();
	return AdjOut.at(u)._map.at(v)._weight;
}

//vertex access

list<Vertex> Graph::vertex_list(void) const {
	list<Vertex> _all;
	for(Vertex v: vertices)
        _all.push_back(v);
	return _all;
}

int Graph::out_degree(Vertex v) const{
	if(vertices.find(v) == vertices.end()) return -1;
	else if(AdjOut.find(v) == AdjOut.end()) return 0;
	else return AdjOut.at(v)._list.size();
}

int Graph::in_degree(Vertex v) const{
	if(!is_directed()) return out_degree(v);
	else if(vertices.find(v) == vertices.end()) return -1;
	else if(AdjIn.find(v) == AdjIn.end()) return 0;
	else return AdjIn.at(v)._list.size();
}

int Graph::degree(Vertex v) const {
	if(vertices.find(v) == vertices.end()) return -1;
	else if(is_directed()) return out_degree(v)+in_degree(v);
	else return out_degree(v);
}

list<Vertex> Graph::out_neighbours(Vertex v) const {
	list<Vertex> out_list;
	if(AdjOut.find(v) != AdjOut.end())
        for(Vertex u : AdjOut.at(v)._list)
            out_list.push_back(u);
	return out_list;
}

list<Vertex> Graph::in_neighbours(Vertex v) const {
	if(!is_directed()) return out_neighbours(v);
	list<Vertex> in_list;
	if(AdjIn.find(v) != AdjIn.end())
        for(Vertex u : AdjIn.at(v)._list)
            in_list.push_back(u);
	return in_list;
}

//dynamic operations

void Graph::insert(Graph::AdjList& adjList, Vertex v, double w){ // insert edge e to vertex v in AdjList
	adjList._list.push_front(v); //insert v in adj list at the begging
	Graph::Edge e; //create edge
	e._weight = w; //set edge weight
	e._pos = adjList._list.begin(); //set position in adj list
	adjList._map[v] = e;
}

void Graph::erase(Graph::AdjList& adjList, Vertex v){
	Graph::Edge e = adjList._map[v]; // find edge with vertex v
	adjList._map.erase(v); // delete vertex and edge from map
	adjList._list.erase(e._pos); // delete vertex
}

void Graph::make_first(Graph::AdjList& adjList, Vertex v){
	adjList._list.push_front(v); // add vertex v to the beginning
	adjList._list.erase(adjList._map[v]._pos); // delete vertex v from old position
	adjList._map[v]._pos = adjList._list.begin(); // set new position to beginning
}

bool Graph::add_edge(Vertex u, Vertex v, double w){
	if(vertices.find(u) == vertices.end() || vertices.find(v) == vertices.end() || is_adjacent(u,v) || u == v)
        return 0;
	if(!is_weighted() && w != 1)
        return 0;

	insert(AdjOut[u], v, w);

	if(is_directed())
        insert(AdjIn[v], u, w);
	else
        insert(AdjOut[v], u, w);

	edges_nr++;
	return 1;
}

bool Graph::remove_edge(Vertex u, Vertex v) {
	if(!is_adjacent(u,v)) return 0;

	erase(AdjOut[u], v);
    if(out_degree(u) == 0) AdjOut.erase(u);

	if(is_directed()){
        erase(AdjIn[v], u);
        if(in_degree(v) == 0)
            AdjIn.erase(v);
	} else {
        erase(AdjOut[v], u);
        if(out_degree(v)==0) AdjOut.erase(v);
	}

	edges_nr--;
	return 1;
}

bool Graph::add_node(Vertex v) {
	if(vertices.find(v) != vertices.end()) return 0;
	vertices.insert(v);
	return 1;
}

bool Graph::remove_node(Vertex v) {
	if(vertices.find(v) == vertices.end()) return 0;
	for(Vertex u: out_neighbours(v))
        remove_edge(v, u);
	for(Vertex u: in_neighbours(v))
        remove_edge(u, v);
	vertices.erase(v);
	dead_vertices.push_back(v);
	return 1;
}

void Graph::reorder(list<Vertex> perm){
	perm.reverse(); // reverse order of the list
	for(Vertex v : perm){

        for(Vertex u: in_neighbours(v))
            make_first(AdjOut[u],v);

		if(is_directed())
            for(Vertex u: out_neighbours(v))
                make_first(AdjIn[u],v);

	}
}

bool Graph::contract(Vertex x, Vertex y, Vertex z) {

	if(!is_adjacent(x,y)) return 0;
	if(z != x && z != y && vertices.find(z) != vertices.end()) return 0; // z alt nod din graf in afara de x sau y

	list<Vertex> out_x = out_neighbours(x), in_x = in_neighbours(x), out_y = out_neighbours(y), in_y = in_neighbours(y);
	remove_node(x); remove_node(y);
	add_node(z);
	for(Vertex v: out_x) add_edge(z,v);
	for(Vertex v: out_y) add_edge(z,v);
	for(Vertex u: in_x) add_edge(u,z);
	for(Vertex u: in_y) add_edge(u,z);
	return z;

}

int Graph::contract(Vertex u, Vertex v){
	Vertex z = (degree(u) <= degree(v)) ? u : v; // if true then u else v
	if(contract(u, v, z)) return z;
	return -1;
}

Vertex Graph::gen_id(void) {
	if(dead_vertices.empty()) return nr_vertices(); // return number of nodes
	Vertex zombie = dead_vertices.front();
	dead_vertices.pop_front();
	return zombie;
}

Vertex Graph::add_node(void){
	Vertex v = gen_id();
	add_node(v);
	return v;
}




void Graph::print_graph(void) {
    cout << "Graph:" << endl;
    for(auto adj : this->AdjOut) {
        cout << adj.first;
        for(auto v : adj.second._list) {
            cout << " -> " << v;
        }
        cout << endl;
    }
}

void Graph::print_edges(void) {
    cout << "Edges:" << endl;
    for(auto adj : this->AdjOut) {
        for(auto v : adj.second._list) {
            if (adj.first < v)
            cout << adj.first << " " << v << endl;
        }
    }
}

void Graph::print_weighted_edges(void) {
    cout << "Edges:" << endl;
    for(auto adj : this->AdjOut) {
        for(auto v : adj.second._list) {
            if (adj.first < v)
            cout << adj.first << " " << v << " " << weight(adj.first, v) << endl;
        }
    }
}


// Read file and create a graph
Graph Graph::read_file(string filename) {

    cout << "File name: " << filename << endl;

    // Open file
    ifstream read_file(filename);
    string text;

    // Get first line with nr_vertices and nr_edges
    getline(read_file, text);
    stringstream text_stream_firs_tline;
    unsigned int nr_vertices;
    unsigned int nr_edges;
    text_stream_firs_tline << text;
    text_stream_firs_tline >> nr_vertices;
    //cout << "nr_vertices = " << nr_vertices << endl;
    text_stream_firs_tline >> nr_edges;
    //cout << "nr_edges = " << nr_edges << endl;
    //cout << nr_vertices << " " << nr_edges << endl;
    Graph G;

    // Get rest of the line with edges
    bool graph_initialized = 0;
    while (getline(read_file, text)) {

        stringstream text_stream;
        unsigned int v;
        unsigned int u;
        double w;
        vector<unsigned int> numbers;
        int number;

        text_stream << text;

        // read a line
        while(!text_stream.eof()) {
            text_stream >> number;
            numbers.push_back(number);
        }

        // call the right constructor
        if (graph_initialized == 0) {
            if (numbers.size() == 2) {
                G = Graph();
            } else if (numbers.size() == 3) {
                G = Graph(false, true);
            }
            graph_initialized = 1;
        }

        // add nodes and edges
        if (numbers.size() == 2) {
            v = numbers[0];
            u = numbers[1];
            G.add_node(v);
            G.add_node(u);
            G.add_edge(v, u);
            //cout << v << " " << u << endl;
        } else if (numbers.size() == 3) {
            v = numbers[0];
            u = numbers[1];
            w = numbers[2];
            G.add_node(v);
            G.add_node(u);
            G.add_edge(v, u, w);
            //cout << v << " " << u << " " << w << endl;
        }

    }
    // Close file
    read_file.close();

    return G;

}
