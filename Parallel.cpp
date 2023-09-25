#include "Graph.hpp"
#include <omp.h>
#include <queue>
#include <iostream>
#include <unordered_set>
#include <float.h>

using namespace std;

void ParallelGraphAlgorithms::make_directed(Graph& G) { // make Graph G directed
	if(!G.is_directed()){ // check if Graph G is undirected

		G.directed = 1; // set Graph G as directed;

		for(unsigned int v : G.vertex_list()) // get every vertex in graph
            G.AdjIn[v]; //create inlist

		for(unsigned int v : G.vertex_list()){ // get every vertex in graph
			if(G.degree(v) > 0){ // if node has edges
				#pragma omp parallel for
				for(unsigned int i = 0 ; i < G.AdjOut[v]._map.bucket_count(); i++){
					for(auto it = G.AdjOut[v]._map.begin(i); it != G.AdjOut[v]._map.end(i); it++){ // begin(i) = first element of bucket i
				 		G.insert(G.AdjIn[it->first], v, 1); // insert vertex in the AdjIn; it->first = vertex; AdjIn[it->first] = AdjList; 1 = unweighted
				 	}
				}
			}
		}
	}
}

void ParallelGraphAlgorithms::isolate(unsigned int v, Graph& G){
	if(G.in_degree(v) > 0){ // checks if the are any edges to V
		Graph::AdjList& Nv = G.AdjIn[v]; //Nv = neighbours of v
		#pragma omp parallel for
		for(unsigned int i = 0 ; i < Nv._map.bucket_count(); i++){
			for(auto it = Nv._map.begin(i); it != Nv._map.end(i); it++){
				G.erase(G.AdjOut[it->first], v);
			}
		}
	}
}

void ParallelGraphAlgorithms::pdfs_rec(unsigned int v, Graph& G){

	cout << v << ' ';

	isolate(v, G); // isolate(v)

	while(G.out_degree(v) > 0){ // while A[v] is nonempty:
		unsigned int u = G.AdjOut[v]._list.front(); // get first element in AdjList of v
		pdfs_rec(u, G); // dfs(A[v].head) //rst unvisited neighbour
	}

}

void ParallelGraphAlgorithms::parallel_dfs(unsigned int s, Graph G){

	make_directed(G);

	cout << "PDFS: ";
	pdfs_rec(s, G); // visit v
	cout << endl << endl;
}

void ParallelGraphAlgorithms::parallel_bfs(unsigned int s, Graph G){

	make_directed(G);

	queue<unsigned int> que; // Q := {}
	que.push(s); // Q.enqueue(s)
	isolate(s, G); // isolate(s)
	cout << "PBFS: ";
	while(!que.empty()){ // while Q is nonempty:
		unsigned int v = que.front();
		que.pop(); // v := Q.dequeue()
		cout << v << ' ';
		for(unsigned int u : G.out_neighbours(v)){ // for every neighbour u of v
			que.push(u); // Q.enqueue(u)
			isolate(u, G); // isolate(u)
		}
	}
	cout << endl << endl;

}

void ParallelGraphAlgorithms::simple_pbfs(unsigned int s, Graph G, ParallelGraphAlgorithms::bfs_t& T){

	// for all v in V in parallel do
	#pragma omp parallel for
	for(unsigned int i = 0; i < G.AdjOut.bucket_count(); i++){
	  for(auto it = G.AdjOut.begin(i); it != G.AdjOut.end(i); it++){ // get all vertices from graph
	  	T._p[it->first] = -1; // set parrent as -1
	   	T._d[it->first] = -1; // set distance to source as -1
	  }
	}

	T._p[s] = s;
	T._d[s] = 0;

	queue<unsigned int> CQ;
	CQ.push(s); // Enqueue r

	unsigned int sz; // size of queue
	while(!CQ.empty()){ // while CQ not empty

		queue<unsigned int> NQ; // NQ empty

		sz = CQ.size(); // used to iterate trough list
		#pragma omp parallel for
		for(unsigned int i = 0; i < sz; i++){

			unsigned int v;
			#pragma omp critical //Dequeue CQ. needs to be atomic
			{
				v = CQ.front();
				CQ.pop();
			}

			#pragma omp parallel for // for each u adjacent to v
			for(unsigned int j = 0; j < G.AdjOut[v]._map.bucket_count(); j++){
				for(auto it = G.AdjOut[v]._map.begin(j); it != G.AdjOut[v]._map.end(j); it++){ // get all neighbors of v
					unsigned int u = it->first;
					#pragma omp critical // if u has no parent. needs to be atomic
					if(T._p[u] == -1){
                        T._p[u] = v; // set v as parrent of u
                        T._d[u] = T._d[v]+1; // set distance to u
                        NQ.push(u); // Enqueue u
					}
				}
			}


		}

		CQ = NQ; // Swap(CQ, NQ);

	}

}


void ParallelGraphAlgorithms::D_1_distributed_memory_BFS(unsigned int s, Graph G, ParallelGraphAlgorithms::bfs_t& T){

    #pragma omp parallel for
	for(unsigned int i = 0; i < G.AdjOut.bucket_count(); i++){
        for(auto it = G.AdjOut.begin(i); it != G.AdjOut.end(i); it++){
            T._p[it->first] = -1;
            T._d[it->first] = -1;
        }
	}

	T._p[s] = s;
	T._d[s] = s;
	unsigned int level = 0;
	unsigned int nr_processors = 4;
	omp_set_num_threads(nr_processors);
	//vector<list<unsigned int>> FS(nr_processors);
	//vector<list<unsigned int>> NS(nr_processors);
	//list<unsigned int> NS_all;

    bool flag(1);
	while(flag) {

        vector<list<unsigned int>> FS(nr_processors);
        vector<list<unsigned int>> NS(nr_processors);
        list<unsigned int> NS_all;

        flag = false;

        unsigned int processor = 0;
        #pragma omp parallel for
        for(unsigned int i = 0; i < G.AdjOut.bucket_count(); i++) {
            for(auto it = G.AdjOut.begin(i); it != G.AdjOut.end(i); it++) {
                if(T._d[it->first] == level) {
                    #pragma omp critical
                    {
                        FS[processor].push_back(it->first);
                        processor++;
                    }

                }
            }
        }

        #pragma omp parallel for
        for(unsigned int i = 0; i < nr_processors; i ++) {
            if(!FS[i].empty()) {
                flag = true;
            }
        }

        /*
        #pragma omp parallel for
        for(unsigned int i = 0; i < nr_processors; i ++) {
            for(auto v : FS[i]) {
                #pragma omp parallel for
                for(unsigned int j = 0; j < G.AdjOut[v]._map.bucket_count(); j++){
                    for(auto it = G.AdjOut[v]._map.begin(j); it != G.AdjOut[v]._map.end(j); it++){
                        //#pragma omp critical
                        //{
                            NS_all.push_back(it->first);
                        //}
                    }
                }
            }
        }
        */

        #pragma omp parallel for
        for(unsigned int i = 0; i < nr_processors; i ++) {
            for(auto v : FS[i]) {
                #pragma omp parallel for
                for(unsigned int j = 0; j < G.AdjOut[v]._map.bucket_count(); j++){
                    for(auto it = G.AdjOut[v]._map.begin(j); it != G.AdjOut[v]._map.end(j); it++){
                        //#pragma omp critical
                        //{
                            NS[i].push_back(it->first);
                        //}
                    }
                }
            }
        }

        #pragma omp parallel for
        for(unsigned int i = 0; i < nr_processors; i ++) {
            for(auto v : NS[i]) {
                NS_all.push_back(v);
            }
        }

        for(auto v : NS_all) {
            if(T._d[v] == -1) {
                T._d[v] = level + 1;
            }
        }

        level++;

	}

}




void ParallelGraphAlgorithms::a_star_hash_distributed(unsigned int s, Graph G, ParallelGraphAlgorithms::bfs_t& T)
{

    // I chose an exact number of processors, to easier initialize the lists open and close, as every processor has his own list
    unsigned int nr_processors = 4;
	omp_set_num_threads(nr_processors);

	// To simulate that every processor has its own list, I created an map of ints and sets, the key of the map representing the processor
    unordered_map<unsigned int, unordered_set<unsigned int>> open(nr_processors);
    unordered_map<unsigned int, unordered_set<unsigned int>> close(nr_processors);

	unordered_map<unsigned int, double> g; // the distance from the source node to the candidate node

	// I initialised every node in grahp
    #pragma omp parallel for
	for(unsigned int i = 0; i < G.AdjOut.bucket_count(); i++){
        for(auto it = G.AdjOut.begin(i); it != G.AdjOut.end(i); it++){
            T._p[it->first] = -1; // Set every the parent of every node as -1, as there is no vertex -1
            T._d[it->first] = UINT_MAX; // Set the distance from source to node as max
            g[it->first] = DBL_MAX; // Set the distance from source to node as max
            // I decided to use both T._d and g for distance, because T._d is initialized with the number of vertices, and later, x = -1 is not a vertex.
            // The algorithm will return only T._p and T._d
        }
	}

	T._p[s] = s; // The parent of the source node would be itself
	T._d[s] = 0; // The distance from source node to itself is 0
	g[s] = 0; // The distance from source node to itself is 0

	open[rand() % nr_processors].insert(s);// Instead of the hash function, I tried to asignate vertices random

	bool flag = true;
	while(flag) {

        flag = false;

        //In original A*, the while loop was executed as long as the open list was not empty
        //Because I have more lists, I use bool varaible. If one list is not empty, the flag becomse true
        #pragma omp parallel for
        for(unsigned int i = 0; i < nr_processors; i ++) {
            if(!open[i].empty()) {
                //cout<<"open[" << i << "] not empty" << endl;
                flag = true;
            }
        }

        //Because every processor has its lists, I thought that every procesor should have its x
        unordered_map<unsigned int, unsigned int> x (nr_processors);

        #pragma omp parallel for
        for(unsigned int i = 0; i < nr_processors; i ++) {
            //Chose an x only for processors who have nodes in the open lists
            if(!open[i].empty()) {

                // First, for every processor I created a variable x set as -1 with max distance
                #pragma omp critical
                {
                    x[i] = -1; // First every x is equal to -1
                    g[x[i]] = DBL_MAX;// Set g for x as max
                    //cout<<"i = "<< i << endl;
                    //cout<<"g["<< x[i] << "] = " << g[x[i]] << endl;
                }
                // Next I compared x with the open nodes and I chose the one with the shortest distance
                for (auto v : open[i]) { // choose vertex with smallest g
                    //cout<<"v = "<< v << endl;
                    //cout<<"g["<< v << "] = " << g[v] << endl;
                    if(g[x[i]] >= g[v]) {
                        #pragma omp critical
                        {
                            x[i] = v;
                        }
                    }
                }
                // Then, as the definition of original A* sais, I moved x from the open list in the close list
                #pragma omp critical
                {
                    //cout<<"i = "<< i << endl;
                    //cout<<"x["<< i << "] = " << x[i] << endl;
                    //cout<<"g["<< x[i] << "] = " << g[x[i]] << endl;
                    open[i].erase(x[i]);
                    close[i].insert(x[i]);
                }

            }

        }
        //cout << endl;

        #pragma omp parallel for
        for(unsigned int i = 0; i < nr_processors; i ++) {
            if(!close[i].empty()) {
                //cout<<"i = "<< i << endl;

                //For every x, the program chooses an y from its neighbours
                #pragma omp parallel for
                for(unsigned int j = 0; j < G.AdjOut[x[i]]._map.bucket_count(); j++) {
                    for(auto it = G.AdjOut[x[i]]._map.begin(j); it != G.AdjOut[x[i]]._map.end(j); it++){

                        unsigned int y = it->first;
                        //cout<<"y = "<< y << endl;
                        //cout << "g[" << y << "] = " << g[y] << endl;

                        //Source node should always be closed
                        if(y == s) {
                            #pragma omp critical
                            {
                                close[i].insert(y);
                            }
                        } else {

                            // If y is opened, but there is a better distance for it, its distance and parent will be updated
                            if(open[i].find(y) != open[i].end() && g[y] > g[x[i]] + G.weight(x[i], y)) {
                                #pragma omp critical // I used a barier because g and T are shared variables/structues
                                {
                                    g[y] = g[x[i]] + G.weight(x[i], y);
                                    T._d[y] = g[x[i]] + G.weight(x[i], y);
                                    T._p[y] = x[i];
                                    //cout << "g[" << y << "] = " << g[y] << endl;
                                    //cout << "d[" << y << "] = " << T._d[y] << endl;
                                    //cout << "p[" << y << "] = " << T._p[y] << endl;
                                }
                            // If y is closed, but there is a better distance for it, its distance and parent will be updated, and it will reopen
                            } else if (close[i].find(y) != close[i].end() && g[y] > g[x[i]] + G.weight(x[i], y)) {
                                #pragma omp critical
                                {
                                    g[y] = g[x[i]] + G.weight(x[i], y);
                                    T._d[y] = g[x[i]] + G.weight(x[i], y);
                                    T._p[y] = x[i];
                                    close[i].erase(y);
                                    open[rand() % nr_processors].insert(y);
                                    //cout << "g[" << y << "] = " << g[y] << endl;
                                    //cout << "d[" << y << "] = " << T._d[y] << endl;
                                    //cout << "p[" << y << "] = " << T._p[y] << endl;
                                }
                            // If y is neither opened or closed, its distance will be computed and its parent set
                            } else if (open[i].find(y) == open[i].end() && close[i].find(y) == close[i].end()) {
                                #pragma omp critical
                                {
                                    g[y] = g[x[i]] + G.weight(x[i], y);
                                    T._d[y] = g[x[i]] + G.weight(x[i], y);
                                    T._p[y] = x[i];
                                    open[rand() % nr_processors].insert(y);
                                    //cout << "g[" << y << "] = " << g[y] << endl;
                                    //cout << "d[" << y << "] = " << T._d[y] << endl;
                                    //cout << "p[" << y << "] = " << T._p[y] << endl;
                                }
                            }
                        }
                    }
                }
            }
        }




	}
	//cout << endl;

}







/*------------------------------------------- Connected components ----------------------------------------------*/

bool ParallelGraphAlgorithms::my_direct_connect(Graph& G, vector<unsigned int>& p, vector<unsigned int>& o){

	bool changes(0);

	// for each edge {v, w} do
	#pragma omp parallel for // first get vertex v
	for(unsigned int i = 0; i < G.AdjOut.bucket_count(); i++){ // get bucket from undorderd_map one by one
		for(auto itv = G.AdjOut.begin(i); itv != G.AdjOut.end(i); itv++){ // get vertex from bucket
			unsigned int v = itv->first; // iterator is a pointer
			unsigned int parent = p[v]; // used later to check for changes

			#pragma omp parallel for // second get neighbors of v vertices w to find edges {v, w}
			for(unsigned int j = 0; j < G.AdjOut[v]._map.bucket_count(); j++){
				for(auto itw = G.AdjOut[v]._map.begin(j); itw != G.AdjOut[v]._map.end(j); itw++){

					unsigned int w = itw->first;

					if(v > w) { // if v > w then
                        #pragma omp critical
                        {
                            p[v] = min(p[v], w); // v.p = min {v.p, w}
                        }
                    } else { // else
                        #pragma omp critical
                        {
                            p[w] = min(p[w], v); // w.p = min {w.p, v}
                        }
                    }

				}
			}

			if(parent != p[v]) {
                changes = 1; // verify if changes for S, A, RA algorithms
            }

		}
	}

	return changes;

}

bool ParallelGraphAlgorithms::my_parent_connect(Graph& G, vector<unsigned int>& p, vector<unsigned int>& o){

	bool changes(0);

	#pragma omp parallel for // for each vertex v do
	for(unsigned int i = 0; i < G.vertices.bucket_count(); i++){ // get buckets from unorderd_map
		for(auto itv = G.vertices.begin(i); itv != G.vertices.end(i); itv++){ // get vertices from buckets
			unsigned int v = *itv; // iterator is a pointer
			o[v] = p[v]; // v.o = v.p
		}
	}

	#pragma omp parallel for // for each edge {v, w} do
	for(unsigned int i = 0; i < G.AdjOut.bucket_count(); i++){
		for(auto itv = G.AdjOut.begin(i); itv != G.AdjOut.end(i); itv++){

			unsigned int v = itv->first;
			unsigned int parent = p[v];

			#pragma omp parallel for
			for(unsigned int j = 0; j < G.AdjOut[v]._map.bucket_count(); j++){
				for(auto itw = G.AdjOut[v]._map.begin(j); itw != G.AdjOut[v]._map.end(j); itw++){

					unsigned int w = itw->first;

					if(o[v] > o[w]) {// if v.o > w.o then
                        #pragma omp critical
                        {
                            p[o[v]] = min(p[o[v]], o[w]); // v.o.p = min{v.o.p, w.o}
                        }
                    } else {
                        #pragma omp critical
                        {
                            p[o[w]] = min(p[o[w]], o[v]); // w.o.p = min{w.o.p, v.o}
                        }
                    }

				}
			}

            if(parent != p[v]) {
                changes = 1;
            }

		}
	}

	return changes;

}

bool ParallelGraphAlgorithms::my_direct_root_connect(Graph& G, vector<unsigned int>& p, vector<unsigned int>& o){

	bool changes(0);

	#pragma omp parallel for // for each vertex v do
	for(unsigned int i = 0; i < G.vertices.bucket_count(); i++){
		for(auto itv = G.vertices.begin(i); itv != G.vertices.end(i); itv++){
			unsigned int v = *itv;
			o[v] = p[v]; // v.o = v.p
		}
	}

	#pragma omp parallel for // for each vertex v do
	for(unsigned int i = 0; i < G.AdjOut.bucket_count(); i++){
		for(auto itv = G.AdjOut.begin(i); itv != G.AdjOut.end(i); itv++){

			unsigned int v = itv->first;
            unsigned int parent = p[v];

            #pragma omp parallel for reduction(min : parent)
            for(unsigned int j = 0; j < G.AdjOut[v]._map.bucket_count(); j++){
                for(auto itw = G.AdjOut[v]._map.begin(j); itw != G.AdjOut[v]._map.end(j); itw++){

                    unsigned int w = itw->first;

                    if(v > w && v == o[v]) {// if v > w && v = v.o then
                        #pragma omp critical
                        {
                            p[v] = min(p[v], w); // v.p = min{v.p, w}
                        }
                    } else if (w == o[w]) { // else if w = w.o then
                        #pragma omp critical
                        {
                            p[w]= min(p[w], v); // w.p = min{w.p, v}
                        }
                    }
                }
            }

            if(parent != p[v]) {
                changes = 1;
            }

		}
	}

	return changes;

}

bool ParallelGraphAlgorithms::my_parent_root_connect(Graph& G, vector<unsigned int>& p, vector<unsigned int>& o){

	bool changes(0);

	#pragma omp parallel for // for each vertex v do
	for(unsigned int i = 0; i < G.vertices.bucket_count(); i++){
		for(auto itv = G.vertices.begin(i); itv != G.vertices.end(i); itv++){
			unsigned int v = *itv;
			o[v] = p[v]; // v.o = v.p
		}
	}

	#pragma omp parallel for // for each edge {v, w} do
	for(unsigned int i = 0; i < G.AdjOut.bucket_count(); i++){
		for(auto itv = G.AdjOut.begin(i); itv != G.AdjOut.end(i); itv++){

			unsigned int v = itv->first;
			unsigned int parent = p[v];

			#pragma omp parallel for
			for(unsigned int j = 0; j < G.AdjOut[v]._map.bucket_count(); j++){
				for(auto itw = G.AdjOut[v]._map.begin(j); itw != G.AdjOut[v]._map.end(j); itw++){

					unsigned int w = itw->first;

					if(o[v] > o[w] && o[o[v]] == o[v]) { // if v.o > w.o and v.o = v.o.o then
                        #pragma omp critical
                        {
                            p[o[v]] = min(p[o[v]], o[w]); // v.o.p = min {v.o.p, w.o}
                        }
                    } else if (o[w] == o[o[w]]) { // else if w.o == w.o.o then
                        #pragma omp critical
                        {
                            p[o[w]] = min(p[o[w]], o[v]); // w.o.p = min {w.o.p, v.o}
                        }
                    }

				}
			}

            if(parent != p[v]) {
                changes = 1;
            }

		}
	}

	return changes;

}




bool ParallelGraphAlgorithms::shortcut(Graph& G, vector<unsigned int>& p, vector<unsigned int>& o){

	bool changes(0);

	#pragma omp parallel for // for each vertex v do
	for(unsigned int i = 0; i < G.vertices.bucket_count(); i++){
		for(auto it = G.vertices.begin(i); it != G.vertices.end(i); it++){
			unsigned int v = *it;
			o[v] = p[v]; // v.o = v.p
		}
	}

	#pragma omp parallel for // for each vertex v do
	for(unsigned int i = 0; i < G.vertices.bucket_count(); i++){
		for(auto it = G.vertices.begin(i); it != G.vertices.end(i); it++){
			unsigned int v = *it;
			p[v] = o[o[v]]; // v.p = v.o.o
			if(p[v] != o[v]) {
                changes = 1;
            }
		}
	}

	return changes;

}

void ParallelGraphAlgorithms::alter(Graph& G, vector<unsigned int>& p){

	#pragma omp parallel for // for each edge {v, w} do
	for(unsigned int i = 0; i < G.AdjOut.bucket_count(); i++){
		for(auto itv = G.AdjOut.begin(i); itv != G.AdjOut.end(i); itv++){

			unsigned int v = itv->first;

			#pragma omp parallel for
			for(unsigned int j = 0; j < G.AdjOut[v]._map.bucket_count(); j++){
				for(auto itw = G.AdjOut[v]._map.begin(j); itw != G.AdjOut[v]._map.end(j); itw++){

					unsigned int w = itw->first;

					if(p[v] != p[w]){ // if v.p = w.p then

					 Graph::Edge e; //dummy

					 #pragma omp critical
					 {
                        G.AdjIn[p[w]]._map[p[v]] = e; // add {w.p, v.p} // G undirected => AdjIn empty
					 }

					}
				}
			}

		}
	}

	G.AdjOut = unordered_map<unsigned int, Graph::AdjList>(); //empty AdjOut
	swap(G.AdjOut, G.AdjIn); // replace {v, w} by {v.p, w.p}

}




void ParallelGraphAlgorithms::my_algorithm_S(Graph G, vector<unsigned int>& p, vector<unsigned int>& o){

	#pragma omp parallel for // LT-framework Initialization: v.p = v for every vertex v
	for(unsigned int i = 0; i < G.AdjOut.bucket_count(); i++){
		for(auto it = G.AdjOut.begin(i); it != G.AdjOut.end(i); it++){
			unsigned int v = it->first;
			p[v] = v;
		}
	}

	bool not_done(1); // no v.p changes
	while(not_done){ // repeat
		not_done = my_parent_connect(G, p, o); // parent-connect
		bool not_star(1); // no v.p changes
		while(not_star) { // repeat until no v.p changes
            not_star = shortcut(G, p, o);
        }
	} // until no v.p changes

}

void ParallelGraphAlgorithms::my_algorithm_RA(Graph G, vector<unsigned int>& p, vector<unsigned int>& o){

	#pragma omp parallel for // LT-framework Initialization: v.p = v for every vertex v
	for(unsigned int i = 0; i < G.AdjOut.bucket_count(); i++){
		for(auto it = G.AdjOut.begin(i); it != G.AdjOut.end(i); it++){
			unsigned int v = it->first;
			p[v] = v;
		}
	}

	bool not_done(1); // check if no v.p changes
	while(not_done){ // repeat
		bool a = my_direct_root_connect(G, p, o); // direct-root-connect
		bool b = shortcut(G,p,o); // shortcut
		not_done = a || b; //
        alter(G,p); // alter
	} // until no v.p changes

}

void ParallelGraphAlgorithms::my_algorithm_A(Graph G, vector<unsigned int>& p, vector<unsigned int>& o){

	#pragma omp parallel for // LT-framework Initialization: v.p = v for every vertex v
	for(unsigned int i = 0; i < G.AdjOut.bucket_count(); i++){
		for(auto it = G.AdjOut.begin(i); it != G.AdjOut.end(i); it++){
			unsigned int v = it->first;
			p[v] = v;
		}
	}

	bool not_done(1); // check if no v.p changes
	while(not_done){ // repeat
		bool a = my_direct_connect(G,p,o); // direct-connect
		bool b = shortcut(G,p,o); // shortcut
		not_done = a || b;
		alter(G,p); // alter
	} // until no v.p changes

}

unordered_map<unsigned int, Graph> ParallelGraphAlgorithms::my_algorithm_R(Graph G, vector<unsigned int>& p, vector<unsigned int>& o){

	#pragma omp parallel for // LT-framework Initialization: v.p = v for every vertex v
	for(unsigned int i = 0; i < G.AdjOut.bucket_count(); i++){
		for(auto it = G.AdjOut.begin(i); it != G.AdjOut.end(i); it++){
			unsigned int v = it->first;
			p[v] = v;
		}
	}

	bool not_done(1); // check if no v.p changes
	while(not_done){ // repeat
		bool a = my_parent_root_connect(G,p,o); // parent-root-connect
		bool b = shortcut(G,p,o); // shortcut
		not_done = a || b;
	} // until no v.p changes

	unordered_set<unsigned int> visited;
    unordered_map<unsigned int, Graph> spanning_trees;

/*
	#pragma omp parallel for
	for(unsigned int i = 0; i < G.vertices.bucket_count(); i++){
		for(auto itv = G.vertices.begin(i); itv != G.vertices.end(i); itv++){

			unsigned int v = *itv;
*/
    #pragma omp parallel for
    for(unsigned int i = 0; i < G.AdjOut.bucket_count(); i++){
        for(auto itv = G.AdjOut.begin(i); itv != G.AdjOut.end(i); itv++){

            unsigned int v = itv->first;

			if(v == p[v]) {
                //#pragma omp critical
                {
                    Graph spanning_tree = Graph(false, false);
                    spanning_trees[v] = spanning_tree;
                    spanning_trees[v].add_node(v);
                    visited.insert(v);
                }
			}

        }
    }

    //#pragma omp parallel for
    for(unsigned int i = 0; i < G.AdjOut.bucket_count(); i++){
        for(auto itv = G.AdjOut.begin(i); itv != G.AdjOut.end(i); itv++){

            unsigned int v = itv->first;

            #pragma omp parallel for
            for(unsigned int j = 0; j < G.AdjOut[v]._map.bucket_count(); j++){
                for(auto itw = G.AdjOut[v]._map.begin(j); itw != G.AdjOut[v]._map.end(j); itw++){

                    unsigned int w = itw->first;

                    if(v < w) {
                        #pragma omp critical
                        {
                            if(visited.find(w) == visited.end()){
                                //cout << "v = " << v << endl;
                                //cout << "w = " << w << endl;
                                //cout << "p[" << v << "] = " << p[v] << endl;
                                //cout << v << " " << w << endl;
                                spanning_trees[p[v]].add_node(w);
                                spanning_trees[p[v]].add_edge(v, w);
                                visited.insert(w);
                            }
                        }
                    }

				}
            }

		}
    }

    return spanning_trees;
}
