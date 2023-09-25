#ifndef SEARCH_HPP_INCLUDED
#define SEARCH_HPP_INCLUDED

#include "Graph.hpp"

#include <unordered_map>
#include <unordered_set>
#include <functional>

list<unsigned int> bfs(unsigned int,const Graph&);

list<unsigned int> my_dfs(unsigned int,const Graph&);

// ------------------------------------------------------------- Weighted Search -------------------------------------------------------------

struct search_output {

	unordered_map<unsigned int, double> _d; //distances from the source
	unordered_map<unsigned int, unsigned int> _p; //pairs (child, parent)

	//search_output(unsigned n): _d(n), _p(n)	{}

};

search_output my_dijkstra_search(unsigned int, const Graph&, search_output&);
search_output my_a_star_search(unsigned int, unsigned int, const Graph&, search_output&);


search_output bi_a_star_search(unsigned int, unsigned int, const Graph&, search_output&);


search_output beam_a_star_search(unsigned int, unsigned int, const Graph&, search_output&, unsigned int);



#endif // SEARCH_HPP_INCLUDED
