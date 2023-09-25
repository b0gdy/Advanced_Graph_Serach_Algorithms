#include "Search.hpp"

#include <queue>
#include <stack>
#include <unordered_set>

#include <iostream>
using namespace std;

/*-------------------------------------------- Graph Search -----------------------------------------*/

list<unsigned int> my_dfs(unsigned int s, const Graph& G){

    list<unsigned int> output;
    unordered_map<unsigned int, bool> visited;
    for(auto i : G.vertex_list()) {
        visited.insert(make_pair(i, false));
    }

    stack<int> stk;
    stk.push(s);

    while(!stk.empty()) {
        int node = stk.top();
        stk.pop();
        if (visited.at(node) == false) {
            output.push_back(node);
            visited.at(node) = true;
        }
        unordered_set<int>::iterator itr;
        for(auto i : G.out_neighbours(node)) {
            if (visited.at(i) == false) {
                stk.push(i);
            }
        }
    }

	return output;

}

list<unsigned int> bfs(unsigned int s, const Graph& G){
	list<unsigned int> output;
	unordered_set<unsigned int> visited;
	queue<unsigned int> que;

	visited.insert(s);
	que.push(s);

	while(!que.empty()){

		unsigned int v = que.front();
		que.pop();
		output.push_back(v);

		for(unsigned int u : G.out_neighbours(v))
		 if(visited.find(u) == visited.end()){
		 	visited.insert(u);
		 	que.push(u);
		 }
	}

	return output;
}



search_output my_dijkstra_search(unsigned int s, const Graph& G, search_output& res){

    for(auto v : G.vertex_list()){
        res._d[v] = UINT_MAX;
        res._p[v] = -1;
    }

    unordered_set<unsigned int> visited;

    res._d[s] = 0;
    res._p[s] = 0;

    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> heap;
    heap.push(make_pair(res._d[s], s));

    while(!heap.empty()) {

        unsigned int v = heap.top().second;
        heap.pop();
        visited.insert(v);

        for(auto u : G.out_neighbours(v)) {
            double _w = G.weight(v, u);

            if(res._d[u] > res._d[v] + _w) {
                res._d[u] = res._d[v] + _w;
                res._p[u] = v;
                if(visited.find(u) == visited.end()) {
                    visited.insert(u);
                    heap.push(make_pair(res._d[u], u));
                }
            }

        }

    }

    return res;

}

search_output my_a_star_search(unsigned int s, unsigned int p, const Graph& G, search_output& res){

    unordered_set<unsigned int> open;
    open.insert(s);
    unordered_set<unsigned int> close;

    unordered_map<unsigned int, double> h; // estimate distance to destination

    for(auto v : G.vertex_list()){
        res._d[v] = UINT_MAX; // in course: g = best known distance to source node
        res._p[v] = -1;
        h[v] = 0;
    }

    res._d[s] = 0;
    res._p[s] = s;

    //unordered_map<unsigned int, double> f; // f = h + g;
    //f[s] = res._d[s] + h[s];

    while(open.size() != 0) {

        unsigned int x = -1;
        res._d[x] = UINT_MAX;

        /*
        f[x] = numeric_limits<unsigned int>::max();
        for (auto i : open) {
            f[i] = res._d[i] + h[i];
        }
        for (auto i : open) { // choose vertex with smallest f
            if(f[x] >= f[i]) {
                x = i;
            }
        }
        */

        for (auto i : open) { // choose vertex with smallest g
            if(res._d[x] >= res._d[i]) {
                x = i;
            }
        }


        open.erase(x);
        close.insert(x);

        for(auto y : G.out_neighbours(x)) {
            if(open.find(y) != open.end() && res._d[y] > res._d[x] + G.weight(x, y)) {
                res._d[y] = res._d[x] + G.weight(x, y);
                res._p[y] = x;
            } else if (close.find(y) != close.end() && res._d[y] > res._d[x] + G.weight(x, y)) {
                res._d[y] = res._d[x] + G.weight(x, y);
                res._p[y] = x;
                close.erase(y);
                open.insert(y);
            } else if (open.find(y) == open.end() && close.find(y) == close.end()) {
                res._d[y] = res._d[x] + G.weight(x, y);
                res._p[y] = x;
                open.insert(y);
            }
        }
    }

    return res;

}



// Hash function
struct hashFunction
{
    size_t operator()(const pair<unsigned int ,unsigned int> &x) const
    {
        return x.first ^ x.second;
    }
};

search_output bi_a_star_search(unsigned int s, unsigned int t, const Graph& G, search_output& res){

    unordered_set<pair<unsigned int, unsigned int>, hashFunction> open;
    open.insert(make_pair(s, t));
    unordered_set<pair<unsigned int, unsigned int>, hashFunction> close;

    unordered_map<pair<unsigned int, unsigned int>, double, hashFunction> g; // distance from source to first node and from target to second node
    unordered_map<pair<unsigned int, unsigned int>, double, hashFunction> h; // estimated distance

    for(auto v : G.vertex_list()){
        res._d[v] = UINT_MAX; // distance from source
        res._p[v] = -1; // parent
        for(auto u : G.vertex_list()) {
            g[make_pair(v, u)] = UINT_MAX;
            h[make_pair(v, u)] = 0;
        }
    }

    res._d[s] = 0;
    res._p[s] = s;

    pair<unsigned int, unsigned int> p0;
    p0 = make_pair(s, t);
    g[p0] = 0;

    while(open.size() != 0) {

        cout << "Open = { ";
        for(auto i : open) {
            cout << "(" << i.first << ", " << i.second << ") ";
        }
        cout << "}" << endl;

        pair<unsigned int, unsigned int> x = make_pair(-1, -1);

        g[x] = numeric_limits<unsigned int>::max();
        //cout << "g[(" << x.first << ", " << x.second << ")] = " << g[x] << endl;
        for (auto i : open) {
            //cout << "g[(" << i.first << ", " << i.second << ")] = " << g[i] << endl;
            if(g[x] >= g[i]) {
                x = i;
            }
        }

        open.erase(x);
        close.insert(x);

        cout << "pair = (" << x.first << ", " << x.second << ")" << endl;
        cout << "g[(" << x.first << ", " << x.second << ")] = " << g[x] << endl;

        for(auto y_s : G.out_neighbours(x.first)) {
            for(auto y_t : G.out_neighbours(x.second)) {
                if(y_s > x.first && y_t < x.second) {

                pair<unsigned int, unsigned int> y = make_pair(y_s, y_t);
                cout << "(" << y.first << ", " << y.second << ")" << endl;
                cout << "g[(" << y.first << ", " << y.second << ")] = " << g[y] << endl;

                if(open.find(y) != open.end() && g[y] > g[x] + G.weight(x.first, y.first) + G.weight(x.second, y.second)) {
                    g[y] = g[x] + G.weight(x.first, y.first) + G.weight(x.second, y.second);
                    res._d[y.first] = res._d[x.first] + G.weight(x.first, y.first);
                    res._p[y.first] = x.first;
                    cout << "g[(" << x.first << ", " << x.second << ")] + G.weight(" << x.first << ", " <<  y.first << ") + G.weight(" << x.second << ", " <<  y.second << ") = " << g[x] + G.weight(x.first, y.first) + G.weight(x.second, y.second) << endl;
                    cout << "p[" << y.first << "] = " << x.first << endl;
                    cout << "d[" << y.first << "] = " << res._d[y.first] << endl;
                } else if (close.find(y) != close.end() && g[y] > g[x] + G.weight(x.first, y.first) + G.weight(x.second, y.second)) {
                    g[y] = g[x] + G.weight(x.first, y.first) + G.weight(x.second, y.second);
                    res._d[y.first] = res._d[x.first] + G.weight(x.first, y.first);
                    res._p[y.first] = x.first;;
                    close.erase(y);
                    open.insert(y);
                    cout << "g[(" << x.first << ", " << x.second << ")] + G.weight(" << x.first << ", " <<  y.first << ") + G.weight(" << x.second << ", " <<  y.second << ") = " << g[x] + G.weight(x.first, y.first) + G.weight(x.second, y.second) << endl;
                    cout << "p[" << y.first << "] = " << x.first << endl;
                    cout << "d[" << y.first << "] = " << res._d[y.first] << endl;
                } else if (open.find(y) == open.end() && close.find(y) == close.end()) {
                    g[y] = g[x] + G.weight(x.first, y.first) + G.weight(x.second, y.second);
                    res._d[y.first] = res._d[x.first] + G.weight(x.first, y.first);
                    res._p[y.first] = x.first;
                    open.insert(y);
                    cout << "g[(" << x.first << ", " << x.second << ")] + G.weight(" << x.first << ", " <<  y.first << ") + G.weight(" << x.second << ", " <<  y.second << ") = " << g[x] + G.weight(x.first, y.first) + G.weight(x.second, y.second) << endl;
                    cout << "p[" << y.first << "] = " << x.first << endl;
                    cout << "d[" << y.first << "] = " << res._d[y.first] << endl;
                }
            }
            }
            res._p[s] = s;
            res._d[s] = 0;
        }

        cout<<endl;

    }
    res._p[s] = s;
    res._d[s] = 0;
    /*
    unsigned int p_aux = 0;
    for(auto v : G.vertex_list()) {
        if(p_aux <)
    }
    for(auto v : G.vertex_list()) {
            if (v != s) {
                cout << "v = " << v << endl;
                cout << "p[" << v << "] = " << res._p[v] << endl;
                cout << "d[" << res._p[v] << "] = " << res._d[res._p[v]] << endl;
                cout << "weight(" << res._p[v] << ", " << v << ") = " << G.weight(res._p[v], v) << endl;
                res._d[v] = res._d[res._p[v]] + G.weight(res._p[v], v);
                cout << "d[" << v << "] = " << res._d[v] << endl;
            }
            cout << endl;
        }
    */

    return res;

}




search_output beam_a_star_search(unsigned int s, unsigned int p, const Graph& G, search_output& res, unsigned int beta){

    unordered_set<unsigned int> open;
    open.insert(s);
    unordered_set<unsigned int> close;

    unordered_map<unsigned int, double> h; // estimate distance to destination

    for(auto v : G.vertex_list()){
        res._d[v] = UINT_MAX; // in course: g = best known distance to source node
        res._p[v] = -1;
        h[v] = 0;
    }

    res._d[s] = 0;
    res._p[s] = s;

    //unordered_map<unsigned int, double> f; // f = h + g;
    //f[s] = res._d[s] + h[s];

    while(open.size() != 0) {

        unsigned int x;


        while (open.size() > beta) {

            x = -1;
            res._d[x] = 0;

            for (auto i : open) {
                if(res._d[x] <= res._d[i]) {
                    x = i;
                }
            }

            open.erase(x);

        }

        x = -1;
        res._d[x] = UINT_MAX;

        /*
        f[x] = numeric_limits<unsigned int>::max();
        for (auto i : open) {
            f[i] = res._d[i] + h[i];
        }
        for (auto i : open) { // choose vertex with smallest f
            if(f[x] >= f[i]) {
                x = i;
            }
        }
        */

        for (auto i : open) { // choose vertex with smallest g
            if(res._d[x] >= res._d[i]) {
                x = i;
            }
        }


        open.erase(x);
        close.insert(x);

        for(auto y : G.out_neighbours(x)) {
            if(open.find(y) != open.end() && res._d[y] > res._d[x] + G.weight(x, y)) {
                res._d[y] = res._d[x] + G.weight(x, y);
                res._p[y] = x;
            } else if (close.find(y) != close.end() && res._d[y] > res._d[x] + G.weight(x, y)) {
                res._d[y] = res._d[x] + G.weight(x, y);
                res._p[y] = x;
                close.erase(y);
                open.insert(y);
            } else if (open.find(y) == open.end() && close.find(y) == close.end()) {
                res._d[y] = res._d[x] + G.weight(x, y);
                res._p[y] = x;
                open.insert(y);
            }
        }
    }

    return res;

}
