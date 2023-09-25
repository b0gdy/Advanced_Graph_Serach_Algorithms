/*
Badea Bogdan-Andrei
Gr. 506
Subiectul 4
*/

#include "Graph.hpp"
#include "Search.hpp"

#include <iostream>

using namespace std;

int main()
{

    string filename;
    cout << "Please enter the name of the to be read: ";
    /*
    Files examples to be entered in the terminal:
    Weighted graph:
    C:\Users\Bogdan\Desktop\Tehnici_avansate_de_programare\Laborator\All_algorithms\Graph.txt
    Unweighted graph:
    C:\Users\Bogdan\Desktop\Tehnici_avansate_de_programare\Laborator\All_algorithms\Graph2.txt
    Unweighted graph for RA:
    C:\Users\Bogdan\Desktop\Tehnici_avansate_de_programare\Laborator\All_algorithms\Graph3.txt
    */
    cin >> filename;
    cout << endl;
    Graph G = G.read_file(filename);
    cout << endl;

    if(G.is_weighted()) {
        G.print_weighted_edges();
    } else {
        G.print_edges();
    }
    cout << endl;


    unsigned int n = G.nr_vertices();
    ParallelGraphAlgorithms::bfs_t T = ParallelGraphAlgorithms::bfs_t(n);
    search_output output;
    vector<unsigned int> p(G.nr_vertices());
    vector<unsigned int> o(G.nr_vertices());

    // In this part we'd have to call the desired algorithm. Uncomment the desired algorithm.



    // Search sequential algorithms
    // Please use file Graph.txt
    //output = my_dijkstra_search(0, G, output);
    //output = my_a_star_search(0, 6, G, output);
    //output = bi_a_star_search(0, 6, G, output);
    //output = beam_a_star_search(0, 6, G, output);

    // Print for weighted graphs
    // The output will be for every node its parent and distance from source
    /*
    for(auto i : G.vertex_list()){
            cout << "d[" << i << "] = " << output._d[i] << endl;
            cout << "p[" << i << "] = " << output._p[i] << endl;
            cout << endl;
        }
    */



    // Parallel algorithms
    // Please use file Graph.txt
    //ParallelGraphAlgorithms::simple_pbfs(0, G, T);
    //ParallelGraphAlgorithms::D_1_distributed_memory_BFS(0, G, T);
    //ParallelGraphAlgorithms::a_star_hash_distributed(0, G, T);
    // Print for paralel graphs
    // The output will be for every node its parent and distance from source
    /*
    for(auto i : G.vertex_list()) {
        cout << "d[" << i << "] = " << T._d[i] << endl;
        cout << "p[" << i << "] = " << T._p[i] << endl;
        cout << endl;
    }
    */




    // Connected Components algorithms
    // Please use file Graph3.txt
    //ParallelGraphAlgorithms::my_algorithm_S(G, p, o);
    //ParallelGraphAlgorithms::my_algorithm_RA(G, p, o);
    //ParallelGraphAlgorithms::my_algorithm_A(G, p, o);
    // Print for paralel graphs
    // The output will be for every node its parent
    /*
    for(auto v: G.vertex_list()) {
        cout << "p[" << v << "] = " << p[v] << endl;
    }
    cout << endl;
    */


    // The last connected-components algorithm was treated separately
    /*
    unordered_map<unsigned int, Graph> spanning_trees;
    spanning_trees = ParallelGraphAlgorithms::my_algorithm_R(G, p, o);
    for(auto v: G.vertex_list()) {
        cout << "p[" << v << "] = " << p[v] << endl;
    }
    cout << endl;
    for(auto g : spanning_trees) {
        g.second.print_edges();
        cout << endl;
    }
    */


    return 0;

}
