// (concept)
// https://www.geeksforgeeks.org/dsa/steiner-tree/ 
// (heuristic algorithm on page 10 and 11)
// https://www.comp.nus.edu.sg/~stevenha/cs4234/lectures/03b.SteinerTree.pdf 

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

// structure to represent an edge
typedef struct {
    int src;
    int dest;
    int weight;
} Edge;

// structure to represent the graph
typedef struct {
    int V;              // total number of vertices
    int E;              // total number of edges
    Edge* edges;        // array of all edges in the graph
    
    // array to flag if a vertex is a mandatory terminal.
    // index is the vertex ID (0 to V-1).
    // true = Terminal, false = Steiner Point.
    bool* is_terminal;  

    // array used during the brute-force step to track which optional 
    // steiner points are currently "turned on" for the current test.
    bool* is_active;    
} Graph;

// structure for the Union-Find cycle detection
typedef struct {
    int parent;
    int rank;
} Subset;

// function to create the graph in memory
Graph* createGraph(int V, int E) {
    Graph* graph = (Graph*)malloc(sizeof(Graph));
    graph->V = V;
    graph->E = E;
    graph->edges = (Edge*)malloc(E * sizeof(Edge));
    
    graph->is_terminal = (bool*)calloc(V, sizeof(bool));
    graph->is_active = (bool*)calloc(V, sizeof(bool)); 
    
    return graph;
}