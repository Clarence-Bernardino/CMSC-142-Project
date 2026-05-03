// (concept)
// https://www.geeksforgeeks.org/dsa/steiner-tree/ 
// (heuristic algorithm on page 10 and 11)
// https://www.comp.nus.edu.sg/~stevenha/cs4234/lectures/03b.SteinerTree.pdf 

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <time.h>

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
    
    // array to flag if a vertex is mandatory.
    // index is the vertex ID (0 to V-1).
    // true = Terminal, false = Steiner Point.
    bool* is_required;  

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
    
    graph->is_required = (bool*)calloc(V, sizeof(bool));
    graph->is_active = (bool*)calloc(V, sizeof(bool)); 
    
    return graph;
}

// function to free the graph and all associated memory
void freeGraph(Graph* graph) {
    if (graph == NULL) {
        return;
    }
    
    // free the edges array
    if (graph->edges != NULL) {
        free(graph->edges);
    }
    
    // free the is_required array
    if (graph->is_required != NULL) {
        free(graph->is_required);
    }
    
    // free the is_active array
    if (graph->is_active != NULL) {
        free(graph->is_active);
    }
    
    // finally free the graph structure itself
    free(graph);
}

// ===== Brute Force Functions Start =====

// union-find: Find the root of the set in which element i belongs
int find(Subset subsets[], int i) {
    if (subsets[i].parent != i) {
        subsets[i].parent = find(subsets, subsets[i].parent); // path compression
    }
    return subsets[i].parent;
}

// union-find: Unite two sets
void Union(Subset subsets[], int x, int y) {
    int xroot = find(subsets, x);
    int yroot = find(subsets, y);

    if (subsets[xroot].rank < subsets[yroot].rank) {
        subsets[xroot].parent = yroot;  // make yroot the parent of xroot
    } else if (subsets[xroot].rank > subsets[yroot].rank) {
        subsets[yroot].parent = xroot;  // make xroot the parent of yroot
    } else {
        subsets[yroot].parent = xroot;  // choose xroot as parent and increment its rank by 1
        subsets[xroot].rank++;
    }
}

// comparison function for qsort to sort edges by weight ascending
int compareEdges(const void* a, const void* b) {
    Edge* edgeA = (Edge*)a;
    Edge* edgeB = (Edge*)b;
    return edgeA->weight - edgeB->weight;
}

// calculates the MST for ONLY the vertices marked as "active"
int activeSubsetMST(Graph* graph) {
    int V = graph->V;
    int edges_added = 0;
    int current_edge_idx = 0;
    int total_weight = 0;

    // count how many vertices are currently active
    int active_nodes = 0;
    for (int v = 0; v < V; v++) {
        if (graph->is_active[v]) {
            active_nodes++;
        }
    }

    // allocate memory for union-find subsets
    Subset* subsets = (Subset*)malloc(V * sizeof(Subset));
    for (int v = 0; v < V; ++v) {
        subsets[v].parent = v;
        subsets[v].rank = 0;
    }

    // process edges in ascending order
    while (edges_added < active_nodes - 1 && current_edge_idx < graph->E) {
        Edge next_edge = graph->edges[current_edge_idx++];

        // skip edge if either endpoint is turned off in this subset
        if (!graph->is_active[next_edge.src] || !graph->is_active[next_edge.dest]) {
            continue;
        }

        int x = find(subsets, next_edge.src);
        int y = find(subsets, next_edge.dest);

        // if including this edge does not cause a cycle
        if (x != y) {
            total_weight += next_edge.weight;
            edges_added++;
            Union(subsets, x, y);
        }
    }

    free(subsets);

    // if we couldn't connect all active nodes, this subset forms a disconnected 
    // forest, not a spanning tree. Return INT_MAX to invalidate this attempt.
    if (active_nodes == 0 || edges_added != active_nodes - 1) {
        return INT_MAX;
    }

    return total_weight;
}

int bruteForceSteinerTree(Graph* graph) {
    // sort edges once globally, rather than re-sorting inside every Kruskal call
    qsort(graph->edges, graph->E, sizeof(Edge), compareEdges);

    int V = graph->V;
    int* steiner_indices = (int*)malloc(V * sizeof(int)); 
    int num_steiner = 0;
    
    // identify which vertices are the optional Steiner points
    for (int v = 0; v < V; v++) {
        if (!graph->is_required[v]) {
            steiner_indices[num_steiner++] = v;
        }
    }

    int min_total_weight = INT_MAX;
    
    // total combinations is 2^num_steiner, use bitwise operations for simpler solutions
    int total_combinations = 1 << num_steiner; 

    // brute force loop generating the power set
    for (int i = 0; i < total_combinations; i++) {
        
        // reset active state: Terminals are ALWAYS on.
        for(int v = 0; v < V; v++) {
            graph->is_active[v] = graph->is_required[v]; 
        }
        
        // use bitwise AND to check the bits of our counter 'i'.
        // if the j-th bit is 1, turn on the j-th Steiner point.
        for (int j = 0; j < num_steiner; j++) {
            if (i & (1 << j)) {
                graph->is_active[steiner_indices[j]] = true;
            }
        }

        // run MST on this specific subset
        int current_weight = activeSubsetMST(graph);

        // update global minimum
        if (current_weight < min_total_weight) {
            min_total_weight = current_weight;
        }
    }

    free(steiner_indices);
    return min_total_weight;
}

// ===== Brute Force Functions End =====

int main() {
    clock_t start_time, end_time;
    double time_spent;

    // test Case 1: The Square Graph (SEE PDF PAGE 6) 
    // vertices 0, 1, 2, 3 are Required corners. Vertex 4 is the Steiner center.
    int V = 5;
    int E = 8;
    Graph* graph = createGraph(V, E);

    // set Required vs Steiner
    graph->is_required[0] = true;
    graph->is_required[1] = true;
    graph->is_required[2] = true;
    graph->is_required[3] = true;
    graph->is_required[4] = false; // the Steiner point in the middle

    // define the perimeter edges of weight 10
    graph->edges[0] = (Edge){0, 1, 10};
    graph->edges[1] = (Edge){1, 2, 10};
    graph->edges[2] = (Edge){2, 3, 10};
    graph->edges[3] = (Edge){3, 0, 10};

    // define the spokes to the center (weight 5)
    graph->edges[4] = (Edge){0, 4, 5};
    graph->edges[5] = (Edge){1, 4, 5};
    graph->edges[6] = (Edge){2, 4, 5};
    graph->edges[7] = (Edge){3, 4, 5};

    printf("--- SANITY CHECK: SQUARE GRAPH ---\n");
    
    start_time = clock();
    int brute_force_result = bruteForceSteinerTree(graph);
    end_time = clock();
    
    time_spent = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Result: %d (Expected: 20)\n", brute_force_result);
    printf("Time: %f seconds\n\n", time_spent);

    freeGraph(graph); // cleanup

    // TEST CASE: The Chain
    // A straight line of nodes where every alternating node is a terminal.
    // T1 -- S1 -- T2 -- S2 -- T3
    int V_chain = 5;
    int E_chain = 4;
    Graph* chainGraph = createGraph(V_chain, E_chain);

    chainGraph->is_required[0] = true;
    chainGraph->is_required[1] = false;
    chainGraph->is_required[2] = true;
    chainGraph->is_required[3] = false;
    chainGraph->is_required[4] = true;

    chainGraph->edges[0] = (Edge){0, 1, 2};
    chainGraph->edges[1] = (Edge){1, 2, 2};
    chainGraph->edges[2] = (Edge){2, 3, 2};
    chainGraph->edges[3] = (Edge){3, 4, 2};
    
    // Expected Optimal Cost: 8

    printf("--- SANITY CHECK: THE CHAIN ---\n");
    
    start_time = clock();
    int brute_force_chain_result = bruteForceSteinerTree(chainGraph);
    end_time = clock();
    
    time_spent = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Result: %d (Expected: 8)\n", brute_force_chain_result);
    printf("Time: %f seconds\n\n", time_spent);

    freeGraph(chainGraph); // cleanup

    // TEST CASE: The Central Hub
    // 3 Terminals on the outside, connected to a web of 10 Steiner points in the center.
    int V_hub = 13; 
    int E_hub = 15;
    Graph* hubGraph = createGraph(V_hub, E_hub);

    hubGraph->is_required[0] = true;
    hubGraph->is_required[1] = true;
    hubGraph->is_required[2] = true;

    for(int i = 3; i < 13; i++) {
        hubGraph->is_required[i] = false;
    }

    hubGraph->edges[0] = (Edge){0, 3, 5};
    hubGraph->edges[1] = (Edge){1, 7, 5};
    hubGraph->edges[2] = (Edge){2, 12, 5};

    int edge_idx = 3;
    for(int i = 3; i < 12; i++) {
        hubGraph->edges[edge_idx++] = (Edge){i, i+1, 10}; 
    }
    hubGraph->edges[edge_idx++] = (Edge){3, 7, 12};
    hubGraph->edges[edge_idx++] = (Edge){7, 12, 15};
    hubGraph->edges[edge_idx++] = (Edge){3, 12, 20};

    printf("--- SANITY CHECK: THE CENTRAL HUB (10 Steiner Points) ---\n");
    
    start_time = clock();
    int brute_force_hub_result = bruteForceSteinerTree(hubGraph);
    end_time = clock();
    
    time_spent = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Result: %d\n", brute_force_hub_result);
    printf("Time: %f seconds\n\n", time_spent);

    freeGraph(hubGraph);

    return 0;
}