//https://networkx.org/documentation/stable/_modules/networkx/algorithms/approximation/steinertree.html#steiner_tree


#include <stdio.h>
#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h> 



// =============================================
// STEINER TREE APPROXIMATION (KMB Algorithm)
// =============================================
//   1. Find shortest paths between terminals (Dijkstra)
//   2. Build a minimum spanning tree on those paths (Prim)
//   3. Expand the MST back into the original graph

#define MAX_NODES 100          
#define INFINITY INT_MAX
int total_graph_operations = 0; // global counter for total graph operations

typedef struct {
    int src;
    int dest;
    int weight;
} Edge;

typedef struct {
    int V;
    int E;
    Edge* edges;
    bool* is_required;  // true = Terminal, false = Steiner Point
} Graph;

Graph* createGraph(int V, int E) {
    Graph* graph = (Graph*)malloc(sizeof(Graph));
    graph->V = V;
    graph->E = E;
    graph->edges = (Edge*)malloc(E * sizeof(Edge));
    graph->is_required = (bool*)calloc(V, sizeof(bool));
    return graph;
}

void freeGraph(Graph* graph) {
    if (graph == NULL) return;
    if (graph->edges != NULL) free(graph->edges);
    if (graph->is_required != NULL) free(graph->is_required);
    free(graph);
}

// ===============================================================================================
int totalNodes;                          // how many nodes are in the graph
int originalGraph[MAX_NODES][MAX_NODES]; // input graph 
int steinerGraph[MAX_NODES][MAX_NODES];  // resulting Steiner tree edges

int terminalNodes[MAX_NODES]; // nodes to connect (required nodes)
int totalTerminals;           // num of terminal nodes

int mstParent[MAX_NODES];     // stores parent of each node in MST

// stores the shortest distance between every pair of terminal nodes
int terminalDistances[MAX_NODES][MAX_NODES];

void loadGraph(Graph* graph) {
    totalNodes = graph->V;

    for (int i = 0; i < totalNodes; i++)
        for (int j = 0; j < totalNodes; j++)
            originalGraph[i][j] = 0;

    for (int i = 0; i < graph->E; i++) {
        int src    = graph->edges[i].src;
        int dest   = graph->edges[i].dest;
        int weight = graph->edges[i].weight;
        originalGraph[src][dest] = weight;
        originalGraph[dest][src] = weight;
    }

    totalTerminals = 0;
    for (int v = 0; v < graph->V; v++)
        if (graph->is_required[v])
            terminalNodes[totalTerminals++] = v;
}



// ============================== DIJKSTRA'S ==============================
// Finds the shortest distance from startnode 
// After running:
//   shortestDist[i] = shortest distance from startNode to node i
//   cameFrom[i]     = came from which node to reach node i

void findShortestPaths(int startNode, int shortestDist[], int cameFrom[]) {
    int alreadyVisited[MAX_NODES] = {0}; // tracks which finalized nodes 

    // initialize all distances, no parent
    for (int i = 0; i < totalNodes; i++) {
        shortestDist[i] = INFINITY;
        cameFrom[i] = -1; // -1 means no parent yet
    }

    shortestDist[startNode] = 0; // distance to itself is 0

    // repeat for every node
    for (int step = 0; step < totalNodes - 1; step++) {

        // pick unvisited node with smallest known distance
        int currentNode = -1;
        for (int j = 0; j < totalNodes; j++) {
            if (!alreadyVisited[j] && (currentNode == -1 || shortestDist[j] < shortestDist[currentNode])) {
                currentNode = j;
            }
        }

        alreadyVisited[currentNode] = 1; // mark as finalized

        // update distances to neighbors of currentNode
        for (int neighbor = 0; neighbor < totalNodes; neighbor++) {
            total_graph_operations++; // track edge evaluations
            int edgeWeight = originalGraph[currentNode][neighbor];
            int newDist = shortestDist[currentNode] + edgeWeight;

            // if theres an edge AND this new path is shorter, update
            if (edgeWeight > 0 &&
                !alreadyVisited[neighbor] &&
                newDist <= shortestDist[neighbor]) {
                shortestDist[neighbor] = newDist;
                cameFrom[neighbor] = currentNode;
            }
        }
    }
}


// ============================== BUILD THE "METRIC GRAPH" ==============================
// Creates a smaller graph where:
//   - nodes are ONLY the terminal nodes
//   - edge weights = shortest path distance between terminal pairs
//

void buildTerminalDistanceGraph() {
    int shortestDist[MAX_NODES];
    int cameFrom[MAX_NODES];

    // for each terminal, run Dijkstra and record distances to other terminals
    for (int i = 0; i < totalTerminals; i++) {
        findShortestPaths(terminalNodes[i], shortestDist, cameFrom);

        for (int j = 0; j < totalTerminals; j++) {
            // store shortest distance from terminal i to terminal j
            terminalDistances[i][j] = shortestDist[terminalNodes[j]];
        }
    }
}


// ============================== PRIM'S MINIMUM SPANNING TREE ==============================
// After running:
//   mstParent[i] = which terminal index is the parent of terminal index i

void buildMinimumSpanningTree(int distanceGraph[MAX_NODES][MAX_NODES], int nodeCount) {
    int cheapestEdge[MAX_NODES]; // cheapest known edge to reach each node
    int inMST[MAX_NODES];        // whether a node is already in the MST

    // start with all nodes unreached
    for (int i = 0; i < nodeCount; i++) {
        cheapestEdge[i] = INFINITY;
        inMST[i] = 0;
    }

    cheapestEdge[0] = 0;    // start from terminal index 0
    mstParent[0] = -1;      // has no parent (root))

    for (int step = 0; step < nodeCount - 1; step++) {

        // pick node not yet in MST with cheapest connection
        int currentNode = -1;
        for (int j = 0; j < nodeCount; j++) {
            if (!inMST[j] && (currentNode == -1 || cheapestEdge[j] < cheapestEdge[currentNode])) {
                currentNode = j;
            }
        }

        inMST[currentNode] = 1; // add to MST

        // update neighbors: if connecting through currentNode is cheaper, update
        for (int neighbor = 0; neighbor < nodeCount; neighbor++) {
            total_graph_operations++; // tracking edge evaluations
            int edgeCost = distanceGraph[currentNode][neighbor];

            if (edgeCost != INFINITY && !inMST[neighbor] && edgeCost < cheapestEdge[neighbor]) {
                mstParent[neighbor] = currentNode; 
                cheapestEdge[neighbor] = edgeCost;
            }
        }
    }
}


// ============================== ADD A PATH TO THE STEINER GRAPH ==============================
// Traces back from 'destinationNode' to 'sourceNode'  using the 'cameFrom' parent array (from Dijkstra), and
//  adds each edge along that path to the Steiner graph

void addShortestPathToSteinerGraph(int sourceNode, int destinationNode, int cameFrom[]) {
    int currentNode = destinationNode;

    // walk backwards from destination to source using parent pointers
    while (currentNode != sourceNode) {
        int parentNode = cameFrom[currentNode];

        // add edge between parentNode and currentNode to Steiner graph
        steinerGraph[parentNode][currentNode] = originalGraph[parentNode][currentNode];
        steinerGraph[currentNode][parentNode] = originalGraph[currentNode][parentNode];

        currentNode = parentNode; // move one step back toward source
    }
}


// ============================== BUILD THE STEINER TREE ==============================
// For each edge in MST (between terminal pairs), find shortest path in  
// original graph and add it to Steiner graph


void buildSteinerTree() {
    int shortestDist[MAX_NODES];
    int cameFrom[MAX_NODES];

    // clear the Steiner graph first
    for (int i = 0; i < totalNodes; i++)
        for (int j = 0; j < totalNodes; j++)
            steinerGraph[i][j] = 0;

    // for each terminal (except the root), connect it to its MST parent
    for (int i = 1; i < totalTerminals; i++) {
        int parentTerminal = terminalNodes[mstParent[i]]; // actual node of parent
        int childTerminal  = terminalNodes[i];            // actual node of child

        // find shortest path from parent terminal to child terminal
        findShortestPaths(parentTerminal, shortestDist, cameFrom);

        // add that paths edges into the Steiner graph
        addShortestPathToSteinerGraph(parentTerminal, childTerminal, cameFrom);
    }
}


// ─────────────────────────────────────────────
// STEP 4: PRUNE UNNECESSARY LEAF NODES
// ─────────────────────────────────────────────
// After building the Steiner graph, some non-terminal
// nodes might be "dangling" (only one connection)
//  add cost but connect nothing useful, so remove
// ─────────────────────────────────────────────
void removeLeaves() {
    int isTerminal[MAX_NODES] = {0};

    // mark which nodes are terminals
    for (int i = 0; i < totalTerminals; i++)
        isTerminal[terminalNodes[i]] = 1;

    int somethingChanged = 1;

    // keep pruning until no more useless leaves exist
    while (somethingChanged) {
        somethingChanged = 0;

        for (int node = 0; node < totalNodes; node++) {
            if (isTerminal[node]) continue; // never remove terminal

            // count how many edges this node has in Steiner graph
            int edgeCount = 0;
            for (int j = 0; j < totalNodes; j++)
                if (steinerGraph[node][j]) edgeCount++;

            // if leaf (only 1 connection) and not terminal, remove
            if (edgeCount == 1) {
                for (int j = 0; j < totalNodes; j++)
                    steinerGraph[node][j] = steinerGraph[j][node] = 0;

                somethingChanged = 1; // removing might create new leaves, check again
            }
        }
    }
}


// ============================== OUTPUT ==============================
void printSteinerTree() {
    int totalCost = 0;

    printf("\n===== Steiner Tree Edges =====\n");

    for (int i = 0; i < totalNodes; i++) {
        for (int j = i + 1; j < totalNodes; j++) {
            if (steinerGraph[i][j]) {
                printf("  Node %d -- Node %d  (weight: %d)\n", i, j, steinerGraph[i][j]);
                totalCost += steinerGraph[i][j];
            }
        }
    }

    printf("==============================\n");
    printf("Total Cost: %d\n", totalCost);
}


void runTestCase(Graph* graph, const char* title) {
    clock_t start_time, end_time;
    double time_spent;

    total_graph_operations = 0; // reset counter per test
    printf("\n--- SANITY CHECK: %s ---\n", title);

    start_time = clock(); 

    loadGraph(graph);
    buildTerminalDistanceGraph();
    buildMinimumSpanningTree(terminalDistances, totalTerminals);
    buildSteinerTree();
    removeLeaves();
    printSteinerTree();

    end_time = clock();    

    time_spent = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("  -> Total Graph Operations (Edge Lookups): %d\n", total_graph_operations);
    printf("Time: %f seconds\n", time_spent);

    freeGraph(graph);
}

int main() {
    // TEST CASE 1: Square Graph
    int V = 5;
    int E = 8;
    Graph* graph = createGraph(V, E);

    graph->is_required[0] = true;
    graph->is_required[1] = true;
    graph->is_required[2] = true;
    graph->is_required[3] = true;
    graph->is_required[4] = false;

    graph->edges[0] = (Edge){0, 1, 10};
    graph->edges[1] = (Edge){1, 2, 10};
    graph->edges[2] = (Edge){2, 3, 10};
    graph->edges[3] = (Edge){3, 0, 10};

    graph->edges[4] = (Edge){0, 4, 5};
    graph->edges[5] = (Edge){1, 4, 5};
    graph->edges[6] = (Edge){2, 4, 5};
    graph->edges[7] = (Edge){3, 4, 5};

    runTestCase(graph, "SQUARE GRAPH");


    // TEST CASE 2: Chain
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

    runTestCase(chainGraph, "THE CHAIN");


    // TEST CASE 3: Central Hub
    int V_hub = 13;
    int E_hub = 15;
    Graph* hubGraph = createGraph(V_hub, E_hub);

    hubGraph->is_required[0] = true;
    hubGraph->is_required[1] = true;
    hubGraph->is_required[2] = true;

    for (int i = 3; i < 13; i++) {
        hubGraph->is_required[i] = false;
    }

    hubGraph->edges[0] = (Edge){0, 3, 5};
    hubGraph->edges[1] = (Edge){1, 7, 5};
    hubGraph->edges[2] = (Edge){2, 12, 5};

    int edge_idx = 3;
    for (int i = 3; i < 12; i++) {
        hubGraph->edges[edge_idx++] = (Edge){i, i + 1, 10};
    }

    hubGraph->edges[edge_idx++] = (Edge){3, 7, 12};
    hubGraph->edges[edge_idx++] = (Edge){7, 12, 15};
    hubGraph->edges[edge_idx++] = (Edge){3, 12, 20};

    runTestCase(hubGraph, "THE CENTRAL HUB (10 Steiner Points)");

    return 0;
}