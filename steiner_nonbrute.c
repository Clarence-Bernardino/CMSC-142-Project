//https://networkx.org/documentation/stable/_modules/networkx/algorithms/approximation/steinertree.html#steiner_tree


#include <stdio.h>

// =============================================
// STEINER TREE APPROXIMATION (KMB Algorithm)
// =============================================
//   1. Find shortest paths between terminals (Dijkstra)
//   2. Build a minimum spanning tree on those paths (Prim)
//   3. Expand the MST back into the original graph

#define MAX_NODES 100          
#define INFINITY  1000000000   

int totalNodes;                          // how many nodes are in the graph
int originalGraph[MAX_NODES][MAX_NODES]; // input graph 
int steinerGraph[MAX_NODES][MAX_NODES];  // resulting Steiner tree edges

int terminalNodes[MAX_NODES]; // nodes to connect (required nodes)
int totalTerminals;           // num of terminal nodes

int mstParent[MAX_NODES];     // stores parent of each node in MST

// stores the shortest distance between every pair of terminal nodes
int terminalDistances[MAX_NODES][MAX_NODES];


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
            int edgeWeight = originalGraph[currentNode][neighbor];
            int newDist = shortestDist[currentNode] + edgeWeight;

            // if theres an edge AND this new path is shorter, update
            if (edgeWeight > 0 && newDist < shortestDist[neighbor]) {
                shortestDist[neighbor] = newDist;
                cameFrom[neighbor] = currentNode; // came from currentNode
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
// nodes might be "dangling" (only one connection).
// They add cost but connect nothing useful, so remove them.
// ─────────────────────────────────────────────
void removeUselessLeaves() {
    int isTerminal[MAX_NODES] = {0};

    // mark which nodes are terminals
    for (int i = 0; i < totalTerminals; i++)
        isTerminal[terminalNodes[i]] = 1;

    int somethingChanged = 1;

    // keep pruning until no more useless leaves exist
    while (somethingChanged) {
        somethingChanged = 0;

        for (int node = 0; node < totalNodes; node++) {
            if (isTerminal[node]) continue; // never remove a terminal

            // count how many edges this node has in the Steiner graph
            int edgeCount = 0;
            for (int j = 0; j < totalNodes; j++)
                if (steinerGraph[node][j]) edgeCount++;

            // if it's a leaf (only 1 connection) and not a terminal, remove it
            if (edgeCount == 1) {
                for (int j = 0; j < totalNodes; j++)
                    steinerGraph[node][j] = steinerGraph[j][node] = 0;

                somethingChanged = 1; // removing it might create new leaves, so check again
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


int main() {
    int totalEdges;

    printf("Enter number of nodes: ");
    scanf("%d", &totalNodes);

    // initialize graph with no edges
    for (int i = 0; i < totalNodes; i++)
        for (int j = 0; j < totalNodes; j++)
            originalGraph[i][j] = 0;

    printf("Enter number of edges: ");
    scanf("%d", &totalEdges);

    printf("Enter each edge as: node1 node2 weight\n");
    for (int i = 0; i < totalEdges; i++) {
        int nodeA, nodeB, weight;
        scanf("%d %d %d", &nodeA, &nodeB, &weight);

        // store in both directions, graph undirected
        originalGraph[nodeA][nodeB] = weight;
        originalGraph[nodeB][nodeA] = weight;
    }

    printf("Enter number of required nodes (nodes that MUST be connected): ");
    scanf("%d", &totalTerminals);

    printf("Enter required node numbers:\n");
    for (int i = 0; i < totalTerminals; i++) {
        scanf("%d", &terminalNodes[i]);
    }

    // ============================== Run algo  ==============================

    buildTerminalDistanceGraph(); // shortest paths between terminals
    buildMinimumSpanningTree(terminalDistances, totalTerminals); // MST on terminals
    buildSteinerTree();           // expand MST back to original graph
    removeUselessLeaves();        //  clean up dangling non-terminal nodes

    printSteinerTree();           // Output

    return 0;
}