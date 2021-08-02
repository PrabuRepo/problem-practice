package com.basic.datastructures.operations;

import java.util.List;

/*
 * Graph Overview:
 * 1.Basic Graph Theory
 * 2.Graph Representations: 
 *      - Adjacency Matrix Representation
 * 		- Adjacency List Representation - Using LL array, Using List of List, Using Map, Using Object
 * 		- Edge List Representation
 * 	    - Incidence Matrix Representation
 *  3.Graph Traversals: DFS/BFS
 *  4.Topological Sort: DFS/BFS(Kahn's Algorithm)
 *  5.Detecting Cycles: DFS/BFS/UF
 *  6.Graph Dynamic Connectivity:
 *      - Union-Find or DisJoint set Data Structure
 *  7.Shortest Path: 
 *      - BFS Algorithm to find the Shortest path
 *  	- Single Source Shortest Paths using Dijikstra's Algorithm
 *      - Single Source Shortest Paths(including negative edge weights) using Bellman-Ford Algorithm
 *      - All-pairs Shortest paths on Dense Graphs(Floyd-Warshall Algorithm)
 *      - All-pairs Shortest paths on Sparse Graphs(Johnson's Algorithm)
 *  8.Minimum Spanning Tree:
 *  	- Prim's Alg
 *      - Kruskal's Alg
 *  9.Connected Components:
 *      - Find Connected Components for Directed and Undirected Graphs
 *  	- Terminologies: Bridges, Articulation Points, Strongly Connected Components
 *      - Algorithms: Tarjan, Kosaraju
 *  10.Bidirectional Search
 *      - Best First
 *      - A* 
 *  11.Graph Color - Bipartite Graph
 *  12.Maximum Flow Problem - Ford Fulkerson
 *  13.Euler Path/Circuit - Fleury’s Algorithm
 *  14.Hamiltonian Path/Circuit - ??
 *  15.Grid Problems: DFS/BFS
 */
public interface GraphOperations {

	// DG - Directed Graph
	// UG - Undirected Graph
	void buildDirectedGraph(int[][] edges);

	void buildUndirectedGraph(int[][] edges);

	void buildWeightedDG(int[][] edges);

	void buildWeightedUG(int[][] edges);

	int findNumberOfNodes(int[][] edges);

	void printGraph();

	/* Note:
	 * Graph Traversals algorithms have 3 states: Unvisited, Visiting, Visited
	 *   Unvisited - Node/Vertex is not visited
	 *   Visiting - Node/Vertex is in recursion flow, not marked as visited 
	 *   Visited - Node/Vertex is visited 
	 */

	//DFS Recursive Approach
	List<Integer> dfs();

	// DFS Iterative Approach using Stack
	List<Integer> dfsIterative();

	// BFS Iterative Approach using Queue
	List<Integer> bfs();

	//Find the no of disconnected sub graphs in the graph - Using both DFS & BFS
	int findDisconnectedGraph();

	/* Topological Sort:
	 * Topological sorting for Directed Acyclic Graph (DAG) is a linear ordering of vertices such that for "every
	 * directed edge uv, vertex u comes before v in the ordering". Topological Sorting for a graph is not possible if the
	 * graph is not a DAG.
	 * Fact: A DAG G has at least one vertex with in-degree 0 and one vertex with out-degree 0.
	 *   Approach1: DFS Algorithm
	 *   Approach2: BFS Algorithm (Kahn's Algorithm)
	 */
	List<Integer> topologicalSort();

	//TODO: Revisit DG & UG graph cycle solutions and merge it as a single approach. Keep only one method here
	/*
	 * TODO: Solve this problem using 3 approaches:
	 * 	1. DFS Algorithm with additional recursion stack array
	 *  2. Using Graph Colors
	 *  3. Using DisjointSet
	 */
	boolean detectCycleInDG();

	boolean detectCycleInUG();

	// MST - Minimum Spanning Tree
	/* A minimum spanning tree (MST) or minimum weight spanning tree for a weighted, connected and undirected graph is a spanning tree with
	 * weight less than or equal to the weight of every other spanning tree. The weight of a spanning tree is the sum of weights given to 
	 * each edge of the spanning tree. A minimum spanning tree has (V - 1) edges where V is the number of vertices in the given graph.
	 */
	void mstPrimsAlg();

	void mstKruskalsAlg();

	// SP - Shortest Path Alg
	void spDijikstraAlg(int source);

	/*
	 * Dijkstra follows Greed Alg; Bellmanford alg follows DP algorithms and it finds shortest paths from src to all vertices in the given graph.
	 * Dijkstra doesn't work for Graphs with negative weight edges, Bellman-Ford works for such graphs. Time complexity of Bellman-Ford is 
	 * O(VE), which is more than Dijkstra. 
	 * 
	 * How does this work? Like other Dynamic Programming Problems, the algorithm calculate shortest paths in bottom-up manner. It first
	 * calculates the shortest distances which have at-most one edge in the path. Then, it calculates shortest paths with at-most 2 edges, 
	 * and so on. After the i-th iteration of outer loop, the shortest paths with at most i edges are calculated. There can be maximum 
	 * |V| - 1 edges in any simple path, that is why the outer loop runs |v|- 1 times.
	 */
	void spBellmanFordAlg(int source);

	/* The Floyd Warshall Algorithm (Dynamic Programming) is for solving the All Pairs Shortest Path problem. The problem is 
	 * to find shortest distances between every pair of vertices in a given edge weighted directed Graph.
	 */
	void spFloydWarshallAlg();

	/* 
	 * Tarjan's strongly connected components algorithm performs a single pass of depth first search. It maintains a stack of vertices that have been
	 * explored by the search but not yet assigned to a component, and calculates "low numbers" of each vertex (an index number of the highest ancestor
	 * reachable in one step from a descendant of the vertex) which it uses to determine when a set of vertices should be popped off the stack into 
	 * a new component.
	 * Time: O(V+E) => O(V+E); Space: O(V+E)
	 */
	void tarjanAlg();

	/*
	 * Kosaraju's algorithm uses two passes of depth first search. The first, in the original graph, is used to choose the order in which the outer 
	 * loop of the second depth first search tests vertices for having been visited already and recursively explores them if not. The second depth 
	 * first search is on the transpose graph of the original graph, and each recursive exploration finds a single new strongly connected component.
	 * Time: 3*O(V+E) => O(V+E); Space: O(V+E)
	 */
	void kosarajuAlg();

	/*
	 * A Bipartite Graph is a graph whose vertices can be divided into two independent sets, U and V such that every edge (u, v) either connects
	 * a vertex from U to V or a vertex from V to U.
	 * It can be solved using DFS & BFS
	 */
	boolean isBipartite();
}