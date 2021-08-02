package com.basic.datastructures;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Queue;
import java.util.Scanner;
import java.util.Set;
import java.util.Stack;

import com.basic.datastructures.operations.GraphOperations;
import com.common.model.EdgeNode;
import com.common.model.GraphNode;
import com.common.utilities.DisjointSet;

/* 
 * Graph DS & Algorithms 
 */
public class Graph {
	public static void main(String[] args) {
		/*GraphOperations graphMatrix = new Graph().new GraphAdjMatrix();
		CommonUtil.testGraph(graphMatrix);
		
		GraphOperations graphAdjList = new Graph().new GraphAdjList();
		CommonUtil.testGraph(graphAdjList);
		*/
		GraphOperations graphAdjMap = new Graph().new GraphAdjListUsingMap();
		CommonUtil.testGraph(graphAdjMap);

		/*GraphOperations graphEdgeList = new Graph().new GraphEdgeList();
		CommonUtil.testGraph(graphEdgeList);*/
	}

	/**
	 * Adjacency Matrix representation: A matrix indicating with each vertex associated with a row and
	 * column.
	 */
	class GraphAdjMatrix implements GraphOperations {
		int N;
		int[][] adjMatrix;

		@Override
		public void buildDirectedGraph(int[][] edges) {
			N = findNumberOfNodes(edges);
			adjMatrix = new int[N][N];
			for (int[] edge : edges) {
				adjMatrix[edge[0]][edge[1]] = 1;
			}
		}

		@Override
		public void buildUndirectedGraph(int[][] edges) {
			N = findNumberOfNodes(edges);
			adjMatrix = new int[N][N];
			for (int[] edge : edges) {
				adjMatrix[edge[0]][edge[1]] = 1;
				adjMatrix[edge[1]][edge[0]] = 1;
			}
		}

		@Override
		public void buildWeightedDG(int[][] edges) {
			N = findNumberOfNodes(edges);
			adjMatrix = new int[N][N];
			for (int[] edge : edges) {
				adjMatrix[edge[0]][edge[1]] = edge[2];
			}
		}

		@Override
		public void buildWeightedUG(int[][] edges) {
			N = findNumberOfNodes(edges);
			adjMatrix = new int[N][N];
			for (int[] edge : edges) {
				adjMatrix[edge[0]][edge[1]] = edge[2];
				adjMatrix[edge[1]][edge[0]] = edge[2];
			}
		}

		// If input is given as direct matrix
		public void buildDirectedGraphUsingMatrix(int[][] matrix) {
			this.N = matrix.length;
			this.adjMatrix = matrix;
		}

		@Override
		public int findNumberOfNodes(int[][] edges) {
			Set<Integer> set = CommonUtil.findNumberOfNodes(edges);
			return set.size();
		}

		@Override
		public void printGraph() {
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
					System.out.print(adjMatrix[i][j] + " ");
				}
				System.out.println();
			}
		}

		@Override
		public List<Integer> dfs() {
			boolean[] visited = new boolean[N];
			List<Integer> result = new ArrayList<>();
			for (int i = 0; i < N; i++) { //Iteration is to handle the disconnected graph
				if (!visited[i]) {
					dfsUtil(adjMatrix, visited, i, result);
				}
			}
			return result;
		}

		private void dfsUtil(int[][] adjMatrix, boolean[] visited, int s, List<Integer> result) {
			int n = adjMatrix.length;
			visited[s] = true;
			result.add(s);
			for (int j = 0; j < n; j++) {
				if (!visited[j] && j != s && adjMatrix[s][j] == 1) {
					dfsUtil(adjMatrix, visited, j, result);
				}
			}
		}

		@Override
		public List<Integer> dfsIterative() {
			return null;
		}

		@Override
		public List<Integer> bfs() {
			boolean[] visited = new boolean[N];
			List<Integer> result = new ArrayList<>();
			for (int i = 0; i < N; i++) { //Iteration is to handle the disconnected graph
				if (!visited[i]) {
					bfsUtil(adjMatrix, visited, i, result);
				}
			}
			result.forEach(k -> System.out.print(k + " - "));
			return result;
		}

		private void bfsUtil(int[][] adjMatrix, boolean[] visited, int s, List<Integer> result) {
			int n = adjMatrix.length;
			Queue<Integer> queue = new LinkedList<>();
			queue.add(s);
			while (!queue.isEmpty()) {
				int top = queue.poll();
				visited[top] = true;
				result.add(top);
				for (int j = 0; j < n; j++) {
					if (!visited[j] && j != top && adjMatrix[top][j] == 1) queue.add(j);
				}
			}
		}

		@Override
		public int findDisconnectedGraph() {
			int count = dfsDisconnectedGraph();
			count = bfsDisconnectedGraph();
			return 0;
		}

		public int dfsDisconnectedGraph() {
			List<Integer> result = new ArrayList<>();
			int groups = 0; //No of disconnected graphs
			boolean[] visited = new boolean[N];
			for (int i = 0; i < N; i++) {
				if (!visited[i]) {
					groups++;
					dfsUtil(adjMatrix, visited, i, result);
				}
			}
			return groups;
		}

		public int bfsDisconnectedGraph() {
			if (adjMatrix.length == 0) return 0;
			List<Integer> result = new ArrayList<>();
			int groups = 0;
			boolean[] visited = new boolean[N];
			for (int i = 0; i < N; i++) {
				if (!visited[i]) {
					groups++;
					bfsUtil(adjMatrix, visited, i, result);
				}
			}
			return groups;
		}

		@Override
		public List<Integer> topologicalSort() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public boolean detectCycleInDG() {
			return false;
		}

		@Override
		public boolean detectCycleInUG() {
			return false;
		}

		@Override
		public void mstPrimsAlg() {
			int[] edgeWeight = new int[N]; // Edge Weight from source vertex
			int[] parent = new int[N]; // Parent used to print the path
			boolean[] visited = new boolean[N];

			// Set max values in edgeWeight array
			Arrays.fill(edgeWeight, Integer.MAX_VALUE);
			edgeWeight[0] = 0;
			parent[0] = -1;

			for (int i = 0; i < N; i++) {
				// Find the minimum value index in the edgeWeight array
				int u = GraphUtil.findMinWeight(edgeWeight, visited);
				visited[u] = true;

				for (int v = 0; v < N; v++) {
					if (adjMatrix[u][v] != 0 && !visited[v] && adjMatrix[u][v] < edgeWeight[v]) {
						parent[v] = u;
						edgeWeight[v] = adjMatrix[u][v];
					}
				}
			}

			GraphUtil.printMST(parent, edgeWeight, N);
		}

		@Override
		public void mstKruskalsAlg() {
			// TODO Auto-generated method stub

		}

		// Time Complexity: O(n^2)
		@Override
		public void spDijikstraAlg(int source) {
			int[] edgeWeight = new int[N]; // Edge Weight from source vertex
			boolean[] visited = new boolean[N];

			// Set max values in edgeWeight array
			Arrays.fill(edgeWeight, Integer.MAX_VALUE);

			edgeWeight[source] = 0;
			// parent[0] = -1;

			for (int i = 0; i < N; i++) {
				int u = GraphUtil.findMinWeight(edgeWeight, visited); // Find the minimum value index in the edgeWeight array
				visited[u] = true;

				for (int v = 0; v < N; v++)
					if (adjMatrix[u][v] != 0 && !visited[v])
						edgeWeight[v] = Math.min(edgeWeight[v], edgeWeight[u] + adjMatrix[u][v]);
			}

			GraphUtil.printSP(edgeWeight, N);
		}

		@Override
		public void spBellmanFordAlg(int source) {

		}

		@Override
		public void spFloydWarshallAlg() {
			int[][] dist = new int[N][N];

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					dist[i][j] = adjMatrix[i][j];

			for (int v = 0; v < N; v++) // Via each vertex/node one by one
				for (int i = 0; i < N; i++)
					for (int j = 0; j < N; j++)
						if (dist[i][v] != Integer.MAX_VALUE && dist[v][j] != Integer.MAX_VALUE)
							dist[i][j] = Math.min(dist[i][j], dist[i][v] + dist[v][j]);

			GraphUtil.printfloydWarshallSP(dist, N);
		}

		@Override
		public void kosarajuAlg() {
			// TODO Auto-generated method stub
		}

		@Override
		public void tarjanAlg() {
			// TODO Auto-generated method stub
		}

		@Override
		public boolean isBipartite() {
			// TODO Auto-generated method stub
			return false;
		}

	}

	/**
	 * Adjacency List representation:Container of vertices, and each vertex has a list of adjacent
	 * vertices
	 */
	class GraphAdjList implements GraphOperations {
		int N;
		LinkedList<Integer>[] adjList;
		LinkedList<GraphNode>[] adjListW;

		// This representation is convenient for non-sequence nodes, char, string etc.
		List<List<Integer>> adjList2 = new ArrayList<>();

		@Override
		public void buildDirectedGraph(int[][] edges) {
			this.N = findNumberOfNodes(edges);
			// Create the instance for each node
			adjList = new LinkedList[N];
			for (int i = 0; i < N; i++)
				adjList[i] = new LinkedList<>();

			// Add edges in the adjList
			for (int[] edge : edges)
				adjList[edge[0]].add(edge[1]);

		}

		@Override
		public void buildUndirectedGraph(int[][] edges) {
			// TODO Auto-generated method stub

		}

		@Override
		public void buildWeightedDG(int[][] edges) {
			this.N = findNumberOfNodes(edges);
			// Create the instance for each node
			adjListW = new LinkedList[N];
			for (int i = 0; i < N; i++)
				adjListW[i] = new LinkedList<>();

			// Add edges in the adjList
			for (int[] edge : edges)
				adjListW[edge[0]].add(new GraphNode(edge[1], edge[2]));
		}

		@Override
		public void buildWeightedUG(int[][] edges) {
			// TODO Auto-generated method stub

		}

		@Override
		public int findNumberOfNodes(int[][] edges) {
			Set<Integer> set = CommonUtil.findNumberOfNodes(edges);
			return set.size();
		}

		@Override
		public void printGraph() {
			for (int i = 0; i < N; i++) {
				System.out.println("\nEdges from Vertex: " + i + "->");
				adjList[i].forEach(v -> System.out.print(v + ", "));
				//or
				//ListIterator<Integer> iterator = adjList[i].listIterator();
				//while (iterator.hasNext()) System.out.print(iterator.next() + ", ");
			}
		}

		@Override
		public List<Integer> dfs() {
			boolean[] visited = new boolean[N];
			List<Integer> result = new ArrayList<>();

			for (int i = 0; i < N; i++) { //Iteration is to handle the disconnected graph
				if (!visited[i]) {
					dfs(i, visited, result);
				}
			}
			result.forEach(k -> System.out.print(k + " - "));
			return result;
		}

		private void dfs(int v, boolean[] visited, List<Integer> result) {
			visited[v] = true;
			result.add(v);
			ListIterator<Integer> listIterator = adjList[v].listIterator();
			while (listIterator.hasNext()) {
				int data = listIterator.next();
				if (!visited[data]) {
					dfs(data, visited, result);
				}
			}
		}

		@Override
		public List<Integer> dfsIterative() {
			boolean[] visited = new boolean[N];
			List<Integer> result = new ArrayList<>();

			for (int i = 0; i < N; i++) { //Iteration is to handle the disconnected graph
				if (!visited[i]) {
					dfsIterative(i, visited, result);
				}
			}
			result.forEach(k -> System.out.print(k + " - "));
			return result;
		}

		private List<Integer> dfsIterative(int v, boolean[] visited, List<Integer> result) {
			Stack<Integer> stack = new Stack<>();
			visited[v] = true;
			stack.push(v);
			while (!stack.isEmpty()) {
				int data = stack.pop();
				result.add(data);
				ListIterator<Integer> listIterator = adjList[data].listIterator();
				while (listIterator.hasNext()) {
					int next = listIterator.next();
					if (!visited[next]) {
						stack.push(next);
						visited[next] = true;
					}
				}
			}

			result.forEach(k -> System.out.print(k + " - "));
			return result;
		}

		@Override
		public List<Integer> bfs() {
			boolean[] visited = new boolean[N];
			List<Integer> result = new ArrayList<>();

			for (int i = 0; i < N; i++) {//Iteration is to handle the disconnected graph
				if (!visited[i]) {
					bfs(i, visited, result);
				}
			}

			result.forEach(k -> System.out.print(k + " - "));
			return result;
		}

		private void bfs(int vertex, boolean[] visited, List<Integer> result) {
			Queue<Integer> queue = new LinkedList<>();
			queue.add(vertex);
			visited[vertex] = true;
			while (!queue.isEmpty()) {
				int data = queue.poll();
				result.add(data);
				ListIterator<Integer> list = adjList[data].listIterator();
				while (list.hasNext()) {
					int next = list.next();
					if (!visited[next]) {
						visited[next] = true;
						queue.add(next);
					}
				}
			}
		}

		public List<List<Integer>> bfsLevelByLevel() {
			boolean[] visited = new boolean[N];
			List<List<Integer>> result = new ArrayList<>();
			for (int i = 0; i < N; i++) {
				if (!visited[i]) {
					bfsLevelByLevel(i, visited, result);
				}
			}
			return result;
		}

		public void bfsLevelByLevel(int source, boolean[] visited, List<List<Integer>> result) {
			Queue<Integer> queue = new LinkedList<>();
			queue.add(source);
			visited[source] = true;

			while (!queue.isEmpty()) {
				int size = queue.size();
				List<Integer> curr = new ArrayList<>();
				while (size-- > 0) {
					int data = queue.poll();
					curr.add(data);
					ListIterator<Integer> list = adjList[data].listIterator();
					while (list.hasNext()) {
						int next = list.next();
						if (!visited[next]) {
							visited[next] = true;
							queue.add(next);
						}
					}
				}
				result.add(new ArrayList<>(curr));
			}

		}

		@Override
		public int findDisconnectedGraph() {
			// TODO Auto-generated method stub
			return 0;
		}

		@Override
		public List<Integer> topologicalSort() {
			List<Integer> result = topologicalSortDfs();
			result = topologicalSortBfs();
			return result;
		}

		public List<Integer> topologicalSortDfs() {
			boolean[] visited = new boolean[N], onStack = new boolean[N];
			LinkedList<Integer> result = new LinkedList<>();

			for (int i = 0; i < N; i++)
				if (!visited[i]) {
					if (topoSortUtil(i, visited, onStack, result)) {
						//If there any cycle return null
						return null;
					}
				}

			return result;
		}

		//This solution is DFS Algorithm + Detect Cycle Algo
		private boolean topoSortUtil(int v, boolean[] visited, boolean[] onStack, LinkedList<Integer> result) {
			// If this condition satisfies, then graph contains cycle
			if (onStack[v]) return true;

			if (visited[v]) return false;

			// Mark vertex as visited and set recursion stack
			visited[v] = true;
			onStack[v] = true;

			ListIterator<Integer> listIterator = adjList[v].listIterator();
			while (listIterator.hasNext()) {
				int next = listIterator.next();
				if (topoSortUtil(next, visited, onStack, result)) return true;
			}

			result.addFirst(v);
			// Reset the recursion stack/Reset after visited the node 
			onStack[v] = false;
			return false;
		}

		public List<Integer> topologicalSortBfs() {
			Queue<Integer> queue = new LinkedList<>();
			List<Integer> linearOrder = new ArrayList<>();
			int[] indegree;
			int count = 0;

			// Step-1: Compute in-degree
			indegree = GraphUtil.indegree(adjList, N);

			// Step-2: Pick all the vertices with in-degree as 0 and add them into a queue
			for (int i = 0; i < N; i++)
				if (indegree[i] == 0) queue.add(i);

			// Step-3:Remove a vertex from the queue
			while (!queue.isEmpty()) {
				int vertex = queue.poll();
				linearOrder.add(vertex);
				// 1.Increment count of visited nodes by 1.
				count++;

				if (adjList[vertex] == null || adjList[vertex].isEmpty()) continue;

				// 2.Decrease in-degree by 1 for all its neighboring nodes
				for (int adjNode : adjList[vertex]) {
					// 3.If in-degree of a neighboring nodes is reduced to zero, then add it to the queue.
					if (--indegree[adjNode] == 0) queue.add(adjNode);
				}
				/*
				ListIterator<Integer> iter = adjList[vertex].listIterator();
				while (iter.hasNext()) {
					int data = iter.next();
					if (--indegree[data] == 0) queue.add(data); // 3.If in-degree of a neighboring nodes is reduced to zero,
																// then add it to the queue.
				}*/
			}
			// Step-4:If count of visited nodes is equal to the number of nodes in the graph then print the topological sort
			if (count == N) {
				linearOrder.forEach(i -> System.out.print(i + " "));
			} else {
				System.out.println("Graph is not an a DAG and also it contains a cycle.");
			}

			return linearOrder;
		}

		/*
		 * Depth First Traversal can be used to detect a cycle in a Graph. 
		 *  - DFS for a connected graph produces a tree. There is a cycle in a graph only if there is a back edge present in the  graph. 
		 *    A back edge is an edge that is from a node to  itself (self-loop) or one of its ancestors in the tree produced by DFS. 
		 *  - For a disconnected graph, Get the DFS forest as output. To detect cycle, check for a cycle in individual trees by checking 
		 *    back edges.
		 *    
		 *  - To detect a back edge, keep track of vertices currently in the recursion stack of function for DFS traversal. If a vertex is
		 *    reached that is already in the recursion stack, then there is a cycle in the tree. The edge that connects the current vertex 
		 *    to the vertex in the recursion stack is a back edge. Use onStack[] array to keep track of vertices in the recursion stack.
		 */
		@Override
		public boolean detectCycleInDG() {
			//onStack maintains 'Visiting' state nodes.
			boolean[] visited = new boolean[N], onStack = new boolean[N];
			for (int i = 0; i < N; i++) {
				if (!visited[i]) {
					if (hasCycle(i, visited, onStack)) return true;
				}
			}
			return false;
		}

		private boolean hasCycle(int vertex, boolean[] visited, boolean[] onStack) {
			// If this condition satisfies, then graph contains cycle
			if (onStack[vertex]) return true;

			// Mark vertex as visited and set recursion stack
			visited[vertex] = true;
			onStack[vertex] = true;

			if (adjList[vertex] != null) {
				ListIterator<Integer> iter = adjList[vertex].listIterator();
				while (iter.hasNext()) {
					int adjVertex = iter.next();
					if (!visited[adjVertex] && hasCycle(adjVertex, visited, onStack)) return true;
				}
			}
			// Reset the recursion stack array
			onStack[vertex] = false;
			return false;
		}

		@Override
		public boolean detectCycleInUG() {
			boolean[] visited = new boolean[N];
			for (int i = 0; i < N; i++) {
				if (!visited[i]) {
					if (hasCycleInUndirectedGraph(i, visited, -1)) return true;
				}
			}
			return false;
		}

		// using DFS
		private boolean hasCycleInUndirectedGraph(int vertex, boolean[] visited, int parent) {
			visited[vertex] = true;
			ListIterator<Integer> iter = adjList[vertex].listIterator();
			while (iter.hasNext()) {
				int adjVertex = iter.next();
				if (!visited[adjVertex]) {
					if (hasCycleInUndirectedGraph(adjVertex, visited, vertex)) return true;
				} else if (adjVertex != parent) return true;
			}
			return false;
		}

		// Using DisjointSet: Union-Find Algorithm can be used to check whether an undirected graph contains cycle or no
		public boolean hasCycleInUndirectedGraph(EdgeNode[] edges, int n, int e) {
			DisjointSet ds = new DisjointSet(n);
			ds.initialize(n);
			for (int i = 0; i < e; i++) {
				if (ds.union(edges[i].src, edges[i].dest)) return true;
			}
			return false;
		}

		@Override
		public void mstPrimsAlg() {
			int[] edgeWeight = new int[N]; // Edge Weight from source vertex
			boolean[] visited = new boolean[N];
			ResultSet[] resultSet = new ResultSet[N];

			// Set max values in edgeWeight array
			Arrays.fill(edgeWeight, Integer.MAX_VALUE);

			edgeWeight[0] = 0;

			java.util.PriorityQueue<GraphNode> queue = new java.util.PriorityQueue<>((a, b) -> a.weight - b.weight);
			queue.add(new GraphNode(0, 0));
			resultSet[0] = new ResultSet();
			resultSet[0].parent = -1;

			while (!queue.isEmpty()) {
				GraphNode src = queue.poll();
				visited[src.vertex] = true;
				Iterator<GraphNode> iter = adjListW[src.vertex].iterator();
				if (iter == null) continue;
				while (iter.hasNext()) {
					GraphNode adjNode = iter.next();
					if (!visited[adjNode.vertex]) {
						edgeWeight[adjNode.vertex] = Math.min(edgeWeight[adjNode.vertex],
								edgeWeight[src.vertex] + adjNode.weight);
						queue.add(new GraphNode(adjNode.vertex, edgeWeight[adjNode.vertex]));
					}
				}
			}

			GraphUtil.printSP(edgeWeight, N);

		}

		@Override
		public void mstKruskalsAlg() {
			// TODO Auto-generated method stub

		}

		// Time Complexity: O(V+E)
		@Override
		public void spDijikstraAlg(int source) {
			int[] edgeWeight = new int[N]; // Edge Weight from source vertex
			boolean[] visited = new boolean[N];

			// Set max values in edgeWeight array
			Arrays.fill(edgeWeight, Integer.MAX_VALUE);

			edgeWeight[source] = 0;

			java.util.PriorityQueue<GraphNode> queue = new java.util.PriorityQueue<>((a, b) -> a.weight - b.weight);
			queue.add(new GraphNode(source, 0));

			while (!queue.isEmpty()) {
				GraphNode src = queue.poll();
				visited[src.vertex] = true;
				Iterator<GraphNode> iter = adjListW[src.vertex].iterator();
				if (iter == null) continue;
				while (iter.hasNext()) {
					GraphNode adjNode = iter.next();
					if (!visited[adjNode.vertex]) {
						edgeWeight[adjNode.vertex] = Math.min(edgeWeight[adjNode.vertex],
								edgeWeight[src.vertex] + adjNode.weight);
						queue.add(new GraphNode(adjNode.vertex, edgeWeight[adjNode.vertex]));
					}
				}
			}

			GraphUtil.printSP(edgeWeight, N);
		}

		@Override
		public void spBellmanFordAlg(int source) {
			// TODO Auto-generated method stub
		}

		@Override
		public void spFloydWarshallAlg() {
			// TODO Auto-generated method stub
		}

		@Override
		public void kosarajuAlg() {
			// TODO Auto-generated method stub
		}

		@Override
		public void tarjanAlg() {
			// TODO Auto-generated method stub

		}

		@Override
		public boolean isBipartite() {
			// TODO Auto-generated method stub
			return false;
		}

	}

	class GraphAdjListUsingMap implements GraphOperations {

		int N;
		/* Map: Key-> Node or Vertex; It can be Integer, Character, String or any custom object
		 * Map: Value -> Value should be Collection. Collection such as List, Set(To avoid duplicate),
		 * 	PriorityQueue(To maintain any order), Head Node of Linked List(Sample: cloneGraph prob)
		 */
		// Map<Integer, List<Integer>> adjMap;
		Map<Integer, Set<Integer>> adjMap;
		// Map<Integer, PriorityQueue<Integer>> adjMap;

		Map<Integer, List<Integer>> adjMapW;

		List<Integer> vertices;

		@Override
		public void buildDirectedGraph(int[][] edges) {
			// Create the instance for each node
			adjMap = new HashMap<>();

			this.N = findNumberOfNodes(edges);

			// Add edges in the adjList
			for (int[] edge : edges) {
				adjMap.putIfAbsent(edge[0], new HashSet<>());
				//or if (!adjMap.containsKey(edge[0])) adjMap.put(edge[0], new HashSet<>());
				adjMap.get(edge[0]).add(edge[1]);
			}
		}

		@Override
		public void buildUndirectedGraph(int[][] edges) {
			// TODO Auto-generated method stub

		}

		@Override
		public void buildWeightedDG(int[][] edges) {
			// TODO Auto-generated method stub

		}

		@Override
		public void buildWeightedUG(int[][] edges) {
			// TODO Auto-generated method stub

		}

		@Override
		public int findNumberOfNodes(int[][] edges) {
			Set<Integer> set = CommonUtil.findNumberOfNodes(edges);
			this.vertices = new ArrayList<Integer>(set);
			return set.size();
		}

		@Override
		public void printGraph() {
			for (int v : vertices) {
				System.out.println("\nEdges from Vertex: " + v + "->");
				if (adjMap.get(v) != null) adjMap.get(v).forEach(k -> System.out.print(k + " "));
			}
		}

		@Override
		public List<Integer> dfs() {
			Set<Integer> visited = new HashSet<>();
			List<Integer> result = new ArrayList<>();
			for (int v : vertices) {
				if (!visited.contains(v)) {
					dfs(v, visited, result);
				}
			}
			return result;
		}

		// Recursive Approach
		private void dfs(int v, Set<Integer> visited, List<Integer> result) {
			visited.add(v);
			result.add(v);
			if (adjMap.get(v) == null) return;
			for (int adjNode : adjMap.get(v)) {
				if (!visited.contains(adjNode)) dfs(adjNode, visited, result);
			}
		}

		@Override
		public List<Integer> dfsIterative() {
			Set<Integer> visited = new HashSet<>();
			List<Integer> result = new ArrayList<>();
			for (int v : vertices) {
				if (!visited.contains(v)) {
					dfsIterative(v, visited, result);
				}
			}
			return result;
		}

		private void dfsIterative(int s, Set<Integer> visited, List<Integer> result) {
			Stack<Integer> stack = new Stack<>();
			visited.add(s);
			stack.push(s);
			while (!stack.isEmpty()) {
				int v = stack.pop();
				result.add(v);
				if (adjMap.get(v) == null) continue;
				for (int adjNode : adjMap.get(v)) {
					if (visited.add(adjNode)) {
						stack.push(adjNode);
					}
				}
			}
		}

		@Override
		public List<Integer> bfs() {
			Set<Integer> visited = new HashSet<>();
			List<Integer> result = new ArrayList<>();
			for (int v : vertices) {
				if (!visited.contains(v)) {
					bfs(v, visited, result);
				}
			}
			return result;
		}

		public void bfs(int source, Set<Integer> visited, List<Integer> result) {
			Queue<Integer> queue = new LinkedList<>();
			queue.add(source);
			visited.add(source);
			while (!queue.isEmpty()) {
				int v = queue.poll();
				result.add(v);
				if (adjMap.get(v) == null) continue;
				for (int adjNode : adjMap.get(v)) {
					if (visited.add(adjNode)) {
						queue.add(adjNode);
					}
				}
			}
		}

		public List<List<Integer>> bfsLevelByLevel() {
			Set<Integer> visited = new HashSet<>();
			List<List<Integer>> result = new ArrayList<>();
			for (int v : vertices) {
				if (!visited.contains(v)) {
					bfsLevelByLevel(v, visited, result);
				}
			}
			return result;
		}

		public void bfsLevelByLevel(int source, Set<Integer> visited, List<List<Integer>> result) {
			Queue<Integer> queue = new LinkedList<>();
			queue.add(source);
			visited.add(source);

			while (!queue.isEmpty()) {
				int size = queue.size();
				List<Integer> curr = new ArrayList<>();
				while (size-- > 0) {
					int v = queue.poll();
					curr.add(v);
					if (adjMap.get(v) == null) continue;
					for (int adjNode : adjMap.get(v)) {
						if (visited.add(adjNode)) {
							queue.add(adjNode);
						}
					}
				}
				result.add(new ArrayList<>(curr));
			}
		}

		@Override
		public int findDisconnectedGraph() {
			// TODO Auto-generated method stub
			return 0;
		}

		@Override
		public List<Integer> topologicalSort() {
			List<Integer> result = topologicalSortDfs();
			result = topologicalSortBfs();
			return result;
		}

		public List<Integer> topologicalSortDfs() {
			//onStack maintains 'Visiting' state nodes.
			Set<Integer> visited = new HashSet<>(), onStack = new HashSet<>();
			LinkedList<Integer> result = new LinkedList<>();
			for (int i = 0; i < N; i++) {
				if (!visited.contains(i)) {
					if (topoSortUtil(i, visited, onStack, result)) return null;
				}
			}
			return result;
		}

		//This solution is DFS Algorithm + Detect Cycle Algo
		private boolean topoSortUtil(int v, Set<Integer> visited, Set<Integer> onStack, LinkedList<Integer> result) {
			// If this condition satisfies, then adjMap contains cycle
			if (onStack.contains(v)) return true;

			// Mark vertex as visited and set recursion stack
			visited.add(v);
			onStack.add(v);

			if (adjMap.get(v) != null) {
				for (int adjVertex : adjMap.get(v)) {
					if (!visited.contains(adjVertex) && topoSortUtil(adjVertex, visited, onStack, result)) return true;
				}
			}
			result.addFirst(v);
			// Reset the recursion stack/Remove the node after visited
			onStack.remove(v);
			return false;
		}

		public List<Integer> topologicalSortBfs() {
			// If any specific ordering is required for Problem, use PriorityQueue
			Queue<Integer> queue = new LinkedList<>();
			List<Integer> linearOrder = new ArrayList<>();
			int[] indegree;
			int count = 0;

			// Step-1: Compute in-degree
			indegree = GraphUtil.indegree(adjMap);

			// Step-2: Pick all the vertices with in-degree as 0 and add them into a queue
			for (int i = 0; i < N; i++)
				if (indegree[i] == 0) queue.add(i);

			// Step-3:Remove a vertex from the queue
			while (!queue.isEmpty()) {
				int vertex = queue.poll();
				linearOrder.add(vertex);
				// 1.Increment count of visited nodes by 1.
				count++;

				// 2.Decrease in-degree by 1 for all its neighboring nodes
				if (adjMap.get(vertex) == null) continue;
				for (int adjNode : adjMap.get(vertex)) {
					// 3.If in-degree of a neighboring nodes is reduced to zero, then add it to the queue.
					if (--indegree[adjNode] == 0) queue.add(adjNode);
				}
			}
			// Step-4:If count of visited nodes is equal to the number of nodes in the graph then print the topological sort
			if (count == N) {
				linearOrder.forEach(i -> System.out.print(i + " "));
			} else {
				System.out.println("Graph is not an a DAG and also it contains a cycle.");
			}
			return linearOrder;
		}

		@Override
		public boolean detectCycleInDG() {
			//onStack maintains 'Visiting' state nodes.
			Set<Integer> visited = new HashSet<>(), onStack = new HashSet<>();
			for (int i = 0; i < N; i++) {
				if (!visited.contains(i)) {
					if (hasCycle(i, visited, onStack)) return false;
				}
			}

			return true;
		}

		private boolean hasCycle(int v, Set<Integer> visited, Set<Integer> onStack) {
			// If this condition satisfies, then adjMap contains cycle
			if (onStack.contains(v)) return true;

			// Mark vertex as visited and set recursion stack
			visited.add(v);
			onStack.add(v);

			if (adjMap.get(v) != null) {
				for (int adjVertex : adjMap.get(v)) {
					if (!visited.contains(adjVertex) && hasCycle(adjVertex, visited, onStack)) return true;
				}
			}
			// Reset the recursion stack/Remove the node after visited
			onStack.remove(v);
			return false;
		}

		@Override
		public boolean detectCycleInUG() {
			// TODO Auto-generated method stub
			return false;
		}

		@Override
		public void mstPrimsAlg() {
			// TODO Auto-generated method stub

		}

		@Override
		public void mstKruskalsAlg() {
			// TODO Auto-generated method stub

		}

		@Override
		public void spDijikstraAlg(int source) {
			// TODO Auto-generated method stub

		}

		@Override
		public void spBellmanFordAlg(int source) {
			// TODO Auto-generated method stub

		}

		@Override
		public void spFloydWarshallAlg() {
			// TODO Auto-generated method stub

		}

		/*
		 *  Solution:
		 *    A bipartite graph is possible if the graph coloring is possible using two colors such that vertices in a set are colored 
		 *    with the same color.
		 *    Note that it is possible to color a cycle graph with even cycle(even no of edges) using two colors.
		 *    It is not possible to color a cycle graph with odd cycle using two colors. 
		 *  There are two approaches to check whether graph can be colored in 2 colors:
		 *  	1. BFS
		 *  	2. DFS
		 *  Time: O(V+E) => O(V+E); Space: O(V+E)
		 */
		public boolean isBipartite() {
			Map<Integer, Integer> vertexColors = new HashMap<>();
			//BFS Approach
			for (int v : vertices) { // Iteration is to handle the disconnected graph
				if (!vertexColors.containsKey(v)) {
					if (!bfsBipartite(v, vertexColors)) return false;
				}
			}

			//DFS Approach
			for (int v : vertices) {
				if (!vertexColors.containsKey(v)) {
					//Two Colors: 0 and 1; Assign 0 to first node
					if (!dfsBipartite(v, vertexColors, 0)) return false;
				}
			}
			return true;
		}

		public boolean bfsBipartite(int v, Map<Integer, Integer> vertexColors) {
			Queue<Integer> queue = new LinkedList<>();

			queue.add(v);
			//Two Colors: 0 and 1; Assign 0 to first node
			vertexColors.put(v, 0);

			while (!queue.isEmpty()) {
				int curr = queue.poll();
				int col = vertexColors.get(curr);

				if (adjMap.get(curr) == null) continue;

				for (int adjNode : adjMap.get(curr)) {
					if (!vertexColors.containsKey(adjNode)) {
						vertexColors.put(adjNode, col == 1 ? 0 : 1); //Assign opposite color from curr/parent node
						queue.add(adjNode);
					} else if (vertexColors.get(adjNode).equals(col)) {
						//Both current and adjacent nodes are same color. So its not Bipartite Graph
						return false;
					}
				}
			}

			return true;
		}

		public boolean dfsBipartite(int v, Map<Integer, Integer> vertexColors, int color) {
			vertexColors.put(v, color);

			if (adjMap.get(v) != null) {
				for (int adjNode : adjMap.get(v)) {
					if (!vertexColors.containsKey(adjNode)) {
						//adjNode should be opposite of current node color, so opp color is passed below
						if (!dfsBipartite(adjNode, vertexColors, color == 1 ? 0 : 1)) return false;
					} else if (vertexColors.get(adjNode).equals(color)) {
						//Both current and adjacent nodes are same color. So its not Bipartite Graph
						return false;
					}
				}
			}

			return true;
		}

		/*
		 * Imp terminologies to understand the Tarjan's Alg:
		 * 	- Visited or Discovery Time: This is the time when a node is visited 1st time while DFS traversal.
		 *  - lowTime: The vertex with minimum discovery or visited time(earliest visited vertex).“Low” value of a node tells the 
		 *    topmost reachable ancestor (with minimum possible Disc value) via the subtree of that node.
		 *  - Tree Edge: It is an edge which is present in the DFS tree obtained after applying DFS on the graph.
		 *  - Back Edge: Back edges point from a node to one of its ancestors in the DFS tree.
		 *  - Forward Edge: Forward edges point from a node to one of its descendants.
		 *  - Cross Edge: Cross edges point from a node to a previously visited node that is neither an ancestor nor a descendant.
		 * Two cases in Tarjan's Algorithm:
		 *  	1.Case1 (Tree Edge): If child v(adj Node) is not visited already, then after DFS of v is complete, then minimum of low[u]
		 *  	  and low[v] will be updated to low[u].
		 *  		  low[u] = min(low[u], low[v]);
		 *      2.Case2 (Back Edge): When child v is already visited, then minimum of low[u] and Disc[v] will be updated to low[u].
		 *      	  low[u] = min(low[u], disc[v]); 
		 *   -  Time Complexity of Tarjan's Algorithm: O(V+E), Space: O(V)
		 */

		int time, sccCount;

		@Override
		public void tarjanAlg() {
			findBridges2();

			findArticulationPoints2();

			//Strongly Connected Components
			findSCC1();
		}

		/* Find Bridges in the Graph:
		 * 	Approach1:  
		 * 		A simple approach is to one by one remove all edges and see if removal of an edge causes disconnected graph. 
		 * 		Following are steps of simple approach for connected graph.
		 * 		    - For every edge (u, v), do following
		 * 				a) Remove (u, v) from graph
		 * 				b) See if the graph remains connected (We can either use BFS or DFS)
		 * 				c) Add (u, v) back to the graph.
		 * 		Time complexity:O(E*(V+E)) for a graph represented using adjacency list.
		 * 
		 *  Approach2: 
		 *  	Using Tarjan's Algorithm
		 */
		public List<int[]> findBridges2() {
			time = 0;
			List<int[]> bridges = new ArrayList<>();

			//Note: visitedTime map also used to track the visited vertex in the graph, so we dont need to have separate hashset to track this.
			Map<Integer, Integer> visitedTime = new HashMap<>();
			Map<Integer, Integer> lowTime = new HashMap<>();

			for (int vertex : adjMap.keySet()) {
				if (!visitedTime.containsKey(vertex)) {
					dfsBridge(vertex, -1, visitedTime, lowTime, bridges);
				}
			}
			return bridges;
		}

		private void dfsBridge(Integer vertex, int parent, Map<Integer, Integer> visitedTime,
				Map<Integer, Integer> lowTime, List<int[]> bridges) {
			//Directly add the current time in both maps
			visitedTime.put(vertex, time);
			lowTime.put(vertex, time);
			time++;

			if (adjMap.get(vertex) != null) {
				for (int adj : adjMap.get(vertex)) {
					/* This condition is required to skip immediate parent node in undirected graph.
					 * In undirected graph, there is a chance to directly returns to immediate parent node.
					 * This should be avoided in this algorithm, otherwise we cant find the bridge in the graph. 
					 */
					if (adj == parent) continue;

					if (!visitedTime.containsKey(adj)) {
						dfsBridge(adj, vertex, visitedTime, lowTime, bridges);

						//Case1: TreeEdge
						lowTime.put(vertex, Math.min(lowTime.get(vertex), lowTime.get(adj)));

						/* Condition to find the bridge in the graph:
						 *  If this condition satisfies, then there should not be back edge between curr vertex and adj vertex, 
						 *  so if this edge is removed, then graph will be disconnected and increases the no of disconnected graphs.
						 */
						if (visitedTime.get(vertex) < lowTime.get(adj)) {
							bridges.add(new int[] { vertex, adj });
						}
					} else {
						// Case2: Backedge
						lowTime.put(vertex, Math.min(lowTime.get(vertex), visitedTime.get(adj)));
					}
				}
			}

		}

		/* Find Articulation Points in the Graph:
		 * 	Approach1:  
		 * 		A simple approach is to one by one remove all vertices and see if removal of a vertex causes disconnected graph. 
		 * 		Following are steps of simple approach for connected graph.
		 * 		    - For every vertex v, do following 
		 * 				a) Remove v from graph
		 * 				b) See if the graph remains connected (We can either use BFS or DFS)
		 * 				c) Add v back to the graph.
		 * 		Time complexity:O(E*(V+E)) for a graph represented using adjacency list.
		 * 
		 *  Approach2: 
		 *  	Using Tarjan's Algorithm
		 */
		public List<Integer> findArticulationPoints2() {
			time = 0;
			List<Integer> articulationPoints = new ArrayList<>();
			Map<Integer, Integer> visitedTime = new HashMap<>();
			Map<Integer, Integer> lowTime = new HashMap<>();
			Map<Integer, Integer> parent = new HashMap<>();

			for (int vertex : adjMap.keySet()) {
				if (!visitedTime.containsKey(vertex)) {
					dfsAp(vertex, parent, visitedTime, lowTime, articulationPoints);
				}
			}
			return articulationPoints;
		}

		private void dfsAp(Integer vertex, Map<Integer, Integer> parent, Map<Integer, Integer> visitedTime,
				Map<Integer, Integer> lowTime, List<Integer> articulationPoints) {
			visitedTime.put(vertex, time);
			lowTime.put(vertex, time);
			time++;
			int childCount = 0;
			for (Integer adj : adjMap.get(vertex)) {
				//This condition is required to skip immediate parent node in undirected graph.
				if (adj.equals(parent.get(vertex))) continue;

				//if adj has not been visited then visit it.
				if (!visitedTime.containsKey(adj)) {
					parent.put(adj, vertex);
					childCount++;
					dfsAp(adj, parent, visitedTime, lowTime, articulationPoints);

					//Case1: TreeEdge
					lowTime.put(vertex, Math.min(lowTime.get(vertex), lowTime.get(adj)));

					// Conditions to find the Articulation Points in the graph:
					// 1. Current vertex is root of DFS tree and has two or more independent children.
					if (parent.get(vertex) == null && childCount >= 2) {
						articulationPoints.add(vertex);
					}
					/* 2.If Current vertex is not root and satisfies below condition, then adj Vertex should not have
					 *   back edge to one of the ancestors (in DFS tree) of current Vertex.
					 *   Note: Condition2 is same as bridge problem above.
					 */
					if (parent.get(vertex) != null && visitedTime.get(vertex) <= lowTime.get(adj)) {
						articulationPoints.add(vertex);
					}

				} else {
					//Case2: BackEdge
					lowTime.put(vertex, Math.min(lowTime.get(vertex), visitedTime.get(adj)));
				}
			}
		}

		/* Find Strongly Connected Components(SCC) in a graph:
		 * A directed graph is strongly connected if there is a path between all pairs of vertices. A graph is said to be strongly connected if every
		 * vertex is reachable from every other vertex. The strongly connected components of an arbitrary directed graph form a partition into subgraphs 
		 * that are themselves strongly connected.
		 * These two algorithms are used to find the SCC of a directed graph:
		 *   - Approach1: Using Kosaraju's Algorithm
		 * 	 - Apporach2: Using Tarjan's Algorithm
		 */
		public int findSCC2() {
			Stack<Integer> stack = new Stack<>();
			Set<Integer> onStack = new HashSet<>();
			//Discovered or Visited Time:
			Map<Integer, Integer> visitedTime = new HashMap<>();
			//Low Time:
			Map<Integer, Integer> lowTime = new HashMap<>();
			time = 0;
			sccCount = 0;

			for (int vertex : adjMap.keySet()) {
				if (!visitedTime.containsKey(vertex)) {
					dfsScc(vertex, visitedTime, lowTime, onStack, stack);
				}
			}
			lowTime.forEach((k, v) -> System.out.println(k + " - " + v));
			return sccCount;
		}

		private void dfsScc(int currVertex, Map<Integer, Integer> visitedTime, Map<Integer, Integer> lowTime,
				Set<Integer> onStack, Stack<Integer> stack) {
			stack.push(currVertex);
			onStack.add(currVertex);
			visitedTime.put(currVertex, time);
			lowTime.put(currVertex, time);
			time++;

			if (adjMap.get(currVertex) != null) {
				for (int adjVertex : adjMap.get(currVertex)) {
					if (!visitedTime.containsKey(adjVertex)) {
						dfsScc(adjVertex, visitedTime, lowTime, onStack, stack);
						//Case1: TreeEdge
						lowTime.put(currVertex, Math.min(lowTime.get(currVertex), lowTime.get(adjVertex)));
					} else if (onStack.contains(adjVertex)) {
						//Note: Here we have additional condition to check the vertex in stack, before updating the lowtime
						//Case2: BackEdge
						lowTime.put(currVertex, Math.min(lowTime.get(currVertex), visitedTime.get(adjVertex)));
					}
				}
			}

			/* Condition to find SCC:
			 * On recursive callback, if we're at the root node (start of SCC). Empty the seen stack until back to root.
			 */
			if (visitedTime.get(currVertex).equals(lowTime.get(currVertex))) {
				while (!stack.isEmpty()) {
					int top = stack.pop();
					onStack.remove(top);
					if (top == currVertex) break;
				}
				sccCount++;
			}
		}

		//TODO: Rewrite this

		@Override
		public void kosarajuAlg() {
			// TODO Auto-generated method stub
		}

		public int findSCC1() {
			/*//it holds vertices by finish time in reverse order.
			Stack<Integer> stack = new Stack<>();
			//holds visited vertices for DFS.
			Set<Integer> visited = new HashSet<>();
			
			//populate stack with vertices with vertex finishing last at the top.
			for (Integer vertex : adjMap.keySet()) {
				if (visited.contains(vertex)) continue;
				DFSUtil(vertex, visited, stack);
			}
			
			//reverse the graph.
			Map<Integer, Set<Integer>> reverseMap = reverseGraph();
			
			//Do a DFS based off vertex finish time in decreasing order on reverse graph..
			visited.clear();
			List<Set<Integer>> result = new ArrayList<>();
			while (!stack.isEmpty()) {
				Integer vertex = reverseGraph.getVertex(stack.poll().getId());
				if (visited.contains(vertex)) {
					continue;
				}
				Set<Integer> set = new HashSet<>();
				DFSUtilForReverseGraph(vertex, visited, set);
				result.add(set);
			}
			return result;
			*/
			return 0;
		}

		private Map<Integer, Set<Integer>> reverseGraph() {
			Map<Integer, Set<Integer>> reverseGraph = new HashMap<>();
			/*
			for (Edge<Integer> edge : graph.getAllEdges()) {
			reverseGraph.addEdge(edge.getVertex2().getId(), edge.getVertex1().getId(), edge.getWeight());
			}
			*/
			return reverseGraph;
		}

		private void DFSUtil(Integer vertex, Set<Integer> visited, Stack<Integer> stack) {
			/*
			visited.add(vertex);
			for (Integer v : vertex.getAdjacentVertexes()) {
			if (visited.contains(v)) {
			continue;
			}
			DFSUtil(v, visited, stack);
			}
			stack.push(vertex);
			*/
		}

		private void DFSUtilForReverseGraph(Integer vertex, Set<Integer> visited, Set<Integer> set) {
			/*
			visited.add(vertex);
			set.add(vertex);
			for (Integer v : vertex.getAdjacentVertexes()) {
			if (visited.contains(v)) {
			continue;
			}
			DFSUtilForReverseGraph(v, visited, set);
			}
			*/
		}

	}

	/**
	 * Edge List representation: Containers for vertices and edges. Vertices contain information only
	 * about the vertex.
	 */
	class GraphEdgeList {
		// For Edge representation
		public int N;
		public int noOfEdges;
		public EdgeNode[] edges;
		public List<Integer> vertices;

		public void buildIncidence(int V, int E, int[][] input) {
			this.N = V;
			this.noOfEdges = E;
			this.edges = new EdgeNode[E];
			for (int i = 0; i < E; i++)
				edges[i] = new EdgeNode(input[i][0], input[i][1], input[i][2]); // src, dest & weight

		}

		public void mstKruskalsAlg() {
			ArrayList<EdgeNode> result = new ArrayList<>();
			Set<Integer> nodes = CommonUtil.findNumberOfNodes(edges);
			int n = nodes.size();
			DisjointSetUsingMap set = new DisjointSetUsingMap();

			for (int i = 0; i < n; i++)
				set.createNode(i);

			Arrays.sort(edges, (u, v) -> u.weight - v.weight);

			for (EdgeNode edge : edges) {
				if (set.unionByRank(edge.src, edge.dest)) {
					result.add(edge);
				}
			}

			for (EdgeNode edge : result)
				System.out.println(edge.src + " - " + edge.dest + "->" + edge.weight);
		}

		//TODO: Revisit and remove the duplicate kruskal algorithm
		public void KruskalsMST(EdgeNode[] edges, int n) {
			ArrayList<EdgeNode> result = new ArrayList<>();
			DisjointSet ds = new DisjointSet(n);

			ds.initialize(n);

			Arrays.sort(edges, (u, v) -> u.weight - v.weight);

			for (EdgeNode edge : edges) {
				if (!ds.union(edge.src, edge.dest)) {
					result.add(edge);
				}
			}

			for (EdgeNode edge : result)
				System.out.println(edge.src + " - " + edge.dest + "->" + edge.weight);
		}

		public void spBellmanFordAlg(EdgeNode[] edges, int n, int e, int source) {
			int[] edgeWeight = new int[n]; // Edge Weight from source vertex
			// int[] parent = new int[n];

			for (int i = 0; i < n; i++)
				edgeWeight[i] = Integer.MAX_VALUE;

			edgeWeight[source] = 0;

			for (int i = 1; i < n; i++) { // outer loop runs |v| ï¿½ 1 times.
				for (int j = 0; j < e; j++) {
					int u = edges[j].src;
					int v = edges[j].dest;
					int w = edges[j].weight;
					if (edgeWeight[u] != Integer.MAX_VALUE) {
						edgeWeight[v] = Math.min(edgeWeight[v], edgeWeight[u] + w);
					}

					// Below logic is used to find the path from src to dest
					/*if (edgeWeight[u] != Integer.MAX_VALUE && edgeWeight[u] + w < edgeWeight[v]) {
						edgeWeight[v] = w + edgeWeight[u]; 
						parent[v] = u;
					}*/
				}
				System.out.println(Arrays.toString(edgeWeight));
			}

			// Try one more time same process, If graph has cycle, then edges weight keep decreases.
			for (int j = 0; j < e; j++) {
				int u = edges[j].src;
				int v = edges[j].dest;
				int w = edges[j].weight;
				if (edgeWeight[u] != Integer.MAX_VALUE && edgeWeight[u] + w < edgeWeight[v])
					System.out.println("Graph contains negative weight cycle");
			}

			printSP(edgeWeight, n);

		}

		private void printSP(int[] weight, int n) {
			System.out.println("Shortest path from source vertex to all the vertex");
			System.out.println(" Vertex " + "  Distance/Weight");
			for (int i = 0; i < n; i++) {
				System.out.println("    " + i + "   -    " + weight[i]);
			}
		}

		public void display() {
			for (EdgeNode edge : edges)
				System.out.println("Src: " + edge.src + " Dest: " + edge.dest + " Weight: " + edge.weight);
		}

	}

	static class MockData {
		public static int[][] mockMatrix1() {
			int[][] matrix = new int[][] { { 0, 2, 0, 6, 0 }, { 2, 0, 3, 8, 5 }, { 0, 3, 0, 0, 7 }, { 6, 8, 0, 0, 9 },
					{ 0, 5, 7, 9, 0 }, };

			return matrix;
		}

		public static int[][] mockMatrix2() {
			int[][] matrix = new int[][] { { 0, 4, 0, 0, 0, 0, 0, 8, 0 }, { 4, 0, 8, 0, 0, 0, 0, 11, 0 },
					{ 0, 8, 0, 7, 0, 4, 0, 0, 2 }, { 0, 0, 7, 0, 9, 14, 0, 0, 0 }, { 0, 0, 0, 9, 0, 10, 0, 0, 0 },
					{ 0, 0, 4, 14, 10, 0, 2, 0, 0 }, { 0, 0, 0, 0, 0, 2, 0, 1, 6 }, { 8, 11, 0, 0, 0, 0, 1, 7, 0 },
					{ 0, 0, 2, 0, 0, 0, 6, 7, 0 }, };
			return matrix;
		}

		public static int[][] mockMatrix3() { // Directed Graph
			int[][] matrix = new int[][] { { 0, 4, 0, 5, 2, 0 }, { 0, 0, 6, 0, 0, 4 }, { 0, 0, 0, 3, 0, 0 },
					{ 0, 0, 0, 0, 0, 1 }, { 0, 0, 0, 0, 0, 3 }, { 0, 0, 0, 0, 0, 0 }, };
			return matrix;
		}

		public static int[][] mockAdjMatrixData4() {
			int INF = Integer.MAX_VALUE;
			int[][] matrix = new int[][] { { 0, 5, INF, 10 }, { INF, 0, 3, INF }, { INF, INF, 0, 1 },
					{ INF, INF, INF, 0 } };
			return matrix;
		}

		// Mock adjacent list with sample data
		public static int[][] mockEdges1() {
			int[][] edges = { { 0, 1 }, { 0, 2 }, { 1, 2 }, { 2, 0 }, { 2, 3 }, { 3, 3 } };
			return edges;
		}

		public static int[][] mockEdges2() {
			int[][] edges = { { 1, 0 }, { 0, 2 }, { 2, 1 }, { 0, 3 }, { 1, 4 } };
			return edges;
		}

		public static int[][] mockEdges3() {
			int[][] edges = { { 5, 2 }, { 5, 0 }, { 4, 0 }, { 4, 1 }, { 2, 3 }, { 3, 1 } };
			return edges;
		}

		public static int[][] mockEdges4() {
			int[][] edges = { { 1, 0 }, { 0, 2 }, { 2, 1 }, { 0, 3 }, { 1, 4 } };
			return edges;
		}

		public static int[][] mockEdges5() {
			int[][] edges = { { 0, 1 }, { 1, 2 }, { 2, 0 }, { 5, 0 }, { 5, 6 }, { 6, 0 }, { 6, 2 }, { 6, 4 }, { 4, 5 },
					{ 3, 4 }, { 7, 5 }, { 3, 7 }, { 7, 3 } };
			return edges;
		}

		public static int[][] mockWeightedEdge1() {
			int[][] edges = { { 0, 1, -1 }, { 0, 2, 4 }, { 1, 2, 3 }, { 1, 3, 2 }, { 1, 4, 2 }, { 3, 2, 5 },
					{ 3, 1, 1 }, { 4, 3, -3 } };
			return edges;
		}
	}

	static class CommonUtil {

		public static Set<Integer> findNumberOfNodes(int[][] edges) {
			HashSet<Integer> set = new HashSet<>();
			for (int[] edge : edges) {
				if (!set.contains(edge[0])) set.add(edge[0]);
				if (!set.contains(edge[1])) set.add(edge[1]);
			}
			return set;
		}

		public static Set<Integer> findNumberOfNodes(EdgeNode[] edges) {
			HashSet<Integer> set = new HashSet<>();
			for (EdgeNode edge : edges) {
				if (!set.contains(edge.src)) set.add(edge.src);
				if (!set.contains(edge.dest)) set.add(edge.dest);
			}
			return set;
		}

		public static void testGraph(GraphOperations graph) {
			Scanner in = new Scanner(System.in);
			char ch;
			int input;
			do {
				System.out.println("Graph Operations:");
				System.out.println("1.Build Graph:");
				System.out.println("2.BFS:");
				System.out.println("3.DFS:");
				System.out.println("4.Topological Sorting:");
				System.out.println("5.Detect Cycle:");
				System.out.println("6.MST Algorithm: Prim's, Kruskal's Alg:");
				System.out.println("7.SP Algorithm: Dijikstras, Bellmanford, Floydwarshall:");
				System.out.println("8.SCC Algorithm: Kosaraju, Tarjan:");
				System.out.print("Enter option:");
				input = in.nextInt();
				switch (input) {
				case 1:
					graph.buildDirectedGraph(MockData.mockEdges5());
					// graph.buildUndirectedGraph(MockData.mockEdges1());
					// graph.buildWeightedDG(MockData.mockEdges1());
					// graph.buildWeightedUG(MockData.mockEdges1());
					System.out.println("\nDisplay:");
					graph.printGraph();
					break;
				case 2:
					System.out.println("Enter starting index:");
					graph.bfs();
					graph.findDisconnectedGraph();
					break;
				case 3:
					System.out.println("Enter starting index:");
					graph.dfs();
					break;
				case 4:
					graph.topologicalSort();
					break;
				case 5:
					graph.detectCycleInDG();
					graph.detectCycleInUG();
					break;
				case 6:
					graph.mstPrimsAlg();
					// graph.mstKruskalsAlg();
					break;
				case 7:
					System.out.println("Enter starting index:");
					graph.spDijikstraAlg(in.nextInt());
					// graph.spBellmanFordAlg(in.nextInt());
					// graph.spFloydWarshallAlg();
					break;
				case 8:
					//graph.sccKosarajuAlg();
					graph.tarjanAlg();
					break;
				default:
					System.out.println("Please enter the valid option!!!");
					break;
				}

				System.out.println("\nDo you want to continue(y/n):");
				ch = in.next().charAt(0);
			} while (ch == 'y' || ch == 'Y');
			System.out.println("****Thank You******");
			in.close();
		}

	}

	class ResultSet {
		int parent;
		int weight;
	}

	static class GraphUtil {

		private static int[] indegree(Map<Integer, Set<Integer>> adjMap) {
			int n = adjMap.size();
			int[] indegree = new int[n];

			for (Map.Entry<Integer, Set<Integer>> entry : adjMap.entrySet()) {
				if (entry.getValue() == null) continue;
				for (int adjNode : entry.getValue())
					indegree[adjNode]++;
			}

			return indegree;
		}

		private static int[] indegree(LinkedList<Integer>[] adjList, int n) {
			int[] indegree = new int[n];

			for (int i = 0; i < n; i++) {
				if (adjList[i].size() > 0) {
					ListIterator<Integer> iterator = adjList[i].listIterator();
					while (iterator.hasNext()) indegree[iterator.next()]++;
				}
			}
			return indegree;
		}

		private static int findMinWeight(int[] weight, boolean[] visited) {
			int index = -1, min = Integer.MAX_VALUE;
			for (int i = 0; i < weight.length; i++) {
				if (weight[i] < min && !visited[i]) {
					min = weight[i];
					index = i;
				}
			}
			return index;
		}

		private static void printSP(int[] weight, int n) {
			System.out.println("Shortest path from source vertex to all the vertex");
			System.out.println(" Vertex " + "  Distance/Weight");
			for (int i = 0; i < n; i++) {
				System.out.println("    " + i + "   -    " + weight[i]);
			}
		}

		private static void printMST(int[] parent, int[] weight, int n) {
			System.out.println(" Edge " + " Weight ");
			for (int i = 1; i < n; i++) {
				System.out.println(i + "-" + parent[i] + "  " + weight[i]); // adjMatrix[i][parent[i]]
			}
		}

		private static void printfloydWarshallSP(int dist[][], int n) {
			System.out.println("Following matrix shows the shortest distances between every pair of vertices");
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < n; ++j) {
					if (dist[i][j] == Integer.MAX_VALUE) System.out.print("INF ");
					else System.out.print(dist[i][j] + "   ");
				}
				System.out.println();
			}
		}

	}
}