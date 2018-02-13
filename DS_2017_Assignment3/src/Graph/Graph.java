package Graph;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;

import GUI.GraphLoader;

/**
 * Represents a non-directional graph where each vertex
 * is a Node object. Connections between nodes are based
 * on the cartesian coordinate system.
 * @author Fabian Tauriello
 */
public class Graph {

	// class variables to help with DFS
	private Map<Node, Node> parents;
	private ArrayList<Node> visitedArray; 			// holds all vertices that've been visited
	private Node [] discoveryOrder;					// the array that contains each vertex in discovery order
	private Node [] finishOrder;					// the array that contains each vertex in finish order
	private int discoverIndex = 0;					// the index that indicates the discovery order
	private int finishIndex = 0;					// the index that indicates the finish orders

	/**
	 * Connects all nodes, building their E, W, S, N, NE, SE, NW, SW edges,
	 * in that order. The nodes should form a 10x10 square grid, and the array
	 * is such that every i'th node % 10 = 9 is a right edge.
	 * See the assignment specification for more information.
	 * @param nodes An array of Node objects to be connected
	 */
	public void connectNodes(Node[] nodes) {
		if (nodes == null) {
			throw new IllegalArgumentException();
		}

		for (int index = 0; index < nodes.length; index++) {
			// determine what row and column the current Node needs to placed in
			int row = index / 10;
			int col = index % 10;

			if (row > 0 && row < 9 && col > 0 && col < 9) { // assign edges to inner nodes
				// add east edge
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index+1]));
				// add west edge
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index-1]));
				// add south edge
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index+10]));
				// add north edge
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index-10]));
				// add north-east edge
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index-9]));
				// add south-east edge
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index+11]));
				// add north-west edge
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index-11]));
				// add south-west edge
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index+9]));
			} else if (row == 0 && col > 0 && col < 9) { // assign edges to top row
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index+1]));
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index-1]));
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index+10]));
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index+11]));
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index+9]));
			} else if (row == 9 && col > 0 && col < 9) { // assign edges to bottom row
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index+1]));
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index-1]));
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index-10]));
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index-9]));
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index-11]));
			} else if (col == 0 && row > 0 && row < 9) { // assign edges to left column
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index+1]));
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index+10]));
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index-10]));
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index-9]));
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index+11]));
			} else if (col == 9 && row > 0 && row < 9) { // assign edges to right column
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index-1]));
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index+10]));
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index-10]));
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index-11]));
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index+9]));
			} else if (row == 0 && col == 0) { // assign edges to top-left corner
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index+1]));
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index+10]));
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index+11]));
			} else if (row == 0 && col == 9) { // assign edges to top-right corner
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index-1]));
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index+10]));
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index+9]));
			} else if (row == 9 && col == 0) { // assign edges to bottom-left corner
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index+1]));
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index-10]));
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index-9]));
			} else if (row == 9 && col == 9) { // assign edges to bottom-right corner
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index-1]));
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index-10]));
				nodes[index].getEdges().add(new Edge(nodes[index], nodes[index-11]));
			}	
		}

	}

	/**
	 * Performs a breadth-first search of the graph and returns the shortest
	 * path from one node to another. This method does not take into account edge weights
	 * @param start The node from which to start searching
	 * @param target The target node to which a path is built
	 * @return An array of Node objects representing the path from start to target, in that order
	 */
	public Node[] breadthFirstSearch(Node start, Node target) {
		if ((start == null || target == null) || (start.equals(target))) {
			throw new IllegalArgumentException();
		}

		// array of identified nodes
		ArrayList<Node> identifiedNodes = new ArrayList<Node>();
		// Queue used to store nodes that are waiting to be visited
		Queue<Node> toVisit = new LinkedList<Node>();
		// map for holding parents - keys must be unique (may not need to specify node again in hashmap generic)
		Map<Node, Node> parentMap = new HashMap<Node, Node>();
		// declare path list to track parents 
		List<Node> path = new ArrayList<Node>();

		// insert the start node into identified array and the queue
		identifiedNodes.add(start);
		toVisit.offer(start);

		// loop until queue is not empty
		while (toVisit.isEmpty() == false) {
			// take a vertex, 'current', out of queue. visit it
			Node current = toVisit.remove();
			// examine each vertex adjacent to 'current' and see if it's identified. 
			for (Edge e : current.getEdges()) {
				Node adjNode = e.getNode2();
				if (!identifiedNodes.contains(adjNode)) {
					identifiedNodes.add(adjNode);
					toVisit.offer(adjNode);
					parentMap.put(adjNode, current);
				}
			}
		}
		// return reversed path
		return reconstructPath(start, target, path, parentMap);
	}

	/**
	 * Performs a depth-first search of the graph and returns the first-found
	 * path from one node to another.
	 * @param start The node from which to start searching
	 * @param target The target node to which a path is built
	 * @return An array of Node objects representing the path from start to target, in that order
	 */
	public Node[] depthFirstSearch(Node start, Node target) {
		if ((start == null || target == null) || (start.equals(target))) {
			throw new IllegalArgumentException();
		}
		
		// Initialize required data structures for DFS
		parents = new HashMap<Node, Node>();
		visitedArray = new ArrayList<Node>();
		discoveryOrder = new Node[100];
		finishOrder = new Node[100];
		discoverIndex = 0;
		finishIndex = 0;

		// call recursive DFS method
		depthFirstSearch(start);

		// declare path list to track parents 
		List<Node> path = new ArrayList<Node>();

		// return reversed path
		return reconstructPath(start, target, path, parents);

	}

	/**
	 * Helper method for DFS above. This a recursive implementation
	 * @param current node by which other nodes are searched. Initially, this is the starting point.
	 */
	public void depthFirstSearch (Node current) {
		// mark the current vertices visited
		visitedArray.add(current);
		discoveryOrder[discoverIndex++] = current; 		// discoverIndex takes is zero and THEN adds 1
		// examine each vertex adjacent to the current vertex
		for (Edge e : current.getEdges()) {
			Node adjacentNode = e.getNode2();
			// process neighbours that have not been visited
			if (!visitedArray.contains(adjacentNode)) {
				// insert current adjacent node into parent map
				parents.put(adjacentNode, current);
				depthFirstSearch(adjacentNode);
			}
		}
		// insert current node into finishOrder array - mark as finished
		finishOrder[finishIndex++] = current;

	}


	/**
	 * Grabs all nodes adjacent to (and beyond) a given node
	 * @param start position to start node search from
	 * @return HashSet of nodes in graph
	 */
	public HashSet<Node> getNodes(Node start) {

		HashSet<Node> identifiedNodes = new HashSet<Node>();
		Queue<Node> toVisit = new LinkedList<Node>();	

		identifiedNodes.add(start);
		toVisit.offer(start);

		// loop until queue is not empty
		while (toVisit.isEmpty() == false) {
			// take a vertex, 'current', out of queue. visit it
			Node current = toVisit.remove();
			// examine each vertex adjacent to 'current' and see if it's identified. 
			for (Edge e : current.getEdges()) {
				Node adjNode = e.getNode2();
				if (!identifiedNodes.contains(adjNode)) {
					identifiedNodes.add(adjNode);
					toVisit.offer(adjNode);
				}
			}
		}
		return identifiedNodes;
	}


	/**
	 * Performs a search of the graph using Dijkstra's algorithm, which takes into
	 * account the edge weight. It should return the least-costly path from one node
	 * to another.
	 * @param start The node from which to start searching
	 * @param target The target node to which a path is built
	 * @return An array of Node objects representing the path from start to target, in that order
	 */
	public Node[] dijkstrasSearch(Node start, Node target) {
		if ((start == null || target == null) || (start.equals(target))) {
			throw new IllegalArgumentException();
		}

		HashSet<Node> closedSet = new HashSet<Node>();								// contains the vertices for which we have computed the shortest distance (S)
		HashSet<Node> openSet = new HashSet<Node>(getNodes(start));					// contains the vertices that we still need to process (V-S)
		HashMap<Node, Double> distancesMap = new HashMap<Node, Double>(); 			// contains the shortest distances between vertices (d[v])
		HashMap<Node, Node> parentMap = new HashMap<Node, Node>(); 					// contains the predecessors of a given node  -- Textbook uses array?
		List<Node> tempPath = new ArrayList<Node>(); 								// backwards representation of path
		
		// initialize closedSet
		closedSet.add(start);

		// adjust openSet to represent all other nodes possible in path
		openSet.remove(start);

		// initialize openSet
		for (Node current : openSet) {
			// set the parent of all nodes to the start position of the search
			parentMap.put(current, start); 
			// see if there is an edge between current edge and start position
			if (getEdge(current, start) != null) {
				// add weight of this edge to distances map
				distancesMap.put(current, getEdge(current, start).getWeight());
			} else {
				// add infinite weight if direct edge doesn't exist (don't know how far it is between nodes current and start)
				distancesMap.put(current, Double.POSITIVE_INFINITY);
			}
		}

		// main loop
		while (openSet.isEmpty() == false) {
			// find the node in openSet with the smallest distance recorded in distances map
			double minDistance = Double.POSITIVE_INFINITY;
			Node nodeToFind = null;
			for (Node current : openSet) {
				if (distancesMap.get(current) < minDistance) {
					minDistance = distancesMap.get(current);
					nodeToFind = current;
				} 
			}
			// remove nodeToFind from openSet
			openSet.remove(nodeToFind);
			// update distances map
			for (Node current : openSet) {
				if (getEdge(current, nodeToFind) != null) {
					double edgeWeight = getEdge(current, nodeToFind).getWeight();
					if (distancesMap.get(nodeToFind) + edgeWeight < distancesMap.get(current)) {
						distancesMap.put(current, distancesMap.get(nodeToFind) + edgeWeight);
						parentMap.put(current, nodeToFind);
					}
				}
			}
		}

		// return reversed path
		return reconstructPath(start, target, tempPath, parentMap);

	}

	/**
	 * Performs a search of the graph using the A* algorithm, which using a heuristic to optimize the
	 * direction of search. It should return the least-costly path from one node to another.
	 * @param start The node from which to start searching
	 * @param target The target node to which a path is built
	 * @return An array of Node objects representing the path from start to target, in that order
	 */
	public Node[] aStarSearch(Node start, Node target) {
		if ((start == null || target == null) || (start.equals(target))) {
			throw new IllegalArgumentException();
		}

		return null;		
	}
	
	/**
	 * Produces h Cost (the heuristic is an estimate of the distance between a node and the goal)
	 * @param node1 
	 * @param node2 
	 * @return the heuristic cost between two nodes
	 */
	private double getHeuristic(Node node1, Node node2) {
		return node1.getPosition().distance(node2.getPosition());
	}
	
	/**
	 * Reconstruct by reversing order and converting it to a Node array
	 * @param start node to look for 
	 * @param target node to begin search from
	 * @param originalPath path to reconstruct
	 * @param parentMap map of all parents in original path
	 * @return
	 */
	private Node[] reconstructPath(Node start, Node target, List<Node> originalPath, Map<Node, Node> parentMap) {
		
		// add parents in map to path list
		Node mover = target;
		originalPath.add(target);
		// track back from target to start node by looking at parents
		while(!parentMap.get(mover).equals(start)) {
			originalPath.add(parentMap.get(mover));
			mover = parentMap.get(mover);
		}
		// add start node to path
		originalPath.add(parentMap.get(mover));

		// reverse path
		Collections.reverse(originalPath);

		// convert path to Node array
		Node[] resultArray = new Node[originalPath.size()];
		resultArray = originalPath.toArray(resultArray);
		
		return resultArray;
	}

	/**
	 * Searches for an edge from the source node to the destination.
	 * @param source The source, or first, node
	 * @param destination The destination, or second, node
	 * @return The edge between the nodes, or null if not found
	 */
	public Edge getEdge(Node source, Node destination) {
		// grab list from source
		ArrayList<Edge> edgesToCheck = source.getEdges();
		// check each edge the source node has.  
		for (Edge e : edgesToCheck) {
			// if one of the edges is pointing to the destination node, return that edge, otherwise return null.
			if (e.getNode2() == destination) {
				return e;
			}
		}
		return null;
	}

	/**
	 * From an array of Node objects, this calculates the total cost of
	 * travelling from the first to last nodes.
	 * @param vertices An array of Node objects representing a path
	 * @return The total cost of travel
	 */
	public double calculateCost(Node[] vertices) {
		double cost = 0;
		// loop for all nodes in array
		for (int index = 0; index < vertices.length-1; index++) {
			// grab the edges for a node
			ArrayList<Edge> nEdges = vertices[index].getEdges();
			// loop for all edges in a node
			for(Edge e : nEdges) {
				// if the destination node is the next node in the given node array (called "vertices"), 
				// add the weight of the edge between nodes to the cost
				if (e.getNode2() == vertices[index+1]) {
					cost += e.getWeight();
				}
			}
		}
		return cost;
	}
}
