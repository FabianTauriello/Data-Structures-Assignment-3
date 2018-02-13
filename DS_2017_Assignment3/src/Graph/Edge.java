package Graph;

import java.io.Serializable;
import java.util.Map;

/**
 * Represents a path, or connection, between two nodes.
 * Created for Data Structures, SP2 2017
 * @author James Baumeister
 * @version 1.1
 */
public class Edge implements Serializable {
	
	private Node node1;
	private Node node2;
	private double weight;
	
	/**
	 * @param node1
	 * @param node2
	 */
	public Edge(Node node1, Node node2) {
		this.node1 = node1;
		this.node2 = node2;
		calculateWeight();
	}
	
	public Node getNode1() {
		return node1;
	}
	
	public void setNode1(Node node1) {
		this.node1 = node1;
	}
	
	public Node getNode2() {
		return node2;
	}
	
	public void setNode2(Node node2) {
		this.node2 = node2;
	}
	
	/**
	 * Calculates the weight of travel from one node to another,
	 * taking into account each node's biome.
	 * The weight should be calculated as:
	 * w(e) = d(p, q) * |(b1 + b2) / 2|
	 * where b1 is the biome value of the source node, and
	 * b2 is the biome value of the target node. d is a
	 * function that calculates the euclidean distance
	 * between two 2D points, p and q.
	 */
	private void calculateWeight() {
		
		// calculate distance between nodes
		double distance = node1.getPosition().distance(node2.getPosition());
		// calculate average weight between nodes 
		double averageWeight = (node1.getBiomeWeight() + node2.getBiomeWeight()) / 2;
		// set weight for this Edge
		weight = distance * averageWeight;		
	}
	
	public double getWeight() {
		return weight;
	}
	
	@Override
	public String toString() {
		return "Edge{" +
				"node1=" + node1 +
				", node2=" + node2 +
				", weight=" + weight +
				'}';
	}
}
