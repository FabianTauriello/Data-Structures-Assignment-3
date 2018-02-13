package Testing;

import Graph.Graph;
import Graph.Node;
import Graph.Edge;
import GUI.GraphLoader;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.ArrayList;

import static junit.framework.TestCase.fail;

/**
 * @author James Baumeister on 4/5/17
 * Most tests assume that connectNodes is working.
 * Complete that method first.
 */
public class GraphTest extends DSUnitTesting {
	
	Graph g;
	Node[] nodes;
	
	@Before
	public void initialise() {
		GraphLoader gl = new GraphLoader();
		nodes = gl.getNodes();
		g = new Graph();
	}
	
	@Test
	public void connectNodes() {
		AssignmentMarker.marks.put("Graph:connectNodes", 10.0f);
		
		for (Node n : nodes) {
			if (n.getEdges().size() != 0) {
				fail("Nodes should not be connected at this point");
			}
		}
		
		g.connectNodes(nodes);
		
		for (Node n : nodes) {
			if (n.getEdges().size() <= 0) {
				fail("Nodes should be connected at this point");
			}
		}
		
		// Test a couple of cases
		Assert.assertEquals("Node 0 should have three edges", 3, nodes[0].getEdges().size());
		Assert.assertEquals("Test Node 0's east edge", nodes[1], nodes[0].getEdges().get(0).getNode2());
		Assert.assertEquals("Test Node 0's south edge", nodes[10], nodes[0].getEdges().get(1).getNode2());
		Assert.assertEquals("Test Node 0's south-east edge", nodes[11], nodes[0].getEdges().get(2).getNode2());
		
		Assert.assertEquals("Node 45 should have 8 edges", 8, nodes[45].getEdges().size());
		Assert.assertEquals("Test Node 45's east edge", nodes[46], nodes[45].getEdges().get(0).getNode2());
		Assert.assertEquals("Test Node 45's west edge", nodes[44], nodes[45].getEdges().get(1).getNode2());
		Assert.assertEquals("Test Node 45's south edge", nodes[55], nodes[45].getEdges().get(2).getNode2());
		Assert.assertEquals("Test Node 45's north edge", nodes[35], nodes[45].getEdges().get(3).getNode2());
		Assert.assertEquals("Test Node 45's north-east edge", nodes[36], nodes[45].getEdges().get(4).getNode2());
		Assert.assertEquals("Test Node 45's south-east edge", nodes[56], nodes[45].getEdges().get(5).getNode2());
		Assert.assertEquals("Test Node 45's north-west edge", nodes[34], nodes[45].getEdges().get(6).getNode2());
		Assert.assertEquals("Test Node 45's south-west edge", nodes[54], nodes[45].getEdges().get(7).getNode2());
	}
	
	@Test
	public void getEdge() {
		AssignmentMarker.marks.put("Graph:getEdge", 5.0f);
		Node n1 = nodes[0];
		Node n2 = nodes[1];
		
		// Some fake edges
		Edge e1 = new Edge(nodes[4], nodes[6]);
		Edge e2 = new Edge(nodes[10], nodes[11]);
		Edge e3 = new Edge(nodes[3], nodes[2]);
		Edge e4 = new Edge(nodes[8], nodes[9]);
		
		// The actual edge
		Edge e = new Edge(nodes[0], nodes[1]);
		
		ArrayList<Edge> edges1 = new ArrayList<Edge>();
		edges1.add(e1);
		edges1.add(e2);
		edges1.add(e);
		edges1.add(e3);
		n1.setEdges(edges1);
		
		ArrayList<Edge> edges2 = new ArrayList<Edge>();
		edges2.add(e4);
		edges2.add(e);
		
		Assert.assertEquals("Test that the correct edge is found", e, g.getEdge(n1, n2));
	}
	
	@Test
	public void calculateCost() {
		AssignmentMarker.marks.put("Graph:calculateCost", 5.0f);
		g.connectNodes(nodes);
		
		Node[] path = {nodes[23], nodes[24], nodes[35], nodes[46], nodes[57]};
		
		
		try {
			Assert.assertEquals("Test that the path cost is correctly calculated", 727.85, g.calculateCost(path),
					0.01);
			
		}
		catch (NullPointerException e) {
			fail("Check that connectNodes is functioning correctly");
		}
	}
	
	@Test
	public void breadthFirstSearch() {
		AssignmentMarker.marks.put("Graph:breadthFirstSearch", 20.0f);
		g.connectNodes(nodes);

		Node[] exp_path = {nodes[0], nodes[10], nodes[20], nodes[30], nodes[40], nodes[50], nodes[60]};
		Node[] actual_path = g.breadthFirstSearch(nodes[0], nodes[60]);
		
		paths_equal(1, exp_path, actual_path);
		
		exp_path = new Node[]{nodes[71], nodes[61], nodes[51], nodes[42], nodes[33], nodes[24], nodes[15]};
		actual_path = g.breadthFirstSearch(nodes[71], nodes[15]);
		
		paths_equal(2, exp_path, actual_path);
		
		exp_path = new Node[]{nodes[3], nodes[4], nodes[5], nodes[6], nodes[17]};
		actual_path = g.breadthFirstSearch(nodes[3], nodes[17]);
//		exp_path = new Node[] {nodes[1], nodes[1]};
//		actual_path = g.breadthFirstSearch(nodes[1], nodes[1]);
		
		paths_equal(3, exp_path, actual_path);
	}
	
	private void paths_equal(int test_num, Node[] exp, Node[] act) {
		if (exp.length != act.length) {
			fail(test_num + ". Paths are of different lengths");
		}
		
		for (int i = 0; i < exp.length; i++) {
			if (exp[i] != act[i]) {
				fail(test_num + ". Actual path is different than expected");
			}
		}
	}
	
	@Test
	public void depthFirstSearch() {
		AssignmentMarker.marks.put("Graph:depthFirstSearch", 15.0f);
		g.connectNodes(nodes);
		
		Node[] exp_path = {nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6], nodes[7],
				nodes[8], nodes[9], nodes[19], nodes[18]}; 
		Node[] actual_path = g.depthFirstSearch(nodes[0], nodes[18]);
		
		paths_equal(1, exp_path, actual_path);
		
		exp_path = new Node[]{nodes[63], nodes[64], nodes[65], nodes[66], nodes[67], nodes[68], nodes[69], nodes[79],
				nodes[78], nodes[77], nodes[76], nodes[75], nodes[74], nodes[73], nodes[72], nodes[71], nodes[70],
				nodes[60], nodes[61], nodes[62], nodes[52]};
		actual_path = g.depthFirstSearch(nodes[63], nodes[52]);
//		exp_path = new Node[] {nodes[1], nodes[1]};
//		actual_path = g.depthFirstSearch(nodes[1], nodes[1]);
		
		paths_equal(2, exp_path, actual_path);
	}
	
	@Test
	public void dijkstrasSearch() {
		AssignmentMarker.marks.put("Graph:dijkstrasSearch", 20.0f);
		g.connectNodes(nodes);
		
		Node[] exp_path = {nodes[0], nodes[1], nodes[2], nodes[13], nodes[23], nodes[32], nodes[42],
				nodes[51], nodes[60]};
		Node[] actual_path = g.dijkstrasSearch(nodes[0], nodes[60]);
		
		paths_equal(1, exp_path, actual_path);
		
		exp_path = new Node[]{nodes[8], nodes[17], nodes[16], nodes[15], nodes[14], nodes[13], nodes[2],
				nodes[1], nodes[0]}; 
		actual_path = g.dijkstrasSearch(nodes[8], nodes[0]);
//		exp_path = new Node[] {nodes[1], nodes[1]};
//		actual_path = g.dijkstrasSearch(nodes[1], nodes[1]);
		
		paths_equal(2, exp_path, actual_path);
	}
	
	@Test
	public void aStarSearch() {
		AssignmentMarker.marks.put("Graph:aStarSearch", 10.0f);
		g.connectNodes(nodes);
		
		Node[] exp_path = {nodes[9], nodes[8], nodes[17], nodes[27], nodes[37], nodes[47], nodes[57],
				nodes[67], nodes[78], nodes[88], nodes[99]};
		Node[] actual_path = g.aStarSearch(nodes[9], nodes[99]);


		
		paths_equal(1, exp_path, actual_path);
	}
}
