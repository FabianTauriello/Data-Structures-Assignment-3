package GUI;

import Graph.Node;

/**
 * @author James Baumeister on 3/11/16
 */
public class Player {
	public static Node position;
	
	public Player(Node startingPosition) {
		position = startingPosition;
	}
	
	@Override
	public String toString() {
		return "Player{" +
				"position=" + position +
				'}';
	}
}
