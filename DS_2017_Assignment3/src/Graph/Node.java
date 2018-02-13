package Graph;

import java.io.Serializable;
import java.util.ArrayList;

/**
 * A node represents one polygon on a terrain map. It has a number of features
 * related to its biome, elevation, position and connections to other nodes.
 * Many of characteristics may not be used.
 * DO NOT MODIFY THIS CLASS. A clean version will be used for marking.
 * Created for Data Structures, SP2 2017
 * @author James Baumeister
 * @version 1.1
 */
public class Node implements Serializable {
	
	private static final long serialVersionUID = -3085575079765859612L;
	
	public enum BIOME {
		BARE(1),
		BEACH(1.1),
		GRASSLAND(1.5),
		LAKE(5),
		OCEAN(10),
		SCORCHED(3),
		SHRUBLAND(2.5),
		SNOW(5),
		SUBTROPICAL_DESERT(4),
		TAIGA(4),
		TEMPERATE_DECIDUOUS_FOREST(5),
		TEMPERATE_DESERT(4.5),
		TEMPERATE_RAIN_FOREST(6),
		TROPICAL_RAIN_FOREST(7),
		TROPICAL_SEASONAL_FOREST(3),
		TUNDRA(2);
		
		private final double weight;
		
		BIOME(double weight) {
			this.weight = weight;
		}
		
		public double getWeight() {
			return weight;
		}
	}
	
	private BIOME biome;
	private boolean coast;
	private double moisture;
	private boolean water;
	private boolean ocean;
	private boolean border;
	private Vector2 position;
	private ArrayList<Edge> edges;
	
	public Node(BIOME biome, boolean coast, double moisture, boolean ocean, boolean water,
	            boolean border, Vector2 position) {
		this.biome = biome;
		this.coast = coast;
		this.moisture = moisture;
		this.ocean = ocean;
		this.water = water;
		this.border = border;
		this.position = position;
		edges = new ArrayList<Edge>();
	}
	
	public BIOME getBiome() {
		return biome;
	}
	
	public void setBiome(BIOME biome) {
		this.biome = biome;
	}
	
	public double getBiomeWeight() { return biome.getWeight(); }
	
	public boolean isCoast() {
		return coast;
	}
	
	public void setCoast(boolean coast) {
		this.coast = coast;
	}
	
	public boolean isBorder() {
		return border;
	}
	
	public void setBorder(boolean border) {
		this.border = border;
	}
	
	public double getMoisture() {
		return moisture;
	}
	
	public void setMoisture(double moisture) {
		this.moisture = moisture;
	}
	
	public boolean isWater() {
		return water;
	}
	
	public void setWater(boolean water) {
		this.water = water;
	}
	
	public boolean isOcean() {
		return ocean;
	}
	
	public void setOcean(boolean ocean) {
		this.ocean = ocean;
	}
	
	public Vector2 getPosition() {
		return position;
	}
	
	public void setPosition(Vector2 position) {
		this.position = position;
	}
	
	public ArrayList<Edge> getEdges() {
		return edges;
	}
	
	public void setEdges(ArrayList<Edge> edges) {
		this.edges = edges;
	}
	
	@Override
	public String toString() {
		return "Node{" +
				"biome=" + biome +
				", coast=" + coast +
				", moisture=" + moisture +
				", water=" + water +
				", ocean=" + ocean +
				", border=" + border +
				", position=" + position +
				'}';
	}
	
	@Override
	public boolean equals(Object o) {
		if (this == o) {
			return true;
		}
		if (o == null || getClass() != o.getClass()) {
			return false;
		}
		
		Node node = (Node) o;
		
		if (coast != node.coast) {
			return false;
		}
		if (Double.compare(node.moisture, moisture) != 0) {
			return false;
		}
		if (water != node.water) {
			return false;
		}
		if (ocean != node.ocean) {
			return false;
		}
		if (border != node.border) {
			return false;
		}
		if (biome != node.biome) {
			return false;
		}
		if (position != null ? !position.equals(node.position) : node.position != null) {
			return false;
		}
		if (edges == null || ((Node) o).getEdges() == null)
			return false;
		
		for (Edge e : edges) {
			if (!((Node) o).getEdges().contains(e)) {
				return false;
			}
		}
		return true;
	}
	
	@Override
	public int hashCode() {
		int result;
		long temp;
		result = biome != null ? biome.hashCode() : 0;
		result = 31 * result + (coast ? 1 : 0);
		temp = Double.doubleToLongBits(moisture);
		result = 31 * result + (int) (temp ^ (temp >>> 32));
		result = 31 * result + (water ? 1 : 0);
		result = 31 * result + (ocean ? 1 : 0);
		result = 31 * result + (border ? 1 : 0);
		result = 31 * result + (position != null ? position.hashCode() : 0);
		return result;
	}
}
