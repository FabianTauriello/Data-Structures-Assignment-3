package Testing;

import Graph.Vector2;
import org.junit.Assert;
import org.junit.Test;

/**
 * @author James Baumeister on 4/5/17
 */
public class VectorTest extends DSUnitTesting {
	
	@Test
	public void distance() {
		AssignmentMarker.marks.put("Vector:distance", 5.0f);
		
		Vector2 v1 = new Vector2(1, 1);
		Vector2 v2 = new Vector2(1, 4);
		
		Assert.assertEquals("1. Checking that the euclidean distance is correctly calculated", 3.0, v1.distance(v2), 0.01);
		Assert.assertEquals("2. Checking that the euclidean distance is correctly calculated", 3.0, v2.distance(v1), 0.01);
		Assert.assertEquals("3. Checking that the euclidean distance is correctly calculated", 0.0, v1.distance(v1), 0.01);
		
		v2.x = 4;
		
		Assert.assertEquals("4. Checking that the euclidean distance is correctly calculated", 4.24, v1.distance(v2), 0.01);
		
		v2.x = -4;
		
		Assert.assertEquals("5. Checking that the euclidean distance is correctly calculated", 5.83, v1.distance(v2), 0.01);
	}
}
