package Testing;

import java.util.ArrayList;
import java.util.HashMap;

import org.junit.runner.JUnitCore;
import org.junit.runner.Result;
import org.junit.runner.notification.Failure;

public class AssignmentMarker {
	public static HashMap<String, Float> marks = new HashMap<String, Float>();
	
	public static ArrayList<String> passed = new ArrayList<String>();
	public static ArrayList<String> failed = new ArrayList<String>();
	
	private static java.util.ArrayList<Failure> failures = new java.util.ArrayList<Failure>();
	
	private static void testrunner(String name, Class c) {
		
		Result test = JUnitCore.runClasses(c);
		failures.addAll(test.getFailures());
		
		for ( int i = 0 ; i < test.getFailures().size() ; ++i ) {
			String testID = test.getFailures().get(i).getDescription().getClassName() + ":" +
					test.getFailures().get(i).getDescription().getMethodName();
			testID = testID.replaceAll("Test", "").replaceAll("test", "");
		}
	}
	
	// Simple test information
	private static void runATest(String name, Class c) {
		System.out.println("\n" + name);
		for ( int i = 0 ; i < name.length() ; ++i )
			System.out.print("-");
		System.out.println();
		
		testrunner(name, c);
	}
	
	public static void main(String[] args) {
		
		failures = new java.util.ArrayList<Failure>();
		passed = new ArrayList<String>();
		failed = new ArrayList<String>();
		marks = new HashMap<String, Float>();
		
		float eP, eL, vP, vL, gP, gL;
		
		System.out.println("Data Structures Assignment #3:\n\tPathfinding\n");
		
		System.out.println("-----------------------------");
		
		runATest("Edge", EdgeTest.class);

		System.out.println("Summary: ");
		{
			float gained = 0.0f;
			for ( int i = 0 ; i < passed.size(); ++i ) {
				if ( passed.get(i).contains("Edge:") )
					gained += marks.get(passed.get(i));
			}

			float lost = 0.0f;
			for ( int i = 0 ; i < failed.size(); ++i ) {
				if ( failed.get(i).contains("Edge:") )
					lost += marks.get(failed.get(i));
			}
			eP = gained;
			eL = lost;
			System.out.println(gained + " marks gained, " + lost + " marks lost.");
		}
		
		runATest("Vector", VectorTest.class);
		
		System.out.println("Summary: ");
		{
			float gained = 0.0f;
			for ( int i = 0 ; i < passed.size(); ++i ) {
				if ( passed.get(i).contains("Vector:") )
					gained += marks.get(passed.get(i));
			}
			
			float lost = 0.0f;
			for ( int i = 0 ; i < failed.size(); ++i ) {
				if ( failed.get(i).contains("Vector:") )
					lost += marks.get(failed.get(i));
			}
			vP = gained;
			vL = lost;
			System.out.println(gained + " marks gained, " + lost + " marks lost.");
		}
		
		runATest("Graph", GraphTest.class);
		
		System.out.println("Summary: ");
		{
			float gained = 0.0f;
			for ( int i = 0 ; i < passed.size(); ++i ) {
				if ( passed.get(i).contains("Graph:") )
					gained += marks.get(passed.get(i));
			}
			
			float lost = 0.0f;
			for ( int i = 0 ; i < failed.size(); ++i ) {
				if ( failed.get(i).contains("Graph:") )
					lost += marks.get(failed.get(i));
			}
			gP = gained;
			gL = lost;
			System.out.println(gained + " marks gained, " + lost + " marks lost.");
		}
		
		System.out.println("-----------------------------");
		System.out.println("\nFailed test details");
		System.out.println("( Test class: test name -> Error details )\n");
		for (Failure failure : failures) {
			String name = failure.getDescription().getClassName().replaceAll("Testing.", "").replaceAll("Test",  "")
					+ ": " +  failure.getDescription().getMethodName();
			System.out.print(name + " -> ");
			if ( failure.getMessage() != null )
				System.out.print(failure.getMessage());
			else
				System.out.print("No failure message present " +
						"(indicates systemic issues somewhere in the codebase)." +
						"\nTrace: " + failure.getTrace());
			System.out.print("\n");
		}
		System.out.println("-----------------------------");
		
		System.out.println("Mark summary:");
		System.out.println("\tEdge:  [gained " + eP + ", lost " + eL + "]");
		System.out.println("\tVector:  [gained " + vP + ", lost " + vL + "]");
		System.out.println("\tGraph:  [gained " + gP + ", lost " + gL + "]");
		System.out.println("Total: [" + (eP + vP + gP) + ", lost " + (eL + vL + gL) +"] (out of: " + (eP + eL + vP + vL + gP + gL) + ")");
	}
}
