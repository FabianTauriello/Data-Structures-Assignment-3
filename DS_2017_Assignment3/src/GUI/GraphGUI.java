package GUI;

import javafx.application.Application;
import javafx.fxml.FXMLLoader;
import javafx.scene.Parent;
import javafx.scene.Scene;
import javafx.scene.canvas.Canvas;
import javafx.stage.Stage;

public class GraphGUI extends Application {
	
	static Canvas canvas;
	static Canvas player_canvas;
	
	public void start(Stage stage) throws Exception {
		stage.setTitle("Data Structures 2017 - Assignment 3");
		stage.setResizable(false);
		
		Parent root = FXMLLoader.load(getClass().getResource("SimpleGraphGUI.fxml"));
		Scene scene = new Scene(root, 720, 600);
		
		canvas = (Canvas)scene.getRoot().getChildrenUnmodifiable().get(1);
		player_canvas = (Canvas)scene.getRoot().getChildrenUnmodifiable().get(2);
		
		new GraphRenderer();
		
		stage.setScene(scene);
		stage.show();
	}
	
	public static void main(String[] args) {
		launch();
	}
}
