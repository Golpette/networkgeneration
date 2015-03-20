//this is used for RandomWalker.java where we need each vertex to have a list of nearest neighbours
import java.io.*;
import java.util.*; 

public class Vertex3{

    int compoundLabel;
    ArrayList<Integer> neighbours = new ArrayList<Integer>();
    ArrayList<Reactions3> reactions = new ArrayList<Reactions3>();

    Vertex3(int cl){
	this.compoundLabel = cl;
	this.neighbours = new ArrayList<Integer>();
	this.reactions = new ArrayList<Reactions3>();
    }


    //method to copy a vertex
    public void setTo(Vertex3 a){
	
	this.compoundLabel = a.compoundLabel;
	int num = a.getNumNeighbours();
	for(int g=0; g<num; g++){
	    this.neighbours.set( g, a.getNeighbour( g )  );
	}
	int numReact = a.getNumReactions();
	for(int r=0; r<numReact; r++){
	    this.reactions.set( r, a.getReaction( r ) );
	}
    }


    // Add a reaction to the nodes reaction list
    public void addReaction(Reactions3 reac){
	this.reactions.add(reac);
    }
    // Get number of possible reactions from node
    public int getNumReactions(){
	return this.reactions.size();
    }
    // Get a specific reaction (Reactions3 object) from list
    public Reactions3 getReaction(int a){
	return this.reactions.get(a);
    }


   //method to add a vertex to the list of neighbouring vertices
    public void addNeighbour(int a_neighbour){
	this.neighbours.add(a_neighbour);
    }
    //method to get size of neighbours list
    public int getNumNeighbours(){
	return this.neighbours.size();
    }
    public int getNeighbour(int a){
	return this.neighbours.get(a);
    }



    //method to remove a neighbour from neighbours list
    //(so we can have self-avoiding walks)
    //a_neighbour is the INDEX in the neighbours array we want to remove
    public void removeNeighbour(int a_neighbour){
	this.neighbours.remove(a_neighbour);
    }



    //get the ArrayList index of a neighbouring vertex
    public int getIndexOfNeighbour(int label){
	return this.neighbours.indexOf(label);
    }



   public  String toString(){
       return "label = " + compoundLabel;
    }


    int getLabel(){
	return compoundLabel;
    }




  
 

}