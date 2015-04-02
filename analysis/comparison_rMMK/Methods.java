import java.util.*;

public class Methods{




    /** Method to get a list of all metabolites in a given network **/
    public static ArrayList<String> getAllMetabs( ArrayList<Reaction> network ){

	ArrayList<String> allMetabs = new ArrayList<String>();

	for( int r=0; r < network.size(); r++ ){

	    // THIS IS OK!?!?
	    ArrayList<String> rxnMetabs = network.get( r ).getReactionMetabolites();

	    for( int m=0; m < rxnMetabs.size(); m++ ){
		if( !allMetabs.contains(  rxnMetabs.get(m)  ) ){
		    String s = rxnMetabs.get(m);
		    allMetabs.add(  s   );
		}
	    }
	}
	return allMetabs;
    }
    // ===============================================================================================






    // Get index of given metabolite in list of Nodes:
    public static int getIndex(ArrayList<Node> nodes, String label){
	int index = 999777999;
	for(int i=0; i<nodes.size(); i++){
	    if( label.equals( nodes.get(i).getLabel() ) ){
		index = i;
	    }
	}
	if( index == 999777999 ){
	    System.out.println("Methods.java -- index of Node label == 999777999; i.e. not in nodes" );
	}
	return index;
    }
    //====================================================================================








}