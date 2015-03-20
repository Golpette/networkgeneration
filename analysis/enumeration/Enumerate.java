import java.io.*;
import java.util.*; 
import java.text.*;

/** 
  Enumerates paths of length L using depth-first search algorithm (DFS)
  Uses:     Compounds.java, Vertex3.java,  Reactions3.java
  Reads in from terminal:    compound_file, reaction_file, desired path length

  FLUX METHOD:  currently does not compute flux, but could add back in the calculation from Heinrich1997 Eq(40), 
   to allow us to allows us to determine if a pathway is at least thermodynamically feasible for a specific 
   set of conditions.
**/


public class Enumerate{


    // SET START AND END POINTS:
    // Metabolites are labelled by an index in the "compounds_..." file. 
    // 182 = Glyceraldehyde-3-phosphate (G3P);   54 = Pyruvate;
    //
    static final int startingNode = 182;  static final int finishingNode = 54; 
    //static final int startingNode = 54;  static final int finishingNode = 182; 



    // Print paths that have this ATP consumption/production
    static final int desiredATP = 2;




    public static void main(String args[])throws FileNotFoundException,IOException{
	
	long startTime = System.currentTimeMillis();


	/* READ IN BARTEK'S COMPOUND AND REACTION FILES  */
	
	//1 --- getcompounds and vertices:--------------------------------------
	BufferedReader in = new BufferedReader( new FileReader(args[0]));
	ArrayList<Compounds> compounds = new ArrayList<Compounds>(); 
	String s = "";
	Scanner scan = new Scanner(s);
	while(true){
	    try{
		s = in.readLine();
		scan = new Scanner(s);
		int index = scan.nextInt();
		double dG = scan.nextDouble();
		String formula = scan.next();
		String name = scan.next();
		int charge = scan.nextInt();
		compounds.add(new Compounds(index, dG, name, formula, charge));
	    }
	    catch(NullPointerException e){break;}
	    
	}

	ArrayList<Vertex3> vertices = new ArrayList<Vertex3>();   //create vertices. Don't initialise list of neighbours though.
	    for(int j=0; j<compounds.size(); j++){
		int compoundLabel = compounds.get(j).getIndex();		
		vertices.add(new Vertex3(compoundLabel));    
	    }






	//2 ------ get reactions and edges of network ------------------------
	    BufferedReader in2 = new BufferedReader( new FileReader(args[1]));
	    //ArrayList<Reactions3> reactions = new ArrayList<Reactions3>();
	    String s2 = "";
	    Scanner scan2 = new Scanner(s2);
	    while(true){
		try{
		    s2 = in2.readLine();
		    scan2 = new Scanner(s2);
		    String type = scan2.next();
		    int substrate = scan2.nextInt();
		    int product = scan2.nextInt();
		    int ATP = scan2.nextInt();
		    int NADH = scan2.nextInt();
		    double deltaG = scan2.nextDouble();
		    // read in all the next string pieces to form the equation
		    // (and its reverse)
		    String ss = "";
		    String ssReverse = "";
		    while(true){
			try{
			    String nextPiece = scan2.next();
			    ss = ss + " " + nextPiece;
			    ssReverse = nextPiece + " " + ssReverse;
			}
			catch(NoSuchElementException e){break;}
		    }
		    String equation = ss;   String rev_equation = ssReverse;

		    // create a new Reactions3 object for each direction and store in "vertices" list of correct metabolite
		    Reactions3 thisReaction = new Reactions3("forward", type, substrate, product, ATP, NADH, deltaG, equation );
		    vertices.get(substrate).addReaction( thisReaction );
		    vertices.get(substrate).addNeighbour( product );
		   
		    Reactions3 revReaction = new Reactions3("backward", type, product, substrate, -ATP, -NADH, -deltaG, rev_equation);
		    vertices.get(product).addReaction( revReaction );
		    vertices.get(product).addNeighbour( substrate );


		}
		catch(NullPointerException e){break;}
		
	    }	    	    
     	    //--------------------------------------------------------------end of REACTIONS	  

















	    //Prints out time
	    double time = System.currentTimeMillis();
	    double timeToReadFiles = (time - startTime)/1000.0;





	    //________________________ENUMERATE ALL PATHS___________________________________________


	    // Read in path length from terminal
	    int L = (int)Double.parseDouble(args[2]);  
	   
	    // Define start and end metabolite
	    int SOURCENODE = startingNode;
	    int FINALNODE = finishingNode;


	    // Current path length (number of reactions in pathway)
	    int l = 0;   

	    // This stores the reactions being used in current pathway
	    ArrayList<Reactions3> reactionsUsed = new ArrayList<Reactions3>();
	    // This stores the (integer) name of the nodes we have visited
	    ArrayList<Integer> listOfTraversedNodes = new ArrayList<Integer>(); 

	    // Add starting point to list of traversed nodes
	    listOfTraversedNodes.add( SOURCENODE ); 




	    // Do the enumeration-----------	    	    
	    int[] numPaths = new int[4];


	    // RECURSIVE DFS ALGORITHM   ( pathways are printed to terminal from INSIDE THIS METHOD )
	    numPaths = takeStep_counter(SOURCENODE, FINALNODE, l, L, listOfTraversedNodes, reactionsUsed, vertices, compounds, numPaths );


	    // numPaths[0] ==  total number of pathways connecting start and end nodes in Bartek's network
	    // numPaths[1] ==  number which have all compounds electrostatically charged
	    // numPaths[2] ==  satisfying ATP restriction; e.g. ATP >= defined amount
	    // numPaths[3] ==  positive flux for the "external metabolite" concentrations defined at top of program (currently not being used)
	    System.out.println();
	    System.out.print("numPaths = ");
	    for(int d=0; d<numPaths.length; d++){
		System.out.print(numPaths[d] + " ");
	    }
	    System.out.println();
	   



	    // Time to read Bartek's files
	    System.out.println();
	    System.out.println("Time to read in Bartek's files = " + timeToReadFiles + " seconds");
	    //Prints out run time
	    long stopTime = System.currentTimeMillis();
	    long elapsedTime = stopTime - startTime;
	    System.out.println("Total runtime = " + elapsedTime/1000.0 + " seconds");
	    System.out.println();




    }
    // end main ------------------------------------------------------------------------------------------------------------
 


















    //###########################################################################################
    /****   Method to enumerate all paths of length L using depth-first search algorithm    
            Modify this method to alter restrictions on which paths to output    ****/
    //###########################################################################################


    public static int[] takeStep_counter(int curNode, int FinalNode, int l, int L, ArrayList<Integer> listOfTraversedNodes, ArrayList<Reactions3> reactionsUsed, ArrayList<Vertex3> vertices, ArrayList<Compounds> compounds, int[] counter){
	
	//  counter[] holds number of paths satisfying:   [ ALLlinear paths in network connecting start and end node (no restrictions),  paths with charged metabolites only, 
	//                                                  ATP restrictions,  flux restrictions (i.e. thermodynamic feasibility)      ]


	// If path is correct length, check to see if it satisfies our restrictions and print it out
    	if( l == L ){


	    // Check it ends at desired node	    
    	    if(  curNode == FinalNode ){
		
    		counter[0]++;		 

		// Find ATP and NADH production of pathway
    		int ATP_path = 0;  int NADH_path = 0;
    		for( int y=0; y<reactionsUsed.size(); y++ ){
    		    ATP_path = ATP_path + (   reactionsUsed.get(y).getATP()   );
    		    NADH_path = NADH_path + (  reactionsUsed.get(y).getNADH()  );
    		}                           


		// Here we want to exclude pathways that contain any compounds with no electrostatic charge
    		ArrayList<Integer> chargeArray = getChargeArray(listOfTraversedNodes, compounds);
		//	if( true ){
		if( !(chargeArray.contains(0)) ){
		    
    		    counter[1]++;
		    
		    // if( true ){          // Use this for gluconeogenesis since we have no ATP restriction
		    // if( ATP_path >= desiredATP ){   
		    if( ATP_path == desiredATP ){
									
    			counter[2]++;
			

			// REMOVED THIS, but could re-add a flux calculation here to check thermodynamic feasibility under a specific
			//               set of external metabolite concentrations
			//
			//	double flux = getFlux(  reactionsUsed  );


			//if( flux > 0 ){
			if(true){
			    
			    counter[3]++;
			    
			    //	Print out path in Bartek's reaction format  --------------------------------------
			    String pathname = "Path_"+( counter[3] );				    
			    for( int xvx=0; xvx < reactionsUsed.size(); xvx++){
				System.out.println( reactionsUsed.get(xvx).toString() );
			    }
			    System.out.println( pathname + " ends here.  ATP produced = " + ATP_path + "." );
			    // -----------------------------------------------------------------------------------
			
			        
    			}//flux>0


		    }//desiredATP
	    
		}//charge

	    }//curnode=finnode

	    
    	    return counter;
    	}
	


	// *****  This small section is the actual DFS algorithm  *****
	
    	for(int j=0; j<vertices.get(curNode).getNumReactions(); j++){

    	    int nextNode = vertices.get(curNode).getReaction(j).getProduct();

    	    if( !(listOfTraversedNodes.contains(nextNode)) ){

    		listOfTraversedNodes.add(nextNode);
    		reactionsUsed.add(  vertices.get(curNode).getReaction(j)  );
    		l++;

		// The DFS search method recursively calls itself
    		takeStep_counter(nextNode, FinalNode, l, L, listOfTraversedNodes, reactionsUsed,  vertices, compounds, counter);

		// Remove last visited node and last reaction, reduce path length by 1
    		listOfTraversedNodes.remove(    (listOfTraversedNodes.size() - 1)    );
    		reactionsUsed.remove(  (reactionsUsed.size()-1)  );
    	       	l--;
    	    }
    	}


    	return counter;
    }

    
    //_____________________________________________________________________________________________________________________









   //Method to obtain compound charge values of a path, store these in an array
    public static ArrayList<Integer> getChargeArray(ArrayList<Integer> route, ArrayList<Compounds> compounds){
	
	int size = route.size();
	ArrayList<Integer> chargeArray = new ArrayList<Integer>();
	
	for(int i=0; i < size; i++){	
	    chargeArray.add( compounds.get( route.get(i)  ).getCharge()  ); 			
	}
	
	return chargeArray;
    }
    //________________________________________________________________________________________________________________________________
    






    
    
    
    //method to round a double to 3 decimal places
    public static double roundToThreeDP(double dbl){
	
        DecimalFormat df3 = new DecimalFormat("#.###");
	return  Double.valueOf( df3.format(dbl) ); 
    }//____________________________________________
    
    
    //method to round a double to 2 decimal places
    public static double roundToTwoDP(double d){
	
        DecimalFormat df = new DecimalFormat("#.##");
	return  Double.valueOf( df.format(d) ); 
    }//____________________________________________
    
    






}//end class
