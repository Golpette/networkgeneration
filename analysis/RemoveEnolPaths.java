import java.io.*;
import java.util.*;


public class RemoveEnolPaths{

    /**
       Class to read in a list of pathways and remove all pathways
       in which an enol compound is present - there are too many of them
       confusing the comparison
     **/

    public static void main(String args[])throws IOException{


	//---- Read in pathways ----
	BufferedReader reactionsIN = new BufferedReader( new FileReader(args[0]) );
	ArrayList<ArrayList<Reactions3>> listOfPaths = new ArrayList<ArrayList<Reactions3>>();
	ArrayList<Reactions3> listOfReactions = new ArrayList<Reactions3>();
	String s2 = "";
	Scanner scan2 = new Scanner(s2);
	while(true){
	    try{
		s2 = reactionsIN.readLine();
		scan2 = new Scanner(s2);
		
		if( s2.contains("Path") ){
		    ArrayList<Reactions3> justholding = new ArrayList<Reactions3>();
		    for( int qwa=0; qwa < listOfReactions.size(); qwa++ ){
			justholding.add( listOfReactions.get( qwa ) );
		    }
		    listOfPaths.add( justholding );
		    listOfReactions.clear();
		}

		else if( !(s2.contains("Path"))  ){
		    String direction = scan2.next();
		    String type = scan2.next();
		    int substrate = scan2.nextInt();
		    int product = scan2.nextInt();
		    int ATP = scan2.nextInt();
		    int NADH = scan2.nextInt();
		    double deltaG = scan2.nextDouble();	    
		    // read in all the next string pieces to form the equation
		    String ss = "";
		    while(true){
			try{
			    ss = ss + " " + scan2.next();
			}
			catch(NoSuchElementException e){break;}
		    }
		    String equation = ss;
		    listOfReactions.add(new Reactions3(direction, type, substrate, product, ATP, NADH, deltaG, equation));
		    
		}
		else{   System.out.println("ERROR: Check format of input file of pathways!");   }
	    }
	    catch(NullPointerException e){break;}
	} // --------- Now have list of paths (where each path is a list of 'Reactions3' reactions)









	// ONLY print pathways that do not contain an enol compound 
	PrintWriter out = new PrintWriter(new FileWriter( "Paths_out_NO_ENOLS.paths" ) );
	int c=1;
	for( int z=0; z<listOfPaths.size(); z++ ){

	    boolean pathContainsEnol = false;

	    for( Reactions3 r: listOfPaths.get(z) ){
		
		String eqn = r.getEquation();

		if( eqn.contains("=CH(OH)") || eqn.contains("CH(OH)=") || eqn.contains("=C(OH)-") || eqn.contains("-C(OH)=")  ){

		    pathContainsEnol = true;
		}

	    }


	    if( ! pathContainsEnol ){
		// print pathway to file
		for( Reactions3 r: listOfPaths.get(z) ){
		    out.println( r.toString() );
		}
		out.println("(Path_"+z+" ends here).  Path_noEnol_"+c+".");
		c++;
	    }

	}



	out.close();








    }



}
