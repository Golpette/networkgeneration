import java.io.*;
import java.util.*;

public class PathChecker{

    /** 
        Class to find pathways with most/all metabolites found in KEGG
    **/



    // Compounds in KEGG  **that were in pathways I was investigating**; i.e. there might be more compounds in
    //   our network that are also in KEGG but I **have not systematically checked this**:
    //
    //    182=G3P, 218=13BPG, 142=3PG, 109=2PG, 172=PEP, 54=PYR
    //    914=OXA, 154=23BPG, 96=GLYCERALDEHYDE, 95=GLYCERATE, 101, dihyroxyacetone, 104=DHAP, LACTATE=43, 
    //    52=ACETONE, 71=ALANINE, 88=3-hydroxypropanoate, 94=glycerol, 97-glycerol-3-phosphate, 116=serine
    //    386=threonine, 902=succinate, 908=malate,926=aspartate, 1068=asparagine, 102=hydroxypyruvate
    /////////////   REMOVED:  ENOLPYRUVATE=168;  

    static Integer[] KEGG_COMPS = new Integer[]{182,218,142,109,172,54,914,154,96,95,104,43,52,71,88,94,97,101,116,386,
					    902, 908,926,1068,102};





    public static void main(String args[])throws IOException{


	List<Integer> KEGG_COMPOUNDS = Arrays.asList( KEGG_COMPS );


	// Read the list of pathways 
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
                    // then end of path reached. add it to listOfPaths.
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







	// Check how many compounds are NOT IN KEGG LIST.
	PrintWriter out = new PrintWriter(new FileWriter( "PathChecker_output.dat" ) );
	int c=1;
	for( int z=0; z<listOfPaths.size(); z++ ){

	    int count=0; 

	    for( Reactions3 r: listOfPaths.get(z) ){
		int sub = r.getSubstrate();
		int prod = r.getProduct();

		if( !KEGG_COMPOUNDS.contains( sub ) || !KEGG_COMPOUNDS.contains( prod ) ){
		    count++;
		}

	    }


	    if( count == 0 ){
		// print pathway to file
		for( Reactions3 r: listOfPaths.get(z) ){
		    out.println( r.toString() );
		}
		out.println("Path_"+z+" ends here. Path_KEGG_"+c+".");
		c++;
	    }



	}



	out.close();












    }



}
