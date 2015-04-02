import java.util.*;
import java.io.*;

/**
   Use "mpmath" library (findroot function, implementing the multidimensional Newton-Raphson
   algorithm to solve for steady state concentrations.
**/



public class PythonScript_ConvKin_NoExternals{



    private final static double RT = 2.5;


    /// Method to keep generating and running python scripts, up to a "MAX_TRIES" number of times, to find
    ///	a steady state. If none found, assume none exist and add a "-100.0" to the concentration list to indicate this fact

    public static void getSolution( ArrayList<Reaction> new_network_simplified, ArrayList<String> sources, ArrayList<String> products, ArrayList<Double> sourcesConc, ArrayList<Double> productsConc,  ArrayList<String> nodeLabels, ArrayList<Double> nodeConcentrations, HashMap<String, Double> extMetabConcentrations, int MAX_TRIES, int m, double conc_min, double conc_max )throws IOException, InterruptedException{
	
	/** It is possible that the starting concentrations chosen will lead to the error "Encountered singular Jacobian". It is often the
	    case that changing just one of these values by an order of magnitude is enough to make system very easily solveable.
            makeScript_findRoot() now chooses random start points
        **/

	boolean foundSoln = false;
	int countTries = 0;

	int cc=0;
	while( !foundSoln ){

	    cc++;
	    
	    if( countTries == MAX_TRIES ){
		/** Then we assume there is no steady state and just reject mutation **/
		/** Easy way to do this, at this stage, is to add a negative concentration to nodeConcentrations
		    (since this results in path being infeasible and thus flux=0) and break this loop  **/

	  	// System.out.println("Have tried " + MAX_TRIES + " initial conditions in python script. Giving up...");		
	
		// Adding this -100.0 to the end of the nodeConcentrations list indicates that no steady state was found
		nodeConcentrations.add( -100.0 );
		break;
	    }
	    countTries++;
	    //System.out.println( "Attempt " + countTries );
	    
	    /** Make mathematica script  */
	    String scriptName_2 = "pyScript_"+Integer.toString( m )+".py";
	    makeScript_findroot( new_network_simplified, sources, products, sourcesConc, productsConc, scriptName_2, extMetabConcentrations  );	
		
	    String pyOutput_2 = "pyOutput"+Integer.toString( m )+".dat";
	    /**  Run python script and waitFor() it to complete  **/
	    String[] command_2 = {  "/bin/sh",  "-c",  " python " + scriptName_2 + " > " + pyOutput_2  };
	    Process process_2 = Runtime.getRuntime().exec(  command_2   ); 
	    process_2.waitFor();
	    
	    // Once the script is executed, remove the file
	    String[] command_rmScript = {  "/bin/sh", "-c", "rm " + scriptName_2 };
	    Process rmScript = Runtime.getRuntime().exec( command_rmScript );
	    rmScript.waitFor();
	    
	    
	    /** Check output file for error message  **/
	    foundSoln = PythonOutputReader.solutionExists(  pyOutput_2  );
	    
	    if( foundSoln ){		
		/** Read in metabolites and their steady state concentrations from mathematica output file. ALSO include sources and products  **/
		PythonOutputReader.getLabelsAndConcentrations( pyOutput_2, nodeLabels, nodeConcentrations );
		addSourcesAndProds_labelsAndConc( sources, products, sourcesConc, productsConc, nodeLabels, nodeConcentrations );



		// NOW CHECK WE HAVE A STEADY-STATE WITH POSITIVE CONCENTRATIONS
		boolean concPos = true;
		//	boolean concOK = true;
		
		concPos = checkConcPos( nodeConcentrations );
		if( !concPos ){
		    // System.out.println("Neg conc in steady state");

		    //------ Keep looking 
		    foundSoln = false;
		}

		    	
	    }
	    
	    /** Whether or not solution found, remove the output file: **/
	    String[] command_rmOutput = {  "/bin/sh", "-c", "rm " + pyOutput_2 };
	    Process rmOutput = Runtime.getRuntime().exec( command_rmOutput );
      	    rmOutput.waitFor();

	    
	}//while loop
	


    }
    //---------------------------------------------------------------------------------------------------







    // Method to generate a mathematica script that can be run from command line (and or JAVA programs). Format is for "FindRoot" for system of eqns at steady state
    //   3/12/13: modified to use random distribution of start points (starting all at 1e-5M was causing occasional problems)
    public static void makeScript_findroot( ArrayList<Reaction> reactions, ArrayList<String> sources, ArrayList<String> products, ArrayList<Double> sourcesConc, ArrayList<Double> prodsConc, String scriptname, HashMap<String, Double> externalMetabConcentrations  )throws IOException{


	if( reactions.size() == 0 ){
	    System.out.println(   "reactions.size() == 0 in makeScript_FindRoot   OMG WTF?!!? "   );
	}


	PrintWriter out = new PrintWriter( new FileWriter( scriptname ) );

	/** Import mpmath */
	out.println("from mpmath import *");
	out.println();

	/** Make list of nodes */
	ArrayList<Node> nodes = makeNodes( reactions );

	/** Generate lists of parameter definitions (as Strings) from each Reaction
            to be defined at beginning of script */
	ArrayList<String> parameters = getParameters( reactions, externalMetabConcentrations  );

	/** Generate ODE for each internal metabolite (Node) */
	buildEqn_ConvenienceKinetics_NoExternals( reactions,  nodes,  externalMetabConcentrations );


	/** Print to pyScript.py */
	/** Sources, products and their fixed concentrations **/
	if( sources.size() != sourcesConc.size() ){ System.out.println("ERROR: Sources and their concentrations don't match;  makeScript_FindRoot()");  }
	for( int i=0; i<sources.size(); i++ ){
	    out.println( "c"+sources.get(i) + " = "  +  convertDoubleString( sourcesConc.get(i) )  +  ";");
	}
	if( products.size() != prodsConc.size() ){ System.out.println("ERROR: Products and their concentrations dont match;  makeScript_FindRoot()");  }
	for( int j=0; j<products.size(); j++ ){
	    out.println( "c"+products.get(j) + " = "  + convertDoubleString( prodsConc.get(j) ) + ";"  );
	}
	/** Print out all external metabolite concentrations **/
	for (Map.Entry<String, Double> entry : externalMetabConcentrations.entrySet()  ) {
	    String em = entry.getKey();
	    em = em.substring( 1,em.length() ).toLowerCase();
	    out.println( em + " = " + convertDoubleString( entry.getValue() ) + ";" ); 
	}

	/** Parameters -  Vf, Ks and Kp values for all reactions **/
	for( int p=0; p<parameters.size(); p++ ){
	    out.println( parameters.get(p) + ";" );
	}
	out.println();
	out.println();


	
	String vars = "";
	String var_names_list = "";
	for( int n=0; n<nodes.size(); n++ ){
	    if( !sources.contains( nodes.get(n).getLabel() )  &&  !products.contains( nodes.get(n).getLabel() ) ){
		vars = vars + "c"+nodes.get(n).getLabel() + ", ";
		var_names_list = var_names_list + "\"c"+nodes.get(n).getLabel() + "\", ";
	    }
	}
      	vars = vars.substring( 0, vars.length()-2 );
      	var_names_list = var_names_list.substring( 0, var_names_list.length()-2 );



	/** Equations - do not print those of source and product metabolites since these are at fixed concentration in flux calculation  **/
	String list_eqs = "";
	for( int n=0; n<nodes.size(); n++ ){
	    if( !sources.contains( nodes.get(n).getLabel() )  &&  !products.contains( nodes.get(n).getLabel() ) ){

		out.println( "def eq"+nodes.get(n).getLabel() + "(" + vars + "):"  );
		out.println( "\t"+"return " + nodes.get(n).getEqn()  );

		list_eqs = list_eqs + "eq"+nodes.get(n).getLabel() +"("+vars+")" + ", ";
	    }
	}
      	list_eqs  = list_eqs.substring( 0, list_eqs.length()-2 );


	out.println();
	out.println("def F("+vars+"):");
	out.println("\t" + "return " + list_eqs );





	// Choose RANDOM starting point for algorithm   -- method 1;  works best. why? 
       	double[] strtpnts = new double[]{1.0E-5, 1.0E-6, 1.0E-7};
	// For some reason this set of start points gives much faster and more consistent determination
	// of the steady state than just selecting randomLy from the logarithmic distribution
	// say from 1e-7 to 1e-3... I don't know why though - the reason I did the random selection was
	// because starting from a "special point", like all concentrations equal, leads to problems.

	String startpts = "";
	for( int n=0;  n<nodes.size(); n++ ){
	    if( !sources.contains( nodes.get(n).getLabel() )  &&  !products.contains( nodes.get(n).getLabel() ) ){

		double startPnt = strtpnts[   (int)( Math.random()*strtpnts.length )   ];
		//	random multiplier from 1-9
		int randMultiplier =    (int)( Math.random()*9 ) + 1;
		String sp = convertDoubleString(   startPnt*randMultiplier   );

		startpts = startpts + " "+sp+", ";
	    }
	}
      	startpts  = startpts.substring( 0, startpts.length()-2 );


	// // Choose RANDOM starting point for algorithm        -- METHOD 2; WHY IS THIS MUCH WORSE?!
	//  --- THIS METHOD was sometimes leading to no steady-state solution being found
	//      in 20 attempts.
	//
	// double strt_min = 1.0E-7;   double strt_max = 1.0E-3;
	//
	// String startpts = "";
	// for( int n=0;  n<nodes.size(); n++ ){
	//     if( !sources.contains( nodes.get(n).getLabel() )  &&  !products.contains( nodes.get(n).getLabel() ) ){
	//
	// 	//Choose concentration uniformly from logarithmic range
	// 	double rand = Math.random() * (  Math.log(strt_max) - Math.log(strt_min)  );
	// 	double ln_c = Math.log( strt_min ) + rand;
	// 	double c = Math.exp( ln_c  );
	//
	// 	String sp = convertDoubleString( c );
	//
	// 	startpts = startpts + " "+sp+", ";
	//     }
	// }
      	// startpts  = startpts.substring( 0, startpts.length()-2 );



	out.println();
	//	out.println("sol = findroot(   F, " + "("+startpts+" )" + "   )");
      	out.println("sol = findroot(   F, " + "("+startpts+" )" + ", " + " solver='mdnewton', maxsteps=30, tol=1e-10 " + "   )");   // chose maxsteps and tol from trial and error

	out.println();
	out.println("variables = [ "+var_names_list+ " ]");


	out.println();
	out.println(  "for i in range( len(sol) ):"  );
	out.println(  "\t" + "print variables[i], sol[i]"  );


	out.close();

    }//================================================================================ end FindRoot generator









    //##############################################################
    //
    // Below are methods required in the script generator:
    //
    //##############################################################









    // Method to convert a double (if in JAVA exponential notation) to the format that mathematica accepts
    /** i.e. 2.345E-6 is NOT correctly recognized by mathematica - it is read as 2.345*(e^-6) where "e" is Euler's constant  
             It needs to be written:  2.345*^-6  

        NO LONGER NEEDED SINCE WE USE "mpmath" python LIBRARY NOW 
    **/
    public static String convertDoubleString( double number ){

    	String num = Double.toString( number );
	//	if( num.contains("E") ){
	//	    /**  replace "E" with "*^"  **/
	// 	    num = num.replace( "E", "*^");
	//  	}

    	return num;
    }//================================================================================




    // Method to generate lists of all Reaction parameters (Vm, Ks, Kp, q and all K_externals)
    public static ArrayList<String> getParameters( ArrayList<Reaction> reactions, HashMap<String, Double> externalConcentrations ){

	ArrayList<String> parameters = new ArrayList<String>();

	for( int i=0; i < reactions.size(); i++ ){

	    Reaction rxn = reactions.get(i);

            /** Default exponential notations of these doubles is not readable by mathematica!  **/
	    String kcat = convertDoubleString( rxn.getkcat() );
	    String Ei = convertDoubleString( rxn.getEi() );
	    String Ks = convertDoubleString( rxn.getKs() );
	    String Kp = convertDoubleString( rxn.getKp() );


	    // We are NOT INCLUDING the kinetics (no Km terms or concentrations of the external metabolites
	    // Instead, they only enter through the modification of the equilibrium constant
	    ArrayList<String> ext_subs = rxn.getExtSubs();
	    ArrayList<String> ext_prods = rxn.getExtProds();
	    /** Check all external metabolites of this Reaction are actually defined  **/
	    for( String s: ext_subs ){
		if(  !externalConcentrations.containsKey( s ) ){
		    System.out.println("ERROR: reaction contains unknown external metabolite in getParameters(). Exiting. ");
		    System.out.println(" --  " + s);
		    System.exit(0);
		}
	    }
	    for( String s: ext_prods ){
		if(  !externalConcentrations.containsKey( s ) ){
		    System.out.println("ERROR: reaction contains unknown external metabolite in getParameters(). Exiting. ");
		    System.out.println(" --  " + s);
		    System.exit(0);
		}
	    }
	    /** Multiply all **external substrate** concentrations  **/
	    double multAllExtSubs = 1.0;
	    for( String s: ext_subs ){
		double em_concentration = externalConcentrations.get( s );
		multAllExtSubs = multAllExtSubs * em_concentration; 
	    }
	    /** Multiply all **external product** concentrations **/
	    double multAllExtProds = 1.0;
	    for( String s: ext_prods ){
		double em_concentration = externalConcentrations.get( s );
		multAllExtProds = multAllExtProds * em_concentration; 
	    }


	    double dG = rxn.getDeltaG();
	    String q = convertDoubleString(  (multAllExtSubs/multAllExtProds)* Math.exp( -(dG/RT) )   );
	    // String q =   (multAllExtSubs/multAllExtProds)* Math.exp( -(dG/RT) );


	    parameters.add("q" + Integer.toString(i) + " = " + q);
	    parameters.add("kcat" +Integer.toString(i) +" = " + kcat );
	    parameters.add("Ei" +Integer.toString(i) +" = " + Ei );
	    //	    parameters.add("Vf" +Integer.toString(i) +" = " + Vf );
	    parameters.add("Ks" +Integer.toString(i) + " = " + Ks);
	    parameters.add("Kp" +Integer.toString(i) + " = " + Kp);



	}

	return parameters;
    }//===================================================================




    // Method to generate Nodes
    public static ArrayList<Node> makeNodes( ArrayList<Reaction> reactions ){

	ArrayList<Node> nodeList = new ArrayList<Node>();	
	/**  Get all metabolites and make a Node object for them  **/
	ArrayList<String> allMetabs = Methods.getAllMetabs( reactions );
	for( String metab: allMetabs){    
	    Node newnode = new Node( metab );
	    nodeList.add(  newnode  );	   
	}	
	
	return nodeList;	
    }//=================================================================






    // Method to build the mathematicaEqn element for each Node (updates the node list) assuming all reactions use "convenience kinetics" 
    //    BUT kinetics of the external metabolites are NOT INCLUDED.                              
    public static void buildEqn_ConvenienceKinetics_NoExternals( ArrayList<Reaction> reactions,  ArrayList<Node> nodes, HashMap<String, Double> externalConcentrations ){

	for( int r=0; r<reactions.size(); r++ ){

	    /** If there are multiple internal substrates in this reaction, they all have the same Michaelis-Menten
	            constant, Ks. The same is true for multiple products, all having Kp.	    
	        NB: since we label all our internal products as numbers, we append a "c" to the front of each label
	            as Mathematica does not allow variable names to be numbers   **/

	    Reaction rxn = reactions.get(r);

	    ArrayList<String> int_subs = rxn.getIntSubs();
	    ArrayList<String> int_prods = rxn.getIntProds();

	    /** Format for convenience kinetics is v = k_forw * ( allSubsMultiplied - allProdsMult/q )
	        where kf = (Vf/denom_1) * ( 1 / ( denom_2a + denom_2b - 1 )  ) and
	        denom_1 as product of all substrate Michaelis-Mentent constants,
	        denom_2a product of (1 + [s_i]/Ks_i) terms and denom_2b the same idea for products  **/
	    // Multiply all substrate (internal and external) Michaelis-Menten constants.
	    String denom_1 = "";
	    for( String s: int_subs ){
		denom_1 = denom_1 + "*Ks"+r;
	    }
	    denom_1 = denom_1.substring( 1, denom_1.length() );

	    // Multply all terms like (1+c123/Ks) for internal and external substrates
	    String denom_2a = "";
	    for( String s: int_subs ){
		String bit = "(1 + " + "c"+s + "/Ks"+r + ")" ;
	        denom_2a = denom_2a + "*" + bit;
	    }
	    denom_2a = denom_2a.substring( 1, denom_2a.length() );
	    
	    // Multiply all terms like (1+c234/Kp) for all internal and external products
	    String denom_2b = "";
	    for( String s: int_prods ){
		String bit = "(1 + " + "c"+s + "/Kp"+r + ")" ;
	        denom_2b = denom_2b + "*" + bit;
	    }
	    denom_2b = denom_2b.substring( 1, denom_2b.length() );


	    // Format for convenience kinetics is: v = k_forw*(allSubs - allProds/q) with k_forw = (Vf/denom_1) * ( 1 / ( denom_2a + denom_2b - 1 )  )
	    //  -- replaced Vf with (k_cat*Ei)
	    //	    String k_forw = "(Vf"+r+"/("+denom_1+"))*(1 / ( " +denom_2a + " + " + denom_2b + " - 1 ) )";         
	    String k_forw = "( kcat"+r+"*Ei"+r+"/("+denom_1+"))*(1 / ( " +denom_2a + " + " + denom_2b + " - 1 ) )";         


	    // Multiply all internal substrates names
	    String all_subs="";
	    for( String s: int_subs ){
		all_subs = all_subs + "*" + "c"+s;
	    }
	    // Multiply all internal product names
	    String all_prods="";
	    for( String s: int_prods ){
		all_prods = all_prods + "*" + "c"+s;
	    }
	    // Get rid of leading asterixes for the simplified equations:
	    String all_subs_2 = all_subs.substring( 1, all_subs.length() );
	    String all_prods_2 = all_prods.substring( 1, all_prods.length() );


	    // Add appropriate equations to all substrates involved
	    for( String s: int_subs ){
	
		int index = Methods.getIndex( nodes, s );
		nodes.get( index ).addToEqn( " - " + k_forw + "*( " + all_subs_2 + " - " + all_prods_2 + "/q"+r + " )" );
	    }
	    // and all products involved
	    for( String s: int_prods ){
	
		int index = Methods.getIndex( nodes, s );
		nodes.get( index ).addToEqn( " + " + k_forw + "*( " + all_subs_2 + " - " + all_prods_2 + "/q"+r + " )" );
	    }
	       
	}

    }// end buildEqn_ConvenienceKinetics_NoExternals() ================================================================================================










    
    
    
    //Method to add sources and products labels and concentrations to respective lists
    public static void addSourcesAndProds_labelsAndConc(ArrayList<String> sources, ArrayList<String> products, ArrayList<Double> sourcesConc, ArrayList<Double> productsConc, ArrayList<String> nodeLabels, ArrayList<Double> nodeConcentrations){
	
	if( sources.size() != sourcesConc.size() ){ System.out.println("sources and sourcesConc don't match. Exiting.");  System.exit(5);   }
	for( int i=0; i<sources.size(); i++){
	    nodeLabels.add( sources.get(i)  );
	    nodeConcentrations.add( sourcesConc.get(i) );
	}
	if( products.size() != productsConc.size()  ){  System.out.println("products and productsConc don't match. Exiting.");  System.exit(5);   }
	for( int j=0; j<products.size(); j++ ){
	    nodeLabels.add(  products.get(j)  );
	    nodeConcentrations.add(  productsConc.get(j)   ); 
	}
    }
    // ====================================================================================================






    // Method to check all concentrations < MAX_CONC ---------------------------------------------------
    public static boolean checkConcentrations( ArrayList<Double> concentrations, double concMin, double concMax ){

	boolean concentrationsOK = true;

	for( int i=0; i < concentrations.size(); i++ ){
	    if(  concentrations.get(i) > concMax  ||  concentrations.get(i) < concMin  ){
		concentrationsOK = false;
	    }
	}
	return concentrationsOK;
    }
    // ==============================================================================================================





    // Method to check all steady-state concentrations are positive ---------------------------------------------------
    public static boolean checkConcPos( ArrayList<Double> concentrations ){

    	boolean concentrationsOK = true;

    	for( int i=0; i < concentrations.size(); i++ ){
    	    if(  concentrations.get(i) < 0  ){
    		concentrationsOK = false;
    	    }
    	}
    	return concentrationsOK;
    }
    // ==============================================================================================================






}//end class
