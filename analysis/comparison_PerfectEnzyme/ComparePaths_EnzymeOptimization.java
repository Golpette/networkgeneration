import java.io.*;
import java.util.*; 
import java.text.*;


/**
  Program to sample the parameter space of external metabolites and calculate the flux of all 
  pathways inputed at each point sampled.

  ==> Using Heinrich1997 'perfect enzyme' calculation and Powell's optimization to find enzyme 
      distribution that maximizes flux

  Output of program is list of parameters followed by list of pathway fluxes.
  These are then analysed in a different program (Anaylse_comparison.java).
**/




public class ComparePaths_EnzymeOptimization{



    // Number of parameter samples
    private static final int NUMBERSETS = 10;    


    //Limits on the intermediate concentrations
    //private static final double conc_min = 1.0E-9;   private static final double conc_max = 0.5;
    private static final double conc_min = 1.0E-7;   private static final double conc_max = 1.0E-2;
    // private static final double conc_min = 0;   private static final double conc_max = 1.0E+20;



 

    public static void main(String args[])throws FileNotFoundException,IOException, InterruptedException{
	
     
	ArrayList<Double> sourceConc = new ArrayList<Double>();
	sourceConc.add(0.0);
	ArrayList<Double> prodConc = new ArrayList<Double>();
	prodConc.add(0.0);


      	long startTime = System.currentTimeMillis();



	// Read the list of reactions corresponding to each pathway and store each list in
	// a list of pathways

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


	}
	//Now have list of paths (where each path is a list of 'Reactions3' reactions




	//###########################################################################################
	//Here are all the --- VALUES WE NEED TO DEFINE --- before sampling the parameter space:
       
	double RT = 2.5;
	
	//Parameter ranges ---------------------------- 

	//Every run, choose the parameter values uniformly randomly from these ranges

       	double SOURCEmin = 1.0E-6;     double SOURCEmax = 1.0E-3;

	double SourceProduct_min = 1.0E-2;  
	double SourceProduct_max = 1.0E+2;

	double ATPvalue = 1.0E-3;        //this value doesn't matter - only ATP/ADP ratio does (and now ATP/AMP). These values WOULD MATTER if we use realistic kinetics.
       	double ATPADP_min = 0.1;    	double ATPADP_max = 100;

	double AMP_min = 0.01E-3;       double AMP_max = 0.1E-3;

	double NADvalue = 1.0E-3;        //this values doesn't matter...
       	double NADNADH_min = 0.1;       double NADNADH_max = 1000;

	double HPO4_min = 1.0E-3;       double HPO4_max = 100.0E-3;
	double CO2_min = 1.0E-6;        double CO2_max = 1.0E-4;

	double PPi_min = 0.1E-3;        double PPi_max = 10.0E-3;

	// ---- NEW EXTERNALS -----------------
	double NH3_min = 1.0E-6;        double NH3_max = 1.0E-3;
	double GLUT_min = 1.0E-3;       double GLUT_max = 100.0E-3;
	double OXO_min = 0.1E-3;        double OXO_max = 10.0E-3;
	//-------------------------------------

	// ############################################################################################



	String filename = args[1];


	PrintWriter out = new PrintWriter( new FileWriter( filename ) );
	out.println("SOURCE" + "\t" + "PRODUCT" + "\t" + "ATP" + "\t" +
		    "ADP" + "\t" + "AMP" + "\t" + "NAD" + "\t" + "NADH" + "\t" + "HPO4" + "\t" +
		    "PPi" + "\t" + "CO2"+ "\t" + "NH3(aq)" + "\t" + "NH2donor(glutamate)" + "\t" + 
                    "NH2acceptor(2-oxo.)" + "\t" +  "fluxes for how ever many paths..." );



	// Use this to keep a running total of flux in each path over all 10,000
	// parameter samples so as to quantify actual performance of paths rather
	// than just the relative measure of how often one performs the best.

	double[] fluxTotals = new double[ listOfPaths.size() ];
	double[] fluxTotals_NORMALISED = new double[ listOfPaths.size() ];
	int[] BESTcount = new int[ listOfPaths.size() ];

	HashMap<String, Double> extMetabConcentrations = new HashMap<String, Double>();
	extMetabConcentrations.put("-atp", 1.0E-3);
	extMetabConcentrations.put("-adp", 1.0E-3);
	extMetabConcentrations.put("-amp", 1.0E-3);
	extMetabConcentrations.put("-nad", 1.0E-3);
	extMetabConcentrations.put("-nadh", 1.0E-3);
	extMetabConcentrations.put("-pi", 1.0E-3);
	extMetabConcentrations.put("-ppi", 1.0E-3);
	extMetabConcentrations.put("-co2", 1.0E-3);
	// new external metabolites required for nitrogen
	// glut = glutamate = NH2donor;    oxo = 2-oxoglutarate = NH2acceptor
	extMetabConcentrations.put("-nh3", 1.0E-3);
	extMetabConcentrations.put("-glut", 1.0E-3);
	extMetabConcentrations.put("-oxo", 1.0E-3);



	// Want to draw random parameter set multiple times 
	for( int sets=0; sets < NUMBERSETS; sets++ ){

	    double SOURCEvalue=0; // ---------<<<<<<<<<< COMMENT OUT IF YOU FIX SOURCEvalue above
	
	    double PRODUCTvalue=0;  double ADPvalue=0;   double AMPvalue=0;
	    double NADHvalue=0;     double HPO4value=0;  double PPivalue=0;   double CO2value=0;
            double NH3value=0;      double GLUTvalue=0;  double OXOvalue=0;



	    //Set the value of the parameters, by semi-selecting ratios this time:----------------------

	    //   Source/Product ratio  and  PRODUCTvalue;

	    // ######  COMMENT THIS OUT IF WE WANT TO FIX [SOURCE] VALUE ###################
	    boolean looking99 = true;
	    while( looking99 ) {
	    	double rand99 = Math.random() * ( Math.log(SOURCEmax) - Math.log(SOURCEmin) );
	    	double LN_Source = Math.log( SOURCEmin ) + rand99;	
	        SOURCEvalue = Math.exp( LN_Source );
	    	if( SOURCEvalue > conc_min  &&  SOURCEvalue < conc_max ){
	    	    looking99 = false;
	    	}
	    }
	    boolean looking = true;
	    while( looking ) {
		double rand1 = Math.random() * ( Math.log(SourceProduct_max) - Math.log(SourceProduct_min) );
		double LN_SourceProduct = Math.log( SourceProduct_min ) + rand1;	
		double SourceProduct = Math.exp( LN_SourceProduct );
		PRODUCTvalue = SOURCEvalue / SourceProduct;
		if( PRODUCTvalue > conc_min  &&  PRODUCTvalue < conc_max ){
		    looking = false;
		}
	    }
	    //  ADPvalue
	    boolean looking2 = true;
	    while( looking2 ) {
		double rand3 = Math.random() * ( Math.log(ATPADP_max) - Math.log(ATPADP_min) );
		double LN_ATPADP = Math.log( ATPADP_min ) + rand3;
		double ATPADP = Math.exp( LN_ATPADP );
		ADPvalue = ATPvalue / ATPADP;
		if( ADPvalue > conc_min  &&  ADPvalue < conc_max ){
		    looking2 = false;
		}
	    }
	    //  NADHvalue ------------------  change this to like HPO4 and CO2?
	    boolean looking3 = true;
	    while( looking3 ) {
		double rand5 = Math.random() * ( Math.log(NADNADH_max) - Math.log(NADNADH_min) );
		double LN_NADNADH = Math.log(NADNADH_min) + rand5;
		double NADNADH = Math.exp( LN_NADNADH );
		NADHvalue = NADvalue / NADNADH;
		if( NADHvalue > conc_min  &&  NADHvalue < conc_max ){
		    looking3 = false;
		}
	    }
	    //  HPO4value
	    double rand7 = Math.random() * ( Math.log(HPO4_max) - Math.log(HPO4_min) );
	    double LN_HPO4 = Math.log(HPO4_min) + rand7;
	    HPO4value = Math.exp( LN_HPO4 );	
	    //  CO2value
	    double rand8 = Math.random() * ( Math.log(CO2_max) - Math.log(CO2_min) );
	    double LN_CO2 = Math.log(CO2_min) + rand8;
	    CO2value = Math.exp( LN_CO2 );


	    //  DO AMP AND PPi  stuff herre --- DECIDE HOW WE SHOULD BE SAMPLING THEM????
	    //  PPivalue
	    double rand9 = Math.random() * ( Math.log(PPi_max) - Math.log(PPi_min) );
	    double LN_PPi = Math.log(PPi_min) + rand9;
	    PPivalue = Math.exp( LN_PPi );	
	    //  AMPvalue
	    double rand10 = Math.random() * ( Math.log(AMP_max) - Math.log(AMP_min) );
	    double LN_AMP = Math.log(AMP_min) + rand10;
	    AMPvalue = Math.exp( LN_AMP );


	    // NEW nitrogen compounds -----
	    double rand11 = Math.random() * ( Math.log(NH3_max) - Math.log(NH3_min) );
	    double LN_NH3 = Math.log(NH3_min) + rand11;
	    NH3value = Math.exp( LN_NH3 );	
	    
	    double rand12 = Math.random() * ( Math.log(GLUT_max) - Math.log(GLUT_min) );
	    double LN_GLUT = Math.log(GLUT_min) + rand12;
	    GLUTvalue = Math.exp( LN_GLUT );

	    double rand13 = Math.random() * ( Math.log(OXO_max) - Math.log(OXO_min) );
	    double LN_OXO = Math.log(OXO_min) + rand13;
	    OXOvalue = Math.exp( LN_OXO );


	    //----------------------------------------------------------------------------


	    // Fix source and prod conc
	    sourceConc.set(0,SOURCEvalue); 
	    prodConc.set(0,PRODUCTvalue);
	    
	    // Set extMetFactors
	    extMetabConcentrations.put(  "-atp", ATPvalue  );
	    extMetabConcentrations.put(  "-adp", ADPvalue  );
	    extMetabConcentrations.put(  "-amp", AMPvalue  );
	    extMetabConcentrations.put(  "-nad", NADvalue  );
	    extMetabConcentrations.put(  "-nadh", NADHvalue  );
	    extMetabConcentrations.put(  "-pi", HPO4value );
	    extMetabConcentrations.put(  "-ppi", PPivalue );
	    extMetabConcentrations.put(  "-co2", CO2value  );
	    //
	    extMetabConcentrations.put(  "-nh3", NH3value );
	    extMetabConcentrations.put(  "-glut", GLUTvalue );
	    extMetabConcentrations.put(  "-oxo", OXOvalue  );


	    //Need an array of fluxes (doubles) corresponding to each path.
	    double[] fluxes = new double[ listOfPaths.size() ];

	    double[] nonOptimized_fluxes = new double[  listOfPaths.size()  ];


	    
	    for( int z = 0;   z < listOfPaths.size();   z++ ){	    
		

		int sizeOfPath = listOfPaths.get(z).size();


		// Convert Reactions3 list to a list of "Reaction" (need this for Eqn Solver)
		//
		// -------------   SHOULD DO THIS AT THE START ONCE ONLY ------------
		ArrayList<Reaction> rxn_list = new ArrayList<Reaction>();
		for( int i=0; i < sizeOfPath; i++ ){
		    Reaction r = Convert_Reactions3_to_Reaction.convert(  listOfPaths.get(z).get(i)   );
		    rxn_list.add( r );
		}
		//  -------------------------------------------------


		// Set up sources and product concentrations in Double lists:  -------  NOT GENERAL!!
		ArrayList<Double> sourcesConc = new ArrayList<Double>();
		ArrayList<Double> productsConc = new ArrayList<Double>();
		sourcesConc.add(  SOURCEvalue  );
		productsConc.add(  PRODUCTvalue  );
		//
		ArrayList<String> nodeLabels = new ArrayList<String>();
		ArrayList<Double> nodeConcentrations = new ArrayList<Double>();


		// -----  STUFF NEEDED FOR PERFECT ENZYME CALCULATION
		//Enzyme distribution

		double totalEnzyme = 1.0E-1;

		double[] Ei = new double[  sizeOfPath   ];
		// flat distribution
		double enz = (totalEnzyme / sizeOfPath);
		for( int e=0; e < Ei.length; e++ ){
		    Ei[ e ] = enz;
		}

		// diffusion limited constant
		double kd = 1.0E8;



		// ========= Optimize enzymes for maximal flux ===============================================================================	
		// Number of variables, nvars, is one less than the number of reactions since we fix [E]_tot
		int nvars = (sizeOfPath - 1);
		// # of starting points in Powell method
		int startPoints = 1;	
		// varmin and varmax are only used to initialize the procedure
		double varmin = enz;   double varmax = enz;     

		fluxes[z] = PowellOptimization_relativeTolerance.find_maximum( nvars, varmin, varmax, rxn_list, extMetabConcentrations, sourcesConc.get(0), productsConc.get(0), totalEnzyme, conc_min, conc_max, kd, RT, startPoints );

		/** This powell optimization is quite problematic. I can manually choose the "settings" to optimize a specific path, under a
		    specific condition, but do not understand the procedure deeply enough to be able to optimize in a general manner...
		    Below is a quick fix which seems to solve the problem for all cases I investigated.  Main issue seems to arise from choice of
		    starting point: starting from a flat profile often finds a (global?) maximum which is hard to improve upon, but this isnt useful if the
		    the pathway is not feasible with a flat profile to begin with. If this is the case, only then do we try the full Powell Method below. 
		    Too big a starting range often makes it impossible to converge on solution; starting from flat avoids this and leads to extremely faster
		    computation times and it appears a robust choice in all test pathways I've investigated.
		 **/
		if( fluxes[z] < 0 ){
		    // Then try optimize from different starting points
		    //System.out.println("1");
		    varmin = totalEnzyme/10000.0;  varmax = totalEnzyme/100.0;
		    startPoints = 100;
		    fluxes[z] = PowellOptimization_relativeTolerance.find_maximum( nvars, varmin, varmax, rxn_list, extMetabConcentrations, sourcesConc.get(0), productsConc.get(0), totalEnzyme, conc_min, conc_max, kd, RT, startPoints );
		}
		if( fluxes[z] < 0 ){
		    // Then try larger range still
		    //System.out.println("2");
		    varmin = totalEnzyme/100000.0;  varmax = totalEnzyme/10.0;
		    startPoints = 100;
		    fluxes[z] = PowellOptimization_relativeTolerance.find_maximum( nvars, varmin, varmax, rxn_list, extMetabConcentrations, sourcesConc.get(0), productsConc.get(0), totalEnzyme, conc_min, conc_max, kd, RT, startPoints );
		}
		/** If still not feasible with this range, just assume it's not of biological interest...  **/
		// ===========================================================================================================================

		
	    }




	    	  	    
	    out.print( SOURCEvalue + "\t" + PRODUCTvalue + "\t" + ATPvalue + "\t" +
		       ADPvalue + "\t" + AMPvalue + "\t" + NADvalue + "\t" + NADHvalue + "\t" + HPO4value + "\t"
		       + PPivalue + "\t"  + CO2value + "\t" + NH3value + "\t" + GLUTvalue + "\t" + OXOvalue );
	    for( int qaz=0; qaz < fluxes.length; qaz++){
		out.print( "\t" + fluxes[qaz] );
	    }
	    out.println();

	    
	}		
	   	   
	out.close();
 



	System.out.println();
	//Prints out run time
	long stopTime = System.currentTimeMillis();
	long elapsedTime = stopTime - startTime;
	System.out.println();
	System.out.println("Runtime = " + elapsedTime/1000.0 + " seconds");
	System.out.println();





    }// end main



















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
    
    




   
}//end class
