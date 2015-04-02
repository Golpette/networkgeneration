import java.util.*;

public class PerfectEnzyme{

    /** 
        Methods needed to calculate flux and steady state concentrations of intermediate metabolites assuming
        Heinrich1997 'perfect enzyme' calculation and a specific enzyme distribution, Ei[].
    **/



    // NOTE:  This method calculates the pathway flux;  NOT the ATP flux!

    public static double getFlux( ArrayList<Reaction> reactionsUsed, HashMap<String, Double> externalMetabs, double concSOURCE, double concPRODUCT, double[] Ei, double kd, double RT ){
    
	double flux=0;

	//Arrays to hold our Qi and Ki values for each reaction
	double[] Q = new double[ reactionsUsed.size() ];
	double[] K = new double[ reactionsUsed.size() ];
	

	//Get all Q values for each step in pathway
	for(int i=0; i < reactionsUsed.size(); i++){	    

	    double reactionG = 0;
	    double q = 0;

	    /** Note: reactions being input to this comparison are already ordered in the
		correct direction  **/	    
	    reactionG =  reactionsUsed.get(i).getDeltaG();
	    q = Math.exp( (-1.0)*reactionG / (RT)   );

	    // correct equilibrium constant to include external metabolite concentrations 
	    double extMetaboliteCorrection = getExtMetFactor( reactionsUsed.get(i), externalMetabs );

	    Q[i] = extMetaboliteCorrection * q ;

	    K[i] = ( kd * Ei[i] * Q[i] ) / ( 1 + Q[i] );	
	
	}// end for loop. Now have array for Q values
	


	//  This is the flux calculation when the enzyme distribution HAS NOT been optimized
    
	//Find D (see notes)  -- this comes from Eqns(5a,b) and Eqns (16a,b) -----------
	double D = 0;   double multQ=1;
	for(int aa = 0; aa < K.length; aa++){
	    multQ = 1;
	
	    for(int bb = aa; bb < Q.length; bb++){
		multQ = multQ * Q[bb];
	    }
	    D = D + ( 1.0/K[aa] ) * multQ;
	}
    
	//Find J (flux)
	double all_Q_multipled = 1;
	for(int cc = 0; cc < Q.length; cc++){
	    all_Q_multipled = all_Q_multipled * Q[cc];
	}
    
	flux = ( 1.0/D ) * (   ( concSOURCE * all_Q_multipled) - concPRODUCT   );

	return flux;	
    }
    //________________________________________________________________________________________________________________________________
    








    // This method is for calculating the concentrations of all metabolites  

    public static ArrayList<Double> getConcentrations( double flux, ArrayList<Reaction> reactionsUsed, HashMap<String, Double> externalMetabs, double concSOURCE, double concPRODUCT, double[] Ei, double kd, double RT ){

   
	//	ArrayList<Double> concentrations = new ArrayList<Double>(); 
	double[] concentrations = new double[ reactionsUsed.size() + 1 ];

	//Arrays to hold our Qi and Ki values for each reaction
	double[] Q = new double[ reactionsUsed.size() ];
	double[] K = new double[ reactionsUsed.size() ]; 
	
	
	//Get all Q AND K values for each step in pathway ....... this is all repeated from above.
	for(int i=0; i < reactionsUsed.size(); i++){	    

	    double reactionG = 0;
	    double q = 0;

	    /** Note: reactions being input to this comparison are already ordered in the
		correct direction  **/	    
	    reactionG =  reactionsUsed.get(i).getDeltaG();
	    q = Math.exp( (-1.0)*reactionG / (RT)   );

	    // correct equilibrium constant to include external metabolite concentrations 
	    double extMetaboliteCorrection = getExtMetFactor( reactionsUsed.get(i), externalMetabs );

	    Q[i] = extMetaboliteCorrection * q ;

	    K[i] = ( kd * Ei[i] * Q[i] ) / ( 1 + Q[i] );

	}


	//--------------Calculate concentrations via Heinrich equation:----------------
	//
	//  [S_j] =  [S_0].Prod_i=1^j{Q_i}  -   J.Sum_i=1^j{ 1/K_i .Prod_m=i^j{ Q_m }  }
	//        =            A            -        B
	//------------------------------------------------------------------------------

        // Source and product concentrations are fixed
	concentrations[0] = concSOURCE;
	concentrations[ (concentrations.length) - 1 ] = concPRODUCT;
	
        // Calculate steady state intermediate concentrations
	for( int j=1; j < (concentrations.length - 1); j++){	  
	    
	    //Product of Q's up to j is:
	    double Q_mult_j = 1;
	    for(int cc = 0; cc < j; cc++){
		Q_mult_j = Q_mult_j * Q[cc];
	    }
	    
	    double A = concSOURCE * Q_mult_j;


	    double d = 0;   double multQ=1;
	    for(int aa = 0;  aa < j;  aa++){       
		multQ=1;		

		for(int bb = aa; bb < j; bb++){
		    multQ = multQ * Q[bb];
		}
		d = d + ( 1.0/K[aa] ) * multQ;
	    }
	    
	    double B = flux * d;

	    concentrations[j] = A - B;
	}
	

	// Convert to ArrayList
	ArrayList<Double> metab_concentrations = new ArrayList<Double>();
	for( int c=0; c < concentrations.length; c++ ){
	    metab_concentrations.add(  concentrations[ c ]  );
	}

	return metab_concentrations;	
    }

    







   // Method to get equilibirum constant correction factor due to external metabolite concentrations 
    public static double getExtMetFactor( Reaction rxn, HashMap<String, Double> externalConcentrations ){

	double factor = 0;
       		
	// We are NOT INCLUDING the kinetics (no Km terms or concentrations of the external metabolites
	// Instead, they only enter through the modification of the equilibrium constant
	ArrayList<String> ext_subs = rxn.getExtSubs();
	ArrayList<String> ext_prods = rxn.getExtProds();
	
	/** Check all external metabolites of this Reaction are actually defined  **/
	for( String s: ext_subs ){
	    if(  !externalConcentrations.containsKey( s ) ){
		System.out.println("ERROR 22: reaction contains unknown external metabolite in PerfectEnzyme.java. Exiting. ");
		System.out.println(" --  " + s);
		System.exit(0);
	    }
	}
	for( String s: ext_prods ){
	    if(  !externalConcentrations.containsKey( s ) ){
		System.out.println("ERROR 22: reaction contains unknown external metabolite in PerfectEnzyme.java. Exiting. ");
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

	
	factor = multAllExtSubs/multAllExtProds;

	return factor;
    }











}
