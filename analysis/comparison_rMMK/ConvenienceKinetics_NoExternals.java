import java.util.*;
import java.lang.Math;



public class ConvenienceKinetics_NoExternals{



    private static final double RT = 2.5;   //(kJ/mol)




    // Method to make Nodes but not bother instantiating listOf_K or listOf_Links 
    public static ArrayList<Node> makeNodes( ArrayList<Reaction> reactions ){
	
	ArrayList<Node> nodes_list = new ArrayList<Node>();	
	
	/**  Get all metabolites and make a Node object for them  **/	
	ArrayList<String> allMetabs = Methods.getAllMetabs( reactions );
	for( String metab: allMetabs){
	    
	    Node newnode = new Node( metab );
	    nodes_list.add(  newnode  );		
	}	
	
	return nodes_list;
    }// End makeNodes() ==================================================================================================== 







    // Method to take list of reactions and all steady state concentrations to calculate flux through each reaction
    public static ArrayList<Double> getAllReactionFluxes( ArrayList<Reaction> network, ArrayList<String> nodeLabels, ArrayList<Double> nodeConcs, HashMap<String, Double> externalConcentrations, double sensitivity ){

	//  If any Reaction flux is < sensitivity (or FLUX_CUTOFF), we set the flux value to zero.

	if( nodeLabels.size() != nodeConcs.size() ){ System.out.println("ERROR 7: in ConvenienceKinetics_NoExternals.getAllReactionFluxes()");  }

	ArrayList<Double> fluxes = new ArrayList<Double>();

	/** For each Reaction in network, calculate the flux, using the **steady-state** concentrations in "NodeConcs"  **/
	for( Reaction rxn: network ){

	    // ( Originally we had "Vf" which is actually k_cat*[Ei]: )
	    double kcat = rxn.getkcat();
	    double Ei = rxn.getEi();
	    double Vf = kcat * Ei;  

	    double Ks = rxn.getKs();
	    double Kp = rxn.getKp();

	    double rxnFlux = Double.NaN;

	    ArrayList<String> int_subs = rxn.getIntSubs();
	    ArrayList<String> ext_subs = rxn.getExtSubs();
	    ArrayList<String> int_prods = rxn.getIntProds();
	    ArrayList<String> ext_prods = rxn.getExtProds();

	    // "Correct" equilibrium const q to account for external metabolites
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


	    // IS THIS CORRECT!??! SHOULD I BE SCALING THE EQUILIBRIUM CONSTANT OR WILL THIS HAPPEN ********************************
	    // AUTOMATICALLY IN THE RIGHT HAND BRACKETS OF THE RATE EQUATION!??!?!                  ******************************** 
	    /** NB: no external metabolites in rate expression - NEED to rescale equilibrium constant **/
	    double dG = rxn.getDeltaG();
	    double q = (multAllExtSubs / multAllExtProds) * Math.exp( -(dG / RT)  );


	    /** Multiply all substrate (internal)  **Michaelis-Menten constants**  **/
	    double denom_1 = 1.0;
	    for( String s: int_subs ){
		denom_1 = denom_1*Ks;
	    }
	    

	    double allSubstratesMultiplied = 1.0;
	    /** Multiply all terms like (1+c123/Ks) for internal and external SUBSTRATES **/
	    double denom_2a = 1.0;
	    for( String s: int_subs ){

		int index = nodeLabels.indexOf( s );
		double conc = nodeConcs.get( index );

		allSubstratesMultiplied = allSubstratesMultiplied * conc; 

		double bit = 1.0 + ( conc / Ks );
	        denom_2a = denom_2a * bit;
	    }

	    double allProductsMultiplied = 1.0;
	    /** Multiply all terms like (1+c234/Kp) for all internal and external PRODUCTS **/
	    double denom_2b = 1.0;
	    for( String s: int_prods ){

		int index = nodeLabels.indexOf( s );
		double conc = nodeConcs.get( index );

		allProductsMultiplied = allProductsMultiplied * conc;

		double bit = 1.0 + ( conc / Kp );

	        denom_2b = denom_2b * bit;
	    }
	    
	    /**   String k_forw = "(Vf"+r+"/("+denom_1+"))*(1 / ( " +denom_2a + " + " + denom_2b + " - 1 ) )";  **/
	    double k_forw = ( Vf/denom_1 )*( 1.0 / ( denom_2a + denom_2b - 1.0  )   );
	    
	    /** Get Reaction flux: **/
	    rxnFlux = k_forw * ( allSubstratesMultiplied - (allProductsMultiplied / q )    );
	    	   

	    /** If absolute size of flux is less than sensitivity, assume it is an equilibrium reaction
	        and set its flux to zero -- this may not be good enough... **/
	    if(  Math.abs(rxnFlux) < sensitivity ){
		rxnFlux = 0;
	    }
	    fluxes.add(  rxnFlux  );

	}
	return fluxes;
	
    } // End of getReactionFluxes_UPDATEDFOREXTERNALS method
    //=========================================================================================================================











    // Get index of given metabolite in list of Nodes:
    public static int getIndex(ArrayList<Node> nodes, String label){
	int index = 999777999;
	for(int i=0; i<nodes.size(); i++){
	    if( label.equals( nodes.get(i).getLabel() ) ){
		index = i;
	    }
	}
	if( index == 999777999 ){
	    System.out.println("MMK.java -- index of Node label == 999777999; i.e. not in nodes" );
	}
	return index;
    }
    //====================================================================================


    // Checks list of Node objects to see if a given metabolite is already present in list
    public static boolean containsNode(ArrayList<Node> nodes, String label){
	boolean contains = false;
	for(int i=0; i<nodes.size(); i++){
	    if(  label.equals( nodes.get(i).getLabel() )  ){
		contains = true;	
	    }
	}
	return contains;
    }
    //====================================================================================




}