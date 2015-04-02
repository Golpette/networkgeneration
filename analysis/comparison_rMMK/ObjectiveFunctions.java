import java.util.*;


public class ObjectiveFunctions{



    // Method to calculate steady state ATP production (consumption)
    public static double getATPproduction(ArrayList<Reaction> reactions,  ArrayList<Double> fluxes){
	// CAREFUL - stupid link; reactions.get(i) corresponds to fluxes.get(i). Fix this.

	double  ATPflux=0;

	/** First need to reorder all Reactions wrt flux  **/
	ArrayList<Reaction> orderedReactions = ReactionReorder.orderReactionsWithFlux( reactions,  fluxes );
	/** and then make all fluxes positive **/
	ArrayList<Double> orderedFluxes = new ArrayList<Double>();
	for( int i=0; i < fluxes.size(); i++ ){
	    orderedFluxes.add(  Math.abs(new Double(fluxes.get(i))  ));
	}


	for( int rxn=0;  rxn < orderedReactions.size();  rxn++ ){
	    /**  See if reaction produces or consumes ATP  **/
	    int numATP = orderedReactions.get(rxn).getATP();
	    ATPflux = ATPflux + ( numATP * orderedFluxes.get( rxn )  );
	}
	return ATPflux;
    }// ================================================================================================







    // Method to calculate steady state NADH production (consumption) -  1 NADH molecule is potentially equivalent to 2.5 ATP molecules! 
    public static double getNADHproduction(ArrayList<Reaction> reactions,  ArrayList<Double> fluxes){
	// CAREFUL - stupid link; reactions.get(i) corresponds to fluxes.get(i). Fix this.

	double  fluxNADH=0;

	/** First need to reorder all Reactions wrt flux  **/
	ArrayList<Reaction> orderedReactions = ReactionReorder.orderReactionsWithFlux( reactions,  fluxes );
	/** and then make all fluxes positive **/
	ArrayList<Double> orderedFluxes = new ArrayList<Double>();
	for( int i=0; i < fluxes.size(); i++ ){
	    orderedFluxes.add(  Math.abs(new Double(fluxes.get(i))  ));
	}

	/** Go through all reactions and sum up number of NADH produced times flux of that reaction **/
	for( int rxn=0;  rxn < orderedReactions.size();  rxn++ ){
	    /**  See if reaction produces or consumes ATP  **/
	    int numNADH = orderedReactions.get(rxn).getNADH();
	    fluxNADH = fluxNADH + ( numNADH * orderedFluxes.get( rxn )  );
	}
	return fluxNADH;
    }// ================================================================================================








    // Do i still use this method?


    //Method to check the flux into each of the 'product' metabolites. (This is how we can tell if any of them have become disconnected - 
    //    their flux will fall below FLUX_CUTOFF). It is likely quicker than checking connectivity toplogically...
    // Can easily alter this to return the flux of all products in some objective function too
    public static boolean checkProductFluxes( ArrayList<Reaction> reactions, ArrayList<Double> fluxes, ArrayList<String> products, double FLUX_CUTOFF ){

	/**  NOTE THAT 'reactions' AND 'fluxes' LIST ELEMENTS MUST MATCH  **/
	if( reactions.size() != fluxes.size() ){ System.out.println( "ERROR: checkProductFluxes();  'reactions' and 'fluxes' do not match" );  }

	boolean allProductsConnected = true;

	/** First need to order the reactions wrt flux  **/
	ArrayList<Reaction> orderedReactions = ReactionReorder.orderReactionsWithFlux( reactions,  fluxes );
	/** and then make all fluxes positive **/
	ArrayList<Double> orderedFluxes = new ArrayList<Double>();
	for( int i=0; i < fluxes.size(); i++ ){
	    orderedFluxes.add(  Math.abs(new Double(fluxes.get(i))  ));
	}

	for( int p=0; p < products.size(); p++ ){

	    String prod = products.get( p );
	    double prodFlux = 0;

	    for( int r=0;  r<orderedReactions.size();  r++  ){

		if( orderedReactions.get( r ).int_subs.contains( prod ) ){
		    prodFlux = prodFlux - orderedFluxes.get( r );
		}
		else if( orderedReactions.get( r ).int_prods.contains( prod ) ){
		    prodFlux = prodFlux + orderedFluxes.get( r );
		}		

	    }
	    if( prodFlux < FLUX_CUTOFF){
		/** Then this product has become disconnected from network  **/   // IS THIS DEFINITELY TRUE?? ALWAYS??
		allProductsConnected = false;
	    }

	}
	return allProductsConnected;
    }
    //================================================================================================










}//end class
