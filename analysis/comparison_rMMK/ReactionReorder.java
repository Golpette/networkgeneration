import java.util.*;
import java.io.*;


// Program to read in Bartek's reaction list and  print out the same list but reordered, so that
// all deltaG's are negative

// Also contains a method to take list of Reaction objects, along with a list of the reaction fluxes,
// and create a new list of Reactions, ordered wrt flux direction!


public class ReactionReorder{


    // Method to "flip" a Reaction object
    public static Reaction flipReaction( Reaction rxn ){

	double dG = (-1) * rxn.getDeltaG();
	String eqn = " ";
	
	int intProds = rxn.int_prods.size(); 
	int extProds = rxn.ext_prods.size(); 
	int intSubs = rxn.int_subs.size(); 
	int extSubs = rxn.ext_subs.size();

	for( int p=0; p < intProds; p++ ){
	    eqn = eqn + " " + rxn.int_prods.get( p );
	}
	for( int ep=0; ep < extProds; ep++ ){
	    eqn = eqn + " " + rxn.ext_prods.get( ep );
	}
	eqn = eqn + " >";
	for( int s=0; s < intSubs; s++ ){
	    eqn = eqn + " " + rxn.int_subs.get( s );
	}
	for( int es=0; es < extSubs; es++ ){
	    eqn = eqn + " " +rxn.ext_subs.get( es );
	}

	String eqn2 = eqn.trim();
	Reaction newReaction = new Reaction( dG, eqn2 );

	return newReaction;
    }
    //========================================================================




    // Method to take a list of Reactions and corresponding fluxes, and flip all reactions which have flux < 0
    public static ArrayList<Reaction> orderReactionsWithFlux(ArrayList<Reaction> network,  ArrayList<Double> fluxes ){
	
	if( network.size() != fluxes.size() ){
	    System.out.println("Error: size of \"network\" and \"fluxes\" should match. Element j of one list should correspond to element j of the other. Something's gone wrong in ReactionReorder.orderReactionsWithFlux() ");
	}

	ArrayList<Reaction> orderedReactions = new ArrayList<Reaction>();

	for( int i=0; i < network.size(); i++ ){

	    Reaction rxn = new Reaction( network.get(i) );
	    
	    if( fluxes.get(i) < 0 ){

		Reaction rxn2 = new Reaction(  flipReaction(rxn)   );
		orderedReactions.add(   rxn2    );
	    }
	    else{
		orderedReactions.add(  rxn  );
	    }
	}
	return orderedReactions;
    }
    //================================================================================






}// end class