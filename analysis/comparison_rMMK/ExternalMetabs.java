import java.util.*;

// Class holds some useful methods

public class ExternalMetabs{
    


    /** Method to get the correction factor to the equilibrium constant due to external metabolite concentrations  **/
    public static double getExtMetFactor( Reaction rxn, double ATP, double ADP, double AMP, double NAD, double NADH, double HPO4, double PPi, double CO2){
	
	ArrayList<String> ext_subs = rxn.getExtSubs();
	ArrayList<String> ext_prods = rxn.getExtProds();

	// System.out.println( rxn );
	// System.out.println( ext_subs.toString() );
	// System.out.println( ext_prods.toString() );
	// System.out.println();
	
	
	//Method to get the multiplicative factor of the external metabolites which alter
	//the equilibrium constant from Bartek's value.
	double extMetFactor = 0;
	double numerator = 1;  //product of all external metabolites USED in reaction
	for( int i=0; i<ext_subs.size(); i++){
	    double mult = Double.NaN;
	    if( ext_subs.get(i).equals("-1") ){ mult = 1.0;  } //H2O=1.0
	    else if( ext_subs.get(i).equals("-2") ){ mult = CO2;  } 
	    else if( ext_subs.get(i).equals("-3") ){ mult = NAD;  } 
	    else if( ext_subs.get(i).equals("-4") ){ mult = NADH;  } 
	    else if( ext_subs.get(i).equals("-5") ){ mult = AMP;  } 
	    else if( ext_subs.get(i).equals("-6") ){ mult = ADP;  } 
	    else if( ext_subs.get(i).equals("-7") ){ mult = ATP;  } 
	    else if( ext_subs.get(i).equals("-8") ){ mult = HPO4;  } 
	    else if( ext_subs.get(i).equals("-9") ){ mult = PPi;  }
	    else{ System.out.println("Unknown external metabolite = " + ext_subs.get(i)  ); System.exit(0); }

      	    numerator = numerator * mult;  
	}
	double denominator = 1;  //product of all external metabolites PRODUCED in reaction
	for( int i=0; i<ext_prods.size(); i++){
	    double mult2 = Double.NaN;
	    if( ext_prods.get(i).equals("-1") ){ mult2 = 1.0;  } //H2O=1.0
	    else if( ext_prods.get(i).equals("-2") ){ mult2 = CO2;  } 
	    else if( ext_prods.get(i).equals("-3") ){ mult2 = NAD;  } 
	    else if( ext_prods.get(i).equals("-4") ){ mult2 = NADH;  } 
	    else if( ext_prods.get(i).equals("-5") ){ mult2 = AMP;  } 
	    else if( ext_prods.get(i).equals("-6") ){ mult2 = ADP;  } 
	    else if( ext_prods.get(i).equals("-7") ){ mult2 = ATP;  } 
	    else if( ext_prods.get(i).equals("-8") ){ mult2 = HPO4;  } 
	    else if( ext_prods.get(i).equals("-9") ){ mult2 = PPi;  }
	    else{ System.out.println("Unknown external metabolite = " + ext_prods.get(i)  ); System.exit(0); }

            denominator = denominator * mult2;  
	}
	extMetFactor = numerator / denominator;
	return extMetFactor;
    }
    //====================================================================================








}