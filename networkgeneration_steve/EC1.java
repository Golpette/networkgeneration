import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;


public class EC1 {
	
	/**
	 *  EC 1
	 *   1.1.1: {CH2(OH), CH(OH)=, -CH(OH)} + nad --> {CHO, CO=, -CO-} + nadh     Implemented
	 *          R-CH(OH)-r'-COOH + nad --> R-CO-r'-H + co2 + nadh                 Implemented        
	 *          ( In EC examples r' is always a single carbon group )
	 *          
	 *   1.2.1: R-CHO + nad + h2o --> R-COOH + nadh                               Implemented
	 *          R-CHO + nad + pi --> R-COp + nadh                                 Implemented
	 *          
	 *   1.3.1: R-CH2-CH2-R' + nad --> R-CH=CH-R' + nadh                          Implemented         
	 *          (should this allow for e.g. -CH(OH)-CH(OH)- --> -C(OH)=C(OH)-  ??)               
	 *          Currently it allows for ALL possible H-atom losses from adjacent carbons.
	 *          This may need to be restricted...
	 */
	
	
	// this method is in EC 5 too - FIX THS
	
	// Method to get Key of Map from given Value (NOTE: assumes mapping is 1-to-1)
	static String getKeyFromValue( Map<String, String> map, String val ){
		String key = "";
		for ( Map.Entry<String, String> entry : map.entrySet()) {
			if ( entry.getValue().equals( val ) ) {
				key = entry.getKey();
			}
		}
		return key;
	}
	
	
	
	
    // Final "H_ADD_MAP" to hold all groups that can have a hydrogen added to them
	// I'm using this in EC 1.3.1 -- Bartek has used a much restricted set of possible 
	// transformations for this class; see his "types" file, think it needs fixed.
    private static final Map<String, String> H_ADD_MAP = createMap();
    private static Map<String, String> createMap() {
       Map<String, String> result = new HashMap<String, String>();
       result.put("CO=", "CHO-");           result.put("=CO", "-CHO");     
       result.put("CH2=", "CH3-");          result.put("=CH2", "-CH3"); 
       result.put("CH(OH)=", "CH2(OH)-");   result.put("=CH(OH)", "-CH2(OH)"); 
       result.put("CHp=", "CH2p-");         result.put("=CHp", "-CH2p"); 
       result.put("-CH=", "-CH2-");         result.put("=CH-", "-CH2-");
       // Note that the next 2 mean that hashmap is not 1-to-1
       result.put("-C(OH)=", "-CH(OH)-");   result.put("=C(OH)-", "-CH(OH)-"); 
       result.put("-Cp=", "-CHp-");         result.put("=Cp-", "-CHp-"); 
       return Collections.unmodifiableMap(result);
   }// (Note this method is specific to this class and double bonds;
    // i.e. it doesn't include "-CH2-" going to "CH3-" for example. See next map:
    
    
    
    
    // Map to lose H from end group making it a molecular "body group"
    private static final Map<String, String> MAKE_BODY = createMap2();
    private static Map<String, String> createMap2() {
       Map<String, String> result = new HashMap<String, String>();
       result.put("-CH3", "-CH2-");          result.put("CH3-", "-CH2-");   
       //
       // can carboxyl be added to any group or only a CH3- , -CH2- etc??  WHAT DOES BAREK DO?
       //
       result.put("-CH2(OH)", "-CH(OH)-");   result.put("CH2(OH)-", "-CH(OH)-");
       result.put("-CHO", "-CO-");           result.put("CHO-", "-CO-");
       result.put("-CH2p", "-CHp-");         result.put("CH2p-", "-CHp-");
       result.put("=CH2", "=CH-");           result.put("CH2=", "-CH=");
       result.put("=CH(OH)", "=C(OH)-");     result.put("CH(OH)=", "-C(OH)=");
       result.put("=CHp", "=Cp-");           result.put("CHp=", "-Cp=");
       return Collections.unmodifiableMap(result);
   }    
	

    
    // Mapping hydroxyl red-ox pars 
    private static final Map<String, String> REDOX_HYDROXY = createMap3();
    private static Map<String, String> createMap3() {
       Map<String, String> result = new HashMap<String, String>();
       result.put("CH2(OH)-", "CHO-");    result.put("-CH2(OH)", "-CHO");
	   result.put("CH(OH)=", "CO=");      result.put("=CH(OH)", "=CO");    
		// above oxidation of enol possible? unlikely to occur in biology anyway... 
	   result.put("-CH(OH)-", "-CO-");
       return Collections.unmodifiableMap(result);
   } 
    
    
    
    



	public static void _1_1( String sub, ArrayList<String> rxns, HashMap<String,Integer> map  ){
		
		
		ArrayList<String> grps = ReactionMechanisms.getGroups( sub );
		
		
		// ======================     C(OH) -> CO    ============================= 
		// (primary alcohol to aldehyde OR secondary alcohol to ketone)
		for( int g=0; g<grps.size(); g++ ){
			String group = grps.get( g );
			
			if( REDOX_HYDROXY.containsKey( group ) ){
				ArrayList<String> new_grps = new ArrayList<String>( grps );
				new_grps.set( g, REDOX_HYDROXY.get( group ) );	
				String tempProd = ReactionMechanisms.formProd( new_grps );  				
				String prod = "";
				// Check product or its palindrome is in map
				if( map.containsKey( tempProd )  ){
					prod = tempProd;
				}
				else if( map.containsKey(  Methods.getMolecularPalindrome( tempProd )   )){
					prod = Methods.getMolecularPalindrome( tempProd );
				}
				//
				if( prod != "" ){
					String rxn_string = sub+" + nad --> "+prod+" + nadh";
					rxns.add( "1.1.1"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
				}
				else{
					System.out.println("Error: 1_1_1 product not in compound list. Sub: "+sub);
				}			
			}
		} // =====================================================================
		
		
		
		
		// *** COUPLED DECAROXYLATION ***
		// ================== R-CH(OH)-r'-COOH + nad --> R-CO-r'-H + nadh ==========================
		//  The  r' is only ever a 1-carbon group in EC classes. Read about mechanism. Is it ever NOTHING?	
		// Coding this as a "carboxylation" reaction was easier

		int numgrps = grps.size();
		String firstgroup = grps.get( 0 );
		String lastgrp = grps.get( numgrps - 1 );
		
		if( grps.size() >= 2 ){  // POSSIBLE BUG HERE??
			
			if( MAKE_BODY.containsKey( firstgroup )  &&  REDOX_HYDROXY.containsValue( grps.get(1) )   ){  
		    // is this sensible getting key from value? 

				//get key for REDOX_HYDROXY replacement
				String redox_rep = getKeyFromValue(  REDOX_HYDROXY, grps.get(1)  );
				
				ArrayList<String> new_grps = new ArrayList<String>( grps );
				// carboxylate
				new_grps.add(0,"COOH-");
				// remove hydrogen from old end group
				new_grps.set(1, MAKE_BODY.get( firstgroup ) );
				// reduce the carbonyl
				new_grps.set(2, redox_rep );
				
				// form new molecule, check its in our list
				String tempProd = ReactionMechanisms.formProd( new_grps );
				
				String prod = "";
				// Check product or its palindrome is in map
				if( map.containsKey( tempProd )  ){
					prod = tempProd;
				}
				else if( map.containsKey(  Methods.getMolecularPalindrome( tempProd )   )){
					prod = Methods.getMolecularPalindrome( tempProd );
				}
				//
				if( prod != "" ){
					// PRINT RXN STRING AS A OXIDATION + DECARBOXYLATION 
					String rxn_string = prod+" + nad ---> "+sub+" + co2(aq) + nadh ";
					rxns.add( "1.1.1"+"\t"+map.get( prod )+"\t"+map.get( sub )+"\t"+rxn_string  );
				}
				else{
					System.out.println("Error_1_1_1: compound not in compound list. "+sub);
				}
				
			}
			
			
			if( MAKE_BODY.containsKey( lastgrp )  &&  REDOX_HYDROXY.containsValue(  grps.get( numgrps-2 )   ) ){   
				
				//get key for REDOX_HYDROXY replacement
				String redox_rep = getKeyFromValue(  REDOX_HYDROXY, grps.get( numgrps-2 )  );
				
				
				ArrayList<String> new_grps = new ArrayList<String>( grps );
				// carboxylate on end
				new_grps.add("-COOH");
				int ng = new_grps.size();
				// remove hydrogen from old end group
				new_grps.set( ng-2, MAKE_BODY.get( lastgrp ) );
				// reduce the carbonyl
				new_grps.set( ng-3, redox_rep );
					
				// form new molecule, check its in our list
				String tempProd = ReactionMechanisms.formProd( new_grps );  
										
				String prod = "";
				// Check product or its palindrome is in map
				if( map.containsKey( tempProd )  ){
					prod = tempProd;
				}
				else if( map.containsKey(  Methods.getMolecularPalindrome( tempProd )   )){
					prod = Methods.getMolecularPalindrome( tempProd );
				}
				//
				if( prod != "" ){
					// PRINT RXN STRING AS A OXIDATION + DECARBOXYLATION 
					String rxn_string = prod+" + nad ---> "+sub+" + co2(aq) + nadh ";
					rxns.add( "1.1.1"+"\t"+map.get( prod )+"\t"+map.get( sub )+"\t"+rxn_string  );
				}
				else{
					System.out.println("ErrorB_1_1_1: compound not in compound list. "+sub);
				}
							
				
			}
		}//=======================================================================================
		

	}
	// ##########################   END EC 1.1.1    ########################


	
	
	
	
	
	
	
	
	

	public static void _2_1( String sub, ArrayList<String> rxns, HashMap<String,Integer> map  ){

		//=======  R-CHO + nad + h2o --> R-COOH + nadh  ==================
		if( sub.contains("CHO-") ){
			String prod = sub.replaceAll( "CHO-", "COOH-" );
			// Check product or its palindrome is in map
			if( map.containsKey(prod) ){
				String rxn_string = sub+" + nad + h2o --> "+prod+" + nadh";
				rxns.add( "1.2.1"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
			}
			else if( map.containsKey( Methods.getMolecularPalindrome(prod)  ) ){
				String prod2 = Methods.getMolecularPalindrome(prod);
				String rxn_string = sub+" + nad + h2o --> "+prod2+" + nadh";
				rxns.add( "1.2.1"+"\t"+map.get(sub)+"\t"+map.get(prod2)+"\t"+rxn_string  );				
			}
			else{
				System.out.println("ErrorC: 1_2_1 product not in compound list: " + prod);
//				System.out.println("Tried to oxidize CHO in "+ sub);
				//System.exit(0);
			}
		}
		// or one the other end:
		if( sub.contains("-CHO") ){
			String prod = sub.replaceAll( "-CHO", "-COOH" );
			if( map.containsKey(prod) ){
				String rxn_string = sub+" + nad + h2o --> "+prod+" + nadh";
				rxns.add( "1.2.1"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
			}
			else if( map.containsKey( Methods.getMolecularPalindrome(prod)  ) ){
				String prod2 = Methods.getMolecularPalindrome(prod);
				String rxn_string = sub+" + nad + h2o --> "+prod2+" + nadh";
				rxns.add( "1.2.1"+"\t"+map.get(sub)+"\t"+map.get(prod2)+"\t"+rxn_string  );				
			}
			else{
				System.out.println("ErrorC2: 1_2_1 product not in compound list: " + prod);
			}
		}// ====================================================================


		//=============  R-CHO + nad + pi --> R-COp + nadh  ===================  stupid: same as above...
		if( sub.contains("CHO-") ){
			String prod = sub.replaceAll( "CHO-", "COp-" );
			// Check product or its palindrome is in map
			if( map.containsKey(prod) ){
				String rxn_string = sub+" + nad + pi --> "+prod+" + nadh";
				rxns.add( "1.2.1"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
			}
			else if( map.containsKey( Methods.getMolecularPalindrome(prod)  ) ){
				String prod2 = Methods.getMolecularPalindrome(prod);
				String rxn_string = sub+" + nad + pi --> "+prod2+" + nadh";
				rxns.add( "1.2.1"+"\t"+map.get(sub)+"\t"+map.get(prod2)+"\t"+rxn_string  );				
			}
			else{
				System.out.println("ErrorD: 1_2_1 product not in compound list: " + prod);
			}

		}
		// or other end:
		if( sub.contains("-CHO") ){
			String prod = sub.replaceAll( "-CHO", "-COp" );
			if( map.containsKey(prod) ){
				String rxn_string = sub+" + nad + pi --> "+prod+" + nadh";
				rxns.add( "1.2.1"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
			}
			else if( map.containsKey( Methods.getMolecularPalindrome(prod)  ) ){
				String prod2 = Methods.getMolecularPalindrome(prod);
				String rxn_string = sub+" + nad + pi --> "+prod2+" + nadh";
				rxns.add( "1.2.1"+"\t"+map.get(sub)+"\t"+map.get(prod2)+"\t"+rxn_string  );				
			}
			else{
				System.out.println("ErrorD2: 1_2_1 product not in compound list: " + prod );
			}
		}// =========================================================================
		
	} //##############   END EC 1.2.1   #########################################################


	
	
	
	


	
	
	
	public static void _3_1( String sub, ArrayList<String> rxns, HashMap<String,Integer> map  ){
		
		// ==============    -R=R- --> -R(H)-R(H)-  ===================
		// Note: we do it this way (REDUCTION) since losing the H would 
		// allow for multiple options, i.e. -CH(OH)- --> -C(OH)=  OR  =C(OH)-
		
		ArrayList<String> grps = ReactionMechanisms.getGroups( sub );
		// check all groups and their n.n.
		for( int g=0; g < grps.size()-1 ; g++ ){
			String group = grps.get( g );
			String nxtgroup = grps.get( g+1 );
			
			if( H_ADD_MAP.containsKey( group ) &&  H_ADD_MAP.containsKey( nxtgroup ) ){
				ArrayList<String> new_grps = new ArrayList<String>( grps );
				// replace both groups with mapping
				new_grps.set( g, H_ADD_MAP.get( group ) );
				new_grps.set( g+1, H_ADD_MAP.get( nxtgroup ) );
				String tempProd = ReactionMechanisms.formProd( new_grps );  				
				String prod = "";
				// Check product or its palindrome is in map
				if( map.containsKey( tempProd )  ){
					prod = tempProd;
				}
				else if( map.containsKey(  Methods.getMolecularPalindrome( tempProd )   )){
					prod = Methods.getMolecularPalindrome( tempProd );
				}
				//
				if( prod != "" ){
					// REVERSE STRING SO THAT ITS AN OXIDATION REACTION (AS IS EC 1)
					String rxn_string = prod+" + nad --> "+sub+" + nadh";
					rxns.add( "1.3.1"+"\t"+map.get(prod)+"\t"+map.get(sub)+"\t"+rxn_string  );
					// ABOVE IS REVERSED
				}
				else{
					System.out.println("Error: 1_3_1 reduction product not in compound list. Sub: "+sub);
				}			
			}
		}// ============================================		
	} // END EC 1.3.1 

	
	


}




