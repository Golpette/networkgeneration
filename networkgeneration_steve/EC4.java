import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;


public class EC4 {
	/**  
	 *     EC 4 = Lyases
	 *      4.1 = Carbon-carbon lyases
	 *    4.1.1 = Carboxy-lyases
	 *              ==>  R-COOH + h2o --> R-H + co2            (should we add further enzymatic restrictions?)
	 *              ==>  COOH-r-CO-R + pi --> r=Cp-R + co2     (h2o on LHS?)
	 * 
	 *      4.2 = Carbon-oxygen lyases
	 *    4.2.1 = Hydro-lyases
	 *      (A)   ==>  R(H)-(OH)R' --> R=R' + h2o         
	 *      (B)   ==>  R(OH)-(OH)R' --> R(=O)-R' + h2o         (& other way round)
	 *      (C)   ==>  R(OH)-COOH  -->  R(=O)-CHO + h2o
	 */
	
	
	
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
       //
       result.put("=CH2", "=CH-");           result.put("CH2=", "-CH=");
       result.put("=CH(OH)", "=C(OH)-");     result.put("CH(OH)=", "-C(OH)=");
       result.put("=CHp", "=Cp-");           result.put("CHp=", "-Cp=");
       return Collections.unmodifiableMap(result);
   }   
    
    
	
	
	
 	public static void _1_1(String sub, ArrayList<String> rxns, HashMap<String,Integer> map ){
 		// Have altered method to be similar to EC1.1.1 coupled decarbox
 		
		ArrayList<String> grps = ReactionMechanisms.getGroups( sub );
 		
		int numgrps = grps.size();
		String firstgroup = grps.get( 0 );
		String lastgrp = grps.get( numgrps - 1 );
		
		if( grps.size() >= 2 ){
			
			if( MAKE_BODY.containsKey( firstgroup )   ){  
				
				ArrayList<String> new_grps = new ArrayList<String>( grps );
				// carboxylate
				new_grps.add(0,"COOH-");
				// remove hydrogen from old end group
				new_grps.set(1, MAKE_BODY.get( firstgroup ) );				
				
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
					// PRINT RXN STRING AS DECARBOXYLATION 
					String rxn_string = prod+" + h2o ---> "+sub+" + co2(aq)";
					rxns.add( "4.1.1"+"\t"+map.get( prod )+"\t"+map.get( sub )+"\t"+rxn_string  );
				}
				else{
					System.out.println("Error_4_1_1: compound not in compound list. "+sub);
				}
			}
			
			
			if( MAKE_BODY.containsKey( lastgrp ) ){   								
				
				ArrayList<String> new_grps = new ArrayList<String>( grps );
				// carboxylate on end
				new_grps.add("-COOH");
				int ng = new_grps.size();
				// remove hydrogen from old end group
				new_grps.set( ng-2, MAKE_BODY.get( lastgrp ) );
					
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
					// PRINT RXN STRING AS A DECARBOXYLATION 
					String rxn_string = prod+" + h2o ---> "+sub+" + co2(aq)";
					rxns.add( "4.1.1"+"\t"+map.get( prod )+"\t"+map.get( sub )+"\t"+rxn_string  );
				}
				else{
					System.out.println("ErrorB_4_1_1: compound not in compound list. "+sub);
				}
							
				
			}
		}//=======================================================================================

 		
 		
 		

		
		//==============  COOH-r-CO-R + pi --> r=Cp-R + co2 ============================  (+ ATP coupling)
		//       Quite ad hoc
		// Not exactly sure of the mechanism here; seems to go through enol intermediate
		// Implementation of this method may be too specific; need it for PEP carboxylase etc.
		
		HashMap<String,String> END_MAP = new HashMap<String,String>();
		END_MAP.put("CH2=", "-CH2-");   	   END_MAP.put("=CH2", "-CH2-");     
		END_MAP.put("CH(OH)=", "-CH(OH)-");    END_MAP.put("=CH(OH)", "-CH(OH)-");
		//END_MAP.put("CO=", "-CO-");            END_MAP.put("=CO", "-CO-");
		// Bartek manually disallowed ketenes.
		END_MAP.put("CHp=", "-CHp-");          END_MAP.put("=CHp", "-CHp-");
		
		// Map for groups being simultaneously reduced and phosphorylated
		HashMap<String, String> REDUC_PHOS = new HashMap<String,String>();
		REDUC_PHOS.put("=Cp-", "-CO-");        REDUC_PHOS.put("-Cp=", "-CO-");
		REDUC_PHOS.put("=CHp", "-CHO");        REDUC_PHOS.put("CHp=", "CHO-");
		

		if( numgrps >= 2 ){
			
			if( REDUC_PHOS.containsKey( grps.get(1) ) &&  END_MAP.containsKey(firstgroup) ){  
//			if( grps.get(1).equals("=Cp-") &&  END_MAP.containsKey(firstgroup) ){  
				
				// this is too specific? eg should it allow CHO- to CHp= ??? BARTEK ALLOWS THIS
				
				ArrayList<String> new_grps = new ArrayList<String>( grps );

				new_grps.set(0, END_MAP.get( firstgroup )  );
				new_grps.set(1, REDUC_PHOS.get( grps.get(1) ) );
//				new_grps.set(1, "-CO-");
				// carboxylate on end
				new_grps.add(0, "COOH-");
				
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
					// PRINT RXN STRING AS A DECARBOXYLATION 
					// couples ATP
					String rxn_string = prod+" + atp ---> "+sub+" + adp + co2(aq)";
					rxns.add( "4.1.1"+"\t"+map.get( prod )+"\t"+map.get( sub )+"\t"+rxn_string  );
					// couples Pi
					String rxn_string222 = prod+" + pi ---> "+sub+" + co2(aq)";
					rxns.add( "4.1.1"+"\t"+map.get( prod )+"\t"+map.get( sub )+"\t"+rxn_string222  );
				}
				else{
					System.out.println("ErrorQ_4_1_1: compound not in compound list. "+sub);
				}				
			}
			// and reversed:
			if( REDUC_PHOS.containsKey( grps.get( numgrps-2 ) ) &&  END_MAP.containsKey( lastgrp ) ){
//			if( grps.get( numgrps-2 ).equals("-Cp=") &&  END_MAP.containsKey( lastgrp ) ){
				
				ArrayList<String> new_grps = new ArrayList<String>( grps );

				new_grps.set(numgrps-1, END_MAP.get( lastgrp )  );
				new_grps.set(numgrps-2, REDUC_PHOS.get( grps.get( numgrps-2 ) ) );
//				new_grps.set(numgrps-2, "-CO-");
				// carboxylate on end
				new_grps.add("-COOH");
				
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
					// PRINT RXN STRING AS A DECARBOXYLATION 
					// couples ATP
					String rxn_string = prod+" + atp ---> "+sub+" + adp + co2(aq)";
					rxns.add( "4.1.1"+"\t"+map.get( prod )+"\t"+map.get( sub )+"\t"+rxn_string  );
					// couples Pi
					String rxn_string222 = prod+" + pi ---> "+sub+" + co2(aq)";
					rxns.add( "4.1.1"+"\t"+map.get( prod )+"\t"+map.get( sub )+"\t"+rxn_string222  );
				}
				else{
					System.out.println("ErrorQ_4_1_1: compound not in compound list. "+sub);
				}				
			}

			
		}
		
		
		
		//============================================================================
			
		
		
		
		
		
		
		
 	}// ####################### END  EC 4.1.1   ####################################

	
 	
 	
 	
 	
 	
 	
 	
 	
 	public static void _2_1(String sub, ArrayList<String> rxns, HashMap<String,Integer> map ){
 		
 		//   (A)   ==>  R(H)-(OH)R' --> R=R' + h2o         
 		//   (B)   ==>  R(OH)-(OH)R' --> R(=O)-R' + h2o       ( & other way round )
 		//   (C)   ==>  R(OH)-COOH  -->  R(=O)-CHO + h2o
 		
 		
		ArrayList<String> grps = ReactionMechanisms.getGroups( sub );

 		
 		//  =============================    (A)    ================================================
		HashMap<String,String> ADD_H = new HashMap<String,String>();
		ADD_H.put("=CH2", "-CH3");            ADD_H.put("CH2=", "CH3-");
		ADD_H.put("=CH(OH)", "-CH2(OH)");     ADD_H.put("CH(OH)=", "CH2(OH)-");
		ADD_H.put("=CO", "-CHO");             ADD_H.put("CO=", "CHO-");
		ADD_H.put("=CHp", "-CH2p");           ADD_H.put("CHp=", "CH2p-");
		ADD_H.put("-CH=", "-CH2-");           ADD_H.put("=CH-", "-CH2-");
		ADD_H.put("-C(OH)=", "-CH(OH)-");     ADD_H.put("=C(OH)-", "-CH(OH)-");
		ADD_H.put("-Cp=", "-CHp-");           ADD_H.put("=Cp-", "-CHp-");
		
		HashMap<String,String> ADD_OH = new HashMap<String,String>();
		ADD_OH.put("=CH2", "-CH2(OH)");       ADD_OH.put("CH2=", "CH2(OH)-");
		ADD_OH.put("=CO", "-COOH");           ADD_OH.put("CO=", "COOH-");          // Should this be included?
		ADD_OH.put("-CH=", "-CH(OH)-");       ADD_OH.put("=CH-", "-CH(OH)-");
		//ADD_OH.put("=CHp", "-CH(OH)p");     ADD_OH.put("CHp=", "CH(OH)p-");      // no geminal-diols
		//ADD_OH.put("=CH(OH)",);             ADD_OH.put("CH(OH)=", );             //     
		//ADD_OH.put("-C(OH)=", );            ADD_OH.put("=C(OH)-", );             //
		//ADD_OH.put("-Cp=", );               ADD_OH.put("=Cp-", );                //		

		for( int g=1; g < grps.size(); g++ ){
			String group = grps.get( g );
			String prev_grp = grps.get( g-1 );
			
			if(  group.charAt(0) == '='  ){
			
				// (i) try add H to left neighbour, OH to this group
				if( ADD_H.containsKey( prev_grp )  &&  ADD_OH.containsKey( group ) ){
					ArrayList<String> new_grps = new ArrayList<String>( grps );
					new_grps.set( (g-1), ADD_H.get( prev_grp ) );
					new_grps.set(  g,    ADD_OH.get( group ) );
					
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
						// WRITE THIS BACKWARDS AS A DEHYDRATION
						String rxn_string = prod+" --> "+sub+" + h2o";
						rxns.add( "4.2.1"+"\t"+map.get( prod )+"\t"+map.get( sub )+"\t"+rxn_string  );
					}
					else{
						System.out.println("Error_4_2_1: dehydration product not in compound list. "+sub);
					}					
				}
				//
				// (ii) try add OH to left neighbour, H to this
				if( ADD_OH.containsKey( prev_grp )  &&  ADD_H.containsKey( group ) ){					
					ArrayList<String> ng = new ArrayList<String>( grps );
					ng.set( (g-1), ADD_OH.get( prev_grp ) );
					ng.set(  g,    ADD_H.get( group ) );
					
					String tempProd = ReactionMechanisms.formProd( ng );  				
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
						// WRITE THIS BACKWARDS AS A DEHYDRATION
						String rxn_string = prod+" --> "+sub+" + h2o";
						rxns.add( "4.2.1"+"\t"+map.get( prod )+"\t"+map.get( sub )+"\t"+rxn_string  );
					}
					else{
						System.out.println("Error-B_4_2_1 dehydration product not in compound list. "+sub);
					}
				}
							
			}
			
		} // ==============================   END (A)  ===========================================================
// 		
 		
 		
 		
 		
 		
 		
 		// ======================  Combine (B) and (C) using :   ===========================
		HashMap<String,String> LOSE_OH_TO_O = new HashMap<String,String>();
		LOSE_OH_TO_O.put("-CH2(OH)", "-CHO");     LOSE_OH_TO_O.put("CH2(OH)-", "CHO-");     
		LOSE_OH_TO_O.put("-CH(OH)-", "-CO-");
		
		HashMap<String,String> LOSE_OH_TO_H = new HashMap<String,String>();
		LOSE_OH_TO_H.put("-CH2(OH)", "-CH3");     LOSE_OH_TO_H.put("CH2(OH)-", "CH3-");     
		LOSE_OH_TO_H.put("-CH(OH)-", "-CH2-");
		LOSE_OH_TO_H.put("-COOH", "-CHO");        LOSE_OH_TO_H.put("COOH-", "CHO-");
		// Bartek did not include these 2 transformations:
		// add "R1=C(OH)-CH(OH)-R2 > R1=CH-CO-R2" to Bartek's code
		LOSE_OH_TO_H.put("=CH(OH)", "=CH2");      LOSE_OH_TO_H.put("CH(OH)=", "CH2=");
		LOSE_OH_TO_H.put("-C(OH)=", "-CH=");      LOSE_OH_TO_H.put("=C(OH)-", "=CH-");

		
		for( int g=0; g < grps.size(); g++ ){
			String group = grps.get( g );
		
			if( LOSE_OH_TO_O.containsKey( group ) ){			
				
				if( g>0 ){
					// check left group
					String prevgroup = grps.get( g-1 );
					
					if( LOSE_OH_TO_H.containsKey( prevgroup ) ){					
						ArrayList<String> new_grps = new ArrayList<String>( grps );
						new_grps.set( (g-1), LOSE_OH_TO_H.get( prevgroup ) );
						// find what g should become:
						new_grps.set( g, LOSE_OH_TO_O.get( group ) ); 
											
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
							String rxn_string = sub+" --> "+prod+" + h2o";
							rxns.add( "4.2.1"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
						}
						else{
							System.out.println("Error: 4_2_1 product not in compound list. "+prod);
						}																		
					}
				}
				//
				//
				if( g < grps.size()-1 ){
					//check right group					
					String nextgroup = grps.get( g+1 );
					if( LOSE_OH_TO_H.containsKey( nextgroup ) ){
						
						ArrayList<String> new_grps = new ArrayList<String>( grps );
						new_grps.set( (g+1), LOSE_OH_TO_H.get( nextgroup ) );
						// find what current group should become:
						new_grps.set(  g, LOSE_OH_TO_O.get( group )  ); 
						
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
							String rxn_string = sub+" --> "+prod+" + h2o";
							rxns.add( "4.2.1"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
						}
						else{
							System.out.println("ErrorB: 4_2_1 product not in compound list. "+prod);
						}																		
					}		
				}	
			}
		}// end group loop
 		
 		// ==========================   End (B) and (C)   ============================================================
 		
 		
 		
 		
 		
 		
 		
 	}// #########   END EC 4.2.1  ################
	
 	
 	
 	
 	
 	

}








