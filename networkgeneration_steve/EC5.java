import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;


public class EC5 {
	/**  EC 5 = Isomerases
	 *    5.3 = Intramolecular oxidoreductases
	 *  5.3.1 = Interconverting aldoses and ketoses  
	 *        ==>   R-CO-CH2(OH) -->  R-CH(OH)-CHO      
	 *        
	 *  5.3.2 = Interconverting keto- and enol-groups (tautomerization) 
	 *        ==> R=C(OH)=R'  -->  RH-CO-R'                                 
	 *        ==> R=CH(OH)    -->  RH-CHO                                   
	 *        
	 *    5.4 = Intramolecular transferases 
	 *  5.4.2 = Phosphotransferases                     :: Implemented generically so that non-end groups included          
	 *        ==> R-CH(OH)-CH2p  -->  R-CHp-CH2(OH) 
	 *        ==> R=C(OH)-CH2p   -->  R=Cp-CH2(OH)                
	 *        ==> R-CHp-COOH  -->  R-CH(OH)-COp                     
	 *        ==> R=Cp-COOH   -->  R=C(OH)-COp                      
	 */
	
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
	
	 // Final "H_MAP" to hold all groups that can have a hydrogen added to them
    private static final Map<String, String> H_MAP = createMap();
    private static Map<String, String> createMap() {
       Map<String, String> result = new HashMap<String, String>();
       result.put("CO=", "CHO-");           result.put("=CO", "-CHO"); 
       result.put("CH2=", "CH3-");          result.put("=CH2", "-CH3"); 
       result.put("CH(OH)=", "CH2(OH)-");   result.put("=CH(OH)", "-CH2(OH)"); 
       result.put("CHp=", "CH2p-");         result.put("=CHp", "-CH2p"); 
       result.put("-CH=", "-CH2-");         result.put("=CH-", "-CH2-");
       result.put("-C(OH)=", "-CH(OH)-");   result.put("=C(OH)-", "-CH(OH)-"); 
       result.put("-Cp=", "-CHp-");         result.put("=Cp-", "-CHp"); 
       return Collections.unmodifiableMap(result);
   }

	 // Final "P_ISO_MAP" to hold all pairs of groups for intramolecular phosphate transfers 
    private static final Map<String, String> P_ISO_MAP = create_P_ISO_MAP();
    private static Map<String, String> create_P_ISO_MAP() {
       Map<String, String> grpMap = new HashMap<String, String>();
		grpMap.put("CH2p-", "CH2(OH)-");    grpMap.put("-CH2p", "-CH2(OH)");
		grpMap.put("CHp=", "CH(OH)=");      grpMap.put("=CHp", "=CH(OH)");
		grpMap.put("-CHp-", "-CH(OH)-");
		grpMap.put("-Cp=", "-C(OH)=");      grpMap.put("=Cp-", "=C(OH)-");
		grpMap.put("COp-", "COOH-");        grpMap.put("-COp", "-COOH");
       return Collections.unmodifiableMap( grpMap );
   }
    
    
    
    
    
	
	
	public static void _3_1(String sub, ArrayList<String> rxns, HashMap<String,Integer> map ){
		
		// Is a non-end group switch allowed? i.e.  R-CO-CH(OH)-R' <-> R-CHOH)-CO-R'  ???
		// Read more about the enzymatic mechanism of 5.3.1; most occur on cyclic compounds.
		
		// =======    R-CO-CH2(OH) -->  R-CH(OH)-CHO   =================
		if( sub.contains("-CO-CH2(OH)")){
			String temp = sub.replace("-CO-CH2(OH)", "-CH(OH)-CHO" );
			String prod="";
			// Check product or its palindrome is in map
			if( map.containsKey( temp )  ){
				prod = temp;
			}
			else if( map.containsKey(  Methods.getMolecularPalindrome( temp )   )){
				prod = Methods.getMolecularPalindrome( temp );
			}
			//
			if( prod != "" ){
				String rxn_string = sub+" --> "+prod;
				rxns.add( "5.3.1"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
			}
			else{
				System.out.println("Error: 5_3_1 product not in compound list. "+prod);
			}			
		}
		if( sub.contains("CH2(OH)-CO-")){
			String temp = sub.replace("CH2(OH)-CO-", "CHO-CH(OH)-" );
			String prod="";
			// Check product or its palindrome is in map
			if( map.containsKey( temp )  ){
				prod = temp;
			}
			else if( map.containsKey(  Methods.getMolecularPalindrome( temp )   )){
				prod = Methods.getMolecularPalindrome( temp );
			}
			//
			if( prod != "" ){
				String rxn_string = sub+" --> "+prod;
				rxns.add( "5.3.1"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
			}
			else{
				System.out.println("ErrorB: 5_3_1 product not in compound list. "+prod);
			}			
		}//=============================================================================
	
	} // ########### END 5.3.1 ###########################################################
	
	
	
	
	

	
	
	public static void _3_2(String sub, ArrayList<String> rxns, HashMap<String,Integer> map ){
	// 5.3.2 is tautomerization

		ArrayList<String> grps = ReactionMechanisms.getGroups( sub );
		for( int g=0; g<grps.size(); g++ ){
			
			String group = grps.get( g );			
			
			// ========  R=CH(OH)-R'  -->  RH-CO-R'  ======================
	 		if( group.equals("=C(OH)-") ){
	 			String lastgroup = grps.get( g-1 );
				if( H_MAP.containsKey( lastgroup ) ){
					ArrayList<String> new_grps = new ArrayList<String>( grps );
					new_grps.set( (g-1), H_MAP.get( lastgroup ) );
					new_grps.set(   g,   "-CO-" );
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
						String rxn_string = sub+" --> "+prod;
						rxns.add( "5.3.2"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
					}
					else{
						System.out.println("ErrorA: 5_3_2 product not in compound list. "+prod);
					}			
				}		
	 		}
	 		if( group.equals("-C(OH)=") ){
	 			String nextgroup = grps.get( g+1 );
				if( H_MAP.containsKey( nextgroup ) ){
					ArrayList<String> new_grps = new ArrayList<String>( grps );
					new_grps.set( (g+1), H_MAP.get( nextgroup ) );
					new_grps.set(   g,   "-CO-" );
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
						String rxn_string = sub+" --> "+prod;
						rxns.add( "5.3.2"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
					}
					else{
						System.out.println("ErrorB: 5_3_2 product not in compound list. "+prod);
					}			
				}		
	 		}// ==================================================================================
			
	 		
	 		// =======   R=CH(OH) --> RH-CHO ================
	 		if( group.equals("=CH(OH)") ){
	 			String lastgroup = grps.get( g-1 );
				if( H_MAP.containsKey( lastgroup ) ){
					ArrayList<String> new_grps = new ArrayList<String>( grps );
					new_grps.set( (g-1), H_MAP.get( lastgroup ) );
					new_grps.set(   g,   "-CHO" );
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
						String rxn_string = sub+" --> "+prod;
						rxns.add( "5.3.2"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
					}
					else{
						System.out.println("ErrorC: 5_3_2 product not in compound list. "+prod);
					}			
				}		
	 		}
	 		if( group.equals("CH(OH)=") ){
	 			String nextgroup = grps.get( g+1 );
				if( H_MAP.containsKey( nextgroup ) ){
					ArrayList<String> new_grps = new ArrayList<String>( grps );
					new_grps.set( (g+1), H_MAP.get( nextgroup ) );
					new_grps.set(   g,   "CHO-" );
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
						String rxn_string = sub+" --> "+prod;
						rxns.add( "5.3.2"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
					}
					else{
						System.out.println("ErrorD: 5_3_2 product not in compound list. "+prod);
					}			
				}		
	 		}
	 	// ==================================================================================
		
		} // end loop through groups
		
	} 		// ##########  END 5.3.2 ###############################################


	
	
	
	
	
	public static void _4_2(String sub, ArrayList<String> rxns, HashMap<String,Integer> map ){
		// ======  All isomerizations switching p group and adjacent OH =========
		ArrayList<String> grps = ReactionMechanisms.getGroups( sub );
		for( int g=0; g < grps.size(); g++ ){
			String group = grps.get( g );
		
			if( group.contains( "OH" ) ){			
				if( g>0 ){
					// check left group
					String prevgroup = grps.get( g-1 );
					if( P_ISO_MAP.containsKey( prevgroup ) ){					
						ArrayList<String> new_grps = new ArrayList<String>( grps );
						new_grps.set( (g-1), P_ISO_MAP.get( prevgroup ) );
						// find what g should become:
						String thisgroup = getKeyFromValue( P_ISO_MAP, group );
						//
						new_grps.set( g, thisgroup ); 	
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
							String rxn_string = sub+" --> "+prod;
							rxns.add( "5.4.2"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
						}
						else{
							System.out.println("Error: 5_4_2 product not in compound list. "+prod);
						}																		
					}
				}
				if( g < grps.size()-1 ){
					//check right group					
					String nextgroup = grps.get( g+1 );
					if( P_ISO_MAP.containsKey( nextgroup ) ){					
						ArrayList<String> new_grps = new ArrayList<String>( grps );
						new_grps.set( (g+1), P_ISO_MAP.get( nextgroup ) );
						// find what current group should become:
						String thisgroup = getKeyFromValue( P_ISO_MAP, group );
						new_grps.set( g, thisgroup ); 
						//
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
							String rxn_string = sub+" --> "+prod;
							rxns.add( "5.4.2"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
						}
						else{
							System.out.println("ErrorB: 5_4_2 product not in compound list. "+prod);
						}																		
					}		
				}	
			}
					
		}// end group loop	
		
	}  // ########################    END EC 5.4.2    ############################
	
	
	
	
	
}// end class



















//// DO WE ALLOW TRANSFER WITH NON_END GROUPS!?    ***NOT CURRENTLY IMPLEMENTED***
//
//
//
//// =======    R-CH(OH)-CH2p -->  R-CHp-CH2(OH)   =================
//if( sub.contains("-CH(OH)-CH2p")){
//	String temp = sub.replace("-CH(OH)-CH2p", "-CHp-CH2(OH)" );
//	String prod="";
//	// Check product or its palindrome is in map
//	if( map.containsKey( temp )  ){
//		prod = temp;
//	}
//	else if( map.containsKey(  Methods.getMolecularPalindrome( temp )   )){
//		prod = Methods.getMolecularPalindrome( temp );
//	}
//	//
//	if( prod != "" ){
//		String rxn_string = sub+" --> "+prod;
//		rxns.add( "5.4.2"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
//	}
//	else{
//		System.out.println("ErrorA: 5_4_2 product not in compound list. "+prod);
//	}			
//}
//if( sub.contains("CH2p-CH(OH)-")){
//	String temp = sub.replace("CH2p-CH(OH)-", "CH2(OH)-CHp-" );
//	String prod="";
//	// Check product or its palindrome is in map
//	if( map.containsKey( temp )  ){
//		prod = temp;
//	}
//	else if( map.containsKey(  Methods.getMolecularPalindrome( temp )   )){
//		prod = Methods.getMolecularPalindrome( temp );
//	}
//	//
//	if( prod != "" ){
//		String rxn_string = sub+" --> "+prod;
//		rxns.add( "5.4.2"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
//	}
//	else{
//		System.out.println("ErrorB: 5_4_2 product not in compound list. "+prod);
//	}			
//}
////=====================================================================
//
//
//
//// =======    R=C(OH)-CH2p -->  R=Cp-CH2(OH)   =================
//if( sub.contains("=C(OH)-CH2p")){
//	String temp = sub.replace("=C(OH)-CH2p", "=Cp-CH2(OH)" );
//	String prod="";
//	// Check product or its palindrome is in map
//	if( map.containsKey( temp )  ){
//		prod = temp;
//	}
//	else if( map.containsKey(  Methods.getMolecularPalindrome( temp )   )){
//		prod = Methods.getMolecularPalindrome( temp );
//	}
//	//
//	if( prod != "" ){
//		String rxn_string = sub+" --> "+prod;
//		rxns.add( "5.4.2"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
//	}
//	else{
//		System.out.println("ErrorC: 5_4_2 product not in compound list. "+prod);
//	}			
//}
//if( sub.contains("CH2p-C(OH)=")){
//	String temp = sub.replace("CH2p-C(OH)=", "CH2(OH)-Cp=" );
//	String prod="";
//	// Check product or its palindrome is in map
//	if( map.containsKey( temp )  ){
//		prod = temp;
//	}
//	else if( map.containsKey(  Methods.getMolecularPalindrome( temp )   )){
//		prod = Methods.getMolecularPalindrome( temp );
//	}
//	//
//	if( prod != "" ){
//		String rxn_string = sub+" --> "+prod;
//		rxns.add( "5.4.2"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
//	}
//	else{
//		System.out.println("ErrorD: 5_4_2 product not in compound list. "+prod);
//	}			
//}
//// ====================================================================================