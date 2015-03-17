import java.util.*;

public class EC2 {	
	
	/**  EC 2 = Transferase
	 *    2.7 = Transferring phosphorus-containing groups
	 *  2.7.1 = with an alcohol group as acceptor
	 *        ==>   R-OH + atp --> R-p + adp                                 
	 *        ==>   R=Cp-R' + adp --> RH-CO-R' + atp  (pyruvate kinase)      
	 *  
	 *  2.7.2 = with a carboxy group as acceptor
	 *        ==>   R-COOH + atp --> R-COp + adp                             
	 *        
	 *  2.7.9 = with paired acceptors  (the dikinase reactions)
	 *        ==>   R=Cp-R' + amp + pi --> RH-CO-R' + atp + h2o              
 	 *        ==>   R=Cp-R' + amp + ppi --> RH-CO-R' + atp + pi
	 */
	
	 // Final "H_MAP" to hold all groups that can have a hydrogen added to them
     private static final Map<String, String> H_MAP = createMap();
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
     // i.e. it doesn't include "-CH2-" going to "CH3-" for example.
     
     
	 // Final "OH_P_MAP" to hold all hydroxyl groups that can be phosphorylated
     private static final Map<String, String> OH_P_MAP = createMap2();
     private static Map<String, String> createMap2() {
        Map<String, String> result = new HashMap<String, String>();
		result.put("CH2(OH)-", "CH2p-");    result.put("-CH2(OH)", "-CH2p");
		result.put("CH(OH)=", "CHp=");      result.put("=CH(OH)", "=CHp");
		result.put("-CH(OH)-", "-CHp-");
		result.put("-C(OH)=", "-Cp=");      result.put("=C(OH)-", "=Cp-");
        return Collections.unmodifiableMap(result);
    }
     
     
     
     
     
     
     
     
 	public static void _7_1(String sub, ArrayList<String> rxns, HashMap<String,Integer> map ){
 		
		ArrayList<String> grps = ReactionMechanisms.getGroups( sub );
		for( int g=0; g<grps.size(); g++ ){
			String group = grps.get( g );
			
			// ========  R-OH + atp --> R-p + adp  ==============
			if( OH_P_MAP.containsKey( group ) ){
				ArrayList<String> new_grps = new ArrayList<String>( grps );
				new_grps.set( g, OH_P_MAP.get( group ) );
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
					String rxn_string = sub+" + atp --> "+prod+" + adp";
					rxns.add( "2.7.1"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
				}
				else{
					System.out.println("ErrorA: 2_7_1 product not in compound list. Sub=" + sub);
				}			
			}
	 		//===================================================
	 				
	 		
	 		// ======= R=Cp-R' + adp --> RH-CO-R' + atp  ========
	 		if( group.equals("=Cp-") ){
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
						String rxn_string = sub+" + adp --> "+prod+" + atp";
						rxns.add( "2.7.1"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
					}
					else{
						System.out.println("ErrorB: 2_7_1 product not in compound list. Sub: " + sub);
					}			
				}		
	 		}
	 		//
	 		if( group.equals("-Cp=") ){
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
						String rxn_string = sub+" + adp --> "+prod+" + atp";
						rxns.add( "2.7.1"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
					}
					else{
						System.out.println("ErrorC: 2_7_1 product not in compound list. Sub: " + sub);
					}			
				}		
	 		}
	 		// ==================================================
	 		
	 		
	 		// ======= CHp=R + adp --> CHO-RH + atp  ========
	 		if( group.equals("CHp=") ){
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
						String rxn_string = sub+" + adp --> "+prod+" + atp";
						rxns.add( "2.7.1"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
					}
					else{
						System.out.println("ErrorC2: 2_7_1 product not in compound list. Sub: " + sub);
					}			
				}		
	 		}
	 		if( group.equals("=CHp") ){
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
						String rxn_string = sub+" + adp --> "+prod+" + atp";
						rxns.add( "2.7.1"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
					}
					else{
						System.out.println("ErrorC3: 2_7_1 product not in compound list. Sub: "+sub);
					}			
				}		
	 		}
	 		// ======================================================================================
	 	
		}
 	}// END  2.7.1
 	
 	
 	
 	
 	
 	
 	

 	public static void _7_2(String sub, ArrayList<String> rxns, HashMap<String,Integer> map ){
 		// ===========  R-COOH + atp --> R-COp + adp  ==============
 		if( sub.contains( "COOH-" ) ){		
    		String temp = sub.replace("COOH-", "COp-" );
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
				String rxn_string = sub+" + atp --> "+prod+" + adp";
				rxns.add( "2.7.2"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
			}
			else{
				System.out.println("ErrorX: 2_7_2 product not in compound list. Sub: "+sub);
			}
 		}
 		//
 		if( sub.contains( "-COOH" ) ){
 	 		String temp = sub.replace("-COOH", "-COp" );
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
				String rxn_string = sub+" + atp --> "+prod+" + adp";
				rxns.add( "2.7.2"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
			}
			else{
				System.out.println("ErrorY: 2_7_2 product not in compound list: sub: "+ sub );
			}
 		}	
 	}
 	// ######## END 2.7.2 ###########
 	
 	
 	
 	
 	
 	
 	public static void _7_9(String sub, ArrayList<String> rxns, HashMap<String,Integer> map ){
 		// .... this is exact same as one of the 2.7.1 transformations...
 		
		ArrayList<String> grps = ReactionMechanisms.getGroups( sub );
		for( int g=0; g<grps.size(); g++ ){
			String group = grps.get( g );
		
	 		//== both dikinase:  R=Cp-R' + (amp+pi)/(amp+ppi) --> RH-CO-R' + (atp+h2o)/(atp+pi) ====
	 		if( group.equals("=Cp-") ){
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
						String rxn_string = sub+" + amp + pi --> "+prod+" + atp + h2o";
						rxns.add( "2.7.9"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
						String rxn_string2 = sub+" + amp + ppi --> "+prod+" + atp + pi";
						rxns.add( "2.7.9"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string2  );
					}
					else{
						System.out.println("ErrorB: 2_7_9 product not in compound list. Sub: "+sub);
					}			
				}		
	 		}
	 		//
	 		if( group.equals("-Cp=") ){
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
						String rxn_string = sub+" + amp + pi --> "+prod+" + atp + h2o";
						rxns.add( "2.7.9"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
						String rxn_string2 = sub+" + amp + ppi --> "+prod+" + atp + pi";
						rxns.add( "2.7.9"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string2  );
					}
					else{
						System.out.println("ErrorC: 2_7_9 product not in compound list. Sub: "+sub);
					}			
				}		
	 		} // =============================================
	 		
	 		
	 		
	 		
	 		// ======= CHp=R + adp --> CHO-RH + atp  ===========
	 		if( group.equals("CHp=") ){
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
						String rxn_string = sub+" + amp + pi --> "+prod+" + atp + h2o";
						rxns.add( "2.7.9"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
						String rxn_string2 = sub+" + amp + ppi --> "+prod+" + atp + pi";
						rxns.add( "2.7.9"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string2  );
					}
					else{
						System.out.println("ErrorD2: 2_7_9 product not in compound list. Sub: "+sub);
					}			
				}		
	 		}
	 		if( group.equals("=CHp") ){
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
						String rxn_string = sub+" + amp + pi --> "+prod+" + atp + h2o";
						rxns.add( "2.7.9"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
						String rxn_string2 = sub+" + amp + ppi --> "+prod+" + atp + pi";
						rxns.add( "2.7.9"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string2  );
					}
					else{
						System.out.println("ErrorD3: 2_7_9 product not in compound list. Sub: "+sub);
					}			
				}		
	 		}
	 		// ======================================================================================
			
			
			
		}// end group loop
 
 	}
 	// ########## END 2.7.9 ##############	

     
     
     
	
	

}
