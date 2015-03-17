import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;


public class EC6 {
	
	/**
	 *  EC 6 = Ligases
	 *   6.4 = Forming Carbon-Carbon Bonds 
	 * 6.4.1 = Only subclass to date
	 *       ==>  RH-CO-COOH + atp + co2 --> COOH-R-CO-COOH + adp + pi  
	 *   	
	 *   Tiny reaction class, only acts on alph-keto acids and acyl-CoA thioesters. [Is this being too restrictive though?]   
	 *   NB:  I don't know if it is the case that the carbon being carboxylated "R" can 
	 *   only be an alkyl group: CH3, CH2, CH; or this is an artefact of metabolism.
	 *      
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
       result.put("=CH2", "=CH-");           result.put("CH2=", "-CH=");
       result.put("=CH(OH)", "=C(OH)-");     result.put("CH(OH)=", "-C(OH)=");
       result.put("=CHp", "=Cp-");           result.put("CHp=", "-Cp=");
       return Collections.unmodifiableMap(result);
   }   
	

    
    
	public static void _4_1(String sub, ArrayList<String> rxns, HashMap<String,Integer> map ){
	/** EC 6.4.1 **/	
		
		ArrayList<String> grps = ReactionMechanisms.getGroups( sub );	
		int numgrps = grps.size();
		
		
		// =========== RH-CO-COOH + atp + co2 --> COOH-R-CO-COOH + adp + pi ============
		//
		// Since we only have linear compounds this is easy:	
		if( sub.contains( "-CO-COOH") && MAKE_BODY.containsKey( grps.get(0) ) ){
			
			ArrayList<String> new_grps = new ArrayList<String>( grps );
			// make end group a body group
			new_grps.set(0, MAKE_BODY.get( grps.get(0) )  );
			// add carboxyl group
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
				String rxn_string = sub+" + atp + co2(aq) --> "+prod+" + adp + pi";
				rxns.add( "6.4.1"+"\t"+map.get( sub )+"\t"+map.get( prod )+"\t"+rxn_string  );
			}
			else{
				System.out.println("Error_6_4_1: compound not in compound list. "+sub);
			}			
		}	
		// or other end:		
		if( sub.contains( "COOH-CO-") &&  MAKE_BODY.containsKey( grps.get(numgrps-1)  )   ){
			
			ArrayList<String> new_grps = new ArrayList<String>( grps );
			// make end group a body group
			new_grps.set( numgrps-1,   MAKE_BODY.get( grps.get( numgrps-1 ) )  );
			// add carboxyl group
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
				String rxn_string = sub+" + atp + co2(aq) --> "+prod+" + adp + pi";
				rxns.add( "6.4.1"+"\t"+map.get( sub )+"\t"+map.get( prod )+"\t"+rxn_string  );
			}
			else{
				System.out.println("Error_6_4_1: compound not in compound list. "+sub);
			}			
			
		}
		// ===========================================================================================	
		
	}
    
    
    
	
	
}
