import java.util.ArrayList;
import java.util.HashMap;

public class EC3 {

	/**  EC 3 = Hydrolase
	 *    3.1 = Acting on ester bonds
	 *  3.1.3 = Phosphoric monoester hydrolases
	 *        ==>   R-p + h20 --> R-OH + pi      
	 *        
	 *    3.6 = Acting on acid anhydrides
	 *  3.6.1 = In phosphorus-containing anhydrides
	 *        ==>   R-COp + h2o --> R-COOH + pi    
	 */
	
	
	public static void _1_3(String sub, ArrayList<String> rxns, HashMap<String,Integer> map ){
		//  =========== R-p + h2o --> R-(OH) + pi ===========
		// Manually createashMap containing all possible R groups and their transform
		//   NEEDS UPDATED IF WE ADD MORE GROUPS, MAYBE RELEVANT NITROGEN ONES
		HashMap<String,String> grpMap = new HashMap<String,String>();
		grpMap.put("CH2p-", "CH2(OH)-");    grpMap.put("-CH2p", "-CH2(OH)");
		grpMap.put("CHp=", "CH(OH)=");      grpMap.put("=CHp", "=CH(OH)");
		grpMap.put("-CHp-", "-CH(OH)-");
		grpMap.put("-Cp=", "-C(OH)=");      grpMap.put("=Cp-", "=C(OH)-");
		//--------------------------------------------------------------- //
		
		ArrayList<String> grps = ReactionMechanisms.getGroups( sub );
		for( int g=0; g<grps.size(); g++ ){
			String group = grps.get( g );
			
			if( grpMap.containsKey( group ) ){
				ArrayList<String> new_grps = new ArrayList<String>( grps );
				new_grps.set( g, grpMap.get( group ) );
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
					String rxn_string = sub+" + h2o --> "+prod+" + pi";
					rxns.add( "3.1.3"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
				}
				else{
					System.out.println("Error: 3_1_3 product not in compound list. "+prod);
				}			
			}
		}
		
	}
	// ---- END  EC 3.1.3
	
	
	
	
	
	public static void _6_1(String sub, ArrayList<String> rxns, HashMap<String,Integer> map ){
		//  ==========   R-COp + h2o  --> R-COOH + pi  ===========
		if( sub.contains("-COp")){
			String temp = sub.replace("-COp", "-COOH" );
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
				String rxn_string = sub+" + h2o --> "+prod+" + pi";
				rxns.add( "3.6.1"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
			}
			else{
				System.out.println("Error: 3_6_1 product not in compound list. "+prod);
			}			
		}
	
		if( sub.contains("COp-") ){
			String temp = sub.replace("COp-", "COOH-" );
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
				String rxn_string = sub+" + h2o --> "+prod+" + pi";
				rxns.add( "3.6.1"+"\t"+map.get(sub)+"\t"+map.get(prod)+"\t"+rxn_string  );
			}
			else{
				System.out.println("Error: 3_6_1 product not in compound list. "+prod);
			}		
		}
		// =================================================================================	
	}
	// End EC 3.6.1
	
	

	
	
}
