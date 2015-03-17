import java.util.regex.Pattern;
import java.util.*;

public class Test {

	static String getKeyFromValue( HashMap<String, String> map, String val ){
		String key = "";
		for ( Map.Entry<String, String> entry : map.entrySet()) {
			if ( entry.getValue().equals( val ) ) {
				key = entry.getKey();
			}
		}
		return key;
	}

	
	public static void main(String[] args){
		
		Rxn cc = new Rxn("5.4.2	1002	1007	CH(OH)=Cp-Cp=CHp --> CHp=C(OH)-Cp=CHp");
		Rxn b = new Rxn("5.4.2	1002	1007	CH(OH)=Cp-Cp=CHp --> CHp=C(OH)-Cp=CHp");
		Rxn a = new Rxn("5.4.2	1007	1002	CHp=C(OH)-Cp=CHp --> CH(OH)=Cp-Cp=CHp");
		Rxn r = new Rxn("5.4.2	1007	1007	CHp=C(OH)-Cp=CHp --> CHp=C(OH)-Cp=CHp");
		
		System.out.println( b );
		
		ArrayList<Rxn> liss = new ArrayList<Rxn>();
		liss.add( cc );
		System.out.println("contains: "+liss.contains(b) );
		
		System.out.println( a==r );
		System.out.println( b.equals(cc) );
		
		
		System.out.print("b==a : "); System.out.println( b==a ); 
		System.out.println("b.equals(a) : " + b.equals(a) );
		System.out.println("b.equivalent(a) : " + b.equivalent( a )  );
		
		System.exit(0);
		
		
		
		
		HashMap<String, Integer> map = new HashMap<String, Integer>();
		map.put("COOH-Cp=CHp", 1);
		map.put("CH(OH)=C(OH)-CH2p", 2);
		map.put("CHp=C(OH)-CH2(OH)", 3);
		map.put("CH2p-CO-COOH", 4);
		map.put("CH3-CH(OH)-CH2(OH)", 5);
		map.put("CH2(OH)-CH(OH)-Cp=CHp", 6);
		map.put("CH2(OH)-CHp-Cp=CH(OH)", 7);

		ArrayList<String> rxns = new ArrayList<String>();

			
		String comp2 = "COOH-CH(OH)-CH(OH)-CH3";
	
		EC4._1_1( comp2, rxns, map  );
		
		//EC5._3_1(  comp2, rxns, map );
//		EC5._3_2(  comp2, rxns, map );
	//	EC5._4_2(  comp2, rxns, map );  
	//	EC2._7_1(  comp2, rxns, map );
	//	EC2._7_2(  comp2, rxns, map );
	//	EC2._7_9(  comp2, rxns, map );
//		EC3._1_3(  comp2, rxns, map );
	//	EC3._6_1(  comp2, rxns, map );
		
		for( String s: rxns ){
			System.out.println( s );				
		}
		System.exit(0);
		
		
		
		
		String c = "CH3-CH(OH)-CH(OH)-COOH";
		
        ArrayList<String> grps = ReactionMechanisms.getGroups( c );	
		System.out.println( grps );
		
		String mol = "";
		for( String g: grps ){
			mol = mol + g;
		}
		System.out.println( mol );
		mol = mol.replace("--", "-");
		System.out.println( mol );
		
		String prod = ReactionMechanisms.formProd( grps );
		System.out.println( "prod = "+prod );
		
		System.exit(0);
		
		
		
		
		
		
		
		
		String mol2 = "CH3-CO-COOH";
		
		System.out.println( Methods.getMolecularPalindrome( mol2 ) );
		
		
		String s1 = "CH2(OH)=CH(OH)-COOH";
		
		if( s1.contains("CH2(OH)-") ){
			System.out.println( "ok" );
		}
		System.exit(0);
		
		String s2 = s1.replaceAll( "CH2\\(OH\\)", "CHO");
		
		//s2 = "I am dan";
		
		System.out.println( s1 + " : " + s2 );
		
		
		
	}
	
}
