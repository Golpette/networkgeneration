import java.util.*;
import java.io.*;

public class ReactionGenerator {

	public static void main( String args[] )throws IOException{
		
		long start_time = System.currentTimeMillis();
				

		// Read desired compound list from file
		ArrayList<Compound> cmpnds = Methods.readCompounds( "CHOP_C4_charged.txt" );


		
		// Store formula<->label mapping in Hash table
		HashMap<String, Integer> map = new HashMap<String, Integer>();
		for( Compound c: cmpnds ){
			map.put( c.getFormula(), c.getLabel() );			
		}
		
		// Store rxns as list of Strings
		ArrayList<String> rxns = new ArrayList<String>();
		
		
		// For each compound in network, check which reactions apply 
		for( Compound c: cmpnds ){
						
		EC1._1_1( c.getFormula(), rxns, map );		
		EC1._2_1( c.getFormula(), rxns, map );
	        EC1._3_1( c.getFormula(), rxns, map ); 
			
		EC2._7_1( c.getFormula(), rxns, map );  
		EC2._7_2( c.getFormula(), rxns, map );  
		EC2._7_9( c.getFormula(), rxns, map );  

		EC3._1_3( c.getFormula(), rxns, map ); 
		EC3._6_1( c.getFormula(), rxns, map ); 
					
		EC4._1_1( c.getFormula(), rxns, map );     		
		EC4._2_1( c.getFormula(), rxns, map ); 
	
		EC5._3_1( c.getFormula(), rxns, map ); 
		EC5._3_2( c.getFormula(), rxns, map ); 
		EC5._4_2( c.getFormula(), rxns, map ); 
	
		EC6._4_1( c.getFormula(), rxns, map );
						
		}
		
		
		// Convert our Strings into Reaction objects
		//   (...should just have done this from the start rather than dealing with Strings)
		ArrayList<Rxn> reactions = new ArrayList<Rxn>();
		for( String r: rxns ){
			reactions.add( new Rxn( r ) );
		}
		// Remove duplicate entries which arise from symmetric molecules
		reactions = removeDuplicateRxns( reactions );
		
		
		
		PrintWriter out = new PrintWriter( new FileWriter( "output_ReactionGenerator.rxns" ) );
		
		for( Rxn r: reactions ){
		out.println( r );
		}
		out.close();
		
		
		long end_time = System.currentTimeMillis();
		double run_time = (double)( end_time - start_time )/1000.0;
		System.out.println("Run time = " + run_time );
		
		
	}
	
	
	
	
	
	
	
	
	static ArrayList< Rxn > removeDuplicateRxns( ArrayList<Rxn> rxns ){
		// Remove exact duplicates
		ArrayList<Rxn> newList = new ArrayList<Rxn>();
		for( Rxn s: rxns ){
			if( ! newList.contains( s ) ){
				newList.add( s );
			}
		}
		
		ArrayList<Rxn> finalList = new ArrayList<Rxn>();
		// Go through new list and remove
		for( Rxn s: newList ){
			boolean keep = true;
			
			//  (a)  x->x which arises from isomerization reactions
			if( s.getSub().equals( s.getProd() )  ){
				keep = false;
			}
			//  (b)  x->y and y->x via same EC class
			for( Rxn rr: finalList ){
				if( rr.equivalent( s )  ){
					keep = false;
				}
			}
			if( keep ){
				finalList.add( s );
			}
					
		}		
		return finalList;
	}
	
	
	
	
	
	
	
	
//	EC1._1_1( c.getFormula(), rxns, map );
//	// I've got 631, Bartek has 670.  Bartek has 70 decarboxylating, I have 31.
//	// Difference is in the r' group; I force it to be a 1-C group, Bartek allows it to
//	// be absent as well. Bartek then disallows the 3 decarbox. rxns that involve ketenes.
//		
//	EC1._2_1( c.getFormula(), rxns, map );                 //agrees with Bartek
//		
//    EC1._3_1( c.getFormula(), rxns, map );                 //Bartek has 145, I have 675 
//    // issue with 1.3.1 s both in Bartek's code and the string in  "types" file. I can alter it to get 675.
//				
//	EC2._7_1( c.getFormula(), rxns, map );                 //agrees with Bartek
//	EC2._7_2( c.getFormula(), rxns, map );                 //agrees with Bartek
//	EC2._7_9( c.getFormula(), rxns, map );                 //agrees with Bartek
//
//	EC3._1_3( c.getFormula(), rxns, map );                 //agrees with Bartek
//	EC3._6_1( c.getFormula(), rxns, map );                 //agrees with Bartek
//	
//	EC4._1_1( c.getFormula(), rxns, map );                 //agrees with Bartek			
//	EC4._2_1( c.getFormula(), rxns, map );                 //agrees with Bartek
//	//(...when Bartek's repeated rxns are removed and extra line is added
//	//to his 4.2.1:  "R1=C(OH)-CH(OH)-R2 > R1=CH-CO-R2" )
//	
//	EC5._3_1( c.getFormula(), rxns, map );                 //agrees with Bartek
//	EC5._3_2( c.getFormula(), rxns, map );                 //agrees with Bartek
//	EC5._4_2( c.getFormula(), rxns, map );                 //agrees with Bartek
//	//(... when add in 4 lines to his code and then remove duplicate reactions) 
//	//(Bartek only implements these transfers when end-groups are involved)
//	
//	EC6._4_1(  c.getFormula(), rxns, map  );               //agrees with Bartek
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
//	static ArrayList<String> removeDuplicateRxns( ArrayList<String> rxns ){
//		// Remove exact duplicates
//		ArrayList<String> newList = new ArrayList<String>();
//		for( String s: rxns ){
//			if( ! newList.contains( s ) ){
//				newList.add( s );
//			}
//		}
//		
//		ArrayList<String> finalList = new ArrayList<String>();
//		// Go through new list and remove
//		for( String s: newList ){
//			boolean keep = true;
//			Scanner scan = new Scanner( s );
//			String EC = scan.next();
//			String sub = scan.next();
//			String prod = scan.next();
//			
//			String equiv = "";
//			
//			if( sub.equals( prod ) ){
//				keep = false;
//			}
//					
//		}
//		//  (a)  x->y and y->x via same EC class
//		//  (b)  x->x which arises from isomerization reactions
//		
//		return finalList;
//	}
	
	
	
	
}
