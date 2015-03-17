import java.util.*;
import java.io.*;

public class CompoundGenerator{
	
	/**
	 * Program to take a set of chemical groups are arrange them into as many linear
	 * molecules as possible. "Biological molecules" are defined as having at least one
	 * hydrogen atom and at least one oxygen atom (either as an O or as part of a phosphate 
	 * group) i.e. disallow hydrocarbons.
	 * 
	 * If you want only compounds that are electrostatically charged in solution (through
	 * a dissociated COOH or phosphate group), uncommment the obvious line of code.
	 * Phosphates contribute -2 to charge, carboxyl group contributes -1.
	 */
	
	// I get different number of CHOPN molecules than Bartek - does he allow molecules that are
	// "hydrocarbons" with an N? (i.e. no O or phosphate group?) 
	
	
	public static void main(String args[])throws IOException{
		
		long start_time = System.currentTimeMillis();
		
		
		PrintWriter out = new PrintWriter(new FileWriter("output_CompoundGenerator.txt")); 
		
		// Program will generate molecules of this length and smaller, down to 2 carbon atoms.
		int maxSize = 4;


		// Note: need to enter the "end groups" in format with bond on RHS: i.e. "CH3-" and NOT "-CH3".
		String[] enders = {"CO=","COp-","CH3-","CH2=","COOH-","CH2(OH)-","CH2p-","CHO-", "CH(OH)=", "CHp="}; 
		String[] body = {"-CH2-","-CH=", "-CH(OH)-", "-C(OH)=","-CHp-","-Cp=","-CO-"};
		// Double hydroxyl or phosphates, -C(OH)(OH)- or -Cpp-, as well as the two groups -CH(OH)p- and -C(OH)p- 
		// are not allowed due to being unstable in solution. No triple bonds, or ester linkages.
		
		// + Nitrogen containing groups
	//	String[] enders = {"CO=","COp-","CH3-","CH2=","COOH-","CH2(OH)-","CH2p-","CHO-", "CH(OH)=", "CHp=", "CH2(NH2)-", "CO(NH2)-"}; 
	//	String[] body = {"-CH2-","-CH=", "-CH(OH)-", "-C(OH)=","-CHp-","-Cp=","-CO-", "-CH(NH2)-"};
		
		
		
		ArrayList<String> compounds = new ArrayList<String>();	

		// Get all compounds, starting with length 2 ones.
		for( int desiredLength = 2; desiredLength <= maxSize;  desiredLength++ ){
    		
			// Get compounds of desired length
			ArrayList<String> cmpnds = new ArrayList<String>();	
			int currentLength = 0;
    	    
			cmpnds =  Methods.generateCompounds( enders,  body, desiredLength, currentLength, cmpnds );
    	    
    	    // Add to full list
    	    for( String s: cmpnds ){
    	    	compounds.add( s );
    	    }
		}
		
		//Counts all possible combinations
		int count=0;
		//Counts all combinations of biological interest
		int countBioCompounds=0;
		
		for( int i=0; i<compounds.size(); i++){
			
			count++;
		//    out.println(  count+ "\t" +compounds.get(i).toString() );

			
			if( compounds.get(i).contains("O")  ||  compounds.get(i).contains("p")  ){
					if( compounds.get(i).contains("H") ){
		                //Then compound is of biological interest
						
						int phosphates = countGroupOccurrences( compounds.get(i), "p" );
						int carboxyl = countGroupOccurrences( compounds.get(i), "COOH");						
						
						int charge = (-2)*phosphates + (-1)*carboxyl;
						
						
						// Restrict to only charged molecules here
						if( charge != 0 ){
				    //	if( true ){	
						    countBioCompounds++;
						    out.println(  countBioCompounds+ "\t" +compounds.get(i).toString() + "\t" + charge );
						}
						
						
					}
			}
		}
		System.out.println();
		System.out.println("number of bio compounds = " + countBioCompounds );
		
		out.close();

		
		long end_time = System.currentTimeMillis();
		double run_time = (double)( end_time - start_time) / 1000.0;
		
		System.out.println("Run_time = " + run_time);
		
		
	}	//end main
	
	
	
	
	
	
	
	
	//Method to count number of occurrences of a String (particular chemical group)
	public static int countGroupOccurrences( String molecule,  String group ){		
		int index = molecule.length()+1;
		int count=0;
		
		while( true ){
			index = molecule.lastIndexOf( group, index-1 );
			if( index >= 0 ){
				count++;
			}
			else{
				break;
			}
		}	
	    return count;
	}
	
	

	
	
}//end class

