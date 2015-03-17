import java.util.*;
import java.io.*;

public class Methods{


	
	public static ArrayList<String> generateCompounds( String[] enders, String[] body, int desiredSize, int currentLength, ArrayList<String> compounds ){
    /** Method to enumerate all possible combinations of groups (for LINEAR MOLECULES ONLY)  **/
		
		for(int i=0; i<enders.length; i++){

			String molecule = enders[i];
			currentLength = 1;
			
			addBody( enders, body, desiredSize, currentLength, molecule, compounds);
		}

		return compounds;
	}

	
	
	
	// This is called in above method
	public static ArrayList<String> addBody( String[] enders, String[] body, int desiredSize, int currentLength, String molecule, ArrayList<String> compounds ){
        
		/** Essentially a depth-first search algorithm for finding all allowed
        *   combinations of "body" and "ender" carbon groups, for molecules of "desiredSize"
        */
		
		// If we are one group short of the desired length, attempt to fit all "ender"
		// groups to molecule.
		if( currentLength == desiredSize-1 ){
			
			char finalBond = molecule.charAt( molecule.length()-1 );

			for( int ee=0; ee<enders.length; ee++){
				if( enders[ee].contains( finalBond+"" ) ){

					//add ender
					if( enders[ee].charAt(0) == finalBond ){//this should never be the case due to form of our enders[]

						molecule = molecule.substring( 0, molecule.length()-1 );
						molecule = molecule + enders[ee];
						currentLength++;

					}
					else if( enders[ee].charAt( enders[ee].length()-1 ) == finalBond ){

						String orientedEnder = enders[ee].substring( 0, enders[ee].length()-1 );
						molecule = molecule + orientedEnder;
						currentLength++;
					}

					// If molecular palindrome has not yet been found, accept it, store it
					boolean palindromeAlreadyFound = false;
					for( int comp=0; comp<compounds.size(); comp++ ){
						// Check if "molecule" has already been found in palindromic form
						if( isMolecularPalindrome( compounds.get(comp), molecule )   ){
							palindromeAlreadyFound = true;
						}
					}
					if( !palindromeAlreadyFound ){
						compounds.add( molecule );
					}

					// Back-track: i.e. remove last group, decrease length by one
					molecule = removeLastGroup( molecule );   
					currentLength--;
				}
			}		    
			return compounds;
		} // end  if(currentLength == desiredSize-1 )

		
		for( int i=0; i<body.length; i++){

			boolean SuitableBodyGroup = true;
			
			//check if body[i] group can be added
			if( molecule.charAt( molecule.length()-1 ) == body[i].charAt( 0 ) ){

				molecule = molecule.substring( 0, molecule.length()-1 );   // - molecule.charAt( molecule.length()-1 );
				molecule = molecule + body[i];
				currentLength++;
			}
			else if(  molecule.charAt( molecule.length()-1 ) == body[i].charAt(  body[i].length()-1  )  ){

				char[] chars = body[i].toCharArray();
				String orientedgroup = "";
				for( int x=1; x<chars.length-1; x++ ){
					orientedgroup = orientedgroup + chars[x];
				}
				orientedgroup = orientedgroup + chars[0];

				molecule = molecule + orientedgroup;
				currentLength++;
			}
			else{
				SuitableBodyGroup = false;				
			}

			if( SuitableBodyGroup ){
				
				// Then keep adding more "body" groups until currentLength = (desiredSize -1).
				addBody( enders, body, desiredSize, currentLength, molecule, compounds);

				molecule = removeLastGroup( molecule );

				currentLength--;
			}
			
		}
		return compounds;

	}  //  End of addBody() method -------------------------
   


	

	
	public static boolean isMolecularPalindrome(String mol1, String mol2){
		/** Method to check if a "palindromic molecule" has already been generated
		 *  i.e.  CH3-CO-CH2(OH) is equivalent to CH2(OH)-CO-CH3.
		 */
		boolean isPalindrome = false;
		int L1 = mol1.length();    int L2 = mol2.length();

		if( L1 == L2 ){		
			String palindrome = getMolecularPalindrome( mol2 );

			if( palindrome.equals( mol1 ) ){
				isPalindrome = true;	
			}
		}
		return isPalindrome;
	} 

	
	
	
	public static String getMolecularPalindrome(String mol1){
		/** Generates molecular palindrome, i.e. CH3-CO-CH2(OH) is 
		 *  transformed to CH2(OH)-CO-CH3.
		 */
		int L1 = mol1.length();    
			
		ArrayList<Character> groupHolder = new ArrayList<Character>();
		String palindrome = "";	    
				
		for(int i=0; i < L1; i++){
			char c = mol1.charAt( i );
			if( c == '-' || c == '=' ){
				String group = "";
				for( int j=0; j<groupHolder.size(); j++){
					group  = group + groupHolder.get(j);
				}
				palindrome = c + group + palindrome;
				groupHolder.clear();
			}
			else if( i == (L1-1) ){
				String group = "";
				for( int j=0; j<groupHolder.size(); j++){
					group  = group + groupHolder.get(j);
				}
				palindrome = group + c + palindrome;
			}
			else{
				groupHolder.add( c );
			}
		}
		return palindrome;
	} 
	
	
	


	
	
	public static String removeLastGroup( String molecule ){
		/** Method to remove the last group from a molecule but leaving the bond identifier **/	
						
		int remFromIndex = 9999;
		for( int i=molecule.length()-2;  i>0; i-- ){
			if( molecule.charAt(i) == '-'  || molecule.charAt(i) == '=' ){
				//Find index of bond we want to leave after removing group
				remFromIndex = i;
				//Leave for loop
				i=0;
			}
		}
		
		// Remove group, leave "dangling bond" on molecule
		for( int i=molecule.length()-1;   i>remFromIndex;   i--){
			molecule = molecule.substring(  0, molecule.length()-1  );   // - molecule.charAt( i ); 
		} 

		return molecule;
	}






	
	
	// Read in compounds file
	
	// What's the correct way to do a readFile method? I need the file to be in 3
	// columns but do not explicitly check this / output any error messages.
	public static ArrayList<Compound> readCompounds( String filename )throws IOException{
		ArrayList<Compound> comps = new ArrayList<Compound>();
		
		BufferedReader in = new BufferedReader(new FileReader( filename ) );
		String s = "";
		Scanner scan = new Scanner( s );
		while( true ){
			try{
				s = in.readLine();
				scan = new Scanner( s );
				int label = scan.nextInt();
				String formula = scan.next();
				int charge = scan.nextInt();
				Compound c = new Compound( label, formula, charge );
				comps.add( c );
				
			}catch( NullPointerException e ){
				break;
			}
//			//  I included this to try and remove the "resource leak" warning but doing so
//			//  actually breaks the program - I don't know why
//			finally{                         
//				scan.close();
//				in.close();
//			}
		}
		
		return comps;
	}
	
	



}