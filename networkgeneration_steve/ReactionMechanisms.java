import java.util.*;

public class ReactionMechanisms {

	
	
	// Separate molecule into groups
	public static ArrayList<String> getGroups( String molecule ){
		ArrayList<String> grps = new ArrayList<String>();
		String g = "";
		for( int i=0; i<molecule.length(); i++ ){
			char c = molecule.charAt( i );
			g = g+c;
			if( c=='='|| c=='-' ){
				grps.add( g );
				g = ""+c;
			}
		}
		// add final group
		grps.add(g);
		return grps;
	}
	
	
	
	
	// Form molecule from groups
	public static String formProd( ArrayList<String> grps  ){
		String prod = "";
		for( String g: grps ){
			prod = prod + g;
		}
		if( prod.contains("-=") || prod.contains("=-")  ){
			System.out.println( "Error: formProd() produced invalid molecule: "+prod );
			//System.exit(0);
		}
		prod = prod.replace( "--", "-");
		prod = prod.replace("==", "=");		
		return prod;
	}

	
	
		
	
}
