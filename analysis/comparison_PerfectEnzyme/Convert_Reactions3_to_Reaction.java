import java.util.*;


public class Convert_Reactions3_to_Reaction{

    /**
       (Needed this new "Reaction" object when I started reversible MMK / common modular rate law)
    **/



    // Convert the "equation" String in old file format (Reactions3) to new file format (Reaction).
    public static String convertRxnString( String rxn_string ){

	String output = "";
	output = rxn_string.replaceAll("H2O","");
	output = output.replaceAll("CO2\\(aq\\)","-co2");
	output = output.replaceAll("NAD_ox","-nad");
	output = output.replaceAll("NAD_red","-nadh");
	output = output.replaceAll("AMP","-amp");
	output = output.replaceAll("ADP","-adp");
	output = output.replaceAll("ATP","-atp");
	output = output.replaceAll("PPi","-ppi");
	output = output.replaceAll("Pi","-pi");

	output = output.replaceAll("NH3\\(aq\\)","-nh3");
	output = output.replaceAll("NH2donor","-glut");
	output = output.replaceAll("NH2acceptor","-oxo");

	output = output.replaceAll("--->", ">");
	output = output.replaceAll("\\+", "");

	return output;
    }




    // MichaelisMenten eqn solver needs list of "Reaction". This method creates "Reaction" object from a "Reactions3" object.
    // Mathematica doesn't accept "CH(OH)-CH3" etc as variable names so we need to replace these               ***  DOES PYTHON?? MIGHT BE BEST TO KEEP INTEGERS THOUGH  ***
    //   with the int label of the substrates
    public static Reaction convert( Reactions3 rxn ){
	
	double dG = rxn.getDeltaG();
	String rxn_string = rxn.getEquation();
	String s = Convert_Reactions3_to_Reaction.convertRxnString( rxn_string );


	// Replace the "CH2(OH)-" names with the integers----------   
	// This only works since we have uni-uni reactions (fine for comparing trunk pathways)  #################### CAREFUL   
	String sub = Integer.toString( rxn.getSubstrate()  );
	String prod = Integer.toString( rxn.getProduct()  );

	String s2 = s;
	Scanner scan2 = new Scanner(s2);	

	boolean subsFound = false;
	String ss = "";
	String reaxn = "";
	while(true){
	    try{
		ss = scan2.next();
		if( ss.contains("C") && ss.contains("H") ){          // ########    ONLY BECAUSE all of our internal metabs contain at least 1 C and 1 H
		    if( !subsFound ){
			reaxn = reaxn + " " + sub;
			subsFound = true;
		    }
		    else{
			reaxn = reaxn + " " + prod;
		    }
		}
		else{
		    reaxn = reaxn + " " + ss;
		}
	    }
	    catch(NoSuchElementException e){break;}
	}
	//--------------------------------------------------------	
	
	Reaction r = new Reaction( dG, reaxn );
	
	return r;
    }






}


