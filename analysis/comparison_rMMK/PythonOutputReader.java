import java.io.*;
import java.util.*;



public class PythonOutputReader{




    // Method to check if output from mathematica script has found a solution (true) or generated error message
    //
    //    So far, the only error I have had is that we "encounter a singular Jacobian" - these all result
    //    in the first line of the output file being an empty String... this is all I test for just now...          **  POSSIBLE BUGS HERE... **
    //
    public static boolean solutionExists(String mathOutput)throws IOException{

	boolean solnExists = true;


	BufferedReader in = new BufferedReader( new FileReader(mathOutput)  );

	//  Python error messages don't seem to be printed to file; so an unsuccessful solve attempt
	//  will just leave a blank file

	if( in.readLine() == null ){
	    solnExists = false;
	}
	else{
	    String s = in.readLine();
	    if(  s.equals("")   ){
		solnExists = false;
	    }
	    else if(  s.contains("Traceback")  ){
		solnExists = false;
	    }
	}


	in.close();


	return solnExists;
    }
    //=====================================================================================================================





    // Method to update nodeLabels and nodeConcentrations lists. 
    public static void getLabelsAndConcentrations(String pyOutput, ArrayList<String> nodeLabels, ArrayList<Double> nodeConcentrations ) throws IOException{

	/** Format of pyOutput file is lines like "c124 1.234567e-4"  **/

	nodeLabels.clear();
	nodeConcentrations.clear();

	BufferedReader in = new BufferedReader( new FileReader( pyOutput )  );
	String s = "";
	Scanner scan = new Scanner(s);
	while(true){
	    try{
		s = in.readLine();
		scan = new Scanner(s);

		String cLabel = scan.next();
		//	String arrow = scan.next();
		String pyConc = scan.next();

		/** Get rid of the leading "c" in the Mathematica metabolite label **/
		String label = cLabel.substring( 1 );
		/** Need to convert Mathematica's double rep back to JAVA'a  **/
		double conc = convertStringDouble( pyConc );

		nodeLabels.add( label );
		nodeConcentrations.add( conc );

	    }
	    catch(NullPointerException e){
		    break;
	    }
	}


    } //end getLabelsAndConcentrations() ==========================================================================================================================






    // Method to take mathematica's exponential format: 1.2345*^-7 and convert this back to JAVA's: 1.2345E-7.
    //
    // NOT AN ISSUE IN PYTHON
    private static double convertStringDouble( String mathematicaExp ){
	// if( mathematicaExp.contains("*^") ){
	//     mathematicaExp = mathematicaExp.replace( "*^", "E");
       	// }

	double conc = Double.parseDouble( mathematicaExp );
	return conc;
    }
    // ==================================================================================================================================================







}//end class