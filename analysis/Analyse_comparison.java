import java.io.*;
import java.util.*; 
import java.text.*;
/*

This code reads in the ParameterSampling_combined.dat data, which is of the form
SOURCE PRODUCT ... FLUX[0] FLUX[1]... containing all fluxes for paths which both do and do not use Phydrolysis reaction.

We will have methods to print out:

(1) Cumulative normalised flux for all paths

(2) Parameters, best path (for the colour phase-space plots)

(3) Counter for number of times each path performs absolute best

(4) Flux normalised by [SOURCE] (out6)

//--------
?? (5) Parameters and best flux (to corroborate Barteks maximisation code)

*/
    
    
    public class Analyse_comparison{
	



	public static void main(String args[])throws IOException{

	    //-------------------------------
	    // We first read in the paths being compared so that we can print them
	    // re-ordered wrt whatever we want.
	    //---------------------------------------------
	    // Read the list of reactions corresponding to each pathway and store each list in
	    // a list of pathways

	    BufferedReader reactionsIN = new BufferedReader( new FileReader(args[0]) );
	    
	    ArrayList<ArrayList<Reactions3>> listOfPaths = new ArrayList<ArrayList<Reactions3>>();
	    ArrayList<Reactions3> listOfReactions = new ArrayList<Reactions3>();
	    
	    String s2 = "";
	    Scanner scan2 = new Scanner(s2);
	    while(true){
		try{
		    s2 = reactionsIN.readLine();
		    scan2 = new Scanner(s2);
		    
		    if( s2.contains("Path") ){
			ArrayList<Reactions3> justholding = new ArrayList<Reactions3>();
			for( int qwa=0; qwa < listOfReactions.size(); qwa++ ){
			    justholding.add( listOfReactions.get( qwa ) );
			}
			listOfPaths.add( justholding );
			listOfReactions.clear();
		    }


		    
		    else if( !(s2.contains("Path"))  ){
			String direction = scan2.next();
			String type = scan2.next();
			int substrate = scan2.nextInt();
			int product = scan2.nextInt();
			int ATP = scan2.nextInt();
			int NADH = scan2.nextInt();
			double deltaG = scan2.nextDouble();	    
			//try to read in all the next string pieces to form the equation
			String ss = "";
			while(true){
			    try{
				ss = ss + " " + scan2.next();
			    }
			    catch(NoSuchElementException e){break;}
			}
			String equation = ss;
			listOfReactions.add(new Reactions3(direction, type, substrate, product, ATP, NADH, deltaG, equation));
		    
		    }
		    else{   System.out.println("ERRORZZ WTF!");   }
		}
		catch(NullPointerException e){break;}
		
		
	    }
	    //Now have list of paths (where each path is a list of 'Reactions3' reactions
	    
	    





	    PrintWriter out1 = new PrintWriter( new FileWriter("an_CumNormFlux.dat") );
	    PrintWriter out2 = new PrintWriter( new FileWriter("an_PhaseSpace.dat") );
	    PrintWriter out3 = new PrintWriter( new FileWriter("an_BestPathCounter.dat") );
	    PrintWriter out4 = new PrintWriter( new FileWriter("an_BestFlux.dat"));
	    PrintWriter out5 = new PrintWriter( new FileWriter("an_FluxTot_notNormalised.dat"));
	    PrintWriter out6 = new PrintWriter( new FileWriter("an_Flux_normalisedBySOURCE.dat") );


	    // Holder variables:
	    ArrayList<Double> cumNormFlux = new ArrayList<Double>();
	    ArrayList<Integer> bestCount = new ArrayList<Integer>();
	    ArrayList<Double> totFlux = new ArrayList<Double>();
	    ArrayList<Double> fluxOverSource = new ArrayList<Double>();



	    BufferedReader data = new BufferedReader( new FileReader( "an_ParameterSampling_x.dat") );
	    String s = "";
	    Scanner scan = new Scanner(s);

	    int line = 0;
	    while(true){

		line++;

		try{
		    s = data.readLine();
		    scan = new Scanner(s);

		    // Read in all the parameters
		    double SOURCE = scan.nextDouble();
		    double PRODUCT = scan.nextDouble();
		    double ATP = scan.nextDouble();
		    double ADP = scan.nextDouble();
		    double AMP = scan.nextDouble();
		    double NAD = scan.nextDouble();
		    double NADH = scan.nextDouble();
		    double HPO4 = scan.nextDouble();
		    double PPi = scan.nextDouble();
		    double CO2 = scan.nextDouble();
		    double NH3 = scan.nextDouble();
		    double NH2donor = scan.nextDouble();
		    double NH2acceptor = scan.nextDouble();

		    //ArrayList containing the rest of the entries
		    //which are the fluxes of all paths
		    ArrayList<Double> fluxes = new ArrayList<Double>();
		    while(true){
			try{
			    fluxes.add( scan.nextDouble() );
			}
			catch(NoSuchElementException e){break;}
		    }


		    // ---  Do comparison here, before reading in next data line ----

		    // ###  1  #######   Do cumulative normalised flux stuff  #######################
		    double bestFlux = getBestFlux( fluxes );
		    if( line == 1 ){
			for( int j=0; j<fluxes.size(); j++){

			    // Here we ignore first line and just use it to create correct sized lists
			    // for cumNormFlux, FluxTot, and bestCount
			    cumNormFlux.add(    0.0     );
			    totFlux.add(    0.0     );
			    bestCount.add(  0  );
			    fluxOverSource.add( 0.0 );
			}
		    }
		    else{
			for( int j=0; j<fluxes.size(); j++ ){

			    // Only consider pathways that have a positive flux ----- IMPORTANT
			    if( fluxes.get(j) > 0 ){

				cumNormFlux.set( j, ( cumNormFlux.get(j) + (fluxes.get(j)/bestFlux) )   );
				totFlux.set(j, (totFlux.get(j) + fluxes.get(j))  );
				fluxOverSource.set( j,  fluxOverSource.get(j) + ( fluxes.get(j)/SOURCE )   );
			    }
			    // ------------------------------------------------------ IMPORTANT

			}
		    }




		    // ###  2  #######  Do phase space stuff here  ##############################
		    int bestPath = getBestPath(  fluxes  );

		    if( bestPath != -1 ) {
			bestCount.set(  bestPath,  (bestCount.get( bestPath ) + 1)   );
		    }

		    out2.println( SOURCE + "\t" + PRODUCT + "\t" + ATP + "\t" +ADP + "\t" + AMP + "\t" + NAD + "\t" + NADH + "\t" + HPO4 +
				  "\t" + PPi + "\t" + CO2 + "\t" + NH3 + "\t" + NH2donor + "\t" + NH2acceptor + "\t" + bestPath);

		    out4.println( SOURCE + "\t" + PRODUCT + "\t" + ATP + "\t" +ADP + "\t" + AMP + "\t" + NAD + "\t" + NADH + "\t" + HPO4 + "\t" + PPi +
				  "\t" + CO2 + "\t" + NH3 + "\t" + NH2donor + "\t" + NH2acceptor + "\t" + bestPath + "\t" + bestFlux);


		}// Stop reading in lines
		catch(NullPointerException e){break;}

	    }


	    for( int j=0; j<cumNormFlux.size(); j++){
		out1.println(  j + " " + ( cumNormFlux.get(j)/(double)line )  );     
		//divide by "line" - # of parameter sets - so that range of cumNormFlux plot is 0->1
	    }
	    for( int j=0; j<totFlux.size(); j++){
		out5.println(  j + " " + ( totFlux.get(j) )  );     
	    }
	    for( int j=0; j<bestCount.size(); j++){
		out3.println( j + " " + bestCount.get(j)  );
	    }
	    for( int j=0; j<fluxOverSource.size(); j++){
		out6.println( j + " " + fluxOverSource.get(j)  );	
	    }


	    out1.close();   out2.close();   out3.close();   out4.close();   out5.close();   out6.close();


	    //================================================================================================

	    //---------------------------------------------------------------------------------------------
	    //   HERE I PRINT OUT ALL THE ABOVE GRAPHS IN THE NEW ORDER (w.r.t. totFlux or CF high to low)
	    //           (don't really need the above figures anymore)
	    //---------------------------------------------------------------------------------------------

	    PrintWriter out8 = new PrintWriter( new FileWriter("an_REORDERED_paths.paths") );
	    // and print out the rest too...
	    PrintWriter out7 = new PrintWriter( new FileWriter("an_FluxTot_notNormalised_REORDERED.dat") );
	    PrintWriter out9 = new PrintWriter( new FileWriter("an_CumNormFlux_REORDERED.dat") );
	    PrintWriter out10 = new PrintWriter( new FileWriter("an_BestPathCounter_REORDERED.dat") );
	    PrintWriter out11 = new PrintWriter( new FileWriter("an_Flux_normalisedBySOURCE_REORDERED.dat") );

	    // First we need to find the new order wrt the desired measure

	    // // This reorders list of paths wrt to **TOTAL FLUX**
	    //++++++++++++++++++++++++++++++++++++++++++++++++++
	    // int numPaths = totFlux.size();
	    // // Make tempArray so 'totFlux' is not being altered.
	    // ArrayList<Double> tempArray = new ArrayList<Double>();
	    // for( int i=0; i<numPaths; i++){
	    // 	tempArray.add(   totFlux.get(i)   );
	    // }
	    // ArrayList<Integer> newOrder = new ArrayList<Integer>();
	    // //Reorder pathways wrt "total flux"
	    // for( int i=0; i < numPaths; i++ ){
	    // 	int biggest = getBiggestFlux(  tempArray );
	    // 	newOrder.add( biggest );
	    // 	tempArray.set( biggest, -2000.0 );
	    // }
	    //++++++++++++++++++++++++++++++++++++++++++++++++++
	    
	    // This reorders list of paths wrt to **CUMULATIVE NORMALISED FLUX**
    	    //++++++++++++++++++++++++++++++++++++++++++++++++++
	    int numPaths = totFlux.size();
	    // Make tempArray so 'cumNormFlux' is not being altered.
	    ArrayList<Double> tempArray = new ArrayList<Double>();
	    for( int i=0; i<numPaths; i++){
	    	tempArray.add(   cumNormFlux.get(i)   );
	    }
	    ArrayList<Integer> newOrder = new ArrayList<Integer>();
	    //Reorder pathways wrt "cumulative normalised flux"
	    for( int i=0; i < numPaths; i++ ){
	    	int biggest = getBiggest_CF(  tempArray );
	    	newOrder.add( biggest );
	    	tempArray.set( biggest, -2000.0 );
	    }
	    //++++++++++++++++++++++++++++++++++++++++++++++++++

	    //    /*********** This chunk allows me to set order by hand so when plotting
	    //                     comparison at diff params / without concentration restrictions,
            //        the CF, flux plots etc will all be ordered the same as in orginal
            //         comparison.

	    // // OK, so now we have newOrder list, ordered wrt to whatever we did above.
	    // // If I want to re-run the comparison with different parameter ranges, for instance
	    // // without restrictions on concentrations, then the paths will be ordered 
	    // // differently than they were with the base set.
	    // // Here I just override the ordering, and set it to what I want, BY HAND.
	    // newOrder.clear();
	    // // GLYCOLYSIS ORDERED WRT CF USING BASE SET:
	    // //	    int[] orderByHand = {18, 19, 7, 0, 6, 4, 5, 22, 17, 10, 21, 20, 3, 2, 1, 16, 12, 13, 14, 15, 25, 9, 11, 8, 23, 24};
	    // // OR GLUCONEOGENESIS:
	    // //	    int[] orderByHand = {14, 9, 21, 19, 18, 8, 17, 4, 11, 33, 22, 2, 32, 34, 0, 39, 5, 3, 1,  37, 30, 35, 40, 29, 24, 31, 23, 38, 36, 27, 25, 28, 26, 6, 10, 12, 7, 16, 15, 13, 20};

	    // // All 106 GLYCOLYTIC paths (with AMP and PPi)
	    // // int[] orderByHand = {81, 85, 82, 17, 46, 48, 0, 86, 14, 49, 51, 10, 94, 11, 79, 36, 16, 15, 55, 54, 18, 20, 91, 88, 58, 21, 12, 96, 23, 7, 57, 56, 95, 89, 92, 4, 1, 60, 9, 8, 37, 13, 39, 84, 2, 5, 87, 32, 33, 93, 30, 31, 6, 47, 90, 3, 50, 19, 59, 38, 80, 83, 22, 72, 64, 66, 68, 73, 69, 70, 67, 103, 62, 63, 41, 42, 71, 65, 104, 26, 43, 27, 44, 105, 24, 25, 28, 29, 34, 35, 40, 45, 52, 53, 61, 74, 75, 76, 77, 78, 97, 98, 99, 100, 101, 102};


	    // // The best 25 Glyc paths, from the above 106
	    // //  int[] orderByHand = {0,   1,  2,   3,  6,  8,  4, 7,  5, 11, 9,  10, 12, 13, 15, 14, 20, 16, 21, 17, 18, 19, 22, 24, 23};
	        
	    // // FINAL:  Glyc, best 25 from all 256/771 paths that had ATP=2  -----------------
	    // //  int[] orderByHand = {0,   1,  2,   3,  4,  6,  7, 5,  8, 9,  10, 11, 12, 13, 14, 16, 15, 18, 17, 19, 20, 21, 22, 23, 24};

	    // //int[] orderByHand = {0, 1, 2, 3, 4, 6, 7, 5, 8, 9, 10, 11, 12, 13, 14, 16, 15, 18, 17, 19, 20, 21, 22, 23, 24};
	    // //    int[] orderByHand = {0, 1, 2, 3, 4, 5, 6, 8, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
	    // //    int[] orderByHand = {0, 1, 2, 4, 3, 8, 5, 6, 7, 9, 10, 13, 11, 12, 16, 14, 15, 17, 18, 20, 22, 19, 21, 24, 23};

	    // // These are the orders from Powell method, glyc and gluc resp.
	    // //int[] orderByHand =  {0, 1, 2, 3, 6, 5, 8, 4, 10, 7, 9, 15, 11, 12, 13, 18, 14, 16, 17, 20, 24, 19, 22, 21, 23};
       	    // int[] orderByHand = {0, 1, 3, 2, 4, 5, 6, 7, 8, 9, 10, 12, 11, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};

	    // for( int kk=0; kk<orderByHand.length; kk++){
	    // 	newOrder.add( orderByHand[kk] );
	    // }
	
	    // ///   *************/

	    //Now take this newOrder list of integers and print out the reactions
	    //and REORDERED path numbers
	    for( int i=0; i<newOrder.size(); i++){
		
	      	//Print out totFluxes REORDERED
		out7.println( i + " " + totFlux.get( newOrder.get(i) )  );
		// and the rest
		out9.println( i + " " + cumNormFlux.get( newOrder.get(i) )/(double)line  );              // ******* THIS SHOULD HAVE BEEN DIVIDED BY THE NUMBER OF SAMPLINGS!!! ********
		out10.println( i + " " + bestCount.get( newOrder.get(i) )  );
		out11.println( i + " " + fluxOverSource.get( newOrder.get(i) )  );


		//Print lines for Gnu scripts  (need to uncomment these one at a time)
		 System.out.println(  "awk '{ if($4=="+newOrder.get(i)+") print $1, $2, $3, \"0x7F7F7F\"; else print $0}' | \\"  );
		//      System.out.println(  "ROW"+i +"(x,y)=(x eq \""+ i + "\") ? y:1/0"  );
		//	System.out.println( "\t"+"\'\' u ($0):(ROW"+i+"(stringcolumn(1),$2)):xtic(1) w boxes lt 1 lc rgb \"#7F7F7F\" ti \""+ i +"\" ,\\" );  // edit this for odd/even lines


		// Print list of paths in new order
		ArrayList<Reactions3> path = listOfPaths.get( newOrder.get(i) );
		for( int h=0; h < path.size(); h++){
		    out8.println( path.get( h ) );
		}
		out8.println( "Path_" + i + " ends here");

	    }

	    out7.close();   out9.close();   out10.close();   out11.close();
	    out8.close();
	    //=========================================================================------------------------------------------


	    System.out.println( newOrder.toString()  );






	}//end main















	//Method to get biggest TOTAL FLUX
	public static int getBiggestFlux( ArrayList<Double> totFlux ){

	    double biggestFlux = -999.0;
	    int place = -3;
	    for(int i=0; i<totFlux.size(); i++){
		if( totFlux.get(i) > biggestFlux ){
		    biggestFlux = totFlux.get(i);
		    place = i;
		}
	    }
	    if( place == -3 ){ System.out.println( "getBiggestFlux() ERROR"); }
	    return place;
	}//-------------------------------------------------------------------

	//Method to get biggest CUMULATIVE NORMALISED FLUX
	public static int getBiggest_CF( ArrayList<Double> cumNormFlux ){

	    double biggestCF = -999.0;
	    int place = -3;
	    for(int i=0; i<cumNormFlux.size(); i++){
		if( cumNormFlux.get(i) > biggestCF ){
		    biggestCF = cumNormFlux.get(i);
		    place = i;
		}
	    }
	    if( place == -3 ){ System.out.println( "getBiggest_CF() ERROR"); }
	    return place;
	}//--------------------------------------------------------------------









	public static double getBestFlux( ArrayList<Double> fluxes ){
	    double bestFlux = -999.0;
	    int bestPath = -444;
	    for( int i=0; i<fluxes.size(); i++){
		if( fluxes.get(i) > bestFlux ){
		    bestPath = i;
		    bestFlux = fluxes.get(i);
		}
	    }

	    if( bestFlux < 0 ){ bestFlux = -1; }

	    return bestFlux;
	}
	//-------------------------------------------------------------------
	public static int getBestPath( ArrayList<Double> fluxes ){
	    double bestFlux = -999.0;
	    int bestPath = -444;
	    for( int i=0; i<fluxes.size(); i++){
		if( fluxes.get(i) > bestFlux ){
		    bestPath = i;
		    bestFlux = fluxes.get(i);
		}
	    }
	    if( bestFlux < 0 ){ bestPath = -1; }

	    return bestPath;
	}
	//---------------------------------------------------------------------









       
    }//end class


