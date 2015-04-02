import java.io.*;
import java.util.*; 
import java.text.*;

/**
Use find_maximum() to compute the optimized flux through a pathway, subject to total enzyme concentration being fixed 
and all internal metabolite concentrations falling within allowed range
 */



public class PowellOptimization_relativeTolerance{


    // Define global cellular condition parameters:
    static double concSOURCE, concPRODUCT, kd, RT; 
    // If we ever change the temperature of the network watch out for these RTs... they're everywhere...

    // Store all external metabolite concentrations here
    static HashMap<String, Double> extMetabConcentrations = new HashMap<String, Double>();

    // Fixed total enzyme concentration:
    static double TOTAL_ENZYME;

    // Number of variables (number of enzymes MINUS ONE)
    static int nvars;   
    
    // Array of variables:
    static double[] vars = new double[nvars];   
    
    // These were defined in original Powell maximization code
    public static int nmax = nvars;
    public static double[]  dir=new double[nmax],   pkt=new double[nmax];

    // Metabolite restrictions:
    static double conc_min, conc_max; 





    // Method to calculate flux, setting it to zero is any conditions are broken
    static double FluxFunction( int nvar, double[] vars,  ArrayList<Reaction> reactionsUsed ){
       
	double flux = 0;

	// Need to calculate the final concentration
	double[] enzymeConcentrations = new double[ nvar + 1 ];  // This is what goes into the flux & conc methods - NOT "vars[]"
	double sumVars = 0;
	for( int i=0; i<vars.length; i++ ){
	    sumVars += vars[i];
	    enzymeConcentrations[i]  =  vars[i];
	}
	double finalConc =  0.99*(TOTAL_ENZYME - sumVars);   // 0.99; debug; trying to find why flat profile beats optimized sometimes
	enzymeConcentrations[ nvar ] = finalConc;


	// Check that no enzyme concentrations are below zero
	boolean enzymesPositive = true;
	for(int ii=0; ii < enzymeConcentrations.length; ii++){   
	    if( enzymeConcentrations[ii] <= 0 ){ 
		enzymesPositive = false;
	    }
	}

	// Check that they are all positive  ------ this is OFTEN not this case... I think this method is very inefficient.
	if( !enzymesPositive ){
	    flux = -1E-20;
	    // System.out.println("Negative [Ei]");
	}
	// Else continue as before
	else{

	    //calculate flux
	    flux = PerfectEnzyme.getFlux( reactionsUsed, extMetabConcentrations, concSOURCE, concPRODUCT, enzymeConcentrations,  kd, RT   );

	    //calculate steady state concentrations, {[Sj]}
	    ArrayList<Double> concentrations = PerfectEnzyme.getConcentrations( flux, reactionsUsed, extMetabConcentrations, concSOURCE, concPRODUCT, enzymeConcentrations, kd, RT  );
	    
	    //Now, if any of the concentrations fall outside allowed range, set flux to zero.
	    for( double conc : concentrations  ){
		if( conc < conc_min || conc > conc_max ){
		    flux = -1E-20;
		    break;
		}
	    }
    
	}
	
	return flux;
    }
    // End FluxFunction() ==============================================================





    static double ff(int n,double x, ArrayList<Reaction> reactionsUsed  ){

	// The ff() function is to update all the parameter values v[i] = dir[i]*x + pkt[i]:
	//     i.e. It's current value + the direction of increase*ratio 
	// and return the value of the objective function at this new set of parameters
	
	double[] v=new double[nmax] ;

	//? Why would this happen?
	if (n>nmax) System.exit (1) ; 

	for (int i=0;i<n;i++){
	    v[i]=dir[i]*x+pkt[i];          // ### CAREFUL ###
	}

	return FluxFunction(n, v, reactionsUsed  ) ;
    }




    static double golden(int nvar,double a,double b,double c, ArrayList<Reaction> reactionsUsed ){

	// This golden() method is the "golden section search", the limit of the Fibonacci search for a large number of
	// function evaluations.  It is a technique for finding the extremum of a strictly unimodal function by successively narrowing the 
	// range of values inside which the extremum is known to exist.

	// When called from linmax() values were hard-coded as a=-2, b=0 and c=2    (why this choice?)

	double fa,fc,fb,d,fd,b0,r=(3.0-Math.sqrt(5))/2.0; // This is not the golden ratio, gr, it is (2-gr), see wikipedia page.

	int n=0 ;
	// ff() returns objective function after updating parameter values via for example:  newValue = oldValue + (dir * a)
	fa=-ff(nvar,a, reactionsUsed ); 
	fb=-ff(nvar,b, reactionsUsed ) ; 
	fc=-ff(nvar,c, reactionsUsed ) ;


	// God only knows - think about this later
	do
	    {
		b0=b ; n++ ;
		if (c-b>b-a)
		    {
			d=b+r*(c-b) ; 
			fd = -ff( nvar,d, reactionsUsed ) ;
			if (fd>fb) { c=d ; fc=fd ; } else { a=b ; b=d ; fa=fb ; fb=fd ; }
		    }
		else
		    {
			d=b-r*(b-a) ;
			fd = -ff( nvar,d, reactionsUsed  ) ;
			if (fd>fb) { a=d ; fa=fd ; } else { c=b ; b=d ; fc=fb ; fb=fd ; }
		    }

		//////////////////////////////////////////////////
	    } while (  Math.abs(b0-d) > 1e-7  ) ;
	/////////////////////////////////////////////////////////////


	// Return d, which is the distance, along the vector of maximum increase, in which to travel
	//   (this will allow us to update all parameters, as in ff(), by setting them to ( pkt[i] + dir[i]*d )
	return (d) ;
    }
    // =============================  golden()  ==========================================



    static double linmax( double p[], double xi[], int n,  ArrayList<Reaction> reactionsUsed  ){

	double x;	int i;

	// Store the parameter values "p[]" and direction vector "xi[]" in the static variables "pkt[]" and "dir[]"
	for (i=0;i<n;i++) pkt[i]=p[i] ;
	for (i=0;i<n;i++) dir[i]=xi[i] ;


	//  I do not understand why smaller brackets perform so much worse???  "1.0" is ridiculously huge.
	//////////////////////////////////////////////////////////////
       	x = golden( n,   -0.1, 0, 0.1,    reactionsUsed   ) ;
	//	x = golden( n,   -0.001, 0, 1.001,    reactionsUsed   ) ;
	//////////////////////////////////////////////////////////////


	for (i=0;i<n;i++) {
	    xi[i]*=x ;   
	    p[i]+=xi[i] ;    // ## CAREFUL ##
	}

	return FluxFunction( n, p, reactionsUsed  ) ;
    }
    // linmax()  ========================================================================




    //  static double powell(double p[],int n,double ftol,    ArrayList<Integer> path,  double[] reactionG, String[] forwBack, String[] reactionType ) {
    static double powell( double p[], int n,double ftol, ArrayList<Reaction> reactionsUsed ) {


	int itmax = 10000;        
	int i,j, ibig, iter;

	double[][] xi=new double[nmax][nmax] ;
	double[] pt=new double[nmax],           ptt=new double[nmax],     xit=new double[nmax] ;
	double fp,del,fptt,t,fret;

	// Initialize "xi[]" as the unit matrix
        for (i=0;i<n;i++) { for (j=0;j<n;j++) xi[i][j]=0 ; xi[i][i]=1 ; }

	// Find objective function using parameters (vars) in "p[]"
	fret=FluxFunction( n,  p,  reactionsUsed  ) ;

	// Initialize "pt[]" as "p[]", i.e. the variables fed into powell method
	for (j=1;j<=n;j++) pt[j-1]=p[j-1] ;

	// Loop until solution is found, or max number of iterations is reached
	iter=0 ;
	for (;;) {
	    iter++ ; 
	    fp=fret ; 
	    ibig=0 ; del=0 ;

	    // Set "xit[]" vector as values in each row of "xi[][]" matrix
	    for (i=1;i<=n;i++) {
		for (j=1;j<=n;j++) {
		    xit[j-1]=xi[j-1][i-1] ;
		}

		// Feed initial parameters "p[]" and "xit[]" into linmax() method
		fret = linmax( p,xit,n, reactionsUsed ) ;

		// If change in objective function is bigger than "del", store "del", and store this row number "i" as "ibig"
		//	if (Math.abs(fp-fret)>del) { del=Math.abs(fp-fret) ; ibig=i ; }
		if ( (Math.abs(fp-fret)/fret)*100 > del ) { del= (Math.abs(fp-fret)/fret)*100;  ibig=i;    }
		//########################### ABOVE HAS BEEN CHANGED TO RELATIVE MEASURE TOO #################################
	    }


	    //	    if ( (2.0*Math.abs(fp-fret) <= ftol )   && (iter>2) ) return fret; #######################################################
	    //        ALTERED THIS FOR A RELATIVE SENSITIVITY SINCE FLUX CAN VARY OVER SEVERAL ORDERS OF MAGNITUDE 
	    double relChange = (Math.abs(fp-fret)/fret)*100;
	    if(  (relChange <= ftol)  &&  (iter>2)  )  return fret;
	    // ####################################################################################################################


	    if (iter>=itmax) return fret; //{ exit (2) ; } // exit if too many iterations    // -- (Bartek commented this out, not me.)

	    // If not, then:
	    for (j=1;j<=n;j++) {
		// Set "ptt[]" parameter values as (2*new - old)  [ NB: linmax() can alter "p[]" ] ??????
		ptt[j-1]=2*p[j-1]-pt[j-1];   // ## CARFEFUL ##

		// Store size of change of each parameter in "xit[]"  (= oldValue - newValue)
		xit[j-1]=p[j-1]-pt[j-1]; 

		// Update "pt[]" param values with those in "p[]"
		pt[j-1]=p[j-1];
	    }

	    fptt = FluxFunction( n, ptt, reactionsUsed ) ;

	    // NOTE: the continue statement means that the code following it will not be executed. Program will jump
	    //   back up and "continue" from the start of the for(;;) loop.

	    if (fptt>=fp) continue ;
	    // So, if latest objective is smaller than the last value,  calculate "t"  using this:  (???)
	    t=2.0*(fp-2.0*fret+fptt)*(fp-fret-del)*(fp-fret-del)-del*(fp-fptt)*(fp-fptt) ;

	    if (t>=0) continue ;
	    // And if t<0, put "p[]" and "xit[]" back into linmax() function 
	    fret = linmax( p,xit,n, reactionsUsed ) ;

	    // Update the 'big column' in "xi[][]" with values in "xit[]" vector ( which were just updated in linmax() )
	    for (j=1;j<=n;j++) xi[j-1][ibig-1]=xit[j-1] ;
	}

    }
    // End of powell() ===========================================================================
    







    static double find_maximum(int _nvars,  double varmin, double varmax, ArrayList<Reaction> reactionsUsed, HashMap<String, Double> _extMetabs, double _concSOURCE, double _concPRODUCT, double _totEnzymeConc, double _metabConc_min, double _metabConc_max, double _kd, double _RT, int _startPoints ){
	/* This is the method we call from another program */


	// Initialize global variables:

	concSOURCE=_concSOURCE;        concPRODUCT=_concPRODUCT;     
	extMetabConcentrations = new HashMap<String, Double>(  _extMetabs  );
	kd=_kd;     RT=_RT;
	conc_min = _metabConc_min;     conc_max = _metabConc_max; 
	TOTAL_ENZYME = _totEnzymeConc;
	nvars=_nvars;
	//-------------------------------
	// RE-INITIALIZE global variables that depend on "nvars":
	// Array of variables:
	vars = new double[nvars];   
	// These were defined in the original optimization code
	nmax = nvars;
	dir=new double[nmax];   pkt=new double[nmax];
	


	int iter ;
	double y,ymax=-1e+12 ;
	double[] var0=new double[nvars] ;

	for (int st=0;   st < _startPoints;   st++) { 
	    for (int i=0; i<nvars; i++){

		// log-uniform sampling
		double qq = Math.log( varmin ) + ( Math.log( varmax ) - Math.log( varmin )  )*Math.random();
	        vars[i] = Math.exp( qq );
	    }


	    y = powell( vars, nvars, 1e-3, reactionsUsed ) ;    // NOTE:


	    if (y>ymax) {
		ymax = y;
		for (int i=0;i<nvars;i++) var0[i]=vars[i] ;
	    }

	}

	for (int i=0;i<nvars;i++) vars[i]=var0[i] ;


	//	System.out.println( "Maximized flux = " + ymax );
	//	System.out.println( "Optimal {[Ei]}  = " + Arrays.toString( vars ) );


	return ymax ;
    } 
    //=============================================================






   
}//end class
