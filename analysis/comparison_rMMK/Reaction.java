import java.util.*;
import java.io.*;

public class Reaction{

    // Reaction: deltaG value and a string like "A > B C"

    // Set these as the default values of enzyme parameters


    /** Define the Reaction parameters individually so we can mutate these in the future */
    /** Vf has been replaced with (k_cat*Ei)  */
    private double kcat;
    /** Enzyme concentration **/
    private double Ei;
    //    private double Vf=1.0E-3;
    /** MM const of (internal) substrates **/
    private double Ks;
    /** MM const of (internal) products **/
    private double Kp;



    //  private double kcat=10;
    //private double Ei=1.0E-4;
    //  private double Ks=1.0E-6;
    // private double Kp=1.0E-3;

    //   /** Diffusion-controlled upper limit **/
    //   private double kd=1.0E8;

    /** NB: having private variables means you need to define a getter method **/
    private double deltaG;     
    private String reaction;


    private final double valueOfAllExternalMMConst = 1.0E-5;
    /** We now include all external metabolites in the kinetics **/
    private HashMap<String, Double> externalMetabs_MM_constants;


    /** Having these lists comes in handy **/
    //    private 
    ArrayList<String> int_subs;
    //    private 
    ArrayList<String> int_prods;
    //    private                        I made this no longer private so that ReactionReorder.java could access them... though I should have written getter methods...
    ArrayList<String> ext_subs;
    //    private 
    ArrayList<String> ext_prods;




    //  Constructor =====================================================
    Reaction( Reaction r ){

	this.deltaG = r.getDeltaG();
	this.reaction = r.getReaction();

	//double vf=r.getVf(); 
	double ks=r.getKs();       double kp=r.getKp();
	double ei=r.getEi();  double kcattt=r.getkcat();  //double kdd=r.getkd();
	this.Ei = ei;
	this.kcat = kcattt;
	//	this.Vf = vf;
	this.Ks = ks;
	this.Kp = kp;

	//	this.kd = kdd;

	// Use reaction string to make seperate lists of int/ext subs/prods:
	this.int_subs = new ArrayList<String>();
	for( int i=0; i < r.int_subs.size(); i++ ){
	    this.int_subs.add(  r.int_subs.get(i)  );
	}

	this.int_prods = new ArrayList<String>();
	for( int i=0; i < r.int_prods.size(); i++ ){
	    this.int_prods.add(  r.int_prods.get(i)  );
	}

	this.ext_subs = new ArrayList<String>();
	for( int i=0; i < r.ext_subs.size(); i++ ){
	    this.ext_subs.add(  r.ext_subs.get(i)  );
	}

	this.ext_prods = new ArrayList<String>();
	for( int i=0; i < r.ext_prods.size(); i++ ){
	    this.ext_prods.add(  r.ext_prods.get(i)  );
	}

	HashMap<String, Double> hmp = r.getExtMetabMMConst();

	this.externalMetabs_MM_constants = new HashMap<String, Double>();
	for( Map.Entry<String, Double> entry:  hmp.entrySet()  ){
	    String key = entry.getKey();
	    double value = entry.getValue();
	    this.externalMetabs_MM_constants.put( key, value );
	}



    }
    //
    // Constructor =======================================================
    Reaction( double deltaG, String rxn ){ 

	this.deltaG = deltaG;
	this.reaction = rxn;


	/**  Now we just store lists of internal/external subs and prods to make some methods I use easier  **/
	this.int_subs = new ArrayList<String>();
	this.int_prods = new ArrayList<String>();
	this.ext_subs = new ArrayList<String>();
	this.ext_prods = new ArrayList<String>();

	Scanner scan = new Scanner( rxn );
	boolean productsreached = false;
	while(true){
	    try{
		String nextbit = scan.next();
		if( !productsreached ){
		    if( nextbit.equals(">") ){
			productsreached = true;
		    }
		    else if( nextbit.contains("-") ){
			this.ext_subs.add( nextbit );
		    }
		    else{
			this.int_subs.add( nextbit );
		    }
		}
		else{
		    if( nextbit.contains("-") ){
			this.ext_prods.add( nextbit );
		    }
		    else{
			this.int_prods.add( nextbit );
		    }
		}		
	    }
	    catch(NoSuchElementException e){break;}
	    
	    /** HashMap contains name and Michaelis-Menten constant of all external metabolites for now, set all to same value **/
	    this.externalMetabs_MM_constants = new HashMap<String, Double>();
	    for( String em: ext_subs ){
		this.externalMetabs_MM_constants.put( em, valueOfAllExternalMMConst );
	    }
	    for( String em: ext_prods ){
		this.externalMetabs_MM_constants.put( em, valueOfAllExternalMMConst );
	    }


	}// end while
      
    }// End constructor=====================================





    // Method to return number of ATP molecules *PRODUCED*
    // ATP is metabolite "-6" in Bartek's network. Currently all reactions only consume/produce 1 ATP maximum.
    public int getATP(){
	int numATP=99999999;
	if(  this.ext_subs.contains("-atp") && this.ext_prods.contains("-atp") ){
	    System.out.println("ERROR: rxn "+this.reaction+ "  has ATP a sub and prod. ");
	}
	else if( this.ext_subs.contains("-atp")  ){
	    numATP = -1;
	}
	else if(  this.ext_prods.contains("-atp")  ){
	    numATP = 1;
	}
	else{
	    numATP = 0;
	}
	return numATP;
    }


    // Method to return number of NADH molecules *PRODUCED*
    // NADH is metabolite "-4" in Bartek's network. Currently all reactions only consume/produce 1 NADH maximum.
    public int getNADH(){
	int numNADH=99999999;
	if(  this.ext_subs.contains("-nadh") && this.ext_prods.contains("-nadh") ){
	    System.out.println("ERROR: rxn "+ this.reaction+ "  has NADH a sub and prod. ");
	}
	else if( this.ext_subs.contains("-nadh")  ){
	    numNADH = -1;
	}
	else if(  this.ext_prods.contains("-nadh")  ){
	    numNADH = 1;
	}
	else{
	    numNADH = 0;
	}
	return numNADH;
    }







    // Getters
    public double getDeltaG(){
	return this.deltaG;
    }
    public String getReaction(){
	return this.reaction;
    }
    // public double getVf(){
    // 	return this.Vf;
    // }
    public double getKs(){
	return this.Ks;
    }
    public double getKp(){
	return this.Kp;
    }
    // public double getkd(){
    // 	return this.kd;
    // }
    public double getEi(){
	return this.Ei;
    }
    public double getkcat(){
	return this.kcat;
    }
    //
    ////////////////////////////////  CHANGED 9/9/2014  ///////////////////////////
    //
    public ArrayList<String> getIntSubs(){
	ArrayList<String> q = new ArrayList<String>();
	for( String s: this.int_subs ){
	    q.add( s );
	}
	return q;
	//
	//	return this.int_subs;
	//
    }
    public ArrayList<String> getIntProds(){
	ArrayList<String> q = new ArrayList<String>();
	for( String s: this.int_prods ){
	    q.add( s );
	}
	return q;
	//
	//	return this.int_prods;
	//
    }
    public ArrayList<String> getExtSubs(){
	ArrayList<String> q = new ArrayList<String>();
	for( String s: this.ext_subs ){
	    q.add( s );
	}
	return q;
	//
	//	return this.ext_subs;
	//
    }
    public ArrayList<String> getExtProds(){
	ArrayList<String> q = new ArrayList<String>();
	for( String s: this.ext_prods ){
	    q.add( s );
	}
	return q;
	//
	//	return this.ext_prods;
	//
    }
    public HashMap<String, Double> getExtMetabMMConst(){
	return this.externalMetabs_MM_constants;
    }


    // Setters
    // public void setVf( double v ){
    // 	this.Vf = v;
    // }
    public void setKs( double k ){
	this.Ks = k;
    }
    public void setKp( double k ){
	this.Kp = k;
    }
    public void setkcat( double k ){
	this.kcat = k;
    }
    public void setEi( double Ei ){
	this.Ei = Ei;
    }




    //Method to check if the Reaction contains a given metabolite
    public boolean containsMetab(String label){

	boolean contains = false;
	for( int i=0; i<this.int_subs.size(); i++){
	    if( this.int_subs.get(i).equals( label ) ){
		contains = true;
	    }
	}
	for( int i=0; i<this.int_prods.size(); i++){
	    if( this.int_prods.get(i).equals( label ) ){
		contains = true;
	    }
	}
	/** Also check the external metabolites **/
	for( int i=0; i<this.ext_subs.size(); i++){
	    if( this.ext_subs.get(i).equals( label ) ){
		contains = true;
	    }
	}
	for( int i=0; i<this.ext_prods.size(); i++){
	    if( this.ext_prods.get(i).equals( label ) ){
		contains = true;
	    }
	}
	return contains;
    }
    //=============================================================



    public String toString(){
	String output = deltaG + "\t" + reaction + "\t" + " Ks="+Ks + ", Kp="+Kp + ", Ei="+Ei + ", kcat="+kcat;
	return output;
    }



    // Method to return a list of internal metabolites used in a reaction
    public ArrayList<String> getReactionMetabolites(){
	ArrayList<String> metabolites = new ArrayList<String>();
	String s = this.reaction;
	Scanner scan = new Scanner(s);
	while(true){
	    try{
		String nextbit = scan.next();		
		if( !nextbit.equals(">")  &&  !nextbit.contains("-") ){
		    // Then we have an internal metabolite
		    metabolites.add( nextbit );
		}
	    }
	    catch(NoSuchElementException e ){break;}
	}
	return metabolites;
    }





    // Overriding the equals method for our Reaction object
    @Override   /** When overriding a method, putting this notation here will cause compile-time error
		    if you do it wrong!  **/ 
    public boolean equals(Object o){

	Reaction a = (Reaction) o;
	boolean equals = true;

	// Want to check if Reaction objects are equal - but note that "-1.0 A B > D E"
	// could be reversed to "1.0 D E > A B"   **or**  "1.0 E D > B A". Want an equals()
	// method that takes this into account.

	ArrayList<String> aIntSubs = a.getIntSubs(); ArrayList<String> aIntProds = a.getIntProds();
        ArrayList<String> aExtSubs = a.getExtSubs(); ArrayList<String> aExtProds = a.getExtProds();
	int numMetabs = aIntSubs.size()+aIntProds.size()+aExtSubs.size()+aExtProds.size();
	int numMetabsTHIS = this.int_subs.size()+this.int_prods.size()+this.ext_subs.size()+this.ext_prods.size();

	if(  numMetabs != numMetabsTHIS ){    equals = false;    }

	for( String ss: this.int_subs ){
	    if(  !aIntSubs.contains( ss ) && !aIntProds.contains(ss)  ){
		equals = false;
	    }
	}
	for( String ss: this.int_prods ){
	    if(  !aIntProds.contains( ss ) && !aIntSubs.contains(ss)   ){
		equals = false;
	    }
	}
	for( String ss: this.ext_prods ){
	    if(  !aExtProds.contains( ss )  &&  !aExtSubs.contains(ss) ){
		equals = false;
	    }
	}
	for( String ss: this.ext_subs ){
	    if(  !aExtSubs.contains( ss )  &&  !aExtProds.contains(ss) ){
		equals = false;
	    }
	}
	// if( this.deltaG == a.deltaG && this.reaction.equals( a.reaction )  ){
	//     equals = true;
	// }
	return equals;
    }
    /**
      In Java, the equals() method that is inherited from Object is:
      public boolean equals(Object other);
      In other words, the parameter must be of type Object.

      The use of the @Override annotation can help a ton with silly mistakes.
      Use it whenever you think you are overriding a super class' or interface's method. 
      That way, if you do it wrong you will get a compile error.
    **/



}
