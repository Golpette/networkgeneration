import java.util.*;


public class Rxn {
	
	/**    Rxn object:  EC class, substrate, product and reaction String	 */
	
	
	String ec;  String sub;  String prod;  String eqn;
	
	
	// Constructor
	public Rxn( String s ){
		
		Scanner scan = new Scanner( s );
		scan = scan.useDelimiter("\\t");
		this.ec = scan.next();
		this.sub = scan.next();
		this.prod = scan.next();
		this.eqn = scan.next();
	
		// to properly close this scanner object, do I need to have a try{
		// and finally{ scan.close(); }  somewhere?
			
	}
	
	
	// Getters
	public String getEC(){
		return this.ec;
	}
	public String getSub(){
		return this.sub;
	}
	public String getProd(){
		return this.prod;
	}
	public String getEqn(){
		return this.eqn;
	}
	
	
	
	@Override
	public boolean equals(Object other){
		//if( other==null ){ return false; }
		//if( other==this ){ return true; }
		//if( !(other instanceof Rxn )){ return false; }
		Rxn r = (Rxn)other;
		if( r.getSub().equals(this.sub)  &&  r.getProd().equals(this.prod)  && 
				r.getEC().equals(this.ec) && r.getEqn().equals(this.eqn) ){
			return true;
		}
		else{  return false;  }
	}
	
	
	
	// Method to check if two Rxns are EQUIVALENT
	public boolean equivalent( Rxn r ){
		boolean equiv = false;
		if( r.getSub().equals(this.prod) && r.getProd().equals(this.sub) && r.getEC().equals(this.ec) ){
			equiv = true;
		}
		return equiv;
	}
	
	
	
	
	public String toString(){
		return this.ec+"\t"+this.sub+"\t"+this.prod+"\t"+this.eqn;
	}
	
	
	
	
}
