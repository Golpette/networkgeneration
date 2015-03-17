
public class Compound {
	
	private int label;
	private String formula;
	private int charge;
	
	// Constructor
	Compound(int l, String f, int c){
		this.label = l;
		this.formula = f;
		this.charge = c;
	}
	
	
	// Getters
	public int getLabel(){
		return this.label;
	}
	public String getFormula(){
		return this.formula;
	}
	public int getCharge(){
		return this.charge;
	}
	
	
	// toString
	public String toString(){
		return "{"+label+","+formula+","+charge+"}";
	}

}
