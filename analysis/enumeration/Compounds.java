
public class Compounds{    
    /**
    Compound object to hold all properties of an internal metabolite
     */

    // fields
    int  index;
    double G;
    String name;
    String formula;
    int charge;
  

    // constructor
    Compounds(int index, double G, String name, String formula, int charge){
	this.index = index;
	this.G = G;
	this.name = name;
	this.formula = formula;
	this.charge = charge;
    }



    // print method
    public String toString() {
	return index + "  " + G + "  " + name + "  " + formula + " " + charge;
    }


    // getters
    int getIndex(){
	return index;
    }
    double getG(){
	return G;
    }
    String getName(){
	return name;
    }
    String getFormula(){
	return formula;
    }
    int getCharge(){
	return charge;
    }



}