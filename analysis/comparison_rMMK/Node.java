import java.util.*;

public class Node{


    private String label;
    private double concentration;
    private double increment;
    //
    //---- CAREFUL: This could be the source of bugs... --------------
       ArrayList<Double> listOf_K;
       ArrayList<String> listOf_Links;
    //----------------------------------------------------------------
    // (These Strings list the labels of all internal metabolites involved separated by a comma.)
    private String mathematicaEqn;



    //Constructor methods
    Node(String lab){
	this.label = lab;
	this.listOf_K = new ArrayList<Double>();
	this.listOf_Links = new ArrayList<String>();
	this.concentration = 1.0E-5;
	this.increment =  0;
	this.mathematicaEqn="";
    }
    Node(String lab, double conc){
	this.label = lab;
	this.listOf_K = new ArrayList<Double>();
	this.listOf_Links = new ArrayList<String>();
	this.concentration = conc;
	this.increment =  0;
	this.mathematicaEqn="";
    }
    Node( Node n ){
	this.label = n.getLabel();
	this.concentration = n.getConcentration();
	this.increment = n.getIncrement();
	this.mathematicaEqn = n.getEqn();

	this.listOf_K = new ArrayList<Double>();
	for( int i=0; i < n.listOf_K.size(); i++ ){
	    this.listOf_K.add(  n.listOf_K.get(i)  );
	}

	this.listOf_Links = new ArrayList<String>();
	for( int i=0; i < n.listOf_Links.size(); i++ ){
	    this.listOf_Links.add(  n.listOf_Links.get(i)   );
	}
    }


    //GETTERS
    public String getLabel(){
	return this.label;
    }
    public double getConcentration(){
	return this.concentration;
    }
    public double getIncrement(){
	return this.increment;
    }
    public String getEqn(){
	return this.mathematicaEqn;
    }


    //SETTERS
    public void setConcentration( double conc ){
	this.concentration = conc;
    }
    public void setIncrement( double inc ){
	this.increment = inc;
    }


    //adders...
    public void ADDK(double k){
	this.listOf_K.add(k);
    }

    public void ADDLINK(String involvedCompounds ){
	this.listOf_Links.add( involvedCompounds );
    }

    // concatenate
    public void addToEqn( String s ){
	this.mathematicaEqn = this.mathematicaEqn + s;
    }
   

    // To String
    public String toString(){
	return this.label + " conc=" + concentration + "  (K:" + listOf_K.toString() + ")  (Metabolites: " + listOf_Links.toString() + ")";
    }




}