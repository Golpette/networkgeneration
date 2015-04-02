import java.io.*;
import java.util.*; 

public class Reactions3{

    String direction;
    String type;
    int substrate;
    int product;
    int ATP;
    int NADH;
    double deltaG;
    String equation;


    public Reactions3(String direction, String type, int substrate, int product, int ATP, int NADH, double deltaG, String equation){
	this.direction = direction;
	this.type = type;
	this.substrate = substrate;
	this.product = product;
	this.ATP = ATP;
	this.NADH = NADH;
	this.deltaG = deltaG;
	this.equation = equation;
    }


    public String toString(){
	return direction + " " + type + " " + substrate + " " + product + " "+ ATP + " " + NADH + " " + deltaG + " " + equation; 
    }


    //Getter and setters
    public String getDirection(){
	return direction;
    }
    public int getSubstrate(){
	return substrate;
    }
    public int getProduct(){
	return product;
    }
    public String getType(){
	return type;
    }
    public int getATP(){
	return ATP;
    }
    public int getNADH(){
	return NADH;
    }
    public double getDeltaG(){
	return deltaG;
    }
    public String getEquation(){
	return equation;
    }


    public void setDeltaG(double G){
	this.deltaG = G;
    }




}
