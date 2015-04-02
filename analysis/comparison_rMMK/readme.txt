This directory contains the code to perform the pathway comparison based on a flux calculation which assumes a simplified reversible Michaelis-Menten kinetics. 
All enzymes are assumed to be equivalent (same kinetic parameters). These parameter values are hard-coded in the main() method of ComparePaths_rMMK.java

Dependencies:  It uses the findroot() function of the "mpmath" python library to solve for steady state 

To run:    "java ComparePaths_rMMK input_file.paths output_file.dat"



- PythonScript_ConvKin_NoExternals.java contains methods to generate the python scripts, makeScript_findroot(), which will be automatically executed from the command line, and  getSolution() will do this for multiple initial conditions until a steady state is found.
- ConvenienceKinetics_NoExternals.java takes steady state metabolite concentrations and calculates the flux through each reaction using the specified reaction kinetics (see NOTE below).
- PythonOutputReader.java checks output generated from the python script to see if a steady state was found and then reads in the steady state concentrations
- ExternalMetabs.java contains method to calculate the modification factor to a reaction's equilibrium constant arising from the fixed values of the external metabolite concentrations, getExtMetFactor().

NOTE: the "ConvKin_NoExternals" label in file name is no longer correct. The kinetics are actually based on the "Common modular rate law" of Wolfram Leibermeister (though this is very similar to the Convenience Kinetics). The "NoExternals" indicates that we do not attribute kinetic parameters for the external metabolites. Instead, these fixed concentrations enter only through the modification of the equilibrium constant of a reaction - see notes. Since in our comparison all reactions are thus uni-uni molecular, all reaction rates reduce to that of the simple reversible Michaelis-Menten kinetics.








==================================================================================


  -  add source and product labels in more sensible way

  -  tidy code? put all parameters/ranges at top of main method, or read them from a file?
                write a class that reads in the list of pathways
                write class to sample parameter values
  
