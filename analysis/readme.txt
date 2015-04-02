Directories contain various code for analysis:

/enumeration = code to enumerate linear pathways between two metabolites 

/comparison_PerfectEnzyme = code to compare flux of pathways using the "perfect enzyme" assumption, as well as the Powell optimization method to calculate the enzyme distribution that maximizes each pathway's flux

/comparison_rMMK = our pathway comparison assuming all enzymes are equivalent (same kinetic parameters) and obey reversible Michaelis Menten kinetics (rMMK)




-  PathChecker.java reads in list of pathways and checks which uses only metabolites
   found in KEGG database (which are entered by hand in this program) 
-  RemoveEnolPaths.java reads in list of pathways and prints out all that do NOT
   contain an enol compound
-  combineOutput_allButFirstAndLast.sh = shell script to combine all pathway comparison data files ommitting first and last lines (first line contains parameter labels, last line may be incomplete if we terminated prgram early)
-  


