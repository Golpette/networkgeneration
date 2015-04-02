This program, Enumerate.java, enumerates linear pathways in our network using the depth-first search (DFS) algorithm.

It requires as input *from the terminal* both the name of the file containing the list of compounds and the file containing the list of reactions in our network, as produced by Bartek's network generation code. It also requires the desired length of pathways to enumerate. So for example, to run this program to enumerate pathways of length 5, compile the files (javac *.java) and type:

"java Enumerate compounds_list__4C_v3_2_2_ext_100.dat reactions__4C_v3_2_2_ext_100.dat 5"

This will print out the L=5 pathways that satisfy the specified restrictions to the screen, along with the number of pathways satisfying particular requirements (see comments in code). To output to a file, add this extra pieace to the command:

"java Enumerate compounds_list__4C_v3_2_2_ext_100.dat reactions__4C_v3_2_2_ext_100.dat 5 > output_filename".

The format of the output file is the format accepted as input by the pathway comparison codes found in the /comparison_xx directories.





The code: -----------------------------------------------------------------------------------

Enumerate.java has some hard-coded parameters and pathway restrictions.

In the main() method, we must set the start and end metabolites (each compound is labelled by a corresponding integer in the compounds file). We also set "desiredATP" if we want to restrict to pathways producing or consuming a specific number of ATP molecules.

The takeStep_counter() method implements the DFS algorithm.  The enumerated pathways are printed out from WITHIN this method.  It also keeps track of the number of pathways satisfying certain restrictions through the counter[] array. The comments explain this and these can be removed/modified as desired.

---------------------------------------------------------------------------------------------











