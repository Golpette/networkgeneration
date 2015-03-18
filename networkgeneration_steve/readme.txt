Java program to generate list of linear molecules composed of C, H, O and P up to specified length as well as all possible reactions between them as defined by action of relevant EC classes.

Two main programs to run:  ---------------------------------------------

CompoundGenerator.java  -- specifies allowed molecular groups and calls the depth first search algorithm (DFS) algorithm from Methods.java to piece these groups together in all possible combinations. Molecules are represented as Strings. In this program we also specify the upper length of molecules to be generated, as well as any chemical restrictions (such as to produce only non-hydrocarbons with a non-zero electrostatic charge).

ReactionGenerator.java -- reads in the list of compounds output by CompoundGenerator.java (specify input file name in this proogram) and checks for all reactions than can act upon each metabolite. Note that the product metabolite must be in the compound list for reaction to be allowed.

-------------------------------------------------------------------------


The reaction mechanisms that ReactionGenerator.java checks for are all defined in their corresponding ECX.java classes. Each method in these files, ECX._Y_Z() contains the possible molecular operations found in the EC class defined by the first 3 numbers "X.Y.Z", and may contain more than one such operation.  Comments in these files, along with the reaction notes, detail the reactions and point to examples in nature.



Some other notes:

- Run time increases greatly with compound length. Generating the compounds and reactions for 4-carbon molecules for instance takes 0.4 and 1.5 seconds respectively. Allowing 5-carbon compounds increases these values to about 4.5 and 50 seconds and allowing 6-carbon molecules sees a huge increase in run time, to around 150 seconds and 50 minutes for compound and reaction generation respectively.

- Compound.java defines a Compound object

- Rxn.java defines a Rxn (reaction) object which was used for removing duplicate reaction entries
