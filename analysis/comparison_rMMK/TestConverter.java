

public class TestConverter{

    public static void main(String[] args){


	Reactions3 r3 = new Reactions3("backward", "1.1.1.oxidation", 54, 43, 0, -1, -24.1077, "NAD_red + CH3-CO-COOH ---> NAD_ox + CH3-CH(OH)-COOH" );

	Reaction r = Convert_Reactions3_to_Reaction.convert( r3  );

	System.out.println( r3 );
	System.out.println( r );


       


    }

}
