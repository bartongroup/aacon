package compbio.correlation;

public class PackTriangular {
	
	public static void main (String args[] ) {
		
		int[][] testArray = { {1,2,3,4,5},
							  {7,8,9,10},
							  {13,14,15},
							  {19,20},
							  {25} } ;
		
		int[] pa = new int[15]; 

		
		for (int j = 0; j < testArray[0].length; j++) {
			
			for (int i = 0; i <= j; i++) {
				
				pa[i + (j-1)*j/2] = testArray[i][j];
				
			}
			
		}
		
		System.out.println(testArray[3][1] +" " + pa[13]);
	
	}
	}


