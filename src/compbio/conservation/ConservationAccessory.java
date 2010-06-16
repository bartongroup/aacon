package compbio.conservation;

import java.util.*;

class ConservationAccessory {

	private static int[][] blosumMatrix() {

	int[][] matrix = { {4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0,-2,-1,0,-4},
			   {-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1-3,-2,-1,-1,-3,-2,-3,-1,0,-1,-4},
			   {-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3,3,0,-1,-4},
			   {-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3,4,1,-1,-4},
			   {0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4},
			   {-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2,0,3,-1,-4},
			   {-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4},
			   {0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3,-1,-2,-1,-4},
			   {-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3,0,0,-1,-4},
			   {-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3,-3,-3,-1,-4},
			   {-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1,-4,-3,-1,-4},
			   {-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2,0,1,-1,-4},
			   {-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1,-3,-1,-1,-4},
			   {-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1,-3,-3,-1,-4},
			   {-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2,-2,-1,-2,-4},
			   {1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2,0,0,0,-4},
			   {0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0,-1,-1,0,-4},
			   {-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3,-4,-3,-2,-4},
			   {-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1,-3,-2,-1,-4},
			   {0,-3,-3,-3,-1,-2,-2-3,-3,3,1,-2,1,-1-2,-2,0,-3,-1,4,-3,-2-1,-4},
			   {-2,-1,3,4,-3,0,1,-1,0,-3,-4,0,-3,-3,-2,0,-1,-4,-3,-3,4,1,-1,-4},
			   {-1,0,0,1,-3,3,4,-2,0,-3,-3,1,-1,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4},
			   {0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2,0,0,-2,-1,-1,-1,-1,-1,-4},
			   {-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,4,-4,-4,-4,-4} };

		       
	return matrix;

	}
			

	static int BlosumPair(char a, char b) {

	int[][] mat = ConservationAccessory.blosumMatrix();

	int hor = 0;

	int ver = 0;

	switch(a) {

	case 'A' : hor = 0;
		   break;
	case 'R' : hor = 1;
		   break;
	case 'N' : hor = 2;
		   break;
	case 'D' : hor = 3;
		   break;
	case 'C' : hor = 4;
		   break;
	case 'Q' : hor = 5;
		   break;
	case 'E' : hor = 6;
		   break;
	case 'G' : hor = 7;
		   break;
	case 'H' : hor = 8;
		   break;
	case 'I' : hor = 9;
		   break;
	case 'L' : hor = 10;
		   break;
	case 'K' : hor = 11;
		   break;
	case 'M' : hor = 12;
		   break;
	case 'F' : hor = 13;
		   break;
	case 'P' : hor = 14;
		   break;
	case 'S' : hor = 15;
		   break;
	case 'T' : hor = 16;
		   break;
	case 'W' : hor = 17;
		   break;
	case 'Y' : hor = 18;
		   break;
	case 'V' : hor = 19;
		   break;
	case 'B' : hor = 20;
		   break;
	case 'Z' : hor = 21;
		   break;
	case 'X' : hor = 22;
		   break;
	case '-' : hor = 23;
		   break;
	default : System.out.println("No such symbol in the matrix");
	
	}

	switch(b) {

	case 'A' : ver = 0;
		   break;
	case 'R' : ver = 1;
		   break;
	case 'N' : ver = 2;
		   break;
	case 'D' : ver = 3;
		   break;
	case 'C' : ver = 4;
		   break;
	case 'Q' : ver = 5;
		   break;
	case 'E' : ver = 6;
		   break;
	case 'G' : ver = 7;
		   break;
	case 'H' : ver = 8;
		   break;
	case 'I' : ver = 9;
		   break;
	case 'L' : ver = 10;
		   break;
	case 'K' : ver = 11;
		   break;
	case 'M' : ver = 12;
		   break;
	case 'F' : ver = 13;
		   break;
	case 'P' : ver = 14;
		   break;
	case 'S' : ver = 15;
		   break;
	case 'T' : ver = 16;
		   break;
	case 'W' : ver = 17;
		   break;
	case 'Y' : ver = 18;
		   break;
	case 'V' : ver = 19;
		   break;
	case 'B' : ver = 20;
		   break;
	case 'Z' : ver = 21;
		   break;
	case 'X' : ver = 22;
		   break;
	case '-' : ver = 23;
		   break;
	default : System.out.println("No such symbol in the matrix");
	
	}

	int pair = mat[hor][ver];

	return pair;

	}

	// Methods created to do operations on points in N dimensional space

	static int[] addPoints(int[] a, int[] b) {

	if(a == null || b == null) {

	throw new IllegalArgumentException("Reference must not be null");

	}

	assert a.length == b.length;

	int[] aPlusB = new int[a.length];

	for ( int i = 0; i < a.length; i++) {

		aPlusB[i] = a[i] + b[i];

	}

	return aPlusB;

	}

	static double[] multPointByScalar(int[] point, double scalar) {

	double[] pointByScalar = new double[point.length];

	for ( int i = 0; i < point.length; i++) {

		pointByScalar[i] = point[i] * scalar;

	}

	return pointByScalar;

	}  

	static double pointDistance(int[] a, double[] b) {


	if(a == null || b == null) {

	throw new IllegalArgumentException("Reference must not be null");

	}

	assert a.length == b.length;

	double sum = 0;

	for ( int i = 0; i < a.length; i++) {

		sum = sum + ((b[i] - a[i]) * (b[i] - a[i]));

	}

	double distance = Math.sqrt(sum);

	return distance;

	}

	// amino acid alphabet char array, basic 20 aa, no gaps included

	static char[] alphabetArray() {

	char[] alp = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'};

	return alp;

	}
	
	static Map<String,HashSet<Character>> mirnySets() {
		
		HashSet<Character> aliphatic = new HashSet<Character>();
		
		aliphatic.add('A');
		aliphatic.add('V');
		aliphatic.add('L');
		aliphatic.add('I');
		aliphatic.add('M');
		aliphatic.add('C');
		
		HashSet<Character> aromatic = new HashSet<Character>();
		
		aromatic.add('F');
		aromatic.add('W');
		aromatic.add('Y');
		aromatic.add('H');
		
		HashSet<Character> polar = new HashSet<Character>();
		
		polar.add('S');
		polar.add('T');
		polar.add('N');
		polar.add('Q');
		
		HashSet<Character> positive = new HashSet<Character>();
		
		positive.add('K');
		positive.add('R');
		
		HashSet<Character> negative = new HashSet<Character>();
		
		negative.add('D');
		negative.add('E');
		
		HashSet<Character> special = new HashSet<Character>();
		
		special.add('G');
		special.add('P');
		
		Map<String,HashSet<Character>> sets = new HashMap<String,HashSet<Character>>();
		
		sets.put("aliphatic", aliphatic);
		sets.put("aromatic", aromatic);
		sets.put("polar", polar);
		sets.put("positive", positive);
		sets.put("negative", negative);
		sets.put("special", special);
		
		return sets;
		
	}
	
	static Map<String,HashSet<Character>> williamsonSets() {
		
		Map<String,HashSet<Character>> sets = new HashMap<String,HashSet<Character>>();
		
		HashSet<Character> set1 = new HashSet<Character>();
		
		set1.add('V');
		set1.add('L');
		set1.add('I');
		set1.add('M');
		
		HashSet<Character> set2 = new HashSet<Character>();
		
		set2.add('F');
		set2.add('W');
		set2.add('Y');
		
		HashSet<Character> set3 = new HashSet<Character>();
		
		set3.add('S');
		set3.add('T');
		
		HashSet<Character> set4 = new HashSet<Character>();
		
		set4.add('N');
		set4.add('Q');
		
		HashSet<Character> set5 = new HashSet<Character>();
		
		set5.add('H');
		set5.add('K');
		set5.add('R');
		
		HashSet<Character> set6 = new HashSet<Character>();
		
		set6.add('D');
		set6.add('E');
		
		HashSet<Character> set7 = new HashSet<Character>();
		
		set7.add('A');
		set7.add('G');
		
		HashSet<Character> set8 = new HashSet<Character>();
		
		set8.add('P');
		
		HashSet<Character> set9 = new HashSet<Character>();
		
		set9.add('C');
		
		sets.put("set1", set1);
		sets.put("set2", set2);
		sets.put("set3", set3);
		sets.put("set4", set4);
		sets.put("set5", set5);
		sets.put("set6", set6);
		sets.put("set7", set7);
		sets.put("set8", set8);
		sets.put("set9", set9);
		
		return sets;
		
		
		}
	
}

