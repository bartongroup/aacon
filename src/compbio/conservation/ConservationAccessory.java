package compbio.conservation;

import java.util.*;

class ConservationAccessory {

	private static double[] blosumMatrix() {

	double[] matrix = { 6,   -2,   -2,   -3,   -1,   -1,   -1,    0,   -2,   -2,   -2,   -1,   -1,   -3,   -1,    2,    0,   -4,   -3,    0,   -2.5, -2.5, -6,   -6,
					   -2,    8,   -1,   -2,   -5,    1,    0,   -3,    0,   -4,   -3,    3,   -2,   -4,   -3,   -1,   -2,   -4,   -3,   -4,   -1.5, -1.5, -6,   -6,
					   -2,   -1,    8,    2,   -4,    0,    0,   -1,    1,   -5,   -5,    0,   -3,   -4,   -3,    1,    0,   -6,   -3,   -4,    5,    5,   -6,   -6,
					   -3,   -2,    2,    9,   -5,    0,    2,   -2,   -2,   -5,   -5,   -1,   -5,   -5,   -2,    0,   -2,   -6,   -5,   -5,    5.5,  5.5, -6,   -6,
					   -1,   -5,   -4,   -5,   13,   -4,   -5,   -4,   -4,   -2,   -2,   -5,   -2,   -4,   -4,   -1,   -1,   -3,   -4,   -1,   -4.5, -4.5, -6,   -6,
					   -1,    1,    0,    0,   -4,    8,    3,   -3,    1,   -4,   -3,    2,   -1,   -5,   -2,    0,   -1,   -3,   -2,   -3,    0,    0,   -6,   -6,
					   -1,    0,    0,    2,   -5,    3,    7,   -3,    0,   -5,   -4,    1,   -3,   -5,   -2,    0,   -1,   -4,   -3,   -4,    1,    1,   -6,   -6,
						0,   -3,   -1,   -2,   -4,   -3,   -3,    8,   -3,   -6,   -5,   -2,   -4,   -5,   -3,    0,   -2,   -4,   -5,   -5,   -1.5, -1.5, -6,   -6,
					   -2,    0,    1,   -2,   -4,    1,    0,   -3,   11,   -5,   -4,   -1,   -2,   -2,   -3,   -1,   -3,   -4,    3,   -5,   -0.5, -0.5, -6,   -6,
					   -2,   -4,   -5,   -5,   -2,   -4,   -5,   -6,   -5,    6,    2,   -4,    2,    0,   -4,   -4,   -1,   -4,   -2,    4,   -5,   -5,   -6,   -6,
					   -2,   -3,   -5,   -5,   -2,   -3,   -4,   -5,   -4,    2,    6,   -4,    3,    1,   -4,   -4,   -2,   -2,   -2,    1,   -5,   -5,   -6,   -6,
					   -1,    3,    0,   -1,   -5,    2,    1,   -2,   -1,   -4,   -4,    7,   -2,   -5,   -2,    0,   -1,   -4,   -3,   -3,   -0.5, -0.5, -6,   -6,
					   -1,   -2,   -3,   -5,   -2,   -1,   -3,   -4,   -2,    2,    3,   -2,    8,    0,   -4,   -2,   -1,   -2,   -1,    1,   -4,   -4,   -6,   -6,
					   -3,   -4,   -4,   -5,   -4,   -5,   -5,   -5,   -2,    0,    1,   -5,    0,    9,   -5,   -4,   -3,    1,    4,   -1,   -4.5, -4.5, -6,   -6,
					   -1,   -3,   -3,   -2,   -4,   -2,   -2,   -3,   -3,   -4,   -4,   -2,   -4,   -5,   11,   -1,   -2,   -5,   -4,   -4,   -2.5, -2.5, -6,   -6,
						2,   -1,    1,    0,   -1,    0,    0,    0,   -1,   -4,   -4,    0,   -2,   -4,   -1,    6,    2,   -4,   -3,   -2,    0.5,  0.5, -6,   -6,
						0,   -2,    0,   -2,   -1,   -1,   -1,   -2,   -3,   -1,   -2,   -1,   -1,   -3,   -2,    2,    7,   -4,   -2,    0,   -1,   -1,   -6,   -6,
					   -4,   -4,   -6,   -6,   -3,   -3,   -4,   -4,   -4,   -4,   -2,   -4,   -2,    1,   -5,   -4,   -4,   16,    3,   -4,   -6,   -6,   -6,   -6,
					   -3,   -3,   -3,   -5,   -4,   -2,   -3,   -5,    3,   -2,   -2,   -3,   -1,    4,   -4,   -3,   -2,    3,   10,   -2,   -4,   -4,   -6,   -6,
						0,   -4,   -4,   -5,   -1,   -3,   -4,   -5,   -5,    4,    1,   -3,    1,   -1,   -4,   -2,    0,   -4,   -2,    6,   -4.5, -4.5, -6,   -6,
					   -2.5, -1.5,  5,    5.5, -4.5,  0,    1,   -1.5, -0.5, -5,   -5,   -0.5, -4,   -4.5, -2.5,  0.5, -1,   -6,   -4,   -4.5,  8.5,  0.5, -6,   -6,
					   -2.5, -1.5,  5,    5.5, -4.5,  0,    1,   -1.5, -0.5, -5,   -5,   -0.5, -4,   -4.5, -2.5,  0.5, -1,   -6,   -4,   -4.5,  0.5,  7.5, -6,   -6,
					   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,
					   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6,   -6  };

	return matrix;

	}
	
	// gonnet matrix yet to be created
			
	private static double[] gonnetMatrix() {
		
		double[] matrix = {  2.4,  -0.6,  -0.3,  -0.3,   0.5,  -0.2,   0.0,   0.5,  -0.8,  -0.8,  -1.2,  -0.4,  -0.7,  -2.3,   0.3,   1.1,   0.6,  -3.6   -2.2,   0.1,  -0.3,  -0.3,  -5.2,  -5.2,
	             			-0.6,   4.7,   0.3,  -0.3,  -2.2,   1.5,   0.4,  -1.0,   0.6,  -2.4,  -2.2,   2.7,  -1.7,  -3.2,  -0.9,  -0.2,  -0.2,  -1.6,  -1.8,  -2.0,   0  ,   0  ,  -5.2,  -5.2,
	             			-0.3,   0.3,   3.8,   2.2,  -1.8,   0.7,   0.9,   0.4,   1.2,  -2.8,  -3.0,   0.8,  -2.2,  -3.1,  -0.9,   0.9,   0.5,  -3.6,  -1.4,  -2.2,   3  ,   3  ,  -5.2,  -5.2,
	             			-0.3,  -0.3,   2.2,   4.7,  -3.2,   0.9,   2.7,   0.1,   0.4,  -3.8,  -4.0,   0.5,  -3.0,  -4.5,  -0.7,   0.5,   0.0,  -5.2,  -2.8,  -2.9,   3.45,  3.45, -5.2,  -5.2,
	             			 0.5,  -2.2,  -1.8,  -3.2,  11.5,  -2.4,  -3.0,  -2.0,  -1.3,  -1.1,  -1.5,  -2.8,  -0.9,  -0.8,  -3.1,   0.1,  -0.5,  -1.0,  -0.5,   0.0,  -2.5,  -2.5,  -5.2,  -5.2,
	             			-0.2,   1.5,   0.7,   0.9,  -2.4,   2.7,   1.7,  -1.0,   1.2,  -1.9,  -1.6,   1.5,  -1.0,  -2.6,  -0.2,   0.2,   0.0,  -2.7,  -1.7,  -1.5,   0.8,   0.8,  -5.2,  -5.2,
	             			 0.0,   0.4,   0.9,   2.7,  -3.0,   1.7,   3.6,  -0.8,   0.4,  -2.7,  -2.8,   1.2,  -2.0,  -3.9,  -0.5,   0.2,  -0.1,  -4.3,  -2.7,  -1.9,   1.8,   1.8,  -5.2,  -5.2,
	             			 0.5,  -1.0,   0.4,   0.1,  -2.0,  -1.0,  -0.8,   6.6,  -1.4,  -4.5,  -4.4,  -1.1,  -3.5,  -5.2,  -1.6,   0.4,  -1.1,  -4.0,  -4.0,  -3.3,   0.25,  0.25, -5.2,  -5.2,
	             			-0.8,   0.6,   1.2,   0.4,  -1.3,   1.2,   0.4,  -1.4,   6.0,  -2.2,  -1.9,   0.6,  -1.3,  -0.1,  -1.1,  -0.2,  -0.3,  -0.8,   2.2,  -2.0,   0.8,   0.8,  -5.2,  -5.2,
	             			-0.8,  -2.4,  -2.8,  -3.8,  -1.1,  -1.9,  -2.7,  -4.5,  -2.2,   4.0,   2.8,  -2.1,   2.5,   1.0,  -2.6,  -1.8,  -0.6,  -1.8,  -0.7,   3.1,  -3.3,  -3.3,  -5.2,  -5.2,
	             			-1.2,  -2.2,  -3.0,  -4.0,  -1.5,  -1.6,  -2.8,  -4.4,  -1.9,   2.8,   4.0,  -2.1,   2.8,   2.0,  -2.3,  -2.1,  -1.3,  -0.7,   0.0,   1.8,  -3.5,  -3.5,  -5.2,  -5.2,
	             			-0.4,   2.7,   0.8,   0.5,  -2.8,   1.5,   1.2,  -1.1,   0.6,  -2.1,  -2.1,   3.2,  -1.4,  -3.3,  -0.6,   0.1,   0.1,  -3.5,  -2.1,  -1.7,   0.65,  0.65, -5.2,  -5.2,
	             			-0.7,  -1.7,  -2.2,  -3.0,  -0.9,  -1.0,  -2.0,  -3.5,  -1.3,   2.5,   2.8,  -1.4,   4.3,   1.6,  -2.4,  -1.4,  -0.6,  -1.0,  -0.2,   1.6,  -2.6,  -2.6,  -5.2,  -5.2,
	             			-2.3,  -3.2,  -3.1,  -4.5,  -0.8,  -2.6,  -3.9,  -5.2,  -0.1,   1.0,   2.0,  -3.3,   1.6,   7.0,  -3.8,  -2.8,  -2.2,   3.6,   5.1,   0.1,  -3.8,  -3.8,  -5.2,  -5.2,
	             			 0.3,  -0.9,  -0.9,  -0.7,  -3.1,  -0.2,  -0.5,  -1.6,  -1.1,  -2.6,  -2.3,  -0.6,  -2.4,  -3.8,   7.6,   0.4,   0.1,  -5.0,  -3.1,  -1.8,  -0.8,  -0.8,  -5.2,  -5.2,
	             			 1.1,  -0.2,   0.9,   0.5,   0.1,   0.2,   0.2,   0.4,  -0.2,  -1.8,  -2.1,   0.1,  -1.4,  -2.8,   0.4,   2.2,   1.5,  -3.3,  -1.9,  -1.0,   0.7,   0.7,  -5.2,  -5.2,
	             			 0.6,  -0.2,   0.5,   0.0,  -0.5,   0.0,  -0.1,  -1.1,  -0.3,  -0.6,  -1.3,   0.1,  -0.6,  -2.2,   0.1,   1.5,   2.5,  -3.5,  -1.9,   0.0,   0.25,  0.25, -5.2,  -5.2,
	             			-3.6,  -1.6,  -3.6,  -5.2,  -1.0,  -2.7,  -4.3,  -4.0,  -0.8,  -1.8,  -0.7,  -3.5,  -1.0,   3.6,  -5.0,  -3.3,  -3.5,  14.2,   4.1,  -2.6,  -4.4,  -4.4,  -5.2,  -5.2,
	             			-2.2,  -1.8,  -1.4,  -2.8,  -0.5,  -1.7,  -2.7,  -4.0,   2.2,  -0.7,   0.0,  -2.1,  -0.2,   5.1,  -3.1,  -1.9,  -1.9,   4.1,   7.8,  -1.1,  -2.1,  -2.1,  -5.2,  -5.2,
	             			 0.1,  -2.0,  -2.2,  -2.9,   0.0,  -1.5,  -1.9,  -3.3,  -2.0,   3.1,   1.8,  -1.7,   1.6,   0.1,  -1.8,  -1.0,   0.0,  -2.6,  -1.1,   3.4,  -2.55, -2.55, -5.2,  -5.2,
	             			-0.3,   0  ,   3  ,   3.45, -2.5,   0.8,   1.8,   0.25,  0.8,  -3.3,  -3.5,   0.65, -2.6,  -3.8,  -0.8,   0.7,   0.25, -4.4,  -2.1,  -2.55,  4.25,  1.3,  -5.2,  -5.2,
	             			-0.3,   0  ,   3  ,   3.45, -2.5,   0.8,   1.8,   0.25,  0.8,  -3.3,  -3.5,   0.65, -2.6,  -3.8,  -0.8,   0.7,   0.25, -4.4,  -2.1,  -2.55,  1.3,   3.15, -5.2,  -5.2,
	             			-5.2,  -5.2,  -5.2,  -5.2 , -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,
	             			-5.2,  -5.2,  -5.2,  -5.2 , -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2,  -5.2 } ;
		
		
		return matrix;
		
	}
	// gonnet pair, yet to be created
	static double gonnetPair(char a, char b) {
		
		double [] mat = ConservationAccessory.gonnetMatrix();
		
		int [] pairIndeces = ConservationAccessory.pairIndex(a, b);
		
		int pairIndex = 24 * pairIndeces[0] + pairIndeces[1];
		
		double pairValue = mat[pairIndex];
		
		return pairValue;
		
	}

	
	// pam250 matrix yet to be created
	
	private static double[] pam250Matrix() {
		
		double[] matrix = { 2,   -2,    0,    0,   -2,    0,    0,    1,   -1,   -1,   -2,   -1,   -1,   -4,    1,    1,    1,   -6,   -3,    0,    0,    0 ,  -8,   -8,
                		   -2,    6,    0,   -1,   -4,    1,   -1,   -3,    2,   -2,   -3,    3,    0,   -4,    0,    0,   -1,    2,   -4,   -2,   -0.5, -0.5, -8,   -8,
                			0,    0,    2,    2,   -4,    1,    1,    0,    2,   -2,   -3,    1,   -2,   -4,   -1,    1,    0,   -4,   -2,   -2,    2,    2,   -8,   -8,
                			0,   -1,    2,    4,   -5,    2,    3,    1,    1,   -2,   -4,    0,   -3,   -6,   -1,    0,    0,   -7,   -4,   -2,    3,    3,   -8,   -8,
                		   -2,   -4,   -4,   -5,   12,   -5,   -5,   -3,   -3,   -2,   -6,   -5,   -5,   -4,   -3,    0,   -2,   -8,    0,   -2,   -4.5, -4.5, -8,   -8,
                            0,    1,    1,    2,   -5,    4,    2,   -1,    3,   -2,   -2,    1,   -1,   -5,    0,   -1,   -1,   -5,   -4,   -2,    1.5,  1.5, -8,   -8,
                            0,   -1,    1,    3,   -5,    2,    4,    0,    1,   -2,   -3,    0,   -2,   -5,   -1,    0,    0,   -7,   -4,   -2,    2,    2,   -8,   -8,
                            1,   -3,    0,    1,   -3,   -1,    0,    5,   -2,   -3,   -4,   -2,   -3,   -5,   -1,    1,    0,   -7,   -5,   -1,    0.5,  0.5, -8,   -8,
                           -1,    2,    2,    1,   -3,    3,    1,   -2,    6,   -2,   -2,    0,   -2,   -2,    0,   -1,   -1,   -3,    0,   -2,    1.5,  1.5, -8,   -8,
                           -1,   -2,   -2,   -2,   -2,   -2,   -2,   -3,   -2,    5,    2,   -2,    2,    1,   -2,   -1,    0,   -5,   -1,    4,   -2,   -2,   -8,   -8,
                           -2,   -3,   -3,   -4,   -6,   -2,   -3,   -4,   -2,    2,    6,   -3,    4,    2,   -3,   -3,   -2,   -2,   -1,    2,   -3.5, -3.5, -8,   -8,
                           -1,    3,    1,    0,   -5,    1,    0,   -2,    0,   -2,   -3,    5,    0,   -5,   -1,    0,    0,   -3,   -4,   -2,    0.5,  0.5, -8,   -8,
                           -1,    0,   -2,   -3,   -5,   -1,   -2,   -3,   -2,    2,    4,    0,    6,    0,   -2,   -2,   -1,   -4,   -2,    2,   -2.5, -2.5, -8,   -8,
                           -4,   -4,   -4,   -6,   -4,   -5,   -5,   -5,   -2,    1,    2,   -5,    0,    9,   -5,   -3,   -3,    0,    7,   -1,   -5,   -5,   -8,   -8,
                            1,    0,   -1,   -1,   -3,    0,   -1,   -1,    0,   -2,   -3,   -1,   -2,   -5,    6,    1,    0,   -6,   -5,   -1,   -1,   -1,   -8,   -8,
                            1,    0,    1,    0,    0,   -1,    0,    1,   -1,   -1,   -3,    0,   -2,   -3,    1,    2,    1,   -2,   -3,   -1,    0.5,  0.5, -8,   -8,
                            1,   -1,    0,    0,   -2,   -1,    0,    0,   -1,    0,   -2,    0,   -1,   -3,    0,    1,    3,   -5,   -3,    0,    0,    0,   -8,   -8,
                           -6,    2,   -4,   -7,   -8,   -5,   -7,   -7,   -3,   -5,   -2,   -3,   -4,    0,   -6,   -2,   -5,   17,    0,   -6,   -5.5, -5.5, -8,   -8,
                           -3,   -4,   -2,   -4,    0,   -4,   -4,   -5,    0,   -1,   -1,   -4,   -2,    7,   -5,   -3,   -3,    0,   10,   -2,   -3,   -3,   -8,   -8,
                            0,   -2,   -2,   -2,   -2,   -2,   -2,   -1,   -2,    4,    2,   -2,    2,   -1,   -1,   -1,    0,   -6,   -2,    4,   -2,   -2,   -8,   -8,
                            0,   -0.5,  2,    3,   -4.5,  1.5,  2,    0.5,  1.5, -2,   -3.5,  0.5, -2.5, -5,   -1,    0.5,  0,   -5.5, -3,   -2,    3,    1.75,-8,   -8,
                            0,   -0.5,  2,    3,   -4.5,  1.5,  2,    0.5,  1.5, -2,   -3.5,  0.5, -2.5, -5,   -1,    0.5,  0,   -5.5, -3,   -2,    1.75,  4,  -8,   -8,
                           -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,
                           -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8,   -8 };
		
		return matrix;
		
	}
	// pam250 pair, yet to be created
	static double pam250Pair(char a, char b) {
		
		double [] mat = ConservationAccessory.pam250Matrix();
		
		int [] pairIndeces = ConservationAccessory.pairIndex(a, b);
		
		int pairIndex = 24 * pairIndeces[0] + pairIndeces[1];
		
		double pairValue = mat[pairIndex];
		
		return pairValue;
		
		
	}
	
// pet91 matrix yet to be created
	
	private static double[] pet91Matrix() {
		
		double[] matrix = {     1.00,   0.36,   0.45,   0.45,   0.36,   0.36,   0.36,   0.55,   0.27,   0.45,   0.36,   0.36,   0.36,   0.18,   0.55,   0.55,   0.64,   0.09,   0.18,   0.55,   0.45 ,  0.45 ,  0.00, 0.00,
                0.36,   1.00,   0.45,   0.36,   0.36,   0.64,   0.45,   0.45,   0.64,   0.18,   0.18,   0.82,   0.27,   0.09,   0.36,   0.36,   0.36,   0.45,   0.27,   0.18,   0.405,  0.405,  0.00, 0.00,
                0.45,   0.45,   1.00,   0.64,   0.36,   0.45,   0.55,   0.45,   0.55,   0.27,   0.18,   0.55,   0.27,   0.18,   0.36,   0.55,   0.55,   0.00,   0.36,   0.27,   0.82,   0.82 ,  0.00, 0.00,
                0.45,   0.36,   0.64,   1.00,   0.18,   0.55,   0.82,   0.55,   0.45,   0.18,   0.09,   0.45,   0.18,   0.00,   0.27,   0.45,   0.36,   0.00,   0.27,   0.27,   0.82,   0.82 ,  0.00, 0.00,
                0.36,   0.36,   0.36,   0.18,   1.00,   0.18,   0.09,   0.36,   0.45,   0.27,   0.18,   0.18,   0.27,   0.45,   0.27,   0.55,   0.36,   0.55,   0.64,   0.27,   0.27,   0.27 ,  0.00, 0.00,
                0.36,   0.64,   0.45,   0.55,   0.18,   1.00,   0.64,   0.36,   0.64,   0.18,   0.27,   0.64,   0.27,   0.09,   0.45,   0.36,   0.36,   0.18,   0.27,   0.18,   0.5 ,   0.5  ,  0.00, 0.00,
                0.36,   0.45,   0.55,   0.82,   0.09,   0.64,   1.00,   0.45,   0.45,   0.18,   0.09,   0.55,   0.18,   0.00,   0.27,   0.36,   0.36,   0.00,   0.09,   0.27,   0.685,  0.685,  0.00, 0.00,
                0.55,   0.45,   0.45,   0.55,   0.36,   0.36,   0.45,   1.00,   0.27,   0.18,   0.09,   0.36,   0.18,   0.00,   0.36,   0.55,   0.36,   0.27,   0.09,   0.27,   0.5 ,   0.5  ,  0.00, 0.00,
                0.27,   0.64,   0.55,   0.45,   0.45,   0.64,   0.45,   0.27,   1.00,   0.18,   0.27,   0.55,   0.27,   0.45,   0.45,   0.36,   0.36,   0.18,   0.82,   0.18,   0.5 ,   0.5  ,  0.00, 0.00,
                0.45,   0.18,   0.27,   0.18,   0.27,   0.18,   0.18,   0.18,   0.18,   1.00,   0.64,   0.18,   0.73,   0.45,   0.27,   0.36,   0.55,   0.09,   0.27,   0.82,   0.225,  0.225,  0.00, 0.00,
                0.36,   0.18,   0.18,   0.09,   0.18,   0.27,   0.09,   0.09,   0.27,   0.64,   1.00,   0.18,   0.73,   0.64,   0.45,   0.27,   0.36,   0.27,   0.36,   0.64,   0.135,  0.135,  0.00, 0.00,
                0.36,   0.82,   0.55,   0.45,   0.18,   0.64,   0.55,   0.36,   0.55,   0.18,   0.18,   1.00,   0.27,   0.00,   0.27,   0.36,   0.36,   0.18,   0.18,   0.18,   0.5 ,   0.5  ,  0.00, 0.00,
                0.36,   0.27,   0.27,   0.18,   0.27,   0.27,   0.18,   0.18,   0.27,   0.73,   0.73,   0.27,   1.00,   0.45,   0.27,   0.36,   0.45,   0.18,   0.27,   0.64,   0.225,  0.225,  0.00, 0.00,
                0.18,   0.09,   0.18,   0.00,   0.45,   0.09,   0.00,   0.00,   0.45,   0.45,   0.64,   0.00,   0.45,   1.00,   0.18,   0.27,   0.27,   0.36,   0.91,   0.45,   0.09,   0.09 ,  0.00, 0.00,
                0.55,   0.36,   0.36,   0.27,   0.27,   0.45,   0.27,   0.36,   0.45,   0.27,   0.45,   0.27,   0.27,   0.18,   1.00,   0.55,   0.55,   0.09,   0.18,   0.36,   0.315,  0.315,  0.00, 0.00,
                0.55,   0.36,   0.55,   0.45,   0.55,   0.36,   0.36,   0.55,   0.36,   0.36,   0.27,   0.36,   0.36,   0.27,   0.55,   1.00,   0.55,   0.18,   0.36,   0.36,   0.5 ,   0.5  ,  0.00, 0.00,
                0.64,   0.36,   0.55,   0.36,   0.36,   0.36,   0.36,   0.36,   0.36,   0.55,   0.36,   0.36,   0.45,   0.27,   0.55,   0.55,   1.00,   0.09,   0.18,   0.45,   0.455,  0.455,  0.00, 0.00,
                0.09,   0.45,   0.00,   0.00,   0.55,   0.18,   0.00,   0.27,   0.18,   0.09,   0.27,   0.18,   0.18,   0.36,   0.09,   0.18,   0.09,   1.00,   0.45,   0.18,   0   ,   0    ,  0.00, 0.00,
                0.18,   0.27,   0.36,   0.27,   0.64,   0.27,   0.09,   0.09,   0.82,   0.27,   0.36,   0.18,   0.27,   0.91,   0.18,   0.36,   0.18,   0.45,   1.00,   0.18,   0.315,  0.315,  0.00, 0.00,
                0.55,   0.18,   0.27,   0.27,   0.27,   0.18,   0.27,   0.27,   0.18,   0.82,   0.64,   0.18,   0.64,   0.45,   0.36,   0.36,   0.45,   0.18,   0.18,   1.00,   0.27,   0.27 ,  0.00, 0.00,
                0.45,   0.405,  0.82,   0.82,   0.27,   0.5 ,   0.685,  0.5 ,   0.5 ,   0.225,  0.135,  0.5 ,   0.225,  0.09,   0.315,  0.5 ,   0.455,  0   ,   0.315,  0.27,   1   ,  0.5925,  0.00, 0.00,
                0.45,   0.405,  0.82,   0.82,   0.27,   0.5 ,   0.685,  0.5 ,   0.5 ,   0.225,  0.135,  0.5 ,   0.225,  0.09,   0.315,  0.5 ,   0.455,  0   ,   0.315,  0.27,   0.5925, 1    ,  0.00, 0.00,
                0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0   ,   0    ,  0.00, 0.00,
                0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0.00,   0   ,   0    ,  0.00, 0.00  } ;

		return matrix;
		
	}
	
	static double pet91Pair(char a, char b) {
		
		double [] mat = ConservationAccessory.pet91Matrix();
		
		int [] pairIndeces = ConservationAccessory.pairIndex(a, b);
		
		int pairIndex = 24 * pairIndeces[0] + pairIndeces[1];
		
		double pairValue = mat[pairIndex];
		
		return pairValue;
		
	}
	
	private static double[] miyataArmonMatrix() {
		
		double[] matrix = { 0, 2.92, 1.78, 2.37, 1.39, 1.92, 2.46, 0.91, 2.17, 2.69, 2.76, 2.96, 2.42, 3.23, 0.06, 0.51, 0.90, 4.23, 3.18, 1.85, 2.075, 2.075, 6, 6,
           2.92, 0, 2.04, 2.34, 3.06, 1.13, 1.45, 3.58, 0.82, 2.49, 2.62, 0.40, 2.29, 2.47, 2.90, 2.74, 2.03, 2.72, 2.02, 2.43, 2.19, 2.19, 6, 6, 
		   1.78, 2.04, 0, 0.65, 2.83, 0.99, 0.85, 1.96, 1.29, 3.37, 3.49, 1.84, 3.08, 3.70, 1.80, 1.31, 1.40, 4.39, 3.42, 2.76, 0.325, 0.325, 6, 6,
		   2.37, 2.34, 0.65, 0, 3.48, 1.47, 0.90, 2.37, 1.72, 3.98, 4.10, 2.05, 3.69, 4.27, 2.40, 1.87, 2.05, 4.88, 3.95, 3.40, 0.325, 0.325, 6, 6,
		   1.39, 3.06, 2.83, 3.48, 0, 2.48, 3.26, 2.22, 2.56, 1.63, 1.65, 3.27, 1.46, 2.24, 1.33, 1.84, 1.45, 3.34, 2.38, 0.86, 3.155, 3.155, 6, 6,
		   1.92, 1.13, 0.99, 1.47, 2.48, 0, 0.84, 2.48, 0.32, 2.57, 2.70, 1.06, 2.30, 2.81, 1.92, 1.65, 1.12, 3.42, 2.48, 2.13, 1.23, 1.23, 6, 6,
		   2.46, 1.45, 0.85, 0.90, 3.26, 0.84, 0, 2.78, 0.96, 3.39, 3.53, 1.14, 3.13, 3.59, 2.48, 2.06, 1.83, 4.08, 3.22, 2.97, 0.875, 0.875, 6, 6,
		   0.91, 3.58, 1.96, 2.37, 2.22, 2.48, 2.78, 0, 2.78, 3.60, 3.67, 3.54, 3.34, 4.14, 0.97, 0.85, 1.70, 5.13, 4.08, 2.76, 2.165, 2.165, 6, 6,
		   2.17, 0.82, 1.29, 1.72, 2.56, 0.32, 0.96, 2.78, 0, 2.45, 2.59, 0.79, 2.19, 2.63, 2.15, 1.94, 1.32, 3.16, 2.27, 2.11, 1.505, 1.505, 6, 6,
           2.69, 2.49, 3.37, 3.98, 1.63, 2.57, 3.39, 3.60, 2.45, 0, 0.14, 2.84, 0.29, 0.61, 2.62, 2.95, 2.14, 1.72, 0.86, 0.85, 3.675, 3.675, 6, 6,
		   2.76, 2.62, 3.49, 4.10, 1.65, 2.70, 3.53, 3.67, 2.59, 0.14, 0, 2.98, 0.41, 0.63, 2.70, 3.04, 2.25, 1.73, 0.94, 0.91, 3.795, 3.795, 6, 6,
		   2.96, 0.40, 1.84, 2.05, 3.27, 1.06, 1.14, 3.54, 0.79, 2.84, 2.98, 0, 2.63, 2.85, 2.94, 2.71, 2.10, 3.11, 2.42, 2.70, 1.945, 1.945, 6, 6,
		   2.42, 2.29, 3.08, 3.69, 1.46, 2.30, 3.13, 3.34, 2.19, 0.29, 0.41, 2.63, 0, 0.82, 2.36, 2.67, 1.86, 1.89, 0.93, 0.62, 3.385, 3.385, 6, 6,
		   3.23, 2.47, 3.70, 4.27, 2.24, 2.81, 3.59, 4.14, 2.63, 0.61, 0.63, 2.85, 0.82, 0, 3.17, 3.45, 2.60, 1.11, 0.48, 1.43, 3.985, 3.985, 6, 6,
		   0.06, 2.90, 1.80, 2.40, 1.33, 1.92, 2.48, 0.97, 2.15, 2.62, 2.70, 2.94, 2.36, 3.17, 0, 0.56, 0.87, 4.17, 3.12, 1.79, 2.1, 2.1, 6, 6,
		   0.51, 2.74, 1.31, 1.87, 1.84, 1.65, 2.06, 0.85, 1.94, 2.95, 3.04, 2.71, 2.67, 3.45, 0.56, 0, 0.89, 4.38, 3.33, 2.15, 1.59, 1.59, 6, 6,
		   0.90, 2.03, 1.40, 2.05, 1.45, 1.12, 1.83, 1.70, 1.32, 2.14, 2.25, 2.10, 1.86, 2.60, 0.87, 0.89, 0, 3.50, 2.45, 1.42, 1.725, 1.725, 6, 6,
		   4.23, 2.72, 4.39, 4.88, 3.34, 3.42, 4.08, 5.13, 3.16, 1.72, 1.73, 3.11, 1.89, 1.11, 4.17, 4.38, 3.50, 0, 1.06, 2.51, 4.635, 4.635, 6, 6,
 		   3.18, 2.02, 3.42, 3.95, 2.38, 2.48, 3.22, 4.08, 2.27, 0.86, 0.94, 2.42, 0.93, 0.48, 3.12, 3.33, 2.45, 1.06, 0, 1.52, 3.685, 3.685, 6, 6,
		   1.85, 2.43, 2.76, 3.40, 0.86, 2.13, 2.97, 2.76, 2.11, 0.85, 0.91, 2.70, 0.62, 1.43, 1.79, 2.15, 1.42, 2.51, 1.52, 0, 3.08, 3.08, 6, 6,
		   2.075, 2.19, 0.325, 0.325, 3.155, 1.23, 0.875, 2.165, 1.505, 3.675, 3.795, 1.945, 3.385, 3.985, 2.1, 1.59, 1.725, 4.635, 3.685, 3.08, 0, 1.0525, 6, 6,
		   2.075, 2.19, 0.325, 0.325, 3.155, 1.23, 0.875, 2.165, 1.505, 3.675, 3.795, 1.945, 3.385, 3.985, 2.1, 1.59, 1.725, 4.635, 3.685, 3.08, 1.0525,0, 6, 6,
		   6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 0.5, 0.5,
           6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 0.5, 0.5 } ;

		return matrix;
	}
	// pam250 pair, yet to be created
	static double miyataArmonPair(char a, char b) {
		
		double [] mat = ConservationAccessory.miyataArmonMatrix();
		
		int [] pairIndeces = ConservationAccessory.pairIndex(a, b);
		
		int pairIndex = 24 * pairIndeces[0] + pairIndeces[1];
		
		double pairValue = mat[pairIndex];
		
		return pairValue;
		
		
	}
	
	static double BlosumPair(char a, char b) {

	double[] mat = ConservationAccessory.blosumMatrix();
	
	int [] pairIndeces = ConservationAccessory.pairIndex(a, b);
	
	int pairIndex = 24 * pairIndeces[0] + pairIndeces[1];
	
	double pairValue = mat[pairIndex];
	
	return pairValue;
	
	}
	
	
	// this method takes charcters and returns their row and column indeces
	// row number is first elemen tof an array an dcolumn number is second
	
	private static int[] pairIndex(char a, char b) {
		
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

	int[] pairIndeces = new int[2];
	
	pairIndeces[0] = hor;
	
	pairIndeces[1] = ver;

	return pairIndeces;

	}

	// Methods created to do operations on points in N dimensional space

	static double[] addPoints(double[] a, double[] b) {

	if(a == null || b == null) {

	throw new IllegalArgumentException("Reference must not be null");

	}

	assert a.length == b.length;

	double[] aPlusB = new double[a.length];

	for ( int i = 0; i < a.length; i++) {

		aPlusB[i] = a[i] + b[i];

	}

	return aPlusB;

	}

	static double[] multPointByScalar(double[] point, double scalar) {

	double[] pointByScalar = new double[point.length];

	for ( int i = 0; i < point.length; i++) {

		pointByScalar[i] = point[i] * scalar;

	}

	return pointByScalar;

	}  

	static double pointDistance(double[] a, double[] b) {


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
	
	// percentage identity method yet to be written
	
	static double percentIdentity(char[] a, char[] b) {
		
		double ident = 0;
		
		return ident;
		
	}
	
	
	static double[] voronoiWeights(AminoAcidMatrix m, int iter) {
		
		AminoAcidMatrix matrix = m;
		
		int iterations = iter;
		
		Random rgen = new Random();
		
		double[] weights = new double[m.numberOfRows()];
		
		char[] randSeq = new char[m.numberOfColumns()];
		
		// satts iterations, doan't really know what does it to, but jon set up the oterations to 1000
		
		for ( int i = 0; i < iterations; i++) {
		
			// generates a random sequence, equal in length to the sequences in the alignment
			
			for (int j = 0; j < matrix.numberOfColumns(); j++) {
				
				int random = rgen.nextInt(matrix.numberOfRows());
				
				randSeq[i] = matrix.getMatrixPosition(random, j);
				
			}
			
			// measure the distance between each sequence and a random sequence generated
			// distance measured as percentage identity
			
			double[] distances = new double[matrix.numberOfRows()];
			
			double closestValue = 0;
			
			for (int a = 0; a < matrix.numberOfRows(); a++) {
				
				distances[a] = 1.0 - ConservationAccessory.percentIdentity(matrix.getRow(a), randSeq);
				
				if(distances[a] < closestValue) {
					
					closestValue = distances[a];
					
				}
				
			}
			
			// collect all the sequences with the closest distance 
			
			List<Integer> closestSeqs = new ArrayList<Integer>();
			
			for (int b = 0; b < distances.length; b++) {
				
				double dis = distances[b];
				
				if ( dis == closestValue) {
					
					closestSeqs.add(b);
				}
			}
			
			// increase by one the weight of the closest sequence
			
			double increase = 1.0 / closestSeqs.size();
			
			for (int c = 0; c < closestSeqs.size(); c++ ) {
				
				int cs = closestSeqs.get(c);
				
				weights[cs] += increase;
				
			}
			
			// repeat iterations times
			
		}
		
		// normalize weights so they sum up to N
		
		double weightSum = 0.0;
		
		for (int d = 0; d < weights.length; d++) {
			
			weightSum += weights[d];
			
		}
		
		double scaleFactor = weightSum / matrix.numberOfRows();
		
		for (int e = 0; e < weights.length; e++) {
			
			weights[e] = weights[e] + scaleFactor;
			
		}
		
		return weights;
		
	}
	
	static double weightOfSequenceVingronArgos (int seqNr, AminoAcidMatrix m) {
		
		double weight = 0.0;
		
		for( int i = 0; i < m.numberOfRows(); i++) {
			
			if (i != seqNr) {
				
				weight += ConservationAccessory.percentIdentity(m.getColumn(seqNr), m.getColumn(i));
				
			}
		
		}
		
		double result  = (1.0 / m.numberOfRows()) * weight ;
		
		return result;
		
	}
	
	static double dissimilarity(char a, char b) {
		
		double dis = (ConservationAccessory.gonnetPair(a,a) - ConservationAccessory.gonnetPair(a,b)) / ConservationAccessory.gonnetPair(a,a);
		
		return dis;
		
		}
	
	static Map<String, HashSet<Character>> taylorSets() {
		
		Map<String, HashSet<Character>> taySets = new HashMap<String, HashSet<Character>>();

		HashSet<Character> positive = new HashSet<Character>();

		positive.add('R');
		positive.add('K');
		positive.add('H');

		HashSet<Character> charged = new HashSet<Character>();

		charged.add('D');
		charged.add('E');
		charged.add('R');
		charged.add('K');
		charged.add('H');

		HashSet<Character> chargedNonH = new HashSet<Character>();

		chargedNonH.add('D');
		chargedNonH.add('E');
		chargedNonH.add('R');
		chargedNonH.add('K');

		HashSet<Character> negative = new HashSet<Character>();

		negative.add('D');
		negative.add('E');

		HashSet<Character> hydrophilicNonPositive = new HashSet<Character>();

		hydrophilicNonPositive.add('S');
		hydrophilicNonPositive.add('N');
		hydrophilicNonPositive.add('D');
		hydrophilicNonPositive.add('E');
		hydrophilicNonPositive.add('Q');

		HashSet<Character> hydrophilic = new HashSet<Character>();


		hydrophilic.add('S');
		hydrophilic.add('N');
		hydrophilic.add('D');
		hydrophilic.add('E');
		hydrophilic.add('Q');
		hydrophilic.add('R');

		HashSet<Character> chargedOrHydrophilic = new HashSet<Character>();

		chargedOrHydrophilic.add('S');
		chargedOrHydrophilic.add('N');
		chargedOrHydrophilic.add('D');
		chargedOrHydrophilic.add('E');
		chargedOrHydrophilic.add('Q');
		chargedOrHydrophilic.add('R');
		chargedOrHydrophilic.add('K');
		chargedOrHydrophilic.add('H');

		HashSet<Character> chargedOrHydrophilicOrP = new HashSet<Character>();


		chargedOrHydrophilicOrP.add('S');
		chargedOrHydrophilicOrP.add('N');
		chargedOrHydrophilicOrP.add('D');
		chargedOrHydrophilicOrP.add('E');
		chargedOrHydrophilicOrP.add('Q');
		chargedOrHydrophilicOrP.add('R');
		chargedOrHydrophilicOrP.add('K');
		chargedOrHydrophilicOrP.add('H');
		chargedOrHydrophilicOrP.add('P');

		HashSet<Character> polarNonAromaticOrChargedOrP = new HashSet<Character>();

		polarNonAromaticOrChargedOrP.add('P');
		polarNonAromaticOrChargedOrP.add('T');
		polarNonAromaticOrChargedOrP.add('S');
		polarNonAromaticOrChargedOrP.add('N');
		polarNonAromaticOrChargedOrP.add('D');
		polarNonAromaticOrChargedOrP.add('E');
		polarNonAromaticOrChargedOrP.add('Q');
		polarNonAromaticOrChargedOrP.add('R');
		polarNonAromaticOrChargedOrP.add('K');
		polarNonAromaticOrChargedOrP.add('H');

		HashSet<Character> polar = new HashSet<Character>();

		polar.add('T');
		polar.add('S');
		polar.add('N');
		polar.add('D');
		polar.add('E');
		polar.add('Q');
		polar.add('R');
		polar.add('K');
		polar.add('H');
		polar.add('W');
		polar.add('Y');

		HashSet<Character> polarOrP = new HashSet<Character>();

		polarOrP.add('T');
		polarOrP.add('S');
		polarOrP.add('N');
		polarOrP.add('D');
		polarOrP.add('E');
		polarOrP.add('Q');
		polarOrP.add('R');
		polarOrP.add('K');
		polarOrP.add('H');
		polarOrP.add('W');
		polarOrP.add('Y');
		polarOrP.add('P');

		HashSet<Character> polarNonAromaticOrCharged = new HashSet<Character>();

		polarNonAromaticOrCharged.add('T');
		polarNonAromaticOrCharged.add('S');
		polarNonAromaticOrCharged.add('N');
		polarNonAromaticOrCharged.add('D');
		polarNonAromaticOrCharged.add('E');
		polarNonAromaticOrCharged.add('Q');
		polarNonAromaticOrCharged.add('R');
		polarNonAromaticOrCharged.add('K');
		polarNonAromaticOrCharged.add('H');

		HashSet<Character> polarNonAromaticNonPositiveOrP = new HashSet<Character>();

		polarNonAromaticNonPositiveOrP.add('T');
		polarNonAromaticNonPositiveOrP.add('S');
		polarNonAromaticNonPositiveOrP.add('N');
		polarNonAromaticNonPositiveOrP.add('D');
		polarNonAromaticNonPositiveOrP.add('E');
		polarNonAromaticNonPositiveOrP.add('Q');
		polarNonAromaticNonPositiveOrP.add('P');

		HashSet<Character> polarNonAromaticNonPositive = new HashSet<Character>();

		polarNonAromaticNonPositive.add('T');
		polarNonAromaticNonPositive.add('S');
		polarNonAromaticNonPositive.add('N');
		polarNonAromaticNonPositive.add('D');
		polarNonAromaticNonPositive.add('E');
		polarNonAromaticNonPositive.add('Q');

		HashSet<Character> smallPolarOrP = new HashSet<Character>();

		smallPolarOrP.add('P');
		smallPolarOrP.add('T');
		smallPolarOrP.add('S');
		smallPolarOrP.add('N');
		smallPolarOrP.add('D');

		HashSet<Character> smallPolar = new HashSet<Character>();

		smallPolar.add('T');
		smallPolar.add('S');
		smallPolar.add('N');
		smallPolar.add('D');

		HashSet<Character> smallHydrophilic = new HashSet<Character>();

		smallHydrophilic.add('S');
		smallHydrophilic.add('N');
		smallHydrophilic.add('D');

		HashSet<Character> tiny = new HashSet<Character>();

		tiny.add('A');
		tiny.add('G');
		tiny.add('S');

		HashSet<Character> tinyOrSmallOrPolar = new HashSet<Character>();

		tinyOrSmallOrPolar.add('A');
		tinyOrSmallOrPolar.add('G');
		tinyOrSmallOrPolar.add('S');
		tinyOrSmallOrPolar.add('T');
		tinyOrSmallOrPolar.add('N');
		tinyOrSmallOrPolar.add('D');

		HashSet<Character> tinyOrSmallOrPolarOrP = new HashSet<Character>();

		tinyOrSmallOrPolarOrP.add('A');
		tinyOrSmallOrPolarOrP.add('G');
		tinyOrSmallOrPolarOrP.add('S');
		tinyOrSmallOrPolarOrP.add('T');
		tinyOrSmallOrPolarOrP.add('N');
		tinyOrSmallOrPolarOrP.add('D');
		tinyOrSmallOrPolarOrP.add('P');

		HashSet<Character> tinyOrNegativeHydrophilicOrT = new HashSet<Character>();

		tinyOrNegativeHydrophilicOrT.add('A');
		tinyOrNegativeHydrophilicOrT.add('G');
		tinyOrNegativeHydrophilicOrT.add('T');
		tinyOrNegativeHydrophilicOrT.add('S');
		tinyOrNegativeHydrophilicOrT.add('N');
		tinyOrNegativeHydrophilicOrT.add('D');
		tinyOrNegativeHydrophilicOrT.add('E');
		tinyOrNegativeHydrophilicOrT.add('Q');

		HashSet<Character> tinyOrNegativeHydrophilicOrTOrP = new HashSet<Character>();

		tinyOrNegativeHydrophilicOrTOrP.add('A');
		tinyOrNegativeHydrophilicOrTOrP.add('G');
		tinyOrNegativeHydrophilicOrTOrP.add('T');
		tinyOrNegativeHydrophilicOrTOrP.add('S');
		tinyOrNegativeHydrophilicOrTOrP.add('N');
		tinyOrNegativeHydrophilicOrTOrP.add('D');
		tinyOrNegativeHydrophilicOrTOrP.add('E');
		tinyOrNegativeHydrophilicOrTOrP.add('Q');
		tinyOrNegativeHydrophilicOrTOrP.add('P');

		HashSet<Character> tinyOrPolarNonAromatic = new HashSet<Character>();

		tinyOrPolarNonAromatic.add('A');
		tinyOrPolarNonAromatic.add('G');
		tinyOrPolarNonAromatic.add('T');
		tinyOrPolarNonAromatic.add('S');
		tinyOrPolarNonAromatic.add('N');
		tinyOrPolarNonAromatic.add('D');
		tinyOrPolarNonAromatic.add('E');
		tinyOrPolarNonAromatic.add('Q');
		tinyOrPolarNonAromatic.add('R');
		tinyOrPolarNonAromatic.add('K');

		HashSet<Character> all = new HashSet<Character>();

		all.add('A');
		all.add('R');
		all.add('N');
		all.add('D');
		all.add('C');
		all.add('Q');
		all.add('E');
		all.add('G');
		all.add('H');
		all.add('I');
		all.add('L');
		all.add('K');
		all.add('M');
		all.add('F');
		all.add('P');
		all.add('S');
		all.add('T');
		all.add('W');
		all.add('Y');
		all.add('V');
		all.add('B');
		all.add('Z');
		all.add('X');
		all.add('-');

		taySets.put("positive", positive);
		taySets.put("charged", charged);
		taySets.put("chargedNonH", chargedNonH);
		taySets.put("negative", negative);
		taySets.put("hydrophilicNonPositive", hydrophilicNonPositive);
		taySets.put("hydrophilic", hydrophilic);
		taySets.put("chargedOrHydrophilic", chargedOrHydrophilic);
		taySets.put("chargedOrHydrophilicOrP", chargedOrHydrophilicOrP);
		taySets.put("polarNonAromaticOrChargedOrP", polarNonAromaticOrChargedOrP);
		taySets.put("polar", polar);
		taySets.put("polarOrP", polarOrP);
		taySets.put("polarNonAromaticOrCharged", polarNonAromaticOrCharged);
		taySets.put("polarNonAromaticNonPositiveOrP", polarNonAromaticNonPositiveOrP);
		taySets.put("polarNoAromaticNonPositive", polarNonAromaticNonPositive);
		taySets.put("smallPolarOrP", smallPolarOrP);
		taySets.put("smallPolar", smallPolar);
		taySets.put("smallHydrophilic", smallHydrophilic);
		taySets.put("tiny", tiny);
		taySets.put("tinyOrSmallOrPolar", tinyOrSmallOrPolar);
		taySets.put("tinyOrSmallOrPolarOrP", tinyOrSmallOrPolarOrP);
		taySets.put("tinyOrNegativeHydrophilicOrT", tinyOrNegativeHydrophilicOrT);
		taySets.put("tinyOrNegativeHydrophilicOrTOrP", tinyOrNegativeHydrophilicOrTOrP);
		taySets.put("tinyOrPolarNonAromatic", tinyOrPolarNonAromatic);
		taySets.put("all", all);

		return taySets;

		}

	}



