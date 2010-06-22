package compbio.conservation;

import java.util.*;

class ConservationAccessory {

	
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

	// percentage identity method yet to be written
	
	static double percentIdentity(char[] a, char[] b) {
		
		if(a == null || b == null) {

			throw new IllegalArgumentException("Reference must not be null");

			}
		
		assert a.length == b.length;
		
		int ident = 0;
		
		for (int i = 0; i < a.length; i++) {
			
			if (a[i] == b[i]) {
				
				ident++;
			}
		
		}
		
		double result = (double) ident / a.length;
		
		return result;
		
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
				
				randSeq[j] = matrix.getMatrixPosition(random, j);
				
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
				
				weight += ConservationAccessory.percentIdentity(m.getRow(seqNr), m.getRow(i));
				
			}
		
		}
		
		double result  = (1.0 / m.numberOfRows()) * weight ;
		
		return result;
		
	}
	
	
}



