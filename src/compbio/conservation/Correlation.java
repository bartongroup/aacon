package compbio.conservation;

import java.util.*;

/**
 * Class has static methods used to create similarity matrices, and correlation matrices.
 * @author agolicz
 *
 */

class Correlation {
	
	private final AminoAcidMatrix alignment;
	
	private final int winWidth;
	
	private final double gapTreshold;
	
	/**
	 * Calculates similarity between two sequences, similarity is calculated as a sum of blosum scores for pairs formed by corresponding amino acids in two sequences.
	 *  
	 * @param seq1 sequence 1
	 * @param seq2 sequence 2
	 * @return  similarity score
	 */
	
	public Correlation(AminoAcidMatrix alignment, int winWidth, double gapTreshold) {
		
		if (alignment == null) {
			
			throw new IllegalArgumentException("Alignment must not be null.");
			
		}
		
		if(winWidth < 0 || winWidth%2 != 1 || winWidth > alignment.numberOfColumns()) {
			
			throw new IllegalArgumentException("ColWidth smaller than zero or an even number or largrt tah the number of columns");
		}
		
		this.alignment = alignment;
		
		this.winWidth = winWidth;
		
		this.gapTreshold = gapTreshold;
		
	}
	
	/**
	 * Calculates sequence similarity as a sum of blosum scores.
	 * 
	 * @param seq1 sequence 1
	 * @param seq2 sequence 2
	 * 
	 * @return similarity
	 */
	
	static int sequenceSimilartyBlosum(char [] seq1, char[] seq2) {
		
		assert seq1.length == seq2.length;
		
		int similarity = 0;
		
		for (int i = 0; i < seq1.length; i++) {
			
			// it's the wrong blosum matrix, need the unmodulated one with integers
			
			similarity += ConservationMatrices.BlosumPair(seq1[i], seq2[i]);
			
		}
		
		return similarity;
		
	}
	
	// end point does not count
	
	
	/**
	 * Calculates local similarity between two sequences. Similarity is a sum of blosum scores for pairs of corresponding amino acids in the two sequences.
	 * 
	 * @param seq1 sequence 1 
	 * @param seq2 sequence 2
	 * @param startPoint index of the starting point of similarity calculation
	 * @param endPoint index of the end point of similarity calculation, end point not included in calculation.
	 * @return array of 2 elements first element is the actual result, is the frame for the new window  
	 */
	
	static int[] localSequenceSimilarityBlosumFirstWindow(char[] seq1, char[] seq2, int endPoint) {
		
		assert seq1.length == seq2.length;
		
		int similarity = 0;
		
		for (int i =0; i < endPoint; i++) {
			
			// it's the wrong blosum matrix, need the unmodulated one with integers
			
			similarity += ConservationMatrices.BlosumPair2(seq1[i], seq2[i]);
			
		}
		
		int similarityNoFirst = similarity - ConservationMatrices.BlosumPair2(seq1[0], seq2[0]);
		
		int[] array = { similarity, similarityNoFirst };
		
		return array;
		
	
	}
	
	
	
	/**
	 * Creates a global similarity matrix.
	 * Always use this method along with local similarity matrix, the internal numbering scheme assures that Pearson correlation will be calculated properly
	 * @param matrix reference to the matrix containing sequences
	 * @return similarity matrix stored as a single array
	 */
	
	private int[] globalSimilarity() {
		
		int[] globalSim = new int[alignment.numberOfRows()* (alignment.numberOfRows() - 1) / 2];
		
		int index = 0;
		
		for (int i = 0; i < alignment.numberOfRows(); i++) {
			
			for (int j = i + 1; j < alignment.numberOfRows(); j++) {
				
				globalSim[index] = Correlation.sequenceSimilartyBlosum(alignment.getRow(i), alignment.getRow(j));
				
				index++;
				
			}
		
		}
		
		return globalSim;
		
		}

	/**
	 * Calculates person coefficients.
	 * 
	 * @return array of Pearson correlation scores indexed by window
	 */
	private double[] calcPearsonCoeff3() {
		
		int[] global = globalSimilarity();
		
		int[][] locals = localSimilarity2();
		
		double[] coeffs = new double[locals.length];
		
		for (int i = 0; i < locals.length; i++) {
			
			coeffs[i] = Correlation.pearson2(locals[i], global);
			
		}
		
		return coeffs;
	}
	
			
	/**
	 * Calculates correlation score for two vectors.
	 * If one of the vectors consists of single value only the correlation value is 0.
	 * @param arr1 vector 1
	 * @param arr2 vector 2
	 * @return score
	 */

	static double pearson2(int[] arr1, int[] arr2){
		
		assert arr1.length == arr2.length;
	
		int arr1Sum = 0;
		
		int arr2Sum = 0;
		
		for (int i = 0; i < arr1.length; i++) {
			
			
			arr1Sum += arr1[i];
			
			arr2Sum += arr2[i];
			
		}
		
		double arr1Ave = (double) arr1Sum / arr1.length;
		
		double arr2Ave = (double) arr2Sum / arr2.length;
		
		double sumProduct = 0;
		
		double a1Sum = 0;
		
		double a2Sum = 0;
		
		for (int i = 0; i < arr1.length; i++) {
			
			double s1 = arr1[i] - arr1Ave;
			
			double s2 = arr2[i] - arr2Ave;
			
			sumProduct += s1 * s2;
			
			a1Sum += s1 * s1;
			
			a2Sum += s2 * s2;
			
		}
		
		if(a1Sum == 0 || a2Sum == 0) {
			
			return 0;
		}
		
		assert Math.sqrt(a1Sum) * Math.sqrt(a2Sum) != 0;
		
		double result = sumProduct / (Math.sqrt(a1Sum) * Math.sqrt(a2Sum));
		
		return result;
		
	}
	
	/**
	 * Calculates local similarity matrices for the alignment.
	 * 
	 * @return 2D array with local similarity matrices
	 */
	
	private int[][] localSimilarity2() {
		

		if(alignment == null) {
			
			throw new IllegalArgumentException("Matrix must not be null");
		}
		
		if (winWidth > alignment.numberOfColumns()) {
			
			throw new ColumnTooWideException("The width of the window is greater than the length of the allignment.");
		}
		
		int nrOfWindows = ((alignment.numberOfColumns() - alignment.numberOfColumns()%winWidth)/winWidth + ((alignment.numberOfColumns() - alignment.numberOfColumns()%winWidth)/winWidth - 1)*(winWidth - 1) + alignment.numberOfColumns()%winWidth);

		int [][] localSim = new int[nrOfWindows][alignment.numberOfRows()*(alignment.numberOfRows() - 1)/2];

		int sum = 0;

		int globalIndex = 0;
		
		int windowNr = 0;

		for(int i = 0; i < alignment.numberOfRows(); i++) {
			
			char[] rowI = alignment.getRow(i);

			for(int j = i + 1; j < alignment.numberOfRows(); j++) {
				
				char[] rowJ = alignment.getRow(j);
				
				windowNr = 0;
				
				sum = 0;
				
				for (int z = 0; z < winWidth; z++) {
					
					int index = 24 * ConservationMatrices.getIndex(rowI[z]) + ConservationMatrices.getIndex(rowJ[z]);
					
					int score = ConservationMatrices.blosum2[index];
					
					sum += score;
					
				}
				
				localSim[windowNr][globalIndex] = sum;
				
				windowNr++;

				for (int k = winWidth; k < alignment.numberOfColumns(); k++) {
						
					int index1 = 24 * ConservationMatrices.getIndex(rowI[k]) + ConservationMatrices.getIndex(rowJ[k]);
						
					int index2 = 24 * ConservationMatrices.getIndex(rowI[k - winWidth]) + ConservationMatrices.getIndex(rowJ[k - winWidth]);
						
					sum = sum - ConservationMatrices.blosum2[index2] + ConservationMatrices.blosum2[index1];
					
					localSim[windowNr][globalIndex] = sum;
					
					windowNr++;
					
				}
				
				globalIndex++;

			}
				
		}

		return localSim;

	}
	
	/** 
	 * Calculates pearson coefficient between global similarity matrix and local similarity matrices.
	 * 
	 * @return array of coefficients indexed by windows
	 */
	
	 double[] calcPearson() {
		
		if(alignment == null) {
			
			throw new IllegalArgumentException("Matrix must not be null");
		}
		
		if (winWidth > alignment.numberOfColumns()) {
			
			throw new ColumnTooWideException("The width of the window is greater than the length of the allignment.");
		}
		
		int[] global = globalSimilarity();

		int nrOfWindows = ((alignment.numberOfColumns() - alignment.numberOfColumns()%winWidth)/winWidth + ((alignment.numberOfColumns() - alignment.numberOfColumns()%winWidth)/winWidth - 1)*(winWidth - 1) + alignment.numberOfColumns()%winWidth);

		double [] coeffs = new double[nrOfWindows];

		int sum = 0;

		int globalIndex = 0;

		int counter = 0;

		int[] winValues = new int[winWidth];
		
		float[][] coeffsRaw = new float[nrOfWindows][5];
		
		for(int i = 0; i < alignment.numberOfRows(); i++) {
			
			char[] rowI = alignment.getRow(i);

			for(int j = i + 1; j < alignment.numberOfRows(); j++) {
				
				char[] rowJ = alignment.getRow(j);

					for (int k = 0; k < nrOfWindows; k++) {

							if( (k == 0 || k%winWidth == 0) && k < alignment.numberOfColumns() - alignment.numberOfColumns()%winWidth) {

								counter = 0;

								sum = 0;
								
								int colRange = (winWidth - 1) / 2;

								int midColumn = k + colRange;
								
								int range = - colRange;
								
								for(int d = 0; d < winWidth; d++) {
									
									int index = 24 * ConservationMatrices.getIndex(rowI[midColumn + range]) + ConservationMatrices.getIndex(rowJ[midColumn + range]);
									
									int score = ConservationMatrices.blosum2[index];
									
									winValues[d] = score;
									
									sum += score;
									
									range++;
									
								}
								
								range = -colRange;

								coeffsRaw[k][0] += sum;
								
								coeffsRaw[k][1] += global[globalIndex];
								
								coeffsRaw[k][2] += sum * sum;
								
								coeffsRaw[k][3] += global[globalIndex] * global[globalIndex];
								
								coeffsRaw[k][4] += global[globalIndex] * sum;

							}

							else {

								sum = sum - winValues[counter] + ConservationMatrices.BlosumPair2(alignment.getRow(i)[k + (winWidth - 1)], alignment.getRow(j)[k + (winWidth - 1)]);
								
								coeffsRaw[k][0] += sum;
								
								coeffsRaw[k][1] += global[globalIndex];
								
								coeffsRaw[k][2] += sum * sum;
								
								coeffsRaw[k][3] += global[globalIndex] * global[globalIndex];
								
								coeffsRaw[k][4] += global[globalIndex] * sum;
								
								counter++;
							}

					}
				
					globalIndex++;
					
			}

		}
		

		for (int n = 0; n < coeffsRaw.length; n++) {
			
			double numerator = coeffsRaw[n][4] - coeffsRaw[n][0] * coeffsRaw[n][1]/global.length;
			
			assert coeffsRaw[n][2] - coeffsRaw[n][0] * coeffsRaw[n][0]/(double)global.length > 0;
			
			assert coeffsRaw[n][3] - coeffsRaw[n][1] * coeffsRaw[n][1]/(double)global.length > 0;
			
			double denominator = Math.sqrt(coeffsRaw[n][2] - coeffsRaw[n][0] * coeffsRaw[n][0]/(double) global.length)  * Math.sqrt(coeffsRaw[n][3] - coeffsRaw[n][1] * coeffsRaw[n][1]/(double)global.length);
			
			coeffs[n] = numerator / denominator;
			
			}
		
		return coeffs;

	}
	
	 /**
	  * Returns correlation score.
	  * @param score
	  * @param normalize
	  * @return
	  */
	 
	 double[] getCorrelationScore(SMERFSColumnScore score, boolean normalize) {
		
		if (alignment == null) {
			
			throw new IllegalArgumentException("Alignment must not be null");
			
		}
		
		double[] results = null;
		
		if (alignment.numberOfRows() < 500) {
			
			results = calcPearsonCoeff3(); 
			
		}
		
		else {
			
			results = calcPearson2();
		}
		
		double[] columnResults;
		
		if (score == SMERFSColumnScore.MAX_SCORE) {
			
			columnResults = giveMaxToColumn(results);
		}
		
		else {
			
			columnResults = giveMidToColumn(results);
		}
		
		rejectOverTreshold(columnResults);
		
		if (normalize == true) {
			
			double [] normalized = ConservationAccessory.normalize01(columnResults);
			
			return normalized;
		}
		
		else {
		
			return columnResults;
		}
		
		
	}
	/**
	 * Finds max within a part of an array, both begin and end delimeters are also considered.
	 * @param scores array
	 * @param begin
	 * @param end
	 * @return max
	 */
	double findMax(double[] scores, int begin, int end) {
		
		if (end < begin) {
			
			throw new IllegalArgumentException("End is smaller than the begin.");
		}
		
		if (begin == end) {
			
			return scores[begin];
			
		}
		
		double max = scores[begin];
		
		for (int i = begin; i < end + 1; i++) {
			
			if(scores[i] > max) {
				
				max = scores[i];
				
			}
		
		}
		
		return max;
	}
	
	/**
	 * Gives a score to the column. The score given is the max score of all the windows it belongs to.
	 * @param windowScores
	 * @return
	 */
	
	private double[] giveMaxToColumn(double[] windowScores) {
		
		double[] scores = new double[alignment.numberOfColumns()];
		
		for (int i = 0; i < winWidth - 1; i++) {
			
			scores[i] = findMax(windowScores, 0, i);
		}
		
		for (int i = winWidth - 1; i < scores.length - (winWidth - 1); i++) {
			
			scores[i] = findMax(windowScores, i - (winWidth - 1), i);
		}
		
		int begin = windowScores.length - 1 -(winWidth - 2);
		
		int end = windowScores.length - 1;
		
		for (int i = scores.length - (winWidth - 1); i < scores.length; i++) {
			
			scores[i] = findMax(windowScores, begin, end);
			
			begin++;
			
		}
		
		return scores;
	}
	
	/**
	 * Calculates local similarity matrices. 
	 * @return 2D array indexed by window
	 */
	
	int[][] localSimilarity() {
		

		if(alignment == null) {
			
			throw new IllegalArgumentException("Matrix must not be null");
		}
		
		if (winWidth > alignment.numberOfColumns()) {
			
			throw new ColumnTooWideException("The width of the window is greater than the length of the allignment.");
		}
		
		int nrOfWindows = ((alignment.numberOfColumns() - alignment.numberOfColumns()%winWidth)/winWidth + ((alignment.numberOfColumns() - alignment.numberOfColumns()%winWidth)/winWidth - 1)*(winWidth - 1) + alignment.numberOfColumns()%winWidth);

		int [][] localSim = new int[nrOfWindows][alignment.numberOfRows()*(alignment.numberOfRows() - 1)/2];

		int sum = 0;

		int globalIndex = 0;

		int counter = 0;

		int[] winValues = new int[winWidth];

		for(int i = 0; i < alignment.numberOfRows(); i++) {
			
			char[] rowI = alignment.getRow(i);

			for(int j = i + 1; j < alignment.numberOfRows(); j++) {
				
				char[] rowJ = alignment.getRow(j);

					for (int k = 0; k < nrOfWindows; k++) {

							if( (k == 0 || k%winWidth == 0) && k < alignment.numberOfColumns() - alignment.numberOfColumns()%winWidth) {

								counter = 0;

								sum = 0;
								
								int colRange = (winWidth - 1) / 2;

								int midColumn = k + colRange;
								
								int range = -colRange; 
								
								for (int d = 0; d < winWidth; d++) {
									
									int index = 24 * ConservationMatrices.getIndex(rowI[midColumn + range]) + ConservationMatrices.getIndex(rowJ[midColumn + range]);
									
									int score = ConservationMatrices.blosum2[index];
									
									winValues[d] = score;
									
									sum += score;
									
									range++;
								}
								
								range = 0;

							}

							else {
								
								int index = 24 * ConservationMatrices.getIndex(rowI[k + (winWidth - 1)]) + ConservationMatrices.getIndex(rowJ[k + (winWidth - 1)]);
								
								int score = ConservationMatrices.blosum2[index];
								
								sum = sum - winValues[counter] + score;
								
								localSim[k][globalIndex] = sum;
								
								counter++;

							}

					}
				
					globalIndex++;
					
			}

		}
	
	return localSim;
		
	}
	
	/**
	 * Calculates pearson coefficient for the alignment.
	 * @return array of scores indexed by window
	 * 
	 */
	private double[] calcPearson2() {
		
		if(alignment == null) {
			
			throw new IllegalArgumentException("Matrix must not be null");
		}
		
		if (winWidth > alignment.numberOfColumns()) {
			
			throw new ColumnTooWideException("The width of the window is greater than the length of the allignment.");
		}
		
		int[] global = globalSimilarity();

		int nrOfWindows = ((alignment.numberOfColumns() - alignment.numberOfColumns()%winWidth)/winWidth + ((alignment.numberOfColumns() - alignment.numberOfColumns()%winWidth)/winWidth - 1)*(winWidth - 1) + alignment.numberOfColumns()%winWidth);

		double [] coeffs = new double[nrOfWindows];

		int sum = 0;
		
		int windowNr = 0;

		int globalIndex = 0;

		float[][] coeffsRaw = new float[nrOfWindows][5];
		
		for(int i = 0; i < alignment.numberOfRows(); i++) {
			
			char[] rowI = alignment.getRow(i);

			for(int j = i + 1; j < alignment.numberOfRows(); j++) {
				
				char[] rowJ = alignment.getRow(j);
				
				windowNr = 0;
				
				sum = 0;
				
				for (int z = 0; z < winWidth; z++) {
					
					int index = 24 * ConservationMatrices.getIndex(rowI[z]) + ConservationMatrices.getIndex(rowJ[z]);
					
					int score = ConservationMatrices.blosum2[index];
					
					sum += score;
					
				}
				
				coeffsRaw[windowNr][0] += sum;
				
				coeffsRaw[windowNr][1] += global[globalIndex];
				
				coeffsRaw[windowNr][2] += sum * sum;
				
				coeffsRaw[windowNr][3] += global[globalIndex] * global[globalIndex];
				
				coeffsRaw[windowNr][4] += global[globalIndex] * sum;
				
				windowNr++;

				for (int k = winWidth; k < alignment.numberOfColumns(); k++) {
						
					int index1 = 24 * ConservationMatrices.getIndex(rowI[k]) + ConservationMatrices.getIndex(rowJ[k]);
						
					int index2 = 24 * ConservationMatrices.getIndex(rowI[k - winWidth]) + ConservationMatrices.getIndex(rowJ[k - winWidth]);
						
					sum = sum - ConservationMatrices.blosum2[index2] + ConservationMatrices.blosum2[index1];
					
					coeffsRaw[windowNr][0] += sum;
					
					coeffsRaw[windowNr][1] += global[globalIndex];
					
					coeffsRaw[windowNr][2] += sum * sum;
					
					coeffsRaw[windowNr][3] += global[globalIndex] * global[globalIndex];
					
					coeffsRaw[windowNr][4] += global[globalIndex] * sum;
					
					windowNr++;
					
				}
				
				globalIndex++;
					
			}

		}
		
		for (int n = 0; n < coeffsRaw.length; n++) {
			
			double numerator = coeffsRaw[n][4] - coeffsRaw[n][0] * coeffsRaw[n][1]/global.length;
			
			assert coeffsRaw[n][2] - coeffsRaw[n][0] * coeffsRaw[n][0]/(double)global.length > 0;
			
			assert coeffsRaw[n][3] - coeffsRaw[n][1] * coeffsRaw[n][1]/(double)global.length > 0;
			
			double denominator = Math.sqrt(coeffsRaw[n][2] - coeffsRaw[n][0] * coeffsRaw[n][0]/(double) global.length)  * Math.sqrt(coeffsRaw[n][3] - coeffsRaw[n][1] * coeffsRaw[n][1]/(double)global.length);
			
			coeffs[n] = numerator / denominator;
			
			}
		
		return coeffs;

		
	}
	
	/**
	 * Gives scores to columns. The middle column gets the window score
	 * @param windowScores
	 * @return
	 */
	
	private double[] giveMidToColumn (double[] windowScores){
		
		double[] columnResults = new double[alignment.numberOfColumns()];
		
		int ends = (winWidth - 1) / 2;
		
		for (int i = 0; i < ends; i ++) {
			
			columnResults[i] = windowScores[0];
			
			columnResults[(columnResults.length - 1) - i] = windowScores[windowScores.length - 1];
		}
		
		for (int i = 0; i < windowScores.length; i++) {
			
			columnResults[i + ends] = windowScores[i];
		}
		
		return columnResults;
	}
	
	void rejectOverTreshold(double[] results ) {
		
		for(int i = 0; i < alignment.numberOfColumns(); i++) {
			
			Map<Character, Integer> colMap = alignment.getTotalAcidsFreqByCol().get(i);
			
			if (colMap.containsKey('-')) {
				
				if (colMap.get('-')/alignment.numberOfRows() > gapTreshold) {
					
					results[i] = 0.0;
				}
			}
		}
	}
	
}
	
		
	
