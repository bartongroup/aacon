package compbio.conservation;

/**
 * Class has static methods used to create similarity matrices, and correaltion matrices.
 * @author agolicz
 *
 */

public class Correlation {
	
	/**
	 * Calculates similarity between two sequences, similarity is calculated as a sum of blosum scores for pairs formed by corresponding amino acids in two sequences.
	 *  
	 * @param seq1 sequence 1
	 * @param seq2 sequence 2
	 * @return  similarity score
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
	
	//static int localSequenceSimilarityNextWindows(int resultFromLastWindow, int lastIndex) {
		
	//	simi
		
		
		
		
	//}
	
	/**
	 * Creates a global similarity matrix.
	 * Always use this method along with local similarity matrix, the internal numbering scheme assures that Pearson correlation will be calculated properly
	 * @param matrix reference to the matrix containing sequences
	 * @return similarity matrix stored as a single array
	 */
	
	static int[] globalSimilarity(AminoAcidMatrix matrix) {
		
		int[] globalSim = new int[matrix.numberOfRows()* (matrix.numberOfRows() - 1) / 2];
		
		int index = 0;
		
		for (int i = 0; i < matrix.numberOfRows(); i++) {
			
			for (int j = i + 1; j < matrix.numberOfRows(); j++) {
				
				globalSim[index] = Correlation.sequenceSimilartyBlosum(matrix.getRow(i), matrix.getRow(j));
				
				index++;
				
			}
		
		}
		
		return globalSim;
		
		}

	/**
	 * Creates a local similarity matrix.
	 * Always use this method along with global similarity matrix, the internal numbering scheme assures that Pearson correlation will be calculated properly
	 * @param matrix reference to the matrix containing sequences
	 * @param startPoint
	 * @param endPoint not included in calculation
	 * @return similarity matrix stored as a single array
	 */
	
	//static int[] localSimilarity(AminoAcidMatrix matrix, int startPoint, int endPoint) {
		
		//	int[] localSim = new int[(1 + (matrix.numberOfRows() - 1)) * (matrix.numberOfRows() - 1) / 2];
		
		//	int index = 0;
		
		//	for (int i = 0; i < matrix.numberOfRows(); i++) {
			
		//		for (int j = i + 1; j < matrix.numberOfRows(); j++) {
				
		//			Correlation.localSequenceSimilarityBlosum(matrix.getRow(i), matrix.getRow(j), startPoint, endPoint);
				
		//			index++;
				
		//	}
		
	//	}
		
	//	return localSim;
		
	//	}

	/**
	 * Calculates Pearson product-moment correlation coefficient between two data sets stored in simple arrays.
	 * 
	 * @param arr1 array 1 
	 * @param arr2 array 2
	 * 
	 * @return Pearson correlation coefficient
	 */
	
	static double getPMCC (int[] arr1, int[] arr2) {
		
		if(arr1.length != arr2.length) {
			
			throw new IllegalArgumentException("Arrays not of equal length");
			
		}
		
		int sum1 = 0;
		
		int sum2 = 0;
		
		int sum1Sq = 0;
		
		int sum2Sq = 0;
		
		int pSum = 0;
		
		int n = arr1.length;
		
		for (int i = 0; i < arr1.length; i++) {
			
			sum1 += arr1[i];
			
			sum2 += arr1[i];
			
			sum1Sq += arr1[i] * arr1[i];
			
			sum2Sq += arr2[i] * arr2[i];
			
			pSum += arr1[i] * arr2[i];
			
			}
		
		double numerator = pSum - sum1 * sum2 / (double) n;
		
		assert sum1Sq - sum1 * sum1 / (double) n > 0;
		
		assert sum2Sq - sum2 * sum2 / (double) n > 0;
		
		double denominator = Math.sqrt((sum1Sq - sum1 * sum1 / (double) n) * (sum2Sq - sum2 * sum2 / (double) n));
		
		return numerator / denominator;
		
		
	}
	
	/**
	 * Values needed for calculating Pearson coefficients ar eindexed as follows
	 * 0 = sum1
	 * 1 = sum2
	 * 2 = sum1Sq
	 * 3 = sum2Sq
	 * 4 = pSum
	 * 
	 * @param matrix
	 * @return
	 */



	static double[] calcPearsonCoeff(AminoAcidMatrix matrix) {
		
		int[] global =  Correlation.globalSimilarity(matrix);
		
		int[][] coeffsRaw = new int[matrix.numberOfColumns()][5]; 
		
		double[] coeffs = new double[coeffsRaw.length];
		
		int globalIndex = 0;
		
		int sum = 0;
		
		int windowCounter = 0;
		 
		int counter = 0;
		 
		assert counter < 6;
		
		int[] winValues = new int[7];
		
		for (int i = 0; i < matrix.numberOfRows(); i++) {
			
			 for (int j = i + 1; j < matrix.numberOfRows(); j++ ) {
				 
				 int k = 0;
				 
				 while (k < matrix.numberOfColumns()) {
					 
					  if (k == 0 || (k % 7 == 0 && k < matrix.numberOfColumns() - (matrix.numberOfColumns()%7))) {
						 
						sum = 0;
						 
						for (int l = 0; l < 7; l++) {
							
						winValues[l] = ConservationMatrices.BlosumPair2(matrix.getRow(i)[k + l], matrix.getRow(j)[k + l]);	 
						
						sum += winValues[l];
						
						}
						
						coeffsRaw[windowCounter][0] += sum;
						
						coeffsRaw[windowCounter][1] += global[globalIndex];
						
						coeffsRaw[windowCounter][2] += sum * sum;
						
						coeffsRaw[windowCounter][3] += global[globalIndex] * global[globalIndex];
						
						coeffsRaw[windowCounter][4] += global[globalIndex] * sum;
						
						k = k + 7;
						
					 }
						
					else {

						sum = sum - winValues[counter] + ConservationMatrices.BlosumPair2(matrix.getRow(i)[k], matrix.getRow(j)[k]);
						
						coeffsRaw[windowCounter][0] += sum;
						
						coeffsRaw[windowCounter][1] += global[globalIndex];
						
						coeffsRaw[windowCounter][2] += sum * sum;
						
						coeffsRaw[windowCounter][3] += global[globalIndex] * global[globalIndex];
						
						coeffsRaw[windowCounter][4] += global[globalIndex] * sum;
						
						k++;
						
						counter++;
						 
						 if (counter == 6) {
							 
							 counter = 0;
							 
						
						}
					 
					 windowCounter++; 
					  
					 
					 }
					 
				 }
				 
				 globalIndex++;
				 
				 windowCounter = 0;
					 
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
	 * Values needed for calculating Pearson coefficients ar eindexed as follows
	 * 0 = sum1Sq
	 * 1 = sum2Sq
	 * 2 = sumProduct
	 * 
	 * @param matrix
	 * @return
	 */

				


	static double[] calcPearsonCoeff2(AminoAcidMatrix matrix) {
	
	int[] global =  Correlation.globalSimilarity(matrix);
	
	int[][] coeffsRaw = new int[matrix.numberOfColumns()][5]; 
	
	double[] coeffs = new double[coeffsRaw.length];
	
	int globalIndex = 0;
	
	int sum = 0;
	
	int windowCounter = 0;
	 
	int counter = 0;
	 
	assert counter < 6;
	
	int[] winValues = new int[7];
	
	float[] mean1 = new float[global.length];
	
	float[] mean2 = new float[global.length];
	
	for (int i = 0; i < matrix.numberOfRows(); i++) {
		
		 for (int j = i + 1; j < matrix.numberOfRows(); j++ ) {
			 
			 for (int k = 0; k < matrix.numberOfColumns(); k++) {
				 
				  if (k == 0 || (k % 7 == 0 && k < matrix.numberOfColumns() - (matrix.numberOfColumns()%7))) {
					 
					sum = 0;
					 
					for (int l = 0; l < 7; l++) {
						
					winValues[l] = ConservationMatrices.BlosumPair2(matrix.getRow(i)[k + l], matrix.getRow(j)[k + l]);	 
					
					sum += winValues[l];
					
					}
					
					if(i == 0 && j == 1) {
						
						mean1[windowCounter] = sum; 
						
						mean2[windowCounter] = global[globalIndex];
					}
					
					else {
						
					float sweep = (globalIndex)/ (float) globalIndex;
					
					float delta1 = sum - mean1[windowCounter];
					
					float delta2 = sum - mean2[windowCounter];
						
					coeffsRaw[windowCounter][0] += delta1 * delta1 * sweep;
					
					coeffsRaw[windowCounter][1] += delta2 * delta2 * sweep;
					
					coeffsRaw[windowCounter][2] += delta1 * delta2 * sweep;
					
					mean1[windowCounter] += delta1 / (float) globalIndex;
					
					mean2[windowCounter] += delta2 / (float) globalIndex;
					
					}
					
					k = k + 6;
					
				 }
					
				else {
					
					sum = sum - winValues[counter] + ConservationMatrices.BlosumPair2(matrix.getRow(i)[k], matrix.getRow(j)[k]);
					
					if(i == 0 && j == 1) {
						
						mean1[windowCounter] = sum;
						
						mean2[windowCounter] = global[globalIndex];
					}
					
					else {

					float sweep = (globalIndex)/ (float) globalIndex;

					float delta1 = sum - mean1[windowCounter];
					
					float delta2 = sum - mean2[windowCounter];
						
					coeffsRaw[windowCounter][0] += delta1 * delta1 * sweep;
					
					coeffsRaw[windowCounter][1] += delta2 * delta2 * sweep;
					
					coeffsRaw[windowCounter][2] += delta1 * delta2 * sweep;
					
					mean1[windowCounter] += delta1 / (float) globalIndex;
					
					mean2[windowCounter] += delta2 / (float) globalIndex;
					
					}
					
					counter++;
					
					if (counter == 6) {
						 
						 counter = 0;
					
					}
				 
				 windowCounter++; 
				 
				 
					 
				 }
				 
			 }
			 
			 globalIndex++;
			 
			 windowCounter = 0;
				 
		 }
			
	}
	
	for (int n = 0; n < coeffsRaw.length; n++) {
		
		float numerator = coeffsRaw[n][2]/ (float) global.length;
		
		assert coeffsRaw[n][0]/(float)global.length > 0;
		
		assert coeffsRaw[n][1]/(float)global.length > 0;
		
		double denominator = Math.sqrt(coeffsRaw[n][0]/(float) global.length)  * Math.sqrt(coeffsRaw[n][1]/(float)global.length);
		
		coeffs[n] = numerator / denominator;
		
		}
	
	return coeffs;
	
}
			

	
		

	static int[][] localSimilarity(AminoAcidMatrix matrix) {
	
		
	int nrOfWindows = ((matrix.numberOfColumns() - matrix.numberOfColumns()%7)/ 7 + ((matrix.numberOfColumns() - matrix.numberOfColumns()%7)/ 7 - 1) * 6 + matrix.numberOfColumns()%7);
	
	int[][] localSim = new int[nrOfWindows][matrix.numberOfRows()*(matrix.numberOfRows() - 1)/ 2];
		
	int globalIndex = 0;
	
	int windowCounter = 0;
	
	int sum = 0;
	
	int counter = 0;
	
	boolean stat = true;
	
	int[] winValues = new int[7];
	
	for (int i = 0; i < matrix.numberOfRows(); i++) {
		
		 for (int j = i + 1; j < matrix.numberOfRows(); j++ ) {
			 
			 int k = 0;
			 
			 while (k < matrix.numberOfColumns()) {
				 
				  if ((k == 0 || (k % 7 == 0 && k < matrix.numberOfColumns() - (matrix.numberOfColumns()%7))) && stat == true) {
					 
					sum = 0;
					
					if ( k == 0) {
						
						for (int l = 0; l < 7; l++) {
							
							winValues[l] = ConservationMatrices.BlosumPair2(matrix.getRow(i)[k + l], matrix.getRow(j)[k + l]);	 
							
							sum += winValues[l];
							
							}
						
					}
					
					else {
					 
						for (int l = 0; l < 7; l++) {
						
							winValues[l] = ConservationMatrices.BlosumPair2(matrix.getRow(i)[k - 7 + l], matrix.getRow(j)[k -7 + l]);	 
					
							sum += winValues[l];
					
						}
						
					}
					
					localSim[windowCounter][globalIndex] = sum;
					
					System.out.println(windowCounter + ":" + sum);
					
					k = k + 7;
					
					stat = false;
					
					windowCounter = windowCounter + 1;
					
				 }
					
				else {
					
					sum = sum - winValues[counter] + ConservationMatrices.BlosumPair2(matrix.getRow(i)[k], matrix.getRow(j)[k]);
					
					System.out.println(windowCounter + ":" + sum);
					
					assert windowCounter < nrOfWindows;
					
					assert globalIndex < matrix.numberOfRows() * (matrix.numberOfRows() - 1) / 2;
					
					localSim[windowCounter][globalIndex] = sum;
					
					counter++;
					
					if (counter == 6) {
						 
						 counter = 0;
						 
						 k++;
						 
					}
			
					k++;
					
					windowCounter = windowCounter + 1;
			
				}
				 
			 }
			 
			 globalIndex++;
			 
			 windowCounter = 0;
				 
		 }
			
		 
	}
	
	return localSim;
	
}
	
	static double[] calcPearsonCoeff3(AminoAcidMatrix matrix, int width) {
		
		int[] global = Correlation.globalSimilarity(matrix);
		
		int[][] locals = Correlation.localSimilarity2(matrix, width);
		
		double[] coeffs = new double[locals.length];
		
		for (int i = 0; i < locals.length; i++) {
			
			coeffs[i] = Correlation.pearson2(locals[i], global);
		}
		
		return coeffs;
	}
	
			


	static double pearson2(int[] arr1, int[] arr2){
		
		assert arr1.length == arr2.length;
	
		int arr1Sum = 0;
		
		int arr2Sum = 0;
		
		for (int i = 0; i < arr1.length; i++) {
			
			arr1Sum += arr1[i];
			
			arr2Sum += arr2[i];
			
		}
		
		float arr1Ave = (float) arr1Sum / arr1.length;
		
		float arr2Ave = (float) arr2Sum / arr2.length;
		
		float sumProduct = 0;
		
		float a1Sum = 0;
		
		float a2Sum = 0;
		
		for (int i = 0; i < arr1.length; i++) {
			
			float s1 = arr1[i] - arr1Ave;
			
			float s2 = arr2[i] - arr2Ave;
			
			sumProduct += s1 * s2;
			
			a1Sum += s1 * s1;
			
			a2Sum += s2 * s2;
			
		}
		
		double result = sumProduct / (Math.sqrt(a1Sum) * Math.sqrt(a2Sum));
		
		return result;
		
	}
	
	
	static int[][] localSimilarity2(AminoAcidMatrix matrix, int colWidth) {

		int nrOfWindows = ((matrix.numberOfColumns() - matrix.numberOfColumns()%colWidth)/colWidth + ((matrix.numberOfColumns() - matrix.numberOfColumns()%colWidth)/colWidth - 1)*(colWidth - 1) + matrix.numberOfColumns()%colWidth);

		int [][] localSim = new int[nrOfWindows][matrix.numberOfRows()*(matrix.numberOfRows() - 1)/2];

		int sum = 0;

		int globalIndex = 0;

		int counter = 0;

		int[] winValues = new int[colWidth];

		for(int i = 0; i < matrix.numberOfRows(); i++) {

			for(int j = i + 1; j < matrix.numberOfRows(); j++) {

					for (int k = 0; k < nrOfWindows; k++) {

							if( (k == 0 || k%colWidth == 0) && k < matrix.numberOfColumns() - matrix.numberOfColumns()%colWidth) {

								counter = 0;

								sum = 0;
								
								int colRange = (colWidth - 1) / 2;

								int midColumn = k + colRange;
								
								int range = -colRange; 
								
								for (int d = 0; d < colWidth; d++) {
									
									int score = ConservationMatrices.BlosumPair2(matrix.getRow(i)[midColumn + range], matrix.getRow(j)[midColumn + range]);
									
									winValues[d] = score;
									
									sum += score;
									
									range++;
								}
								
								range = 0;

								//int sum0 =	ConservationMatrices.BlosumPair2(matrix.getRow(i)[midColumn - 3], matrix.getRow(j)[midColumn - 3]);

								//int sum1 =  ConservationMatrices.BlosumPair2(matrix.getRow(i)[midColumn - 2], matrix.getRow(j)[midColumn - 2]);
									
								//int sum2 =  ConservationMatrices.BlosumPair2(matrix.getRow(i)[midColumn - 1], matrix.getRow(j)[midColumn - 1]);

								//int sum3 =  ConservationMatrices.BlosumPair2(matrix.getRow(i)[midColumn], matrix.getRow(j)[midColumn]);

								//int sum4 =  ConservationMatrices.BlosumPair2(matrix.getRow(i)[midColumn + 1], matrix.getRow(j)[midColumn + 1]);

								//int sum5 =  ConservationMatrices.BlosumPair2(matrix.getRow(i)[midColumn + 2], matrix.getRow(j)[midColumn + 2]);

								//int sum6 =  ConservationMatrices.BlosumPair2(matrix.getRow(i)[midColumn + 3], matrix.getRow(j)[midColumn + 3]);

								//winValues[0] = sum0;

								//winValues[1] = sum1;

								//winValues[2] = sum2;

								//winValues[3] = sum3;

								//winValues[4] = sum4;

								//winValues[5] = sum5;

								//winValues[6] = sum6;

								//sum = sum0 + sum1 + sum2 + sum3 + sum4 + sum5 + sum6;

								localSim[k][globalIndex] = sum;
								
								//System.out.println( k + ":" + globalIndex);

							}

							else {

								sum = sum - winValues[counter] + ConservationMatrices.BlosumPair2(matrix.getRow(i)[k + (colWidth - 1)], matrix.getRow(j)[k + (colWidth - 1)]);
								
								localSim[k][globalIndex] = sum;
								
								//System.out.println(k + ":" + globalIndex);
								
								counter++;

							}

					}
				
					globalIndex++;
					
			}

		}
		
		return localSim;

	}
	
	static double[] calcPearson3(AminoAcidMatrix matrix) {
		
		int[] global = Correlation.globalSimilarity(matrix);

		int nrOfWindows = ((matrix.numberOfColumns() - matrix.numberOfColumns()%7)/7 + ((matrix.numberOfColumns() - matrix.numberOfColumns()%7)/7 - 1)*6 + matrix.numberOfColumns()%7);

		double [] coeffs = new double[nrOfWindows];

		int sum = 0;

		int globalIndex = 0;

		int counter = 0;

		int[] winValues = new int[7];
		
		float[][] coeffsRaw = new float[nrOfWindows][3];
		
		float[] mean1 = new float[nrOfWindows];
		
		float[] mean2 = new float[nrOfWindows];

		for(int i = 0; i < matrix.numberOfRows(); i++) {

			for(int j = i + 1; j < matrix.numberOfRows(); j++) {

					for (int k = 0; k < nrOfWindows; k++) {

							if( (k == 0 || k%7 == 0) && k < matrix.numberOfColumns() - matrix.numberOfColumns()%7) {

								counter = 0;

								sum = 0;

								int midColumn = k + 3;

								int sum0 =	ConservationMatrices.BlosumPair2(matrix.getRow(i)[midColumn - 3], matrix.getRow(j)[midColumn - 3]);

								int sum1 =  ConservationMatrices.BlosumPair2(matrix.getRow(i)[midColumn - 2], matrix.getRow(j)[midColumn - 2]);
									
								int sum2 =  ConservationMatrices.BlosumPair2(matrix.getRow(i)[midColumn - 1], matrix.getRow(j)[midColumn - 1]);

								int sum3 =  ConservationMatrices.BlosumPair2(matrix.getRow(i)[midColumn], matrix.getRow(j)[midColumn]);

								int sum4 =  ConservationMatrices.BlosumPair2(matrix.getRow(i)[midColumn + 1], matrix.getRow(j)[midColumn + 1]);

								int sum5 =  ConservationMatrices.BlosumPair2(matrix.getRow(i)[midColumn + 2], matrix.getRow(j)[midColumn + 2]);

								int sum6 =  ConservationMatrices.BlosumPair2(matrix.getRow(i)[midColumn + 3], matrix.getRow(j)[midColumn + 3]);

								winValues[0] = sum0;

								winValues[1] = sum1;

								winValues[2] = sum2;

								winValues[3] = sum3;

								winValues[4] = sum4;

								winValues[5] = sum5;

								winValues[6] = sum6;

								sum = sum0 + sum1 + sum2 + sum3 + sum4 + sum5 + sum6;
								
								if(i == 0 && j == 1) {
									
									mean1[k] = sum; 
									
									mean2[k] = global[globalIndex];
								}
								
								else {
									
								float sweep = (globalIndex)/ (float) globalIndex;
								
								float delta1 = sum - mean1[k];
								
								float delta2 = sum - mean2[k];
									
								coeffsRaw[k][0] += delta1 * delta1 * sweep;
								
								coeffsRaw[k][1] += delta2 * delta2 * sweep;
								
								coeffsRaw[k][2] += delta1 * delta2 * sweep;
								
								mean1[k] += delta1 / (float) globalIndex;
								
								mean2[k] += delta2 / (float) globalIndex;
								
								}
								
								
								//System.out.println( k + ":" + globalIndex);

							}

							else {

								sum = sum - winValues[counter] + ConservationMatrices.BlosumPair2(matrix.getRow(i)[k + 6], matrix.getRow(j)[k + 6]);
								
								//System.out.println(k + ":" + globalIndex);
								
								counter++;
								
								if(i == 0 && j == 1) {
									
									mean1[k] = sum;
									
									mean2[k] = global[globalIndex];
								}
								
								else {

								float sweep = (globalIndex)/ (float) globalIndex;

								float delta1 = sum - mean1[k];
								
								float delta2 = sum - mean2[k];
									
								coeffsRaw[k][0] += delta1 * delta1 * sweep;
								
								coeffsRaw[k][1] += delta2 * delta2 * sweep;
								
								coeffsRaw[k][2] += delta1 * delta2 * sweep;
								
								mean1[k] += delta1 / (float) globalIndex;
								
								mean2[k] += delta2 / (float) globalIndex;
								
								}
								

							}

					}
				
					globalIndex++;
					
			}

		}
		
		for (int n = 0; n < coeffsRaw.length; n++) {
			
			float numerator = coeffsRaw[n][2]/ (float) global.length;
			
			assert coeffsRaw[n][0]/(float)global.length > 0;
			
			assert coeffsRaw[n][1]/(float)global.length > 0;
			
			double denominator = Math.sqrt(coeffsRaw[n][0]/(float) global.length)  * Math.sqrt(coeffsRaw[n][1]/(float)global.length);
			
			coeffs[n] = numerator / denominator;
			
			}
		
		return coeffs;

	}
	

	
	static double[] calcPearson4(AminoAcidMatrix matrix, int colWidth) {
		
		int[] global = Correlation.globalSimilarity(matrix);

		int nrOfWindows = ((matrix.numberOfColumns() - matrix.numberOfColumns()%colWidth)/colWidth + ((matrix.numberOfColumns() - matrix.numberOfColumns()%colWidth)/colWidth - 1)*(colWidth - 1) + matrix.numberOfColumns()%colWidth);

		double [] coeffs = new double[nrOfWindows];

		int sum = 0;

		int globalIndex = 0;

		int counter = 0;

		int[] winValues = new int[colWidth];
		
		float[][] coeffsRaw = new float[nrOfWindows][5];
		
		for(int i = 0; i < matrix.numberOfRows(); i++) {

			for(int j = i + 1; j < matrix.numberOfRows(); j++) {

					for (int k = 0; k < nrOfWindows; k++) {

							if( (k == 0 || k%colWidth == 0) && k < matrix.numberOfColumns() - matrix.numberOfColumns()%colWidth) {

								counter = 0;

								sum = 0;
								
								int colRange = (colWidth - 1) / 2;

								int midColumn = k + colRange;
								
								int range = - colRange;
								
								for(int d = 0; d < colWidth; d++) {
									
									int score = ConservationMatrices.BlosumPair2(matrix.getRow(i)[midColumn + range], matrix.getRow(j)[midColumn + range]);
									
									winValues[d] = score;
									
									sum += score;
									
									range++;
									
								}
								
								range = -colRange;

								//int sum0 =	ConservationMatrices.BlosumPair2(matrix.getRow(i)[midColumn - 3], matrix.getRow(j)[midColumn - 3]);

								//int sum1 =  ConservationMatrices.BlosumPair2(matrix.getRow(i)[midColumn - 2], matrix.getRow(j)[midColumn - 2]);
									
								//int sum2 =  ConservationMatrices.BlosumPair2(matrix.getRow(i)[midColumn - 1], matrix.getRow(j)[midColumn - 1]);

								//int sum3 =  ConservationMatrices.BlosumPair2(matrix.getRow(i)[midColumn], matrix.getRow(j)[midColumn]);

								//int sum4 =  ConservationMatrices.BlosumPair2(matrix.getRow(i)[midColumn + 1], matrix.getRow(j)[midColumn + 1]);

								//int sum5 =  ConservationMatrices.BlosumPair2(matrix.getRow(i)[midColumn + 2], matrix.getRow(j)[midColumn + 2]);

								//int sum6 =  ConservationMatrices.BlosumPair2(matrix.getRow(i)[midColumn + 3], matrix.getRow(j)[midColumn + 3]);

								//winValues[0] = sum0;

								//winValues[1] = sum1;

								//winValues[2] = sum2;

								//winValues[3] = sum3;

								//winValues[4] = sum4;

								//winValues[5] = sum5;

								//winValues[6] = sum6;

								//sum = sum0 + sum1 + sum2 + sum3 + sum4 + sum5 + sum6;
								
								coeffsRaw[k][0] += sum;
								
								coeffsRaw[k][1] += global[globalIndex];
								
								coeffsRaw[k][2] += sum * sum;
								
								coeffsRaw[k][3] += global[globalIndex] * global[globalIndex];
								
								coeffsRaw[k][4] += global[globalIndex] * sum;
								
								
								//System.out.println( k + ":" + globalIndex);

							}

							else {

								sum = sum - winValues[counter] + ConservationMatrices.BlosumPair2(matrix.getRow(i)[k + (colWidth - 1)], matrix.getRow(j)[k + (colWidth - 1)]);
								
								coeffsRaw[k][0] += sum;
								
								coeffsRaw[k][1] += global[globalIndex];
								
								coeffsRaw[k][2] += sum * sum;
								
								coeffsRaw[k][3] += global[globalIndex] * global[globalIndex];
								
								coeffsRaw[k][4] += global[globalIndex] * sum;
								
								
								//System.out.println(k + ":" + globalIndex);
								
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
	
	
}

	
		
	
