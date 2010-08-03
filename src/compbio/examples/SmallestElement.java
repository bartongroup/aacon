package compbio.examples;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.*;



import compbio.conservation.AminoAcidMatrix;
import compbio.conservation.ConservationAccessory;
import compbio.conservation.Correlation;
import compbio.util.FastaSequence;
import compbio.util.SequenceUtil;
public class SmallestElement {
	
	static void printCoeffs(double[] arr1, double[] arr2) {
		
		assert arr1.length == arr2.length;
			
			for (int i = 0; i < arr1.length; i++) {
				
				System.out.println(i + ": " + arr1[i] + "-" + arr2[i]);
			}
		}
	
	static void print2DArray(int[][] arr1, int[][] arr2) {
		
		for (int i = 0; i < arr1.length; i++) {
			
			for (int j = 0; j < arr1[j].length; j++) {
				
				System.out.println(arr1[i][j] + "-" + arr2[i][j]);
			}
		}
	}
	
	
	public static void main (String[] args) {
		
		Integer one = new Integer(1);
		Integer two = new Integer(2);
		Integer three = new Integer(3);
		
		Map<String,Integer> map = new HashMap<String, Integer>();
		
		map.put("one", one);
		map.put("two",two );
		map.put("three", three);
		
		Collection<Integer> c = map.values();
		
		int min = Collections.min(c);
		
		System.out.println(min);
			
			String filePath  = "/homes/agolicz/alignments/alignment1";
			
			InputStream inStr = null;
			
			List<FastaSequence> fastaSeqs = null;
			
			try {
				
				inStr = new FileInputStream(filePath);
				
			}
			
			catch (IOException e) {
				
				System.out.println("Can not find file");
			
			}
			
			try {
				
				fastaSeqs = SequenceUtil.readFasta(inStr);
			}
			
			catch (IOException e) {
				
				System.out.println("Sth wrong with reading the file");
			}
				
				
			AminoAcidMatrix matrix = new AminoAcidMatrix(fastaSeqs);
			
			Correlation corr = new Correlation(matrix, 5);
			
			//double[] coeffs = corr.calcPearson4(matrix, 5);
			
			//double[] coeffs2 = Correlation.calcPearsonCoeff3(matrix, 5);
			
			//int[][] localSim = Correlation.localSimilarity2(matrix);
			
			//int[] globalSim = Correlation.globalSimilarity(matrix);
			
			//printCoeffs(coeffs, coeffs2);
			
			//double[] scores = corr.giveMaxToColumn(coeffs);
			
			//double[] scores2 = Correlation.getCorrelationScore(matrix, 5, false);
			
			//printCoeffs(scores, scores2);
			
			//System.out.println(matrix.numberOfColumns());
			
			long time1 = System.currentTimeMillis();
			
			//int[][] sim = Correlation.localSimilarity2(matrix, 5);
			
			long time2 = System.currentTimeMillis();
			
			long diff = time2 - time1; 
			
			System.out.println("Operation took: " + diff);
			
			long time3 = System.currentTimeMillis();
			
			//int[][] sim2 = Correlation.localSimilarity2_2(matrix, 5);
			
			long time4 = System.currentTimeMillis();
			
			long diff2 = time4 - time3; 
			
			System.out.println("Operation took: " + diff2);
			
			//print2DArray(sim, sim2);
			
			//int[] global = Correlation.globalSimilarity(matrix);
			
			//int[][] locals = Correlation.localSimilarity2_2(matrix, 5);
			
			//double[] coeffs3 = new double[locals.length];
			
			//for (int i = 0; i < locals.length; i++) {
				
			//	coeffs3[i] = Correlation.pearson2(locals[i], global);
			//}
			
			//double[] columnResults = new double[matrix.numberOfColumns()];
			
			//int ends = (5 - 1) / 2;
			
			//for (int i = 0; i < ends; i ++) {
				
			//	columnResults[i] = coeffs3[0];
				
			//	columnResults[(columnResults.length - 1) - i] = coeffs3[coeffs3.length - 1];
			//}
			
			//for (int i = 0; i < coeffs3.length; i++) {
				
			//	columnResults[i + ends] = coeffs3[i];
			//}
			
			
			//double[] normalized = ConservationAccessory.normalize01(columnResults);
			
			//double[] helper = new double[normalized.length];
			
			//printCoeffs(normalized, helper);
			
			}

		
	}


