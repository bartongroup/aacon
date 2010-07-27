package compbio.conservation;

import java.util.*;

public class SubFamiliesConservation {
	
	final List<char[][]> subGroups;
	
	SubFamiliesConservation(AminoAcidMatrix alignment, int[][] groups) {
		
		subGroups = alignment.splitAlignment(groups);
		
	}
	
	String[][] groupsConsStat;
	
	String[][][] groupPairsConsStat;
	
	private static Map<Method, Double> treshholds = new EnumMap<Method, Double>(Method.class);
	
	static {
		
		treshholds.put(Method.ARMON_SCORE, 0.75);
		
		treshholds.put(Method.GERSTEIN_SCORE, 0.75);
		
		treshholds.put(Method.JORES_SCORE, 0.75);
		
		treshholds.put(Method.KABAT_SCORE, 0.75);
		
		treshholds.put(Method.KARLIN_SCORE, 0.75);
		
		treshholds.put(Method.LANDGRAF_SCORE, 0.75);
		
		treshholds.put(Method.MIRNY_SCORE, 0.75);
		
		treshholds.put(Method.NOT_LANCET_SCORE, 0.75);
		
		treshholds.put(Method.SANDER_SCORE, 0.75);
		
		treshholds.put(Method.SCHNEIDER_SCORE, 0.75);
		
		treshholds.put(Method.SHENKIN_SCORE, 0.75);
		
		treshholds.put(Method.TAYLOR_SCORE_GAPS, 0.75);
		
		treshholds.put(Method.TAYLOR_SCORE_NO_GAPS, 0.75);
		
		treshholds.put(Method.THOMPSON_SCORE, 0.75);
		
		treshholds.put(Method.VALDAR_SCORE, 0.75);
		
		treshholds.put(Method.WILLIAMSON_SCORE, 0.75);
		
		treshholds.put(Method.ZVELIBIL_SCORE, 0.75);
		
		
		
	}
	
	static double getTreshhold( Method method) {
		
		double tresh = treshholds.get(method);
		
		return tresh;
		
	}
	/**
	 * Merges two subfamilies into one. A result is an alignment in which the sequences are put in such an order as if the two alignments were put together by placing the second beneth the first one. 
	 * @param sub1 2D array with first alignment
	 * @param sub2 2D array with second alignment
	 * @return matrix representing merged alignments
	 */

	static char[][] mergeSubFamilies(char[][] sub1, char[][] sub2) {
		
		for (int i = 0; i < sub1.length; i++) {
			
			if (sub1[0].length != sub1[i].length) {
				
				throw new SequencesNotEquallyLongException("Sequence nr: " + i + " is of different length than teh first sequence.");
			}
		}
		
		for (int i = 0; i < sub2.length; i++) {
			
			if (sub2[0].length != sub2[i].length) {
				
				throw new SequencesNotEquallyLongException("Sequence nr: " + i + " is of different length than teh first sequence.");
			}
		}
		
		if (sub1[0].length != sub2[0].length) {
			
			throw new SequencesNotEquallyLongException("Sequences in the provided subfamilies are not of equal length.");
		}
		
		char[][] merged = new char[sub1.length + sub2.length][sub1[0].length];
		
		for (int i = 0; i < sub1.length; i++) {
		
			merged[i] = sub1[i];
		}
		
		for (int i = 0; i < sub2.length; i++) {
			
			merged[sub1.length + i] = sub2[i];
		}
		
		return merged;
	}
	
	/**
	 * Calculates conservation in subfamilies.
	 * @param method method name 
	 * @param subGroups list of subfamily alignments
	 * @param normalize if result is to be  normalized set true
	 * @return
	 */
	
	double[][] subgrupsConservation(Method method, List<char[][]> subGroups, boolean normalize) {
		
		if (method == null) {
			
			throw new IllegalArgumentException("Method must not be null.");
		}
		
		if (subGroups == null) {
			
			throw new IllegalArgumentException("List of subgroups/subclasses must not be null.");
			
		}
		
		double[][] conservation = new double[subGroups.size()][subGroups.get(0)[0].length];
		
		ConservationScores2 scores = null;
		
		AminoAcidMatrix sub = null;
		
		String[] names = null;
		
		for ( int i = 0; i < subGroups.size(); i++) {
			
			sub = new AminoAcidMatrix(subGroups.get(i),names );
			
			scores = new ConservationScores2(sub);
			
			conservation[i] = scores.calculateScore(method, normalize);
			
		}
		
		return conservation;
		
	}
	
	/**
	 * Matrix that stores newly formed subfamily pairs.
	 * It also contains info about which subfamilies did the pairs come form.
	 * Example : subPair[x][y] x is the number of the first sequence in the pair(numbering starts form 0).
	 *                         the number of the second sequence is x + 1 + y 
	 * @param subGroups
	 * @return
	 */
	
	
	char[][][][] subFamilyPairs(List<char[][]> subGroups) {
		
		if (subGroups == null) {
			
			throw new IllegalArgumentException("SubGroups list must not be null");
		}
		
		char[][][][] pairs = new char[subGroups.size()][][][];
		
		int startLength = subGroups.size() - 1;
		
		for (int i = 0; i < subGroups.size(); i++) {
			
			pairs = new char[i][startLength][][]; 
			
			for (int j = 0; j < startLength; j++) {
				
				pairs[i][j] = mergeSubFamilies(subGroups.get(i), subGroups.get(j));
			}
			
			startLength--;
		}
		
		return pairs;
	}
	
	/**
	 * Matrix that stores conservation scores for newly formed subfamily pairs.
	 * It also contains info about which subfamilies did the pairs come form.
	 * Example : subPair[x][y] x is the number of the first sequence in the pair(numbering starts form 0).
	 *                         the number of the second sequence is x + 1 + y 
	 * @param subGroups
	 * @return
	 */
	
	double[][][] subFamilyPairsConservation(Method method, List<char[][]> subGroups, boolean normalize) {
		
		if (subGroups == null) {
			
			throw new IllegalArgumentException("SubGroups list must not be null");
		}
		
		double[][][] pairs = new double[subGroups.size()][][];
		
		ConservationScores2 scores = null;
		
		String[] names = null;
		
		char[][] merged = null;
		
		int startLength = subGroups.size() - 1;
		
		for (int i = 0; i < subGroups.size(); i++) {
			
			pairs[i] = new double[i][startLength]; 
			
			for (int j = 0; j < startLength; j++) {
				
				merged = mergeSubFamilies(subGroups.get(i), subGroups.get(j));
				
				scores = new ConservationScores2(new AminoAcidMatrix(merged, names));
				
				pairs[i][j] = scores.calculateScore(method, normalize);
				
			}
			
			startLength--;
		}
		
		return pairs;
	}
	
	void subFamilyResults(double[][] subResults) {
		
		
	}
}
