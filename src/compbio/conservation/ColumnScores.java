package compbio.conservation;

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class ColumnScores {

	// private enum Method { kabatScore, joresScore, schneiderScore,
	// shenkinScore, gersteinScore, SmallestTaylorSetGaps,
	// SmallestTaylorSetNoGaps, zvelibilScore, karlinScore, armonscore,
	// thompsonScore, notLancetScore, mirnyScore, williamsonScore,
	// landgrafScore, sanderScore, valdarScore };
	/**
	 * Stores a map of scores and method names. Scores stored are not
	 * normalized.
	 */
	// private final Map<Method, double[]> scores;
	/**
	 * Stores reference to amino acid matrix
	 */
	// private AminoAcidMatrix matrix;
	/**
	 * Constructor
	 * 
	 * @param matrix
	 *            AminoAcidMatrix based on the alignment
	 */
	// public ColumnScores(AminoAcidMatrix matrix) {
	// this.matrix = matrix;
	// scores = new EnumMap<Method, double[]>(Method.class);
	// }
	/**
	 * Checks if all but one residues in the column are gaps.
	 * 
	 * @return true if all but one residues are gaps, false if not
	 */
	private static boolean allButOneGaps(int columnNr, AminoAcidMatrix matrix) {

		if (columnNr > matrix.numberOfColumns() - 1) {
			throw new IllegalArgumentException(
					"Column number greater than number of columns in teh matrix.");
		}
		Map<Character, Integer> acidsInt = matrix.getTotalAcidsFreqByCol().get(columnNr);
		if (acidsInt.containsKey('-')
				&& acidsInt.get('-') == matrix.getInverseMatrix()[columnNr].length - 1) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * Checks whether there is only one residue type in the column. Gap is not
	 * counted as a residue type.
	 * 
	 * @return returns true if there is only one residue type in the column,
	 *         false if not
	 */
	private static boolean oneResidueTypeNoGaps(int columnNr, AminoAcidMatrix matrix) {

		if (columnNr > matrix.numberOfColumns() - 1) {
			throw new IllegalArgumentException(
					"Column number greater than number of columns in teh matrix.");
		}
		Map<Character, Integer> acidsInt = matrix.getTotalAcidsFreqByCol().get(columnNr);
		if (acidsInt.size() == 1 && acidsInt.containsKey('-') == false) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * Checks whether column contains gaps.
	 * 
	 * @return true if column contains gaps
	 */
	private static boolean containsGaps(int columnNr, AminoAcidMatrix matrix) {

		if (columnNr > matrix.numberOfColumns() - 1) {
			throw new IllegalArgumentException(
					"Column number greater than number of columns in teh matrix.");
		}
		if (matrix.getTotalAcidsFreqByCol().get(columnNr).containsKey('-')) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * Counts the number of different amino acids in the column. Gap is counted
	 * as a 21 amino acid.
	 * 
	 * @return the number of different amino acids
	 */
	private static int numberOfAcidsWithGap(int columnNr, AminoAcidMatrix matrix) {

		if (columnNr > matrix.numberOfColumns() - 1) {
			throw new IllegalArgumentException(
					"Column number greater than number of columns in teh matrix.");
		}
		return matrix.getTotalAcidsFreqByCol().get(columnNr).size();
	}

	/**
	 * Counts the number of different amino acids in the column. Gap not counted
	 * as a 21 amino acid.
	 * 
	 * @return the number of different amino acids
	 */
	private static int numberOfAcidsNoGap(int columnNr, AminoAcidMatrix matrix) {

		if (columnNr > matrix.numberOfColumns() - 1) {
			throw new IllegalArgumentException(
					"Column number greater than number of columns in teh matrix.");
		}
		if (containsGaps(columnNr, matrix) == true) {
			return matrix.getTotalAcidsFreqByCol().get(columnNr).size() - 1;
		} else {
			return matrix.getTotalAcidsFreqByCol().get(columnNr).size();
		}
	}

	/**
	 * Calculates the number of times the most common amino acid in the column
	 * occurs. Does not count gap as 21 amino acid.
	 * 
	 * @return times the most common amino acid occurs
	 */
	private static int mostCommonNumber(int columnNr, AminoAcidMatrix matrix) {

		if (columnNr > matrix.numberOfColumns() - 1) {
			throw new IllegalArgumentException(
					"Column number greater than number of columns in teh matrix.");
		}
		int max = 0;
		Map<Character, Integer> acidsInt = matrix.getTotalAcidsFreqByCol().get(columnNr);
		Set<Character> keys = acidsInt.keySet();
		Iterator<Character> itr = keys.iterator();
		while (itr.hasNext()) {
			Character key = itr.next();
			int value = acidsInt.get(key);
			if (key != '-' && value > max) {
				max = value;
			}
		}
		assert max != 0 : "Zero in the most Common Number";
		return max;
	}

	// Map<Character,Integer> getAcidsIntMap() {
	// return acidsIntMap;
	// }
	/**
	 * Gets length of the column.
	 * 
	 * @return length of the column
	 */
	// int length() {
	// return columnArr.length;
	// }
	/**
	 * Calculates Kabat score for the alignment.
	 * 
	 * @param normalize
	 *            to be set true if the results are to be normalized, false
	 *            otherwise
	 * @return array of results, index corresponding to the column index,
	 *         indexing starts with 0
	 */
	static double kabatScore(AminoAcidMatrix matrix, int columnNumber) {

		if (columnNumber > matrix.numberOfColumns() - 1) {
			throw new IllegalArgumentException(
					"Column number greater than number of columns in teh matrix.");
		}
		assert mostCommonNumber(columnNumber, matrix) > 0
				&& mostCommonNumber(columnNumber, matrix) < matrix.getInverseMatrix()[columnNumber].length + 1;
		double result = matrix.getInverseMatrix()[columnNumber].length
				* (double) numberOfAcidsNoGap(columnNumber, matrix)
				/ mostCommonNumber(columnNumber, matrix);
		return result;
	}

	/**
	 * Calculates Jores score for the alignment. When calculating the number of
	 * distinct pair formed calculates a combination of two out of the amino
	 * acid types present in the column. Than adds the self pairs. However when
	 * there is only one of a type in the column no self pair is formed.
	 * 
	 * @param normalize
	 *            is to be set true if data is to be normalized false otherwise
	 * @return array of Jores scores ,index correspond to the column index,
	 *         indexing starts with 0
	 */
	static double joresScore(AminoAcidMatrix matrix, int colNr) {

		if (colNr > matrix.numberOfColumns() - 1) {
			throw new IllegalArgumentException(
					"Column number greater than number of columns in teh matrix.");
		}
		double result = 0.0;
		// special case #1 one residue type only
		// special case #2 all but one are gaps
		boolean oneRes = oneResidueTypeNoGaps(colNr, matrix);
		boolean allButOne = allButOneGaps(colNr, matrix);
		if (oneRes == true || allButOne == true) {
			if (oneRes == true) {
				result = 1.0;
				return result;
			} else {
				int length = matrix.getInverseMatrix().length;
				result = length * (length - 1) / 2;
				return result;
			}
		} else {
			int samePairs = 0;
			int differentPairs = 0;
			Map<Character, Integer> acidsIntMap = matrix.getTotalAcidsFreqByCol().get(colNr);
			Map<Character, Integer> acidsIntMapCopy = new HashMap<Character, Integer>(acidsIntMap);
			Set<Character> keys = acidsIntMapCopy.keySet();
			keys.remove('-');
			int types = keys.size();
			differentPairs = types * (types - 1) / 2;
			Iterator<Integer> itr = acidsIntMapCopy.values().iterator();
			while (itr.hasNext()) {
				if (itr.next() > 1) {
					samePairs++;
				}
			}
			int totalPairs = samePairs + differentPairs;
			Iterator<Character> itr2 = acidsIntMapCopy.keySet().iterator();
			int max1 = 0;
			Character maxKey = null;
			while (itr2.hasNext()) {
				Character key = itr2.next();
				int value = acidsIntMapCopy.get(key);
				if (value > max1) {
					maxKey = key;
					max1 = value;
				}
			}
			acidsIntMapCopy.remove(maxKey);
			Iterator<Integer> itr3 = acidsIntMapCopy.values().iterator();
			int max2 = 0;
			while (itr3.hasNext()) {
				int value = itr3.next();
				if (value > max2) {
					max2 = value;
				}
			}
			int mostFreqNr = 0;
			if (max2 == 0) {
				mostFreqNr = (max1) * (max1 - 1) / 2;
			} else {
				if (max1 == max2) {
					mostFreqNr = max1 * max2;
				} else {
					int same = (max1) * (max1 - 1) / 2;
					int diff = max1 * max2;
					if (same > diff) {
						mostFreqNr = same;
					} else {
						mostFreqNr = diff;
					}
				}
			}
			int length = matrix.getInverseMatrix()[colNr].length;
			result = ((double) totalPairs / (double) mostFreqNr) * (length * (length - 1) / 2);
			return result;
		}
	}

	// Symbol Enthropy Scores
	/**
	 * Calculates Schneider score for the alignment.
	 * 
	 * @param normalize
	 *            is to be set true if data is to be normalized false otherwise
	 * @return array of Schneider scores ,index correspond to the column index,
	 *         indexing starts with 0
	 */
	static double schneiderScore(AminoAcidMatrix matrix, int colNr) {

		if (colNr > matrix.numberOfColumns() - 1) {
			throw new IllegalArgumentException(
					"Column number greater than number of columns in teh matrix.");
		}
		double result = 0.0;
		double normal = 1.0 / Math.log(20.0);
		result = ShannonEnthropy.ShannonLn(matrix.getTotalAcidsFreqByCol().get(colNr), matrix
				.getInverseMatrix()[colNr].length)
				* normal;
		assert result >= 0 && result <= 1;
		return result;
	}

	/**
	 * Calculates Shenkin score for the alignment.
	 * 
	 * @param normalize
	 *            is to be set true if data is to be normalized false otherwise
	 * @return array of Shenkin scores ,index correspond to the column index,
	 *         indexing starts with 0
	 */
	static double shenkinScore(AminoAcidMatrix matrix, int colNr) {

		if (colNr > matrix.numberOfColumns() - 1) {
			throw new IllegalArgumentException(
					"Column number greater than number of columns in teh matrix.");
		}
		double result = 0.0;
		result = Math.pow(2.0, ShannonEnthropy.ShannonLog2(matrix.getTotalAcidsFreqByCol().get(
				colNr), matrix.getInverseMatrix()[colNr].length)) * 6.0;
		assert result >= 6 && result <= 120;
		return result;
	}

	/**
	 * Calculates Gerstein score.
	 * 
	 * @param normalize
	 *            is to be set true if data is to be normalized false otherwise
	 * @return array of Gerstein scores ,index correspond to the column index,
	 *         indexing starts with 0
	 */
	static double gersteinScore(AminoAcidMatrix matrix, int colNr) {

		if (colNr > matrix.numberOfColumns() - 1) {
			throw new IllegalArgumentException(
					"Column number greater than number of columns in teh matrix.");
		}
		double result = 0.0;
		result = -(ShannonEnthropy.ShannonLn(matrix.totalAcidsFrequency(), matrix.numberOfColumns()
				* matrix.numberOfRows()))
				- (-ShannonEnthropy.ShannonLn(matrix.getTotalAcidsFreqByCol().get(colNr), matrix
						.getInverseMatrix()[colNr].length));
		return result;
	}

	/**
	 * Calculates Taylor score. Returns the number of elements in the smallest
	 * set covering the whole column. Gaps sore is included in the largest set.
	 * Therefore a column containing gaps is automatically given the lowest
	 * possibel score.
	 * 
	 * @param normalize
	 *            is to be set true if data is to be normalized false otherwise
	 * @return array of Taylor scores ,index correspond to the column index,
	 *         indexing starts with 0
	 */
	static double taylorScoreGaps(AminoAcidMatrix matrix, int colNr) {

		if (colNr > matrix.numberOfColumns() - 1) {
			throw new IllegalArgumentException(
					"Column number greater than number of columns in teh matrix.");
		}
		Map<String, HashSet<Character>> setMap = ConservationSets.taylorSets();
		double smallestSetSize = 0.0;
		Map<String, Integer> repSets = new HashMap<String, Integer>();
		Set<String> setMapKeys = setMap.keySet();
		Iterator<String> itr = setMapKeys.iterator();
		while (itr.hasNext()) {
			String key = itr.next();
			if (setMap.get(key).containsAll(matrix.getTotalAcidsFreqByCol().get(colNr).keySet())) {
				repSets.put(key, new Integer(setMap.get(key).size()));
			}
		}
		smallestSetSize = Collections.min(repSets.values());
		assert smallestSetSize > 0 && smallestSetSize < 25;
		return smallestSetSize;
	}

	/**
	 * Does a very similar thing to SmallestTaylorSetGaps but does not take gaps
	 * into account at all.
	 * 
	 * @param normalize
	 *            is to be set true if data is to be normalized false otherwise
	 * @return array of Taylor scores ,index correspond to the column index,
	 *         indexing starts with 0
	 */
	static double taylorScoreNoGaps(AminoAcidMatrix matrix, int colNr) {

		if (colNr > matrix.numberOfColumns() - 1) {
			throw new IllegalArgumentException(
					"Column number greater than number of columns in teh matrix.");
		}
		Map<String, HashSet<Character>> setMap = ConservationSets.taylorSets();
		Set<String> setMapKeys = setMap.keySet();
		double smallestSetSize = 0.0;
		Map<Character, Integer> acidsMapNoGaps = new HashMap<Character, Integer>(matrix
				.getTotalAcidsFreqByCol().get(colNr));
		if (acidsMapNoGaps.containsKey('-')) {
			acidsMapNoGaps.remove('-');
		}
		Map<String, Integer> repSets = new HashMap<String, Integer>();
		Iterator<String> itr = setMapKeys.iterator();
		while (itr.hasNext()) {
			String key = itr.next();
			if (setMap.get(key).containsAll(acidsMapNoGaps.keySet())) {
				repSets.put(key, new Integer(setMap.get(key).size()));
			}
		}
		smallestSetSize = Collections.min(repSets.values());
		assert smallestSetSize > 0 && smallestSetSize < 25;
		return smallestSetSize;
	}

	/**
	 * Gives a score of one to ten based on whether all the amino acids though
	 * the column maintain or fail to maintain a certain trait.
	 * 
	 * @param normalize
	 *            is to be set true if data is to be normalized false otherwise
	 * @return array of Zvelibil scores ,index correspond to the column index,
	 *         indexing starts with 0
	 */
	// this score does not need to normalized, it is normalizrd by dafault, I
	// leave the normalization for clarity but it is not needed
	static double zvelibilScore(AminoAcidMatrix matrix, int colNr) {

		if (colNr > matrix.numberOfColumns() - 1) {
			throw new IllegalArgumentException(
					"Column number greater than number of columns in teh matrix.");
		}
		double result = 0.0;
		double finalResult = 0.0;
		Map<String, HashSet<Character>> setMap = ConservationSets.zvelibilSets();
		Set<String> keys = setMap.keySet();
		Iterator<String> itr = keys.iterator();
		while (itr.hasNext()) {
			if (setMap.get(itr.next()).containsAll(
					matrix.getTotalAcidsFreqByCol().get(colNr).keySet())) {
				result++;
			}
		}
		assert result >= 0 && result < 11;
		finalResult = result * 0.1;
		return finalResult;
	}

	/**
	 * Calculates Karlin score. done exactly like in jon's code: gaps not
	 * counted, if all but 1 gaps conservation score 0 although it looks to me
	 * that the formula suggests that if one doesn't count gaps than the score
	 * of all but one gaps would be undefined.
	 * 
	 * @param normalize
	 *            is to be set true if data is to be normalized false otherwise
	 * @return array of Karlin scores ,index correspond to the column index,
	 *         indexing starts with 0
	 */
	static double karlinScore(AminoAcidMatrix matrix, int colNr) {

		if (colNr > matrix.numberOfColumns() - 1) {
			throw new IllegalArgumentException(
					"Column number greater than number of columns in the matrix.");
		}
		double finalSum;
		double blosumSum = 0.0;
		char[] curColumn = matrix.getInverseMatrix()[colNr];
		for (int a = 0; a < curColumn.length; a++) {
			if (curColumn[a] != '-') {
				int idxA = ConservationMatrices.getIndex(curColumn[a]);
				int pairAIndex = 24 * idxA + idxA;
				double aSelf = ConservationMatrices.blosum[pairAIndex];
				assert aSelf > 0;
				for (int b = a + 1; b < curColumn.length; b++) {
					if (curColumn[b] != '-') {
						int idxB = ConservationMatrices.getIndex(curColumn[b]);
						int pairABIndex = 24 * idxA + idxB;
						double pairScore = ConservationMatrices.blosum[pairABIndex];
						int pairBIndex = 24 * idxB + idxB;
						double bSelf = ConservationMatrices.blosum[pairBIndex];
						assert bSelf > 0;
						blosumSum = blosumSum + ((pairScore) / (Math.sqrt(aSelf * bSelf)));
					}
				}
			}
		}
		finalSum = blosumSum
				* (2.0 / (matrix.getInverseMatrix()[colNr].length * (matrix.getInverseMatrix()[colNr].length - 1)));
		assert finalSum >= -1 && finalSum <= 1;
		return finalSum;
	}

	/**
	 * Calculates Armon score.
	 * 
	 * @param normalize
	 *            is to be set true if data is to be normalized false otherwise
	 * @return array of Armon scores ,index correspond to the column index,
	 *         indexing starts with 0
	 */
	// creates an array containing all the amino acids and gaps present in the
	// column
	// iterates through that array twice(nested loops), finds all the possible
	// pairs
	// that can be formed by aa present
	// gap is considered the 21 aminoacid
	static double armonScore(AminoAcidMatrix matrix, int colNr) {

		if (colNr > matrix.numberOfColumns() - 1) {
			throw new IllegalArgumentException(
					"Column number greater than number of columns in teh matrix.");
		}
		double scoreSum = 0.0;
		int arrayLength = matrix.getTotalAcidsFreqByCol().get(colNr).keySet().size();
		Character[] acidsPresent = new Character[arrayLength];
		int arrayIndex = 0;
		Set<Character> keys = matrix.getTotalAcidsFreqByCol().get(colNr).keySet();
		Iterator<Character> itr = keys.iterator();
		while (itr.hasNext()) {
			acidsPresent[arrayIndex] = itr.next();
			arrayIndex++;
		}
		for (int a = 0; a < acidsPresent.length; a++) {
			char charA = acidsPresent[a];
			int idxA = ConservationMatrices.getIndex(charA);
			for (int b = a + 1; b < acidsPresent.length; b++) {
				char charB = acidsPresent[b];
				int idxB = ConservationMatrices.getIndex(charB);
				int pairABIndex = 24 * idxA + idxB;
				scoreSum = scoreSum + ConservationMatrices.miyataArmon[pairABIndex];
			}
		}
		return scoreSum;
	}

	/**
	 * Calculates Thompson score. Gaps are accounted for.
	 * 
	 * @param normalize
	 *            is to be set true if data is to be normalized false otherwise
	 * @return array of Thompson scores ,index correspond to the column index,
	 *         indexing starts with 0
	 */
	// amino acids are viewed as points in k dimensional space
	// an average point is calculated
	// score is the distance between the av point and the actual point
	static double thompsonScore(AminoAcidMatrix matrix, int colNr) {

		if (colNr > matrix.numberOfColumns() - 1) {
			throw new IllegalArgumentException(
					"Column number greater than number of columns in teh matrix.");
		}
		double[] sum = null;
		double[] meanPoint = null;
		double distance = 0.0;
		double nonGapsFraction = 0.0;
		char[] alp = Alphabet.alphabetArray();
		assert alp != null && alp.length != 0;
		double result = 0.0;
		double[][] points = new double[matrix.getInverseMatrix()[colNr].length][alp.length];
		char[] curColumn = matrix.getInverseMatrix()[colNr];
		for (int a = 0; a < curColumn.length; a++) {
			int idxA = ConservationMatrices.getIndex(curColumn[a]);
			for (int b = 0; b < alp.length; b++) {
				int idxB = ConservationMatrices.getIndex(alp[b]);
				int pairABIndex = 24 * idxA + idxB;
				points[a][b] = ConservationMatrices.blosum[pairABIndex];
			}
			if (sum == null) {
				sum = points[a];
			} else {
				sum = ConservationAccessory.addPoints(sum, points[a]);
			}
		}
		assert sum != null;
		meanPoint = ConservationAccessory.multPointByScalar(sum, 1.0 / curColumn.length);
		for (int c = 0; c < curColumn.length; c++) {
			distance = distance + ConservationAccessory.pointDistance(points[c], meanPoint);
		}
		if (matrix.getTotalAcidsFreqByCol().get(colNr).keySet().contains('-')) {
			nonGapsFraction = (double) (curColumn.length - matrix.getTotalAcidsFreqByCol().get(
					colNr).get('-'))
					/ (double) curColumn.length;
		} else {
			nonGapsFraction = 1.0;
		}
		result = nonGapsFraction * 1.0 / matrix.getInverseMatrix()[colNr].length * distance;
		return result;
	}

	// causes some math problem because denominator can be 0, that's a formula
	// flaw
	// nothing can be done about it
	/**
	 * Calculates NotArmon score. Valdar's version of Armon score that does not
	 * deal with dividing by zero.
	 * 
	 * @param normalize
	 *            is to be set true if data is to be normalized false otherwise
	 * @return array of NotLancet scores ,index correspond to the column index,
	 *         indexing starts with 0
	 */
	static double notLancetScore(AminoAcidMatrix matrix, int colNr) {

		if (colNr > matrix.numberOfColumns() - 1) {
			throw new IllegalArgumentException(
					"Column number greater than number of columns in teh matrix.");
		}
		double result = 0.0;
		char[] curColumn = matrix.getInverseMatrix()[colNr];
		Map<Character, Integer> acidsInt = matrix.getTotalAcidsFreqByCol().get(colNr);
		Set<Character> keys = acidsInt.keySet();
		Iterator<Character> itr1 = keys.iterator();
		while (itr1.hasNext()) {
			char key1 = itr1.next();
			int idx1 = ConservationMatrices.getIndex(key1);
			Iterator<Character> itr2 = keys.iterator();
			while (itr2.hasNext()) {
				char key2 = itr2.next();
				int idx2 = ConservationMatrices.getIndex(key2);
				int pair12Index = 24 * idx1 + idx2;
				double blosum = ConservationMatrices.blosum[pair12Index];
				result = result
						+ ((((double) acidsInt.get(key1) / (double) curColumn.length
								* (double) acidsInt.get(key2) / curColumn.length)) * blosum);
			}
		}
		return result;
	}

	/**
	 * Calculates Mirny Score.
	 * 
	 * @param normalize
	 *            is to be set true if data is to be normalized false otherwise
	 * @return array of Mirny scores ,index correspond to the column index,
	 *         indexing starts with 0
	 */
	// calculates shannon enthropy but based on the sets, has to calculate
	// number of amino acids tat belong
	// a particular set, stores them in a hashmap
	// reads in a hashmap with the sets needed and creates a hashmap with set
	// names as keys and number of aa belonging to set as value
	static double mirnyScore(AminoAcidMatrix matrix, int colNr) {

		if (colNr > matrix.numberOfColumns() - 1) {
			throw new IllegalArgumentException(
					"Column number greater than number of columns in teh matrix.");
		}
		double mirnySum = 0.0;
		Map<String, HashSet<Character>> mirnySets = ConservationSets.mirnySets();
		assert mirnySets != null && !mirnySets.isEmpty();
		Set<String> mirnyKeys = mirnySets.keySet();
		assert !mirnyKeys.isEmpty();
		Iterator<String> mirnyKeysItr = mirnyKeys.iterator();
		Set<Character> acInKeys = matrix.getTotalAcidsFreqByCol().get(colNr).keySet();
		assert !acInKeys.isEmpty();
		Map<String, Integer> setsFreq = new HashMap<String, Integer>();
		while (mirnyKeysItr.hasNext()) {
			String mirnyKey = mirnyKeysItr.next();
			Iterator<Character> acInKeysItr = acInKeys.iterator();
			while (acInKeysItr.hasNext()) {
				Character acInKey = acInKeysItr.next();
				if (mirnySets.get(mirnyKey).contains(acInKey)) {
					Integer count = setsFreq.get(mirnyKey);
					if (count == null) {
						setsFreq.put(mirnyKey, matrix.getTotalAcidsFreqByCol().get(colNr).get(
								acInKey));
					} else {
						setsFreq.put(mirnyKey, count
								+ matrix.getTotalAcidsFreqByCol().get(colNr).get(acInKey));
					}
				}
			}
		}
		// this assertion will not work if we feed an empty column
		assert !setsFreq.isEmpty();
		// setsFreqKeys = setsFreq.keySet();
		Iterator<String> setsFreqKeysItr = setsFreq.keySet().iterator();
		while (setsFreqKeysItr.hasNext()) {
			String setFreqKey = setsFreqKeysItr.next();
			double pI = (double) setsFreq.get(setFreqKey)
					/ (double) matrix.getInverseMatrix()[colNr].length;
			mirnySum = mirnySum + (pI * Math.log(pI));
		}
		return mirnySum;
	}

	/**
	 * Calculates Williamson score.
	 * 
	 * @param normalize
	 *            is to be set true if data is to be normalized false otherwise
	 * @return array of Williamson scores ,index correspond to the column index,
	 *         indexing starts with 0
	 */
	static double williamsonScore(AminoAcidMatrix matrix, int colNr) {

		if (colNr > matrix.numberOfColumns() - 1) {
			throw new IllegalArgumentException(
					"Column number greater than number of columns in teh matrix.");
		}
		double willSum = 0.0;
		Map<String, HashSet<Character>> willSets = ConservationSets.williamsonSets();
		assert willSets != null && !willSets.isEmpty();
		Set<String> willKeys = willSets.keySet();
		assert !willKeys.isEmpty();
		Iterator<String> willKeysItr = willKeys.iterator();
		Set<Character> acInKeys = matrix.getTotalAcidsFreqByCol().get(colNr).keySet();
		assert !acInKeys.isEmpty();
		Map<String, Integer> setsFreq = new HashMap<String, Integer>();
		while (willKeysItr.hasNext()) {
			String willKey = willKeysItr.next();
			Iterator<Character> acInKeysItr = acInKeys.iterator();
			while (acInKeysItr.hasNext()) {
				Character acInKey = acInKeysItr.next();
				if (willSets.get(willKey).contains(acInKey)) {
					Integer count = setsFreq.get(willKey);
					if (count == null) {
						setsFreq.put(willKey, matrix.getTotalAcidsFreqByCol().get(colNr).get(
								acInKey));
					} else {
						setsFreq.put(willKey, count
								+ matrix.getTotalAcidsFreqByCol().get(colNr).get(acInKey));
					}
				}
			}
		}
		// this assertion will not work if we feed an empty column
		assert !setsFreq.isEmpty();
		Set<String> setsFreqKeys = setsFreq.keySet();
		Iterator<String> setsFreqKeysItr = setsFreqKeys.iterator();
		while (setsFreqKeysItr.hasNext()) {
			String setFreqKey = setsFreqKeysItr.next();
			// FIXME Pi in the logarithm needs to be divided by average pi, do
			// it once u get classes to be nested
			double pI = (double) setsFreq.get(setFreqKey)
					/ (double) matrix.getInverseMatrix()[colNr].length;
			double piAve = (double) matrix.totalAcidsWillSets().get(setFreqKey)
					/ matrix.numberOfRows();
			willSum = willSum + (pI * Math.log(pI / piAve));
		}
		return willSum;
	}

	/**
	 * Calculates Landgraf score.
	 * 
	 * @param normalize
	 *            is to be set true if data is to be normalized false otherwise
	 * @return array of Landgraf scores ,index correspond to the column index,
	 *         indexing starts with 0
	 */
	static double landgrafScore(AminoAcidMatrix matrix, int colNr) {

		if (colNr > matrix.numberOfColumns() - 1) {
			throw new IllegalArgumentException(
					"Column number greater than number of columns in teh matrix.");
		}
		double sum = 0.0;
		double result = 0.0;
		char[] curColumn = matrix.getInverseMatrix()[colNr];
		double[] voronoiWeights = matrix.getVoronoiWeights(1000);
		for (int a = 0; a < curColumn.length; a++) {
			double voronoiA = voronoiWeights[a];
			for (int b = a + 1; b < curColumn.length; b++) {
				double disIJ = ConservationMatrices.dissimilarity(curColumn[a], curColumn[b]);
				double disJI = ConservationMatrices.dissimilarity(curColumn[b], curColumn[a]);
				sum = sum + voronoiA * disIJ + voronoiWeights[b] * disJI;
			}
		}
		result = sum / curColumn.length;
		return result;
	}

	/**
	 * Calculates sander score.
	 * 
	 * @param normalize
	 *            is to be set true if data is to be normalized false otherwise
	 * @return array of Sander scores ,index correspond to the column index,
	 *         indexing starts with 0
	 */
	static double sanderScore(AminoAcidMatrix matrix, int colNr) {

		if (colNr > matrix.numberOfColumns() - 1) {
			throw new IllegalArgumentException(
					"Column number greater than number of columns in teh matrix.");
		}
		double sum = 0.0;
		double moderator = 0.0;
		double result = 0.0;
		// char[][] inverseMatrix = matrix.getInverseMatrix();
		char[] curColumn = matrix.getInverseMatrix()[colNr];
		double[][] percent_identity = matrix.getPercentIdentity();
		for (int a = 0; a < curColumn.length; a++) {
			int aIdx = 24 * ConservationMatrices.getIndex(curColumn[a]);
			for (int b = a + 1; b < curColumn.length; b++) {
				double identity = 1 - percent_identity[a][b];
				int pairIndex = aIdx + ConservationMatrices.getIndex(curColumn[b]);
				sum += identity * ConservationMatrices.pam250[pairIndex];
				moderator += identity;
			}
		}
		result = sum * moderator;
		return result;
	}

	/**
	 * Calculates Valdar score.
	 * 
	 * @param normalize
	 *            is to be set true if data is to be normalized false otherwise
	 * @return array of Sander scores ,index correspond to the column index,
	 *         indexing starts with 0
	 */
	static double valdarScore(AminoAcidMatrix matrix, int colNr) {

		if (colNr > matrix.numberOfColumns() - 1) {
			throw new IllegalArgumentException(
					"Column number greater than number of columns in teh matrix.");
		}
		double sum = 0.0;
		double moderator = 0.0;
		double result = 0.0;
		char[] curColumn = matrix.getInverseMatrix()[colNr];
		double[] vingronArgosWeights = matrix.vingronArgosWeights();
		for (int a = 0; a < curColumn.length; a++) {
			int aIdx = 24 * ConservationMatrices.getIndex(curColumn[a]);
			for (int b = a + 1; b < curColumn.length; b++) {
				double mod = vingronArgosWeights[a] * vingronArgosWeights[b];
				int pairIndex = aIdx + ConservationMatrices.getIndex(curColumn[b]);
				sum += mod * ConservationMatrices.pet91[pairIndex];
				moderator += mod;
			}
		}
		result = sum * moderator;
		return result;
	}
}
