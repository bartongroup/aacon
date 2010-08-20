package compbio.conservation;

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

/**
 * TODO remove? This class does not seems to be used!
 * 
 * @author pvtroshin
 */
public class _ConservationScores {

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
	private final AminoAcidMatrix matrix;

	/**
	 * Constructor
	 * 
	 * @param matrix
	 *            AminoAcidMatrix based on the alignment
	 */
	public _ConservationScores(AminoAcidMatrix matrix) {

		this.matrix = matrix;
		// scores = new EnumMap<Method, double[]>(Method.class);
	}

	/**
	 * Checks if all but one residues in the column are gaps.
	 * 
	 * @return true if all but one residues are gaps, false if not
	 */
	private boolean allButOneGaps(int columnNr) {

		if (matrix.getTotalAcidsFreqByCol().get(columnNr).containsKey('-')
				&& matrix.getTotalAcidsFreqByCol().get(columnNr).get('-') == matrix
						.getInverseMatrix()[columnNr].length - 1) {
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
	private boolean oneResidueTypeNoGaps(int columnNr) {

		if (matrix.getTotalAcidsFreqByCol().get(columnNr).size() == 1
				&& matrix.getTotalAcidsFreqByCol().get(columnNr).containsKey('-') == false) {
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
	private boolean containsGaps(int columnNr) {

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
	private int numberOfAcidsWithGap(int columnNr) {

		return matrix.getTotalAcidsFreqByCol().get(columnNr).size();
	}

	/**
	 * Counts the number of different amino acids in the column. Gap not counted
	 * as a 21 amino acid.
	 * 
	 * @return the number of different amino acids
	 */
	private int numberOfAcidsNoGap(int columnNr) {

		if (this.containsGaps(columnNr) == true) {
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
	private int mostCommonNumber(int columnNr) {

		int max = 0;
		Set<Character> keys = matrix.getTotalAcidsFreqByCol().get(columnNr).keySet();
		Iterator<Character> itr = keys.iterator();
		while (itr.hasNext()) {
			Character key = itr.next();
			if (key != '-' && matrix.getTotalAcidsFreqByCol().get(columnNr).get(key) > max) {
				max = matrix.getTotalAcidsFreqByCol().get(columnNr).get(key);
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
	double[] kabatScore(boolean normalize) {

		double[] result = new double[matrix.getInverseMatrix().length];
		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
			assert mostCommonNumber(i) > 0
					&& mostCommonNumber(i) < matrix.getInverseMatrix()[i].length + 1;
			result[i] = matrix.getInverseMatrix()[i].length * (double) numberOfAcidsNoGap(i)
					/ mostCommonNumber(i);
		}
		// scores.put(Method.kabatScore, result);
		if (normalize) {
			double[] normalized = ConservationAccessory.inversedNormalize01(result);
			return normalized;
		} else {
			return result;
		}
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
	double[] joresScore(boolean normalize) {

		double[] result = new double[matrix.getInverseMatrix().length];
		Map<Character, Integer> acidsIntMap = null;
		Map<Character, Integer> acidsIntMapCopy = null;
		Set<Character> keys = null;
		Iterator<Integer> itr = null;
		Iterator<Character> itr2 = null;
		Iterator<Integer> itr3 = null;
		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
			// special case #1 one residue type only
			// special case #2 all but one are gaps
			if (oneResidueTypeNoGaps(i) == true || allButOneGaps(i) == true) {
				if (oneResidueTypeNoGaps(i) == true) {
					result[i] = 1.0;
				} else {
					result[i] = matrix.getInverseMatrix().length
							* (matrix.getInverseMatrix().length - 1) / 2;
				}
			} else {
				int samePairs = 0;
				int differentPairs = 0;
				acidsIntMap = matrix.getTotalAcidsFreqByCol().get(i);
				acidsIntMapCopy = new HashMap<Character, Integer>(acidsIntMap);
				keys = acidsIntMapCopy.keySet();
				keys.remove('-');
				int types = keys.size();
				differentPairs = types * (types - 1) / 2;
				itr = acidsIntMapCopy.values().iterator();
				while (itr.hasNext()) {
					if (itr.next() > 1) {
						samePairs++;
					}
				}
				int totalPairs = samePairs + differentPairs;
				itr2 = acidsIntMapCopy.keySet().iterator();
				int max1 = 0;
				Character maxKey = null;
				while (itr2.hasNext()) {
					Character key = itr2.next();
					if (acidsIntMapCopy.get(key) > max1) {
						maxKey = key;
						max1 = acidsIntMapCopy.get(key);
					}
				}
				acidsIntMapCopy.remove(maxKey);
				itr3 = acidsIntMapCopy.values().iterator();
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
				result[i] = ((double) totalPairs / (double) mostFreqNr)
						* ((matrix.getInverseMatrix()[i].length)
								* (matrix.getInverseMatrix()[i].length - 1) / 2);
			}
		}
		// scores.put(Method.joresScore, result);
		if (normalize) {
			double[] normalized = ConservationAccessory.inversedNormalize01(result);
			return normalized;
		} else {
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
	double[] schneiderScore(boolean normalize) {

		double[] result = new double[matrix.getInverseMatrix().length];
		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
			double normal = 1.0 / Math.log(20.0);
			result[i] = ShannonEnthropy.ShannonLn(matrix.getTotalAcidsFreqByCol().get(i), matrix
					.getInverseMatrix()[i].length)
					* normal;
			assert result[i] >= 0 && result[i] <= 1;
		}
		// scores.put(Method.schneiderScore, result);
		if (normalize) {
			double[] normalized = ConservationAccessory.inversedNormalize01(result);
			return normalized;
		} else {
			return result;
		}
	}

	/**
	 * Calculates Shenkin score for the alignment.
	 * 
	 * @param normalize
	 *            is to be set true if data is to be normalized false otherwise
	 * @return array of Shenkin scores ,index correspond to the column index,
	 *         indexing starts with 0
	 */
	double[] shenkinScore(boolean normalize) {

		double[] result = new double[matrix.getInverseMatrix().length];
		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
			result[i] = Math.pow(2.0, ShannonEnthropy.ShannonLog2(matrix.getTotalAcidsFreqByCol()
					.get(i), matrix.getInverseMatrix()[i].length)) * 6.0;
			assert result[i] >= 6 && result[i] <= 120;
		}
		// scores.put(Method.shenkinScore, result);
		if (normalize) {
			double[] normalized = ConservationAccessory.inversedNormalize01(result);
			return normalized;
		} else {
			return result;
		}
	}

	/**
	 * Calculates Gerstein score.
	 * 
	 * @param normalize
	 *            is to be set true if data is to be normalized false otherwise
	 * @return array of Gerstein scores ,index correspond to the column index,
	 *         indexing starts with 0
	 */
	double[] gersteinScore(boolean normalize) {

		double[] result = new double[matrix.getInverseMatrix().length];
		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
			result[i] = -(ShannonEnthropy.ShannonLn(matrix.totalAcidsFrequency(), matrix
					.numberOfColumns()
					* matrix.numberOfRows()))
					- (-ShannonEnthropy.ShannonLn(matrix.getTotalAcidsFreqByCol().get(i), matrix
							.getInverseMatrix()[i].length));
		}
		// scores.put(Method.gersteinScore, result);
		if (normalize) {
			double[] normalized = ConservationAccessory.inversedNormalize01(result);
			return normalized;
		} else {
			return result;
		}
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
	double[] SmallestTaylorSetGaps(boolean normalize) {

		Map<String, HashSet<Character>> setMap = ConservationSets.taylorSets();
		double[] smallestSetSize = new double[matrix.getInverseMatrix().length];
		Map<String, Integer> repSets = null;
		Set<String> setMapKeys = null;
		Iterator<String> itr = null;
		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
			repSets = new HashMap<String, Integer>();
			setMapKeys = setMap.keySet();
			itr = setMapKeys.iterator();
			while (itr.hasNext()) {
				String key = itr.next();
				if (setMap.get(key).containsAll(matrix.getTotalAcidsFreqByCol().get(i).keySet())) {
					repSets.put(key, new Integer(setMap.get(key).size()));
				}
			}
			smallestSetSize[i] = Collections.min(repSets.values());
			assert smallestSetSize[i] > 0 && smallestSetSize[i] < 25;
		}
		// scores.put(Method.SmallestTaylorSetGaps, smallestSetSize);
		if (normalize) {
			double[] normalized = ConservationAccessory.inversedNormalize01(smallestSetSize);
			return normalized;
		} else {
			return smallestSetSize;
		}
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
	double[] SmallestTaylorSetNoGaps(boolean normalize) {

		Map<String, HashSet<Character>> setMap = ConservationSets.taylorSets();
		Set<String> setMapKeys = setMap.keySet();
		double[] smallestSetSize = new double[matrix.getInverseMatrix().length];
		Map<Character, Integer> acidsMapNoGaps = null;
		Map<String, Integer> repSets = null;
		Iterator<String> itr = null;
		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
			acidsMapNoGaps = new HashMap<Character, Integer>(matrix.getTotalAcidsFreqByCol().get(i));
			if (acidsMapNoGaps.containsKey('-')) {
				acidsMapNoGaps.remove('-');
			}
			repSets = new HashMap<String, Integer>();
			itr = setMapKeys.iterator();
			while (itr.hasNext()) {
				String key = itr.next();
				if (setMap.get(key).containsAll(acidsMapNoGaps.keySet())) {
					repSets.put(key, new Integer(setMap.get(key).size()));
				}
			}
			smallestSetSize[i] = Collections.min(repSets.values());
			assert smallestSetSize[i] > 0 && smallestSetSize[i] < 25;
		}
		// scores.put(Method.SmallestTaylorSetNoGaps, smallestSetSize);
		if (normalize) {
			double[] normalized = ConservationAccessory.inversedNormalize01(smallestSetSize);
			return normalized;
		} else {
			return smallestSetSize;
		}
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
	double[] zvelibilScore(boolean normalize) {

		double[] result = new double[matrix.getInverseMatrix().length];
		double[] finalResult = new double[matrix.getInverseMatrix().length];
		Map<String, HashSet<Character>> setMap = ConservationSets.zvelibilSets();
		Set<String> keys = setMap.keySet();
		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
			Iterator<String> itr = keys.iterator();
			while (itr.hasNext()) {
				if (setMap.get(itr.next()).containsAll(
						matrix.getTotalAcidsFreqByCol().get(i).keySet())) {
					result[i]++;
				}
			}
			assert result[i] >= 0 && result[i] < 11;
			finalResult[i] = result[i] * 0.1;
		}
		// scores.put(Method.zvelibilScore, result);
		if (normalize) {
			double[] normalized = ConservationAccessory.normalize01(finalResult);
			return normalized;
		} else {
			return finalResult;
		}
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
	double[] karlinScore(boolean normalize) {

		double[] finalSum = new double[matrix.getInverseMatrix().length];
		double[] blosumSum = new double[matrix.getInverseMatrix().length];
		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
			for (int a = 0; a < matrix.getInverseMatrix()[i].length; a++) {
				if (matrix.getInverseMatrix()[i][a] != '-') {
					for (int b = a + 1; b < matrix.getInverseMatrix()[i].length; b++) {
						if (matrix.getInverseMatrix()[i][b] != '-') {
							double pairScore = ConservationMatrices.BlosumPair(matrix
									.getInverseMatrix()[i][a], matrix.getInverseMatrix()[i][b]);
							double aSelf = ConservationMatrices.BlosumPair(matrix
									.getInverseMatrix()[i][a], matrix.getInverseMatrix()[i][a]);
							assert aSelf > 0;
							double bSelf = ConservationMatrices.BlosumPair(matrix
									.getInverseMatrix()[i][b], matrix.getInverseMatrix()[i][b]);
							assert bSelf > 0;
							blosumSum[i] = blosumSum[i]
									+ ((pairScore) / (Math.sqrt(aSelf * bSelf)));
						}
					}
				}
			}
			finalSum[i] = blosumSum[i]
					* (2.0 / (matrix.getInverseMatrix()[i].length * (matrix.getInverseMatrix()[i].length - 1)));
			assert finalSum[i] >= -1 && finalSum[i] <= 1;
		}
		// scores.put(Method.karlinScore, finalSum);
		if (normalize) {
			double[] normalized = ConservationAccessory.normalize01(finalSum);
			return normalized;
		} else {
			return finalSum;
		}
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
	double[] armonScore(boolean normalize) {

		double[] scoreSum = new double[matrix.getInverseMatrix().length];
		Set<Character> keys = null;
		Iterator<Character> itr = null;
		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
			int arrayLength = matrix.getTotalAcidsFreqByCol().get(i).keySet().size();
			Character[] acidsPresent = new Character[arrayLength];
			int arrayIndex = 0;
			keys = matrix.getTotalAcidsFreqByCol().get(i).keySet();
			itr = keys.iterator();
			while (itr.hasNext()) {
				acidsPresent[arrayIndex] = itr.next();
				arrayIndex++;
			}
			for (int a = 0; a < acidsPresent.length; a++) {
				char charA = acidsPresent[a];
				for (int b = a + 1; b < acidsPresent.length; b++) {
					char charB = acidsPresent[b];
					scoreSum[i] = scoreSum[i] + ConservationMatrices.miyataArmonPair(charA, charB);
				}
			}
		}
		// scores.put(Method.armonscore, scoreSum);
		if (normalize) {
			double[] normalized = ConservationAccessory.inversedNormalize01(scoreSum);
			return normalized;
		} else {
			return scoreSum;
		}
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
	double[] thompsonScore(boolean normalize) {

		double[] sum = null;
		double[] meanPoint = null;
		double distance = 0.0;
		double nonGapsFraction = 0.0;
		char[] alp = Alphabet.alphabetArray();
		assert alp != null && alp.length != 0;
		double[] result = new double[matrix.getInverseMatrix().length];
		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
			double[][] points = new double[matrix.getInverseMatrix()[i].length][alp.length];
			for (int a = 0; a < matrix.getInverseMatrix()[i].length; a++) {
				for (int b = 0; b < alp.length; b++) {
					points[a][b] = ConservationMatrices.BlosumPair(matrix.getInverseMatrix()[i][a],
							alp[b]);
				}
				if (sum == null) {
					sum = points[a];
				} else {
					sum = ConservationAccessory.addPoints(sum, points[a]);
				}
			}
			assert sum != null;
			meanPoint = ConservationAccessory.multPointByScalar(sum, 1.0 / matrix
					.getInverseMatrix()[i].length);
			for (int c = 0; c < matrix.getInverseMatrix()[i].length; c++) {
				distance = distance + ConservationAccessory.pointDistance(points[c], meanPoint);
			}
			if (matrix.getTotalAcidsFreqByCol().get(i).keySet().contains('-')) {
				nonGapsFraction = (double) (matrix.getInverseMatrix()[i].length - matrix
						.getTotalAcidsFreqByCol().get(i).get('-'))
						/ (double) matrix.getInverseMatrix()[i].length;
			} else {
				nonGapsFraction = 1.0;
			}
			result[i] = nonGapsFraction * 1.0 / matrix.getInverseMatrix()[i].length * distance;
			sum = null;
			meanPoint = null;
			distance = 0.0;
			nonGapsFraction = 0.0;
		}
		// scores.put(Method.thompsonScore, result);
		if (normalize) {
			double[] normalized = ConservationAccessory.inversedNormalize01(result);
			return normalized;
		} else {
			return result;
		}
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
	double[] notLancetScore(boolean normalize) {

		double[] result = new double[matrix.getInverseMatrix().length];
		Set<Character> keys = null;
		Iterator<Character> itr1 = null;
		Iterator<Character> itr2 = null;
		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
			keys = matrix.getTotalAcidsFreqByCol().get(i).keySet();
			itr1 = keys.iterator();
			while (itr1.hasNext()) {
				char key1 = itr1.next();
				itr2 = keys.iterator();
				while (itr2.hasNext()) {
					char key2 = itr2.next();
					double blosum = ConservationMatrices.BlosumPair(key1, key2);
					result[i] = result[i]
							+ ((((double) matrix.getTotalAcidsFreqByCol().get(i).get(key1)
									/ (double) matrix.getInverseMatrix()[i].length
									* (double) matrix.getTotalAcidsFreqByCol().get(i).get(key2) / matrix
									.getInverseMatrix()[i].length)) * blosum);
				}
			}
		}
		// scores.put(Method.notLancetScore, result);
		if (normalize) {
			double[] normalized = ConservationAccessory.normalize01(result);
			return normalized;
		} else {
			return result;
		}
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
	double[] mirnyScore(boolean normalize) {

		double[] mirnySum = new double[matrix.getInverseMatrix().length];
		Map<String, HashSet<Character>> mirnySets = ConservationSets.mirnySets();
		assert mirnySets != null && !mirnySets.isEmpty();
		Set<String> mirnyKeys = mirnySets.keySet();
		assert !mirnyKeys.isEmpty();
		// Iterator<String> mirnyKeysItr = null;
		// Set<Character> acInKeys = null;
		// Map<String,Integer> setsFreq = null;
		// Iterator<Character> acInKeysItr = null;
		// Set<String> setsFreqKeys = null;
		// Iterator<String> setsFreqKeysItr = null;
		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
			Iterator<String> mirnyKeysItr = mirnyKeys.iterator();
			Set<Character> acInKeys = matrix.getTotalAcidsFreqByCol().get(i).keySet();
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
							setsFreq.put(mirnyKey, matrix.getTotalAcidsFreqByCol().get(i).get(
									acInKey));
						} else {
							setsFreq.put(mirnyKey, count
									+ matrix.getTotalAcidsFreqByCol().get(i).get(acInKey));
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
				System.out.println(i + " " + "Sets tester1 " + " key " + setFreqKey + " value "
						+ setsFreq.get(setFreqKey));
				double pI = (double) setsFreq.get(setFreqKey)
						/ (double) matrix.getInverseMatrix()[i].length;
				mirnySum[i] = mirnySum[i] + (pI * Math.log(pI));
			}
		}
		// scores.put(Method.mirnyScore, mirnySum);
		if (normalize) {
			double[] normalized = ConservationAccessory.normalize01(mirnySum);
			return normalized;
		} else {
			return mirnySum;
		}
	}

	// double[] mirnyScore2(boolean normalize) {
	// double[] mirnySum = new double[matrix.getInverseMatrix().length];
	// for(int i = 0; i < matrix.getInverseMatrix().length; i++) {
	// mirnySum[i] = Column.mirnyScore(matrix.getTotalAcidsFreqByCol().get(i),
	// matrix.getInverseMatrix()[i]);
	// }
	// if(normalize) {
	// double[] normalized = ConservationAccessory.normalize01(mirnySum);
	// return normalized;
	// }
	// else {
	// return mirnySum;
	// }
	// }
	/**
	 * Calculates Williamson score.
	 * 
	 * @param normalize
	 *            is to be set true if data is to be normalized false otherwise
	 * @return array of Williamson scores ,index correspond to the column index,
	 *         indexing starts with 0
	 */
	double[] williamsonScore(boolean normalize) {

		double[] willSum = new double[matrix.getInverseMatrix().length];
		Map<String, HashSet<Character>> willSets = ConservationSets.williamsonSets();
		assert willSets != null && !willSets.isEmpty();
		Set<String> willKeys = willSets.keySet();
		assert !willKeys.isEmpty();
		Iterator<String> willKeysItr = null;
		Set<Character> acInKeys = null;
		Map<String, Integer> setsFreq = null;
		Iterator<Character> acInKeysItr = null;
		Set<String> setsFreqKeys = null;
		Iterator<String> setsFreqKeysItr = null;
		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
			willKeysItr = willKeys.iterator();
			acInKeys = matrix.getTotalAcidsFreqByCol().get(i).keySet();
			assert !acInKeys.isEmpty();
			setsFreq = new HashMap<String, Integer>();
			while (willKeysItr.hasNext()) {
				String willKey = willKeysItr.next();
				acInKeysItr = acInKeys.iterator();
				while (acInKeysItr.hasNext()) {
					Character acInKey = acInKeysItr.next();
					if (willSets.get(willKey).contains(acInKey)) {
						Integer count = setsFreq.get(willKey);
						if (count == null) {
							setsFreq.put(willKey, matrix.getTotalAcidsFreqByCol().get(i).get(
									acInKey));
						} else {
							setsFreq.put(willKey, count
									+ matrix.getTotalAcidsFreqByCol().get(i).get(acInKey));
						}
					}
				}
			}
			// this assertion will not work if we feed an empty column
			assert !setsFreq.isEmpty();
			setsFreqKeys = setsFreq.keySet();
			setsFreqKeysItr = setsFreqKeys.iterator();
			while (setsFreqKeysItr.hasNext()) {
				String setFreqKey = setsFreqKeysItr.next();
				System.out.println("Sets tester " + " key " + setFreqKey + " value "
						+ setsFreq.get(setFreqKey));
				// FIXME Pi in the logarithm needs to be divided by average pi,
				// do it once u get classes to be nested
				double pI = (double) setsFreq.get(setFreqKey)
						/ (double) matrix.getInverseMatrix()[i].length;
				double piAve = (double) matrix.totalAcidsWillSets().get(setFreqKey)
						/ matrix.numberOfRows();
				willSum[i] = willSum[i] + (pI * Math.log(pI / piAve));
			}
		}
		// scores.put(Method.williamsonScore, willSum);
		if (normalize) {
			double[] normalized = ConservationAccessory.normalize01(willSum);
			return normalized;
		} else {
			return willSum;
		}
	}

	/**
	 * Calculates Landgraf score.
	 * 
	 * @param normalize
	 *            is to be set true if data is to be normalized false otherwise
	 * @return array of Landgraf scores ,index correspond to the column index,
	 *         indexing starts with 0
	 */
	double[] landgrafScore(boolean normalize) {

		double sum = 0;
		double[] result = new double[matrix.getInverseMatrix().length];
		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
			for (int a = 0; a < matrix.getInverseMatrix()[i].length; a++) {
				for (int b = a + 1; b < matrix.getInverseMatrix()[i].length; b++) {
					double disIJ = ConservationMatrices.dissimilarity(
							matrix.getInverseMatrix()[i][a], matrix.getInverseMatrix()[i][b]);
					double disJI = ConservationMatrices.dissimilarity(
							matrix.getInverseMatrix()[i][b], matrix.getInverseMatrix()[i][a]);
					sum = sum + (matrix.getVoronoiWeights(1000)[a] * disIJ)
							+ (matrix.getVoronoiWeights(1000)[b] * disJI);
				}
			}
			result[i] = sum / matrix.getInverseMatrix()[i].length;
			sum = 0;
		}
		// scores.put(Method.landgrafScore, result);
		if (normalize) {
			double[] normalized = ConservationAccessory.inversedNormalize01(result);
			return normalized;
		} else {
			return result;
		}
	}

	/**
	 * Calculates sander score.
	 * 
	 * @param normalize
	 *            is to be set true if data is to be normalized false otherwise
	 * @return array of Sander scores ,index correspond to the column index,
	 *         indexing starts with 0
	 */
	double[] sanderScore(boolean normalize) {

		double sum = 0.0;
		double moderator = 0.0;
		double[] result = new double[matrix.getInverseMatrix().length];
		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
			for (int a = 0; a < matrix.getInverseMatrix()[i].length; a++) {
				for (int b = a + 1; b < matrix.getInverseMatrix()[i].length; b++) {
					sum = sum
							+ (1 - matrix.getPercentIdentity()[a][b])
							* ConservationMatrices.pam250Pair(matrix.getInverseMatrix()[i][a],
									matrix.getInverseMatrix()[i][b]);
				}
			}
			for (int a = 0; a < matrix.getInverseMatrix()[i].length; a++) {
				for (int b = a + 1; b < matrix.getInverseMatrix()[i].length; b++) {
					moderator = moderator + (1 - matrix.getPercentIdentity()[a][b]);
				}
			}
			result[i] = sum * moderator;
			sum = 0.0;
			moderator = 0.0;
		}
		// scores.put(Method.sanderScore, result);
		if (normalize) {
			double[] normalized = ConservationAccessory.normalize01(result);
			return normalized;
		} else {
			return result;
		}
	}

	/**
	 * Calculates Valdar score.
	 * 
	 * @param normalize
	 *            is to be set true if data is to be normalized false otherwise
	 * @return array of Sander scores ,index correspond to the column index,
	 *         indexing starts with 0
	 */
	double[] valdarScore(boolean normalize) {

		double sum = 0.0;
		double moderator = 0.0;
		double[] result = new double[matrix.getInverseMatrix().length];
		for (int i = 0; i < matrix.getInverseMatrix().length; i++) {
			for (int a = 0; a < matrix.getInverseMatrix()[i].length; a++) {
				for (int b = a + 1; b < matrix.getInverseMatrix()[i].length; b++) {
					sum = sum
							+ matrix.vingronArgosWeights()[a]
							* matrix.vingronArgosWeights()[b]
							* ConservationMatrices.pet91Pair(matrix.getInverseMatrix()[i][a],
									matrix.getInverseMatrix()[i][b]);
				}
			}
			for (int a = 0; a < matrix.getInverseMatrix()[i].length; a++) {
				for (int b = a + 1; b < matrix.getInverseMatrix()[i].length; b++) {
					moderator = moderator + matrix.vingronArgosWeights()[a]
							* matrix.vingronArgosWeights()[b];
				}
			}
			result[i] = sum * moderator;
			sum = 0.0;
			moderator = 0.0;
		}
		// scores.put(Method.valdarScore, result);
		if (normalize) {
			double[] normalized = ConservationAccessory.normalize01(result);
			return normalized;
		} else {
			return result;
		}
	}
}
