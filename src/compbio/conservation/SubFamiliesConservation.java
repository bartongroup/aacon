package compbio.conservation;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.EnumMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import compbio.common.SequencesNotEquallyLongException;
import compbio.data.sequence.FastaSequence;
import compbio.data.sequence.SequenceUtil;

/**
 * Calculates conservation in subfamilies provided. Kabat and Jores methods can
 * not be used here. Otherwise an exception may be thrown. Better tresholds need
 * to be set.
 * 
 * @author agolicz
 */
public class SubFamiliesConservation {

	final List<char[][]> subGroups;
	final AminoAcidMatrix alignment;

	SubFamiliesConservation(AminoAcidMatrix alignment, int[][] groups) {

		this.alignment = alignment;
		subGroups = alignment.splitAlignment(groups);
	}

	ConservationStatus[][] groupsConsStat;
	ConservationStatus[][][] groupPairsConsStat;
	String[][] groupsProperties;
	String[][][] groupPairsProperties;
	private static Map<Method, Double> treshholds = new EnumMap<Method, Double>(Method.class);
	static {
		treshholds.put(Method.ARMON, 0.75);
		treshholds.put(Method.GERSTEIN, 0.75);
		treshholds.put(Method.JORES, 0.75);
		treshholds.put(Method.KABAT, 0.75);
		treshholds.put(Method.KARLIN, 0.75);
		treshholds.put(Method.LANDGRAF, 0.75);
		treshholds.put(Method.MIRNY, 0.75);
		treshholds.put(Method.NOT_LANCET, 0.75);
		treshholds.put(Method.SANDER, 0.75);
		treshholds.put(Method.SCHNEIDER, 0.75);
		treshholds.put(Method.SHENKIN, 0.55);
		treshholds.put(Method.TAYLOR_GAPS, 0.75);
		treshholds.put(Method.TAYLOR_NO_GAPS, 0.75);
		treshholds.put(Method.THOMPSON, 0.75);
		treshholds.put(Method.VALDAR, 0.75);
		treshholds.put(Method.WILLIAMSON, 0.75);
		treshholds.put(Method.ZVELIBIL, 0.70);
	}

	static double getTreshhold(Method method) {

		double tresh = treshholds.get(method);
		return tresh;
	}

	/**
	 * Merges two subfamilies into one. A result is an alignment in which the
	 * sequences are put in such an order as if the two alignments were put
	 * together by placing the second beneth the first one.
	 * 
	 * @param sub1
	 *            2D array with first alignment
	 * @param sub2
	 *            2D array with second alignment
	 * @return matrix representing merged alignments
	 */
	static char[][] mergeSubFamilies(char[][] sub1, char[][] sub2) {

		for (int i = 0; i < sub1.length; i++) {
			if (sub1[0].length != sub1[i].length) {
				throw new SequencesNotEquallyLongException("Sequence nr: " + i
						+ " is of different length than teh first sequence.");
			}
		}
		for (int i = 0; i < sub2.length; i++) {
			if (sub2[0].length != sub2[i].length) {
				throw new SequencesNotEquallyLongException("Sequence nr: " + i
						+ " is of different length than teh first sequence.");
			}
		}
		if (sub1[0].length != sub2[0].length) {
			throw new SequencesNotEquallyLongException(
					"Sequences in the provided subfamilies are not of equal length.");
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
	 * 
	 * @param method
	 *            method name
	 * @param subGroups
	 *            list of subfamily alignments
	 * @param normalize
	 *            if result is to be normalized set true
	 * @return
	 */
	double[][] subgrupsConservation(Method method, boolean normalize) {

		if (method == null) {
			throw new IllegalArgumentException("Method must not be null.");
		}
		if (subGroups == null) {
			throw new IllegalArgumentException("List of subgroups/subclasses must not be null.");
		}
		double[][] conservation = new double[subGroups.size()][subGroups.get(0)[0].length];
		groupsProperties = new String[subGroups.size()][subGroups.get(0)[0].length];
		ConservationScores2 scores = null;
		AminoAcidMatrix sub = null;
		String[] names = null;
		for (int i = 0; i < subGroups.size(); i++) {
			sub = new AminoAcidMatrix(subGroups.get(i), names, null);
			scores = new ConservationScores2(sub);
			conservation[i] = scores.calculateScore(method, normalize);
			groupsProperties[i] = zvelibilPropertiesAlignment(sub);
		}
		return conservation;
	}

	ColumnInfo[][] subgrupsConservation3(int width, double gapTreshold,
			SMERFSColumnScore scoreType, boolean normalize) {

		// if (method == null) {
		// throw new IllegalArgumentException("Method must not be null.");
		// }
		if (subGroups == null) {
			throw new IllegalArgumentException("List of subgroups/subclasses must not be null.");
		}
		ColumnInfo[][] columns = new ColumnInfo[this.alignment.numberOfColumns()][subGroups.size()];
		// double[][] conservation = new
		// double[subGroups.size()][subGroups.get(0)[0].length];
		// groupsProperties = new
		// String[subGroups.size()][subGroups.get(0)[0].length];
		ConservationScores2 scores = null;
		AminoAcidMatrix sub = null;
		String[] names = null;
		for (int i = 0; i < subGroups.size(); i++) {
			sub = new AminoAcidMatrix(subGroups.get(i), names, null);
			scores = new ConservationScores2(sub);
			double[] conservation = ConservationClient.getSMERFS(alignment, width, scoreType,
					gapTreshold, normalize);
			String[] prop = zvelibilPropertiesAlignment(sub);
			for (int j = 0; j < columns.length; j++) {
				String group = i + "";
				columns[j][i] = new ColumnInfo(group, conservation[j], prop[j]);
			}
		}
		return columns;
	}

	/**
	 * Calculates conservation in subfamilies
	 * 
	 * @param method
	 * @param normalize
	 * @return
	 */
	ColumnInfo[][] subgrupsConservation2(Method method, boolean normalize) {

		if (method == null) {
			throw new IllegalArgumentException("Method must not be null.");
		}
		if (subGroups == null) {
			throw new IllegalArgumentException("List of subgroups/subclasses must not be null.");
		}
		ColumnInfo[][] columns = new ColumnInfo[this.alignment.numberOfColumns()][subGroups.size()];
		// double[][] conservation = new
		// double[subGroups.size()][subGroups.get(0)[0].length];
		// groupsProperties = new
		// String[subGroups.size()][subGroups.get(0)[0].length];
		ConservationScores2 scores = null;
		AminoAcidMatrix sub = null;
		String[] names = null;
		for (int i = 0; i < subGroups.size(); i++) {
			sub = new AminoAcidMatrix(subGroups.get(i), names, null);
			scores = new ConservationScores2(sub);
			double[] conservation = scores.calculateScore(method, normalize);
			String[] prop = zvelibilPropertiesAlignment(sub);
			for (int j = 0; j < columns.length; j++) {
				String group = i + "";
				columns[j][i] = new ColumnInfo(group, conservation[j], prop[j]);
			}
		}
		return columns;
	}

	/**
	 * Matrix that stores newly formed subfamily pairs. It also contains info
	 * about which subfamilies did the pairs come form. Example : subPair[x][y]
	 * x is the number of the first sequence in the pair(numbering starts form
	 * 0). the number of the second sequence is x + 1 + y
	 * 
	 * @param subGroups
	 * @return
	 */
	char[][][][] subFamilyPairs() {

		if (subGroups == null) {
			throw new IllegalArgumentException("SubGroups list must not be null");
		}
		char[][][][] pairs = new char[subGroups.size()][][][];
		int startLength = subGroups.size() - 1;
		for (int i = 0; i < subGroups.size(); i++) {
			pairs = new char[i][startLength][][];
			int index = 0;
			for (int j = i + 1; j < this.subGroups.size(); j++) {
				pairs[i][index] = mergeSubFamilies(subGroups.get(i), subGroups.get(j));
				index++;
			}
			startLength--;
		}
		return pairs;
	}

	/**
	 * Matrix that stores conservation scores for newly formed subfamily pairs.
	 * It also contains info about which subfamilies did the pairs come form.
	 * Example : subPair[x][y] x is the number of the first sequence in the
	 * pair(numbering starts form 0). the number of the second sequence is x + 1
	 * + y
	 * 
	 * @param subGroups
	 * @return
	 */
	double[][][] subFamilyPairsConservation(Method method, boolean normalize) {

		if (subGroups == null) {
			throw new IllegalArgumentException("SubGroups list must not be null");
		}
		double[][][] pairs = new double[subGroups.size()][][];
		groupPairsProperties = new String[subGroups.size()][][];
		ConservationScores2 scores = null;
		String[] names = null;
		char[][] merged = null;
		int startLength = subGroups.size() - 1;
		for (int i = 0; i < subGroups.size(); i++) {
			pairs[i] = new double[startLength][];
			groupPairsProperties[i] = new String[startLength][];
			int index = 0;
			for (int j = i + 1; j < this.subGroups.size(); j++) {
				merged = mergeSubFamilies(subGroups.get(i), subGroups.get(j));
				AminoAcidMatrix mat = new AminoAcidMatrix(merged, names, null);
				scores = new ConservationScores2(mat);
				pairs[i][index] = scores.calculateScore(method, normalize);
				groupPairsProperties[i][index] = zvelibilPropertiesAlignment(mat);
			}
			startLength--;
		}
		return pairs;
	}

	/**
	 * Calculates SMERFS for subfamilies
	 * 
	 * @param width
	 * @param scoreType
	 * @param gapTreshold
	 * @param normalize
	 * @return
	 */
	ColumnInfo[][] subFamilyPairsConservation3(int width, SMERFSColumnScore scoreType,
			double gapTreshold, boolean normalize) {

		if (subGroups == null) {
			throw new IllegalArgumentException("SubGroups list must not be null");
		}
		// double[][][] pairs = new double[subGroups.size()][][];
		// groupPairsProperties = new String[subGroups.size()][][];
		ConservationScores2 scores = null;
		String[] names = null;
		char[][] merged = null;
		int index = 0;
		ColumnInfo[][] info = new ColumnInfo[this.alignment.numberOfColumns()][this.subGroups
				.size()
				* (this.subGroups.size() - 1) / 2];
		// int startLength = subGroups.size() - 1;
		for (int i = 0; i < subGroups.size(); i++) {
			// pairs[i] = new double[startLength][];
			// groupPairsProperties[i] = new String[startLength][];
			// int index = 0;
			for (int j = i + 1; j < this.subGroups.size(); j++) {
				merged = mergeSubFamilies(subGroups.get(i), subGroups.get(j));
				AminoAcidMatrix mat = new AminoAcidMatrix(merged, names, null);
				// scores = new ConservationScores2(mat);
				// pairs[i][index] = scores.calculateScore(method, normalize);
				// groupPairsProperties[i][index] =
				// zvelibilPropertiesAlignment(mat);
				double[] score = ConservationClient.getSMERFS(alignment, width, scoreType,
						gapTreshold, normalize);
				String[] props = zvelibilPropertiesAlignment(mat);
				for (int a = 0; a < score.length; a++) {
					String group = i + "-" + j;
					info[a][index] = new ColumnInfo(group, score[a], props[a]);
				}
				index++;
			}
			// startLength--;
		}
		return info;
	}

	/**
	 * Calculates conservation for pairs of subfamilies.
	 * 
	 * @param method
	 * @param normalize
	 * @return
	 */
	ColumnInfo[][] subFamilyPairsConservation2(Method method, boolean normalize) {

		if (subGroups == null) {
			throw new IllegalArgumentException("SubGroups list must not be null");
		}
		// double[][][] pairs = new double[subGroups.size()][][];
		// groupPairsProperties = new String[subGroups.size()][][];
		ConservationScores2 scores = null;
		String[] names = null;
		char[][] merged = null;
		int index = 0;
		ColumnInfo[][] info = new ColumnInfo[this.alignment.numberOfColumns()][this.subGroups
				.size()
				* (this.subGroups.size() - 1) / 2];
		// int startLength = subGroups.size() - 1;
		for (int i = 0; i < subGroups.size(); i++) {
			// pairs[i] = new double[startLength][];
			// groupPairsProperties[i] = new String[startLength][];
			// int index = 0;
			for (int j = i + 1; j < this.subGroups.size(); j++) {
				merged = mergeSubFamilies(subGroups.get(i), subGroups.get(j));
				AminoAcidMatrix mat = new AminoAcidMatrix(merged, names, null);
				scores = new ConservationScores2(mat);
				// pairs[i][index] = scores.calculateScore(method, normalize);
				// groupPairsProperties[i][index] =
				// zvelibilPropertiesAlignment(mat);
				double[] score = scores.calculateScore(method, normalize);
				String[] props = zvelibilPropertiesAlignment(mat);
				for (int a = 0; a < score.length; a++) {
					String group = i + "-" + j;
					info[a][index] = new ColumnInfo(group, score[a], props[a]);
				}
				index++;
			}
			// startLength--;
		}
		return info;
	}

	/**
	 * Gives conservation category based on the results
	 * 
	 * @param subResults
	 * @param method
	 */
	void subFamilyResults(double[][] subResults, Method method) {

		groupsConsStat = new ConservationStatus[subResults.length][subResults[0].length];
		for (int i = 0; i < subResults.length; i++) {
			for (int j = 0; j < subResults[0].length; j++) {
				groupsConsStat[i][j] = ConservationStatus.getStatus(subResults[i][j], method);
			}
		}
	}

	/**
	 * Gives conservation category based on the results.
	 * 
	 * @param subPairs
	 * @param method
	 */
	void subPairsResults(double[][][] subPairs, Method method) {

		groupPairsConsStat = new ConservationStatus[subPairs.length][][];
		for (int i = 0; i < subPairs.length; i++) {
			groupPairsConsStat[i] = new ConservationStatus[subPairs[i].length][];
			for (int j = 0; j < groupPairsConsStat[i].length; j++) {
				groupPairsConsStat[i][j] = new ConservationStatus[subPairs[i][j].length];
				for (int k = 0; k < groupPairsConsStat[i][j].length; k++) {
					groupPairsConsStat[i][j][k] = ConservationStatus.getStatus(subPairs[i][j][k],
							method);
				}
			}
		}
	}

	void printResults(PrintWriter print) {

		for (int i = 0; i < groupsConsStat.length; i++) {
			for (int j = 0; j < groupsConsStat[i].length; j++) {
				print.println("Group: " + i + " Column: " + j + " Status: "
						+ ConservationStatus.stringReps(groupsConsStat[i][j]));
			}
			print.println();
		}
		for (int i = 0; i < groupPairsConsStat.length; i++) {
			for (int j = 0; j < groupPairsConsStat[i].length; j++) {
				for (int k = 0; k < groupPairsConsStat[i][j].length; k++) {
					int gr1 = i;
					int gr2 = i + j + 1;
					print.println("GroupPair: " + gr1 + "-" + gr2 + " Column: " + k + " Status: "
							+ ConservationStatus.stringReps(groupPairsConsStat[i][j][k]));
				}
				print.println();
			}
			print.println();
		}
	}

	public static void main(String[] args) {

		String filePath = "/homes/agolicz/pFamilies/Family6/PF12407_seed.txt";
		InputStream inStr = null;
		List<FastaSequence> fastaSeqs = null;
		try {
			inStr = new FileInputStream(filePath);
		} catch (IOException e) {
			System.out.println("Can not find file");
		}
		try {
			fastaSeqs = SequenceUtil.readFasta(inStr);
		} catch (IOException e) {
			System.out.println("Sth wrong with reading the file");
		}
		AminoAcidMatrix matrix = new AminoAcidMatrix(fastaSeqs, null);
		int[][] groups = { { 0, 1 }, { 2 }, { 3, 4 } };
		SubFamiliesConservation cons = new SubFamiliesConservation(matrix, groups);
		ColumnInfo[][] cons1 = cons.subgrupsConservation2(Method.ZVELIBIL, true);
		ColumnInfo[][] cons2 = cons.subFamilyPairsConservation2(Method.ZVELIBIL, true);
		// cons.subFamilyResults(cons1, Method.ZVELIBIL);
		// cons.subPairsResults(cons2, Method.ZVELIBIL);
		String fileName = "familyRespF6_ZVELIBIL_0.7_new.txt";
		PrintWriter print = null;
		try {
			print = new PrintWriter(new BufferedWriter(new FileWriter(fileName)));
		} catch (IOException ex) {
			System.out.println("Problem writing" + fileName);
		}
		int width1 = 20;
		int width2 = 20;
		String format1 = "%-" + width1 + "s";
		String format2 = "%-" + width2 + "s";
		print.printf(format1, "");
		print.printf(format1, "Group name");
		print.printf(format2, "Conservation");
		print.println("Properties");
		for (int i = 0; i < cons1.length; i++) {
			print.println("Column nr: " + i);
			for (int j = 0; j < cons1[i].length; j++) {
				print.printf(format1, "");
				cons1[i][j].printInfo(Method.ZVELIBIL, print);
			}
			for (int j = 0; j < cons2[i].length; j++) {
				print.printf(format1, "");
				cons2[i][j].printInfo(Method.ZVELIBIL, print);
			}
		}
		print.close();
	}

	String zvelibilProperties(AminoAcidMatrix matrix, int colNr) {

		String properties = "";
		Map<String, HashSet<Character>> setMap = ConservationSets.zvelibilSets();
		Set<String> keys = setMap.keySet();
		Iterator<String> itr = keys.iterator();
		while (itr.hasNext()) {
			String key = itr.next();
			if (setMap.get(key).containsAll(matrix.getTotalAcidsFreqByCol().get(colNr).keySet())) {
				if (key.contains("Compl") != true) {
					properties += key + " ";
				}
			}
		}
		return properties;
	}

	String[] zvelibilPropertiesAlignment(AminoAcidMatrix matrix) {

		String[] properties = new String[matrix.numberOfColumns()];
		for (int i = 0; i < matrix.numberOfColumns(); i++) {
			properties[i] = zvelibilProperties(matrix, i);
		}
		return properties;
	}
}
