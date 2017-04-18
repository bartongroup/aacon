/*
 * Copyright (c) 2010 Agnieszka Golicz & Peter Troshin 
 * 
 * Amino Acid Conservation @version: 1.0 
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the Apache License version 2 as published by the
 * Apache Software Foundation This library is distributed in the hope that it
 * will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Apache
 * License for more details. A copy of the license is in apache_license.txt. It
 * is also available here: http://www.apache.org/licenses/LICENSE-2.0.txt 
 * Any republication or derived work distributed in source code form must 
 * include this copyright and license notice.
 * 
 */
package compbio.conservation;

// creates a matrix of aa in multiple alignment
// gets fasta sequences and puts them into matrix
// might have to check if all sequences are equal length(ask)
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import compbio.common.NotAnAminoAcidException;
import compbio.common.SequencesNotEquallyLongException;
import compbio.data.sequence.Alignment;
import compbio.data.sequence.FastaSequence;

/**
 * This class provides representation of an alignment as a matrix implemented as
 * 2D array. Rows correspond to the sequences. Columns correspond to the
 * vertical columns in the alignment consisting of amino acids with the same
 * index in all the sequences. The only condition is that all the the sequences
 * in the fasta file are of the same length. This class creates a matrix which
 * has rows of equal length. If constructor is fed an alignment with sequences
 * of not equal length an exception is thrown.
 * 
 * @author Agnieszka Golicz & Peter Troshin
 */
final class AminoAcidMatrix {

	/**
	 * Stores the matrix.
	 */
	private final char[][] matrix;
	/**
	 * Stores inverse matrix
	 */
	private final char[][] inverseMatrix;
	/**
	 * Holds names of the sequences in the alignment.
	 */
	private String[] sequenceNames;
	/**
	 * Stores occurances of amino acids in columns, columns indexed starting
	 * fromm 0
	 */
	private List<Map<Character, Integer>> acidsIntMap;
	/**
	 * Holds the in the indices of Xs changed to gaps. The row number is a key
	 * and the column number is value
	 */
	private char[] gaps;
	// private List<HashMap<Integer, Integer>> xToGapSubs;
	/**
	 * The total occurrence of amino acids in the whole alignment.
	 */
	private Map<Character, Integer> totalFrequency;
	/**
	 * The total number of amino acids in the whole alignment belonging to each
	 * Williamson Set.
	 */
	private Map<String, Integer> willSetsTotal;
	/**
	 * Vingron Argos weights of the the sequences.
	 */
	private double[] vingronArgosWeights;
	/**
	 * Percent identity.
	 */
	private double[][] percentIdentity;
	/**
	 * Weights according to Voronoi.
	 */
	private double[] voronoiWeighths;

	/**
	 * This constructor take a 2D array of characters as parameter and a string
	 * of sequence names. Indices of names must be corresponding to the indices
	 * of sequences in the 2D array. If no names are provided sequnece numbers
	 * starting with 0 serve as names.
	 * 
	 * @param alignment
	 *            2D array with alignment
	 * @param names
	 *            string if sequence names/IDs
	 * 
	 *            public AminoAcidMatrix(char[][] alignment, String[] names,
	 *            Character[] gapChars) {
	 * 
	 *            for (int i = 0; i < alignment.length; i++) { if
	 *            (alignment[0].length != alignment[i].length) { int number = i
	 *            + 1; throw new
	 *            SequencesNotEquallyLongException("Sequence number: " + number
	 *            + "ID: " + this.sequenceNames[i] +
	 *            " is of differnet length than the first sequence"); } } if
	 *            (names != null) { if (alignment.length != names.length) {
	 *            throw new IllegalArgumentException(
	 *            "Number sequence IDs not equall to the number of sequences");
	 *            } } if (names == null) { this.sequenceNames = new
	 *            String[alignment.length]; for (int i = 0; i <
	 *            alignment.length; i++) { sequenceNames[i] = i + "."; } } else
	 *            { this.sequenceNames = names; } if (gapChars != null) {
	 *            this.gaps = gapChars; } this.matrix = new
	 *            char[alignment.length][alignment[0].length]; inverseMatrix =
	 *            new char[alignment[0].length][alignment.length]; for (int i =
	 *            0; i < alignment.length; i++) { for (int j = 0; j <
	 *            alignment[i].length; j++) { char ch = alignment[i][j]; if
	 *            (gapChars == null) { if (ch == '.' || ch == '*' || ch == ' '
	 *            || ch == 'X') { ch = '-'; } } else { for (int d = 0; d <
	 *            gaps.length; d++) { if (ch == gaps[d]) { ch = '-'; } } }
	 *            matrix[i][j] = ch; inverseMatrix[j][i] = ch; } } }
	 */

	/**
	 * This constructor constructor allows manual creation of only one column.
	 * Might be of help if somebody wants to check the the functionality of the
	 * class without feeding it the whole alignment.
	 * 
	 * @param column
	 *            array that represents one column in the alignment
	 */
	AminoAcidMatrix(char[] column) {

		Set<Character> alp = Alphabet.alphabet();
		matrix = new char[column.length][1];
		inverseMatrix = new char[1][column.length];
		for (int i = 0; i < column.length; i++) {
			if (alp.contains(column[i]) == false) {
				throw new NotAnAminoAcidException(
						"Illegal chracter in the column");
			}
			matrix[i][0] = column[i];
			inverseMatrix[0][i] = column[i];
		}
		this.calTotalAcidsFreqByCol();
	}

	/**
	 * Takes characters given and creates a a matrix 3x3. Used for tests only.
	 * 
	 * @param p1
	 *            first row first char
	 * @param p2
	 *            first row second char
	 * @param p3
	 *            first row third char
	 * @param p4
	 *            second row first char
	 * @param p5
	 *            second row second char
	 * @param p6
	 *            second row third char
	 * @param p7
	 *            third row first char
	 * @param p8
	 *            third row second char
	 * @param p9
	 *            third row third char
	 */
	AminoAcidMatrix(char p1, char p2, char p3, char p4, char p5, char p6,
			char p7, char p8, char p9) {

		Set<Character> alp = Alphabet.alphabet();
		if (alp.contains(p1) == false) {
			throw new NotAnAminoAcidException("Illegal chracter in the column");
		}
		if (alp.contains(p2) == false) {
			throw new NotAnAminoAcidException("Illegal chracter in the column");
		}
		if (alp.contains(p3) == false) {
			throw new NotAnAminoAcidException("Illegal chracter in the column");
		}
		if (alp.contains(p4) == false) {
			throw new NotAnAminoAcidException("Illegal chracter in the column");
		}
		if (alp.contains(p5) == false) {
			throw new NotAnAminoAcidException("Illegal chracter in the column");
		}
		if (alp.contains(p6) == false) {
			throw new NotAnAminoAcidException("Illegal chracter in the column");
		}
		if (alp.contains(p7) == false) {
			throw new NotAnAminoAcidException("Illegal chracter in the column");
		}
		if (alp.contains(p8) == false) {
			throw new NotAnAminoAcidException("Illegal chracter in the column");
		}
		if (alp.contains(p9) == false) {
			throw new NotAnAminoAcidException("Illegal chracter in the column");
		}
		matrix = new char[3][3];
		matrix[0][0] = p1;
		matrix[0][1] = p2;
		matrix[0][2] = p3;
		matrix[1][0] = p4;
		matrix[1][1] = p5;
		matrix[1][2] = p6;
		matrix[2][0] = p7;
		matrix[2][1] = p8;
		matrix[2][2] = p9;
		inverseMatrix = new char[3][3];
		inverseMatrix[0][0] = p1;
		inverseMatrix[0][1] = p4;
		inverseMatrix[0][2] = p7;
		inverseMatrix[1][0] = p2;
		inverseMatrix[1][1] = p5;
		inverseMatrix[1][2] = p8;
		inverseMatrix[2][0] = p3;
		inverseMatrix[2][1] = p6;
		inverseMatrix[2][2] = p9;
	}

	/**
	 * Make a AminoAcidMatrix from Alignment
	 * 
	 * @param alignment
	 *            the alignment
	 */
	AminoAcidMatrix(Alignment alignment) {

		this(alignment.getSequences(), new char[] { alignment.getMetadata()
				.getGapchar() });
	}

	/**
	 * Takes characters given and creates a a matrix 5x3. Used for tests only.
	 */
	AminoAcidMatrix(char p1, char p2, char p3, char p4, char p5, char p6,
			char p7, char p8, char p9, char p10, char p11, char p12, 
			char p13, char p14, char p15) {

		Set<Character> alp = Alphabet.alphabet();
		if (alp.contains(p1) == false) {
			throw new NotAnAminoAcidException("Illegal chracter in the column");
		}
		if (alp.contains(p2) == false) {
			throw new NotAnAminoAcidException("Illegal chracter in the column");
		}
		if (alp.contains(p3) == false) {
			throw new NotAnAminoAcidException("Illegal chracter in the column");
		}
		if (alp.contains(p4) == false) {
			throw new NotAnAminoAcidException("Illegal chracter in the column");
		}
		if (alp.contains(p5) == false) {
			throw new NotAnAminoAcidException("Illegal chracter in the column");
		}
		if (alp.contains(p6) == false) {
			throw new NotAnAminoAcidException("Illegal chracter in the column");
		}
		if (alp.contains(p7) == false) {
			throw new NotAnAminoAcidException("Illegal chracter in the column");
		}
		if (alp.contains(p8) == false) {
			throw new NotAnAminoAcidException("Illegal chracter in the column");
		}
		if (alp.contains(p9) == false) {
			throw new NotAnAminoAcidException("Illegal chracter in the column");
		}
		if (alp.contains(p10) == false) {
			throw new NotAnAminoAcidException("Illegal chracter in the column");
		}
		if (alp.contains(p11) == false) {
			throw new NotAnAminoAcidException("Illegal chracter in the column");
		}
		if (alp.contains(p12) == false) {
			throw new NotAnAminoAcidException("Illegal chracter in the column");
		}
		if (alp.contains(p13) == false) {
			throw new NotAnAminoAcidException("Illegal chracter in the column");
		}
		if (alp.contains(p14) == false) {
			throw new NotAnAminoAcidException("Illegal chracter in the column");
		}
		if (alp.contains(p15) == false) {
			throw new NotAnAminoAcidException("Illegal chracter in the column");
		}
		matrix = new char[3][5];
		matrix[0][0] = p1;
		matrix[0][1] = p2;
		matrix[0][2] = p3;
		matrix[0][3] = p4;
		matrix[0][4] = p5;
		matrix[1][0] = p6;
		matrix[1][1] = p7;
		matrix[1][2] = p8;
		matrix[1][3] = p9;
		matrix[1][4] = p10;
		matrix[2][0] = p11;
		matrix[2][1] = p12;
		matrix[2][2] = p13;
		matrix[2][3] = p14;
		matrix[2][4] = p15;
		inverseMatrix = new char[5][3];
		inverseMatrix[0][0] = p1;
		inverseMatrix[0][1] = p6;
		inverseMatrix[0][2] = p11;
		inverseMatrix[1][0] = p2;
		inverseMatrix[1][1] = p7;
		inverseMatrix[1][2] = p12;
		inverseMatrix[2][0] = p3;
		inverseMatrix[2][1] = p8;
		inverseMatrix[2][2] = p13;
		inverseMatrix[3][0] = p4;
		inverseMatrix[3][1] = p9;
		inverseMatrix[3][2] = p14;
		inverseMatrix[4][0] = p5;
		inverseMatrix[4][1] = p10;
		inverseMatrix[4][2] = p15;
		
		assert inverseMatrix != null && inverseMatrix[0] != null;
		for (char[] column : inverseMatrix) {
			char[] gapchars;
			gapchars = DEFAULT_GAP_CHARS;
			boolean gaponly = true;
			Arrays.sort(gapchars);
			for (char charr : column) {
				if (Arrays.binarySearch(gapchars, charr) < 0) {
					gaponly = false;
					break;
				}
			}
			// Even if only one column contains gaps the whole alignment is
			// invalid
			if (gaponly) {
				String message = "Input has badly aligned sequences with columns "
						+ "containing nothing but the gaps. "
						+ "Conservation methods cannot be calculated for "
						+ "such an alignment !" + " \nGap characters are : "
						+ Arrays.toString(gapchars);
				System.err.println(message);
				throw new IllegalArgumentException(message);
			}
		}
		
	}
	
	
	private static final char[] DEFAULT_GAP_CHARS = { '.', '*', ' ', 'X', '-' };

	/**
	 * Constructor that reads in a fasta file and creates an amino acid matrix.
	 * Gaps indicated by any sign are now indicated by a dash. The unknown amino
	 * acid X is replaced by a dash.
	 * 
	 * @param seqs
	 *            the alignment
	 * @param gapChars
	 *            the set of gap characters
	 * @throws IllegalArgumentException
	 *             if alignment contains gap only columns
	 * @throws SequencesNotEquallyLongException
	 *             if the length of the sequences in the alignment is different
	 * @throws NotAnAminoAcidException
	 *             if sequences contains letters different from AA alphabet
	 *             {@link Alphabet}
	 */
	AminoAcidMatrix(List<FastaSequence> seqs, char[] gapChars) {

		if (gapChars != null) {
			this.gaps = gapChars;
		} else {
			this.gaps = DEFAULT_GAP_CHARS;
		}
		Set<Character> alph = Alphabet.alphabet();
		int sequenceNr = seqs.size();
		FastaSequence seq = seqs.get(0);
		String firstSequence = seq.getSequence();
		int sequenceLength = firstSequence.length();
		matrix = new char[sequenceNr][sequenceLength];
		inverseMatrix = new char[sequenceLength][sequenceNr];
		sequenceNames = new String[sequenceNr];
		for (int i = 0; i < sequenceNr; i++) {
			FastaSequence s = seqs.get(i);
			char[] sequenceChars = s.getSequence().toCharArray();
			sequenceNames[i] = s.getId();
			// System.out.println(sequenceNames[i]);
			if (sequenceChars.length != sequenceLength) {
				int seqNr = i + 1;
				String message = "Sequence number " + seqNr + "(id: "
						+ sequenceNames[i] + ")"
						+ " is of differen length than previous sequences.";
				throw new SequencesNotEquallyLongException(message);
			}
			for (int j = 0; j < sequenceLength; j++) {
				char ch = sequenceChars[j];
				for (int d = 0; d < gaps.length; d++) {
					if (ch == gaps[d]) {
						ch = '-';
					}
				}
				// if (ch == 'X') {
				// this.xToGapSubs.add(new HashMap<Integer, Integer>());
				// this.xToGapSubs.get(this.xToGapSubs.size() - 1).put(i, j);
				// ch = '-';
				// }
				if (!alph.contains(ch)) {
					String legals = Alphabet.legalCharacterstoString();
					int seqNr = i + 1;
					String message = "Illegal character in sequence number "
							+ seqNr + "(sequence ID: " + sequenceNames[i]
							+ "). Illegal character: " + ch
							+ " is at position: " + j + ". "
							+ "List of legal characters: " + legals + ". ";
					throw new NotAnAminoAcidException(message);
				}
				// assert alph.contains(ch) : "Illegal character in the matrix";
				// System.out.println(ch);
				matrix[i][j] = ch;
				inverseMatrix[j][i] = ch;
			}
		}
		this.calTotalAcidsFreqByCol();
		this.gapOnlyColumnsCheck();
	}

	private void gapOnlyColumnsCheck() {
		assert inverseMatrix != null && inverseMatrix[0] != null;
		for (char[] column : inverseMatrix) {
			gapOnlyColumnCheck(column, this.gaps);
		}
	}

	/**
	 * 
	 * @param column
	 * @param gapchars
	 * @throws IllegalArgumentException
	 *             if column contains gapchars only
	 */
	static void gapOnlyColumnCheck(char[] column, char[] gapchars) {
		boolean gaponly = true;
		Arrays.sort(gapchars);
		for (char charr : column) {
			if (Arrays.binarySearch(gapchars, charr) < 0) {
				gaponly = false;
				break;
			}
		}
		// Even if only one column contains gaps the whole alignment is
		// invalid
		if (gaponly) {
			String message = "Input has badly aligned sequences with columns "
					+ "containing nothing but the gaps. "
					+ "Conservation methods cannot be calculated for "
					+ "such an alignment !" + " \nGap characters are : "
					+ Arrays.toString(gapchars);
			System.err.println(message);
			throw new IllegalArgumentException(message);
		}
	}

	/**
	 * Gets the number of columns.
	 * 
	 * @return number of columns
	 */
	int numberOfColumns() {
		return matrix[0].length;
	}

	/**
	 * Gets the number of rows. E.g. number of sequences in the alignment
	 * 
	 * @return number of rows
	 */
	int numberOfRows() {
		return matrix.length;
	}

	/**
	 * Gets the whole matrix
	 * 
	 * @return matrix
	 */
	char[][] getMatrix() {

		return matrix;
	}

	/**
	 * Gets the inverse matrix
	 * 
	 * @return inverseMatrix
	 */
	char[][] getInverseMatrix() {

		return inverseMatrix;
	}

	/**
	 * Returns character at a specified position.
	 * 
	 * @param row
	 * @param column
	 * @return character
	 */
	char getMatrixPosition(int row, int column) {

		char position = matrix[row][column];
		return position;
	}

	/**
	 * Gets a column specified by given number. Indexing starts with 0.
	 * 
	 * @param number
	 * @return column as an array of characters
	 */
	char[] getColumn(int number) {

		return inverseMatrix[number];
	}

	/**
	 * Gets row specified by a given number. Indexing starts with 0.
	 * 
	 * @param number
	 * @return row as an array of characters
	 */
	char[] getRow(int number) {
		assert number < this.numberOfRows();
		return matrix[number];
	}

	private void calTotalAcidsFreqByCol() {
		acidsIntMap = new ArrayList<Map<Character, Integer>>(
				this.numberOfColumns());
		for (int i = 0; i < this.numberOfColumns(); i++) {
			acidsIntMap.add(Alphabet.calculateOccurance(this.inverseMatrix[i]));
		}
	}

	synchronized List<Map<Character, Integer>> getTotalAcidsFreqByCol() {
		if (acidsIntMap == null) {
			this.calTotalAcidsFreqByCol();
		}
		return acidsIntMap;
	}

	/**
	 * Gets the total amino acid occurrence in the whole alignment.
	 * 
	 * @return map with characters as keys and occurrence as values.
	 */
	synchronized Map<Character, Integer> totalAcidsFrequency() {
		if (totalFrequency == null) {
			this.calTotalAcidsFrequency();
		}
		return totalFrequency;
	}

	/**
	 * Gets the total number of amino acids belonging to particular Williamson
	 * sets.
	 * 
	 * @return map with set names as keys and number of amino acids belonging to
	 *         the set as value.
	 */
	synchronized Map<String, Integer> totalAcidsWillSets() {
		if (willSetsTotal == null) {
			this.calTotalAcidsWillSets();
		}
		return willSetsTotal;
	}

	/**
	 * Calculates the total amino acid occurrence in the whole alignment.
	 */
	private void calTotalAcidsFrequency() {

		Map<Character, Integer> totalFreq = new HashMap<Character, Integer>();
		Set<Character> alph = Alphabet.alphabet();
		for (int i = 0; i < this.numberOfRows(); i++) {
			for (int j = 0; j < this.numberOfColumns(); j++) {
				Character ch = matrix[i][j];
				if (ch == '.' || ch == '*' || ch == ' ' || ch == 'X') {
					ch = '-';
				}
				assert alph.contains(ch) : "Illegal character in the matrix";
				Integer count = totalFreq.get(ch);
				if (count == null) {
					totalFreq.put(ch, 1);
				} else {
					totalFreq.put(ch, count + 1);
				}
			}
		}
		totalFrequency = totalFreq;
	}

	/**
	 * Calculates the total number of amino acids belonging to particular
	 * Williamson sets.
	 */
	private void calTotalAcidsWillSets() {

		Map<String, HashSet<Character>> sets = ConservationSets
				.williamsonSets();
		Map<String, Integer> setsFreq = new HashMap<String, Integer>();
		Set<String> setsKeys = sets.keySet();
		Iterator<String> setsKeysItr = setsKeys.iterator();

		totalAcidsFrequency();

		Map<Character, Integer> totalFreq = totalFrequency;
		Set<Character> totalFreqKeys = totalFreq.keySet();
		while (setsKeysItr.hasNext()) {
			String setsKey = setsKeysItr.next();
			Iterator<Character> totalFreqItr = totalFreqKeys.iterator();
			while (totalFreqItr.hasNext()) {
				Character totalFreqKey = totalFreqItr.next();
				if (sets.get(setsKey).contains(totalFreqKey)) {
					Integer count = setsFreq.get(setsKey);
					if (count == null) {
						setsFreq.put(setsKey, totalFreq.get(totalFreqKey));
					} else {
						setsFreq.put(setsKey,
								count + totalFreq.get(totalFreqKey));
					}
				}
			}
		}
		willSetsTotal = setsFreq;
	}

	/**
	 * Calculates a weight of particular sequence according to Vingron-Argos
	 * model.
	 * 
	 * @param seqNr
	 * @return sequence weight
	 */
	private double weightOfSequenceVingronArgos(int seqNr) {

		double weight = 0.0;
		for (int i = 0; i < this.numberOfRows(); i++) {
			if (i != seqNr) {
				weight += ConservationAccessory.percentIdentity(
						this.getRow(seqNr), this.getRow(i));
			}
		}
		double result = (1.0 / this.numberOfRows()) * weight;
		return result;
	}

	/**
	 * Calculates the weight for all the sequences in the alignment. Weight
	 * calculated according to Vingron-Argos model
	 */
	private void weightOfSequencesVingronArgos() {

		vingronArgosWeights = new double[this.numberOfRows()];
		for (int i = 0; i < this.numberOfRows(); i++) {
			vingronArgosWeights[i] = this.weightOfSequenceVingronArgos(i);
		}
	}

	/**
	 * Gets the values of Vingron-Argos weights for the whole alignment. Indices
	 * in the weights array correspond to the indices of the sequences in the
	 * matrix.
	 * 
	 * @return an array of weights, indices correspond to sequence numbers
	 */
	synchronized double[] vingronArgosWeights() {
		if (vingronArgosWeights == null) {
			this.weightOfSequencesVingronArgos();
		}
		return vingronArgosWeights;
	}

	/**
	 * Calculates percent identity for all the sequences in the alignment.
	 * Stores calculated values in a 2D array. The sequence number index the
	 * values in the array. For example percentage identity of sequence 0 and 5
	 * perIden[0][5]
	 */
	private void calPercentIdentity() {

		percentIdentity = new double[this.numberOfRows()][this.numberOfRows()];
		int ident = 0;
		for (int i = 0; i < this.numberOfRows(); i++) {
			for (int j = i + 1; j < this.numberOfRows(); j++) {
				for (int a = 0; a < this.numberOfColumns(); a++) {
					if (this.getRow(i) == this.getRow(j)) {
						ident++;
					}
				}
				double result = (double) ident / this.numberOfColumns();
				percentIdentity[i][j] = result;
			}
		}
	}

	/**
	 * Gets percentage identity for the whole alignment.
	 * 
	 * @return 2D array; the sequences' numbers index the values in the array.
	 */
	synchronized double[][] getPercentIdentity() {
		if (percentIdentity == null) {
			this.calPercentIdentity();
		}
		return percentIdentity;
	}

	/**
	 * Scary, scary method. Calculates weights according to Voronoi scheme.
	 * 
	 * @param iter
	 */
	void voronoiWeights(int iterNr) {

		int iterations = iterNr;
		Random rgen = new Random();
		double[] weights = new double[numberOfRows()];
		char[] randSeq = new char[numberOfColumns()];
		// sets iterations, don't really know what does it to, but jon set up
		// the iterations to 1000
		for (int i = 0; i < iterations; i++) {
			// generates a random sequence, equal in length to the sequences in
			// the alignment
			for (int j = 0; j < numberOfColumns(); j++) {
				int random = rgen.nextInt(numberOfRows());
				randSeq[j] = getMatrixPosition(random, j);
			}
			// measure the distance between each sequence and a random sequence
			// generated
			// distance measured as percentage identity
			double[] distances = new double[numberOfRows()];
			double closestValue = 1000;
			for (int a = 0; a < numberOfRows(); a++) {
				distances[a] = 1.0 - ConservationAccessory.percentIdentity(
						getRow(a), randSeq);
				if (distances[a] < closestValue) {
					closestValue = distances[a];
				}
			}
			// collect all the sequences with the closest distance
			List<Integer> closestSeqs = new ArrayList<Integer>();
			for (int b = 0; b < distances.length; b++) {
				double dis = distances[b];
				if (dis == closestValue) {
					closestSeqs.add(b);
				}
			}
			// increase by one the weight of the closest sequence
			double increase = 1.0 / closestSeqs.size();
			for (int c = 0; c < closestSeqs.size(); c++) {
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
		double scaleFactor = weightSum / numberOfRows();
		for (int e = 0; e < weights.length; e++) {
			weights[e] = weights[e] + scaleFactor;
		}
		voronoiWeighths = weights;
	}

	/**
	 * Gets Voronoi weights.
	 * 
	 * @param iterNr
	 * @return array containing voronoi weights, indices correspond to sequence
	 *         indices in matrix
	 */
	synchronized double[] getVoronoiWeights(int iterNr) {
		if (voronoiWeighths == null) {
			this.voronoiWeights(iterNr);
		}
		return voronoiWeighths;
	}

	/**
	 * Prints alignment into a file.
	 * 
	 * @param tagWidth
	 *            width of the name tag field
	 * @param resultWidth
	 *            width of the result field
	 * @param outputFile
	 *            output file path
	 * @throws IOException
	 */
	List<FastaSequence> getAlignment() throws IOException {
		List<FastaSequence> fastaseqs = new ArrayList<FastaSequence>();
		for (int i = 0; i < this.numberOfRows(); i++) {
			StringBuffer sb = new StringBuffer();
			for (int j = 0; j < this.getRow(i).length; j++) {
				sb.append(this.getRow(i)[j]);
			}
			fastaseqs.add(new FastaSequence(sequenceNames[i], sb.toString()));
		}
		return fastaseqs;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result
				+ ((acidsIntMap == null) ? 0 : acidsIntMap.hashCode());
		result = prime * result + Arrays.hashCode(matrix);
		result = prime * result + Arrays.hashCode(sequenceNames);
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		AminoAcidMatrix other = (AminoAcidMatrix) obj;
		if (acidsIntMap == null) {
			if (other.acidsIntMap != null)
				return false;
		} else if (!acidsIntMap.equals(other.acidsIntMap))
			return false;
		if (!Arrays.equals(matrix, other.matrix))
			return false;
		if (!Arrays.equals(sequenceNames, other.sequenceNames))
			return false;
		return true;
	}

	/**
	 * Splits alignment into provided groups. Groups provided by giving sequence
	 * index(number) in the alignment. Indexing starts with 0.
	 * 
	 * @param groups
	 * @return
	 * 
	 *         List<char[][]> splitAlignment(int[][] groups) {
	 * 
	 *         List<char[][]> groupsList = new ArrayList<char[][]>(); char[][]
	 *         group = null; for (int i = 0; i < groups.length; i++) { group =
	 *         new char[groups[i].length][]; for (int j = 0; j <
	 *         groups[i].length; j++) { group[j] = matrix[groups[i][j]]; }
	 *         groupsList.add(group); } return groupsList; }
	 */

}
