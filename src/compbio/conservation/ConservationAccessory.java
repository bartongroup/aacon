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

import java.io.PrintWriter;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import compbio.data.sequence.ConservationMethod;

final class ConservationAccessory {

	/**
	 * Adds two points in N dimmentional space.
	 * 
	 * @param a
	 * @param b
	 * @return sum of points
	 */
	static double[] addPoints(double[] a, double[] b) {

		if (a == null || b == null) {
			throw new IllegalArgumentException("Reference must not be null");
		}
		assert a.length == b.length;
		double[] aPlusB = new double[a.length];
		for (int i = 0; i < a.length; i++) {
			aPlusB[i] = a[i] + b[i];
		}
		return aPlusB;
	}

	/**
	 * Multiplies point by scalar.
	 * 
	 * @param point
	 * @param scalar
	 * @return
	 */
	static double[] multPointByScalar(double[] point, double scalar) {

		double[] pointByScalar = new double[point.length];
		for (int i = 0; i < point.length; i++) {
			pointByScalar[i] = point[i] * scalar;
		}
		return pointByScalar;
	}

	/**
	 * Calculates Euclidean distance between two points in N dimmensional space.
	 * 
	 * @param a
	 * @param b
	 * @return distance
	 */
	static double pointDistance(double[] a, double[] b) {

		if (a == null || b == null) {
			throw new IllegalArgumentException("Reference must not be null");
		}
		assert a.length == b.length;
		double sum = 0;
		for (int i = 0; i < a.length; i++) {
			sum = sum + ((b[i] - a[i]) * (b[i] - a[i]));
		}
		double distance = Math.sqrt(sum);
		return distance;
	}

	/**
	 * Calculates percentage identity of two sequences. If two amino acids are
	 * identical they get a score of 1, if not a score of 0. The sum of scores
	 * divided by the length of the sequences. Sequences have to be the same
	 * length.
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	static double percentIdentity(char[] a, char[] b) {

		if (a == null || b == null) {
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

	/**
	 * Calculates Voronoi weights for all the sequences in the alignment.
	 * 
	 * @param m
	 * @param iter
	 * @return
	 */
	static double[] voronoiWeights(AminoAcidMatrix m, int iter) {

		AminoAcidMatrix matrix = m;
		int iterations = iter;
		Random rgen = new Random();
		double[] weights = new double[m.numberOfRows()];
		char[] randSeq = new char[m.numberOfColumns()];
		// satts iterations, doan't really know what does it to, but jon set up
		// the oterations to 1000
		for (int i = 0; i < iterations; i++) {
			// generates a random sequence, equal in length to the sequences in
			// the alignment
			for (int j = 0; j < matrix.numberOfColumns(); j++) {
				int random = rgen.nextInt(matrix.numberOfRows());
				randSeq[j] = matrix.getMatrixPosition(random, j);
			}
			// measure the distance between each sequence and a random sequence
			// generated
			// distance measured as percentage identity
			double[] distances = new double[matrix.numberOfRows()];
			double closestValue = 0;
			for (int a = 0; a < matrix.numberOfRows(); a++) {
				distances[a] = 1.0 - ConservationAccessory.percentIdentity(
						matrix.getRow(a), randSeq);
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
		double scaleFactor = weightSum / matrix.numberOfRows();
		for (int e = 0; e < weights.length; e++) {
			weights[e] = weights[e] + scaleFactor;
		}
		return weights;
	}

	/**
	 * Calculates Vingron and Argos sequence weight for a particular sequence.
	 * 
	 * @param seqNr
	 * @param m
	 * @return sequence weight
	 */
	static double weightOfSequenceVingronArgos(int seqNr, AminoAcidMatrix m) {

		double weight = 0.0;
		for (int i = 0; i < m.numberOfRows(); i++) {
			if (i != seqNr) {
				weight += ConservationAccessory.percentIdentity(
						m.getRow(seqNr), m.getRow(i));
			}
		}
		double result = (1.0 / m.numberOfRows()) * weight;
		return result;
	}

	/**
	 * Returns scores normalized form 0 to 1.
	 * 
	 * @param scores
	 *            array of scores to be normalized
	 * @return array of normalized scores
	 */
	/**
	 * Returns scores normalized form 0 to 1.
	 * 
	 * There is one problem with this method. If it is fed a vector where all
	 * the numbers are of equal value they are given a value of 1. However, this
	 * value does not realistically represent the conservation score. It is
	 * there to avoid method giving NaN values. If such a situation happens a
	 * message will be written to the command window.
	 * 
	 * @param scores
	 *            array of scores to be normalized
	 * @return array of normalized scores
	 */

	static double[] normalize01(double[] scores, ConservationMethod method) {

		double[] normalized = new double[scores.length];
		double max = scores[0];
		double min = scores[0];
		for (int i = 0; i < scores.length; i++) {
			if (scores[i] > max) {
				max = scores[i];
			}
			if (scores[i] < min) {
				min = scores[i];
			}
		}

		if (max == min) {
			System.err.println("Scores in method: " + method.toString()
					+ " could not have been normalized properly, "
					+ "not normalized scores are returned.");
			return null;
		}

		if (min < 0) {
			double minAbs = min * -1;
			double[] shifted = new double[scores.length];
			for (int i = 0; i < shifted.length; i++) {
				shifted[i] = scores[i] + minAbs;
			}
			max = max + minAbs;
			min = min + minAbs;
			for (int i = 0; i < normalized.length; i++) {
				normalized[i] = round(((shifted[i] - min) / (max - min)), 4);
			}
		} else {
			for (int i = 0; i < scores.length; i++) {
				normalized[i] = round((scores[i] - min) / (max - min), 4);
			}
		}

		return normalized;
	}

	/**
	 * Returns inversed version of the normalized scores. for each score gives a
	 * score equal to 1 - the original score.
	 * 
	 * @param scores
	 *            the array of scores supplied.
	 * @return the array of numbers equal 1 - normalized score.
	 */
	/**
	 * Returns inversed version of the normalized scores. for each score gives a
	 * score equal to 1 - the original score.
	 * 
	 * There is one problem with this method. If it is fed a vector where all
	 * the numbers are of equal value they are given a value of 1. However, this
	 * value does not realistically represent the conservation score. It is
	 * there to avoid method giving NaN values. If such a situation happens a
	 * message will be written to the command window.
	 * 
	 * @param scores
	 *            the array of scores supplied.
	 * 
	 * @return the array of numbers equal 1 - normalized score.
	 */
	static double[] inversedNormalize01(double[] scores,
			ConservationMethod method) {
		double[] normalized = normalize01(scores, method);
		// Can't normalize
		if (normalized == null) {
			return null;
		}
		double[] inversed = new double[normalized.length];
		for (int i = 0; i < inversed.length; i++) {
			inversed[i] = round(1 - normalized[i], 4);
		}
		return inversed;
	}

	static void printArrayOfDouble(double[] arr1, PrintWriter print,
			int precision) {

		ConservationFormatter.NUMBER_FORMAT.setMaximumFractionDigits(precision);

		assert arr1 != null : "Nothing to print!";
		PrintWriter writer = print;
		if (print == null) {
			writer = new PrintWriter(System.out);
		}

		for (int i = 0; i < arr1.length; i++) {
			writer.printf(ConservationFormatter.NUMBER_FORMAT.format(arr1[i]));
			writer.print(" ");
		}
		writer.println();
		writer.flush();
	}
	
	/**
	 * Rounds doubles to the n decimal place
	 * One of the preferred solutions at:
	 * http://stackoverflow.com/questions/153724/how-to-round-a-number-to-n-decimal-places-in-java
	 * 
	 * @return rounded number to the n decimal place
	 */
	public static double round(double d, int decimalPlace){
		
	    BigDecimal bd = new BigDecimal(Double.toString(d));
	    bd = bd.setScale(decimalPlace,BigDecimal.ROUND_HALF_UP);
	    return bd.doubleValue();
	}

}
