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

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.List;
import java.util.concurrent.ExecutionException;

import org.testng.annotations.Test;

import compbio.data.sequence.FastaSequence;
import compbio.data.sequence.SequenceUtil;

public class CorrelationTester {

	void printCoeffs(double[] arr1) {

		// assert arr1.length == arr2.length;

		for (int i = 0; i < arr1.length; i++) {

			System.out.println(i + ": " + arr1[i]);
		}
	}

	@Test
	void Corr1Tester() throws InterruptedException, ExecutionException {

		String filePath = "/homes/agolicz/alignments/alignment1";

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

		AminoAcidMatrix matrix = new AminoAcidMatrix(fastaSeqs, null);

		// Correlation corr = new Correlation(matrix, 5, 0.1);

		// double[] coeffs = corr.calcPearson();

		// double[] coeffs2 = Correlation.calcPearsonCoeff3(matrix, 5);

		// int[][] localSim = Correlation.localSimilarity2(matrix);

		// int[] globalSim = Correlation.globalSimilarity(matrix);

		// printCoeffs(coeffs, coeffs2);

		// double[] scores = corr.giveMaxToColumn(coeffs);

		// printCoeffs(coeffs);

	}

}
