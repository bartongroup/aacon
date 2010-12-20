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

import static org.testng.Assert.fail;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import compbio.data.sequence.ConservationMethod;
import compbio.data.sequence.FastaSequence;
import compbio.data.sequence.SMERFSConstraints;
import compbio.data.sequence.SequenceUtil;
import compbio.data.sequence.UnknownFileFormatException;
import compbio.util.Timer;

public class SlowMethodTester {

	public final static String DATA_PATH = "test/data";

	static final String TINY_AL = "small.align";
	static final String SMALL_AL = "TO1296.fasta.align";
	static final String AVG_AL = "avg.aln.fa";
	static final String AVG_AL2 = "avg2.aln.fa";
	static final String AVG_AL4 = "avg4.aln.fa";
	static final String LARGE_AL = "1000x3000Dna.aln.fa";
	static ExecutorFactory efactory;

	@BeforeClass
	public void init() {
		ExecutorFactory.initExecutor();
	}

	@Test(invocationCount = 5)
	public void testSadler() {
		try {
			Timer timer = new Timer(TimeUnit.MILLISECONDS);
			List<FastaSequence> sequences = SequenceUtil
					.readFasta(new FileInputStream(new File(DATA_PATH
							+ File.separator + AVG_AL4)));
			System.out.println("Loading sequences: " + timer.getStepTime());

			AminoAcidMatrix alignment = new AminoAcidMatrix(sequences, null);
			System.out.println("Converting to Matrix: " + timer.getStepTime());

			Conservation scores = new Conservation(alignment, false,
					ExecutorFactory.getExecutor());
			System.out.println("Constructing conservation scores: "
					+ timer.getStepTime());

			double[] result = scores.calculateScore(ConservationMethod.VALDAR);
			System.out.println("Calculating sadler scores: "
					+ timer.getStepTime());

			// Conservation.printResults(result, Method.SANDER);
			scores.outputResults(new File("results.txt" + Math.random()),
					Format.RESULT_NO_ALIGNMENT);
			System.out.println("Total: " + timer.getTotalTime());

		} catch (FileNotFoundException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		} catch (IOException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		}
	}

	@Test(enabled = false)
	public void testLandgraf() {
		try {
			Timer timer = new Timer(TimeUnit.MILLISECONDS);
			List<FastaSequence> sequences = SequenceUtil
					.readFasta(new FileInputStream(new File(DATA_PATH
							+ File.separator + AVG_AL)));
			System.out.println("Loading sequences: " + timer.getStepTime());

			AminoAcidMatrix alignment = new AminoAcidMatrix(sequences, null);
			System.out.println("Converting to Matrix: " + timer.getStepTime());

			Conservation scores = new Conservation(alignment, true,
					ExecutorFactory.getExecutor());
			System.out.println("Constructing conservation scores: "
					+ timer.getStepTime());

			double[] result = scores
					.calculateScore(ConservationMethod.LANDGRAF);
			System.out.println("Calculating sadler scores: "
					+ timer.getStepTime());

			System.out.println("#LADGRAF " + Arrays.toString(result));
			System.out.println("Total: " + timer.getTotalTime());

		} catch (FileNotFoundException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		} catch (IOException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		}
	}

	@Test(enabled = false)
	public void testSMERFS() {
		try {
			Timer timer = new Timer(TimeUnit.MILLISECONDS);
			// 675706 ms
			List<FastaSequence> sequences = SequenceUtil
					.readFasta(new FileInputStream(new File(DATA_PATH
							+ File.separator + AVG_AL)));
			System.out.println("Loading sequences: " + timer.getStepTime());

			AminoAcidMatrix alignment = new AminoAcidMatrix(sequences, null);
			System.out.println("Converting to Matrix: " + timer.getStepTime());

			Conservation scores = new Conservation(alignment, true,
					ExecutorFactory.getExecutor());
			System.out.println("Constructing conservation scores: "
					+ timer.getStepTime());

			double[] result = scores.getSMERFS(11, SMERFSConstraints.MID_SCORE,
					0.2);
			System.out.println("Calculating SMERFS scores: "
					+ timer.getStepTime());

			// this is a wrong call!
			System.out.println(Arrays.toString(result));
			System.out.println("Total: " + timer.getTotalTime());

		} catch (FileNotFoundException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		} catch (IOException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		}
	}

	@Test(expectedExceptions = IllegalArgumentException.class)
	public void reproduceIssue1() {
		try {
			Timer timer = new Timer(TimeUnit.MILLISECONDS);
			List<FastaSequence> sequences = SequenceUtil
					.readFasta(new FileInputStream(new File(DATA_PATH
							+ File.separator + AVG_AL)));
			System.out.println("Loading sequences: " + timer.getStepTime());

			AminoAcidMatrix alignment = new AminoAcidMatrix(sequences, null);
			System.out.println("! "
					+ Arrays.toString(alignment.getInverseMatrix()[17]));
			System.out.println("Converting to Matrix: " + timer.getStepTime());

			Conservation scores = new Conservation(alignment, true,
					ExecutorFactory.getExecutor());
			System.out.println("Constructing conservation scores: "
					+ timer.getStepTime());

			double[] result = scores.calculateScore(ConservationMethod.JORES);
			System.out.println("Calculating sadler scores: "
					+ timer.getStepTime());

			System.out.println("#KABAT " + Arrays.toString(result));
			System.out.println("Total: " + timer.getTotalTime());

		} catch (FileNotFoundException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		} catch (IOException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		}
	}

	@Test(expectedExceptions = IllegalArgumentException.class)
	public void reproduceIssue4() {
		try {
			Timer timer = new Timer(TimeUnit.MILLISECONDS);
			List<FastaSequence> sequences = SequenceUtil
					.readFasta(new FileInputStream(new File(DATA_PATH
							+ File.separator + AVG_AL)));
			System.out.println("Loading sequences: " + timer.getStepTime());

			AminoAcidMatrix alignment = new AminoAcidMatrix(sequences, null);
			System.out.println("Converting to Matrix: " + timer.getStepTime());

			Conservation scores = new Conservation(alignment, true,
					ExecutorFactory.getExecutor());
			System.out.println("Constructing conservation scores: "
					+ timer.getStepTime());

			double[] result = scores.calculateScore(ConservationMethod.JORES);
			System.out.println("Calculating sadler scores: "
					+ timer.getStepTime());

			System.out.println("#JORES " + Arrays.toString(result));
			System.out.println("Total: " + timer.getTotalTime());

		} catch (FileNotFoundException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		} catch (IOException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		}
	}

	@Test()
	public void testAPIgetConservation() {
		try {
			Timer timer = new Timer(TimeUnit.MILLISECONDS);
			File input = new File(DATA_PATH + File.separator + SMALL_AL);

			List<FastaSequence> sequences = SequenceUtil
					.readFasta(new FileInputStream(input));
			System.out.println("Loading sequences: " + timer.getStepTime());

			Conservation cons = Conservation.getConservation(input, true,
					ExecutorFactory.getExecutor());

			double[] results = cons.calculateScore(ConservationMethod.ARMON);

			AminoAcidMatrix alignment = new AminoAcidMatrix(sequences, null);
			System.out.println("Converting to Matrix: " + timer.getStepTime());

			Conservation scores = new Conservation(alignment, true,
					ExecutorFactory.getExecutor());
			System.out.println("Constructing conservation scores: "
					+ timer.getStepTime());

			double[] result = scores.calculateScore(ConservationMethod.ARMON);
			Assert.assertTrue(Arrays.equals(results, result));

			System.out.println("Calculating sadler scores: "
					+ timer.getStepTime());

			System.out.println("#JORES " + Arrays.toString(result));
			System.out.println("Total: " + timer.getTotalTime());

		} catch (FileNotFoundException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		} catch (IOException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		} catch (UnknownFileFormatException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		}
	}

	public static void main(String[] args) {
		ExecutorFactory.initExecutor();

		try {
			Timer timer = new Timer(TimeUnit.MILLISECONDS);

			List<FastaSequence> sequences = SequenceUtil
					.readFasta(new FileInputStream(new File(DATA_PATH
							+ File.separator + AVG_AL)));
			System.out.println("Loading sequences: " + timer.getStepTime());

			AminoAcidMatrix alignment = new AminoAcidMatrix(sequences, null);
			System.out.println("Converting to Matrix: " + timer.getStepTime());

			Conservation scores = new Conservation(alignment, true,
					ExecutorFactory.getExecutor());
			System.out.println("Constructing conservation scores: "
					+ timer.getStepTime());

			double[] result = scores.calculateScore(ConservationMethod.JORES);
			System.out.println("Calculating sadler scores: "
					+ timer.getStepTime());

			System.out.println("#KABAT " + Arrays.toString(result));
			System.out.println("Total: " + timer.getTotalTime());

		} catch (FileNotFoundException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		} catch (IOException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		}
	}
}
