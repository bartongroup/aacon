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

import static org.testng.Assert.assertNotNull;
import static org.testng.Assert.fail;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;

import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import compbio.data.sequence.ConservationMethod;
import compbio.data.sequence.FastaSequence;
import compbio.data.sequence.SMERFSConstraints;
import compbio.data.sequence.SequenceUtil;
import compbio.data.sequence.UnknownFileFormatException;

public class ConservationTester {

	static ExecutorFactory efactory;
	private Map<ConservationMethod, double[]> norm_results = null;
	private Map<ConservationMethod, double[]> results = null;

	private static File input = new File(SlowMethodTester.DATA_PATH
			+ File.separator + "avg4.aln.fa");

	@BeforeClass
	public void init() {
		// Make serial calculation
		ExecutorFactory.initExecutor(1);

		List<FastaSequence> sequences = null;
		try {
			sequences = SequenceUtil.readFasta(new FileInputStream(input));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			fail(e.getMessage());
		} catch (IOException e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
		assertNotNull(sequences);

		Conservation scores = Conservation.getConservation(sequences, true,
				ExecutorFactory.getExecutor());
		EnumSet<ConservationMethod> set = EnumSet
				.allOf(ConservationMethod.class);

		norm_results = scores.calculateScores(set);
		scores = Conservation.getConservation(sequences, false,
				ExecutorFactory.getExecutor());
		results = scores.calculateScores(set);
		// Shutdown the executor and complete submitted tasks
		ExecutorFactory.getExecutor().shutdown();
		try {
			ExecutorFactory.getExecutor().awaitTermination(5, TimeUnit.SECONDS);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	}

	@Test()
	public void testSingleMethods() {
		// Re-initiate executor with multiple threads
		ExecutorFactory.initExecutor(0);
		Conservation cons;
		try {
			cons = Conservation.getConservation(input, false,
					ExecutorFactory.getExecutor());
			double[] result = cons.calculateScore(ConservationMethod.VALDAR);

			cons = Conservation.getConservation(input, false,
					ExecutorFactory.getExecutor());
			double[] result2 = cons.calculateScore(ConservationMethod.VALDAR);

			Assert.assertEquals(results.get(ConservationMethod.VALDAR).length,
					result.length);

			System.out.println(Arrays.toString(result));
			// System.out.println(Arrays.toString(result2));
			// System.out.println(Arrays.toString(results.get(Method.VALDAR)));

			Assert.assertTrue(Arrays.equals(result2, result));

		} catch (FileNotFoundException e) {
			e.printStackTrace();
			fail(e.getMessage());
		} catch (IOException e) {
			e.printStackTrace();
			fail(e.getMessage());
		} catch (UnknownFileFormatException e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
	}

	@Test()
	public void testMethods() {
		try {
			// Re-initiate executor with multiple threads
			ExecutorFactory.initExecutor(0);

			Conservation c = Conservation.getConservation(input, false,
					ExecutorFactory.getExecutor());
			Map<ConservationMethod, double[]> apiresults = c
					.calculateScores(EnumSet.allOf(ConservationMethod.class));
			assertNotNull(results);
			for (ConservationMethod method : apiresults.keySet()) {
				if (method == ConservationMethod.LANDGRAF) {
					// Landgrap results never repeats as they have random
					// element
					Assert.assertNotNull(apiresults.get(method));
					continue;
				}
				double[] result = results.get(method);
				double[] apiresult = apiresults.get(method);
				Assert.assertNotNull(apiresult);
				System.out.println(Arrays.toString(result));
				System.out.println(Arrays.toString(apiresult));

				Assert.assertTrue(Arrays.equals(apiresult, result),
						"Methods results: " + method.toString()
								+ " is not equal!");
			}
			// c.printResults(Format.RESULT_NO_ALIGNMENT);
			// System.out.println(results);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			fail(e.getMessage());
		} catch (IOException e) {
			e.printStackTrace();
			fail(e.getMessage());
		} catch (UnknownFileFormatException e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
	}

	@Test()
	public void testLoadClustalAl() {
		try {
			// Re-initiate executor with multiple threads
			ExecutorFactory.initExecutor(0);

			input = new File(SlowMethodTester.DATA_PATH + File.separator
					+ "avg4.clustal.fa");

			Conservation c = Conservation.getConservation(input, false,
					ExecutorFactory.getExecutor());
			Map<ConservationMethod, double[]> apiresults = c
					.calculateScores(EnumSet.allOf(ConservationMethod.class));
			assertNotNull(results);
			for (ConservationMethod method : apiresults.keySet()) {
				if (method == ConservationMethod.LANDGRAF) {
					// Landgrap results never repeats as they have random
					// element
					Assert.assertNotNull(apiresults.get(method));
					continue;
				}
				double[] result = results.get(method);
				double[] apiresult = apiresults.get(method);
				Assert.assertNotNull(apiresult);

				Assert.assertTrue(Arrays.equals(apiresult, result),
						"Methods results: " + method.toString()
								+ " is not equal!");
			}
			// c.printResults(Format.RESULT_NO_ALIGNMENT);
			// System.out.println(results);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			fail(e.getMessage());
		} catch (IOException e) {
			e.printStackTrace();
			fail(e.getMessage());
		} catch (UnknownFileFormatException e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
	}

	@Test()
	public void testCustomSMERFS() {
		try {
			// Re-initiate executor with multiple threads
			ExecutorFactory.initExecutor(0);

			Conservation c = Conservation.getConservation(input, false,
					ExecutorFactory.getExecutor());
			double[] apiresults = c.getSMERFS(11, SMERFSConstraints.MID_SCORE,
					0.23);
			assertNotNull(results);
			assertNotNull(apiresults);
			Assert.assertEquals(apiresults.length,
					results.get(ConservationMethod.SMERFS).length);
			Assert.assertFalse(Arrays.equals(apiresults,
					results.get(ConservationMethod.SMERFS)));

		} catch (FileNotFoundException e) {
			e.printStackTrace();
			fail(e.getMessage());
		} catch (IOException e) {
			e.printStackTrace();
			fail(e.getMessage());
		} catch (UnknownFileFormatException e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
	}

}
