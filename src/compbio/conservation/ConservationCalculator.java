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

import java.util.ArrayList;
import java.util.EnumMap;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import compbio.data.sequence.Alignment;
import compbio.data.sequence.FastaSequence;
import compbio.data.sequence.Method;

/**
 * A public API for conservation calculation methods.
 * 
 * Input in all methods are
 * 
 * 1) The list of FastaSequence objects holding the aligned sequences. All
 * sequences must be of the same length (as this is an alignment). The input can
 * be conveniently loaded from either the Clustal formatted alignment file or
 * the file containing a list of FASTA formatted sequences using the following
 * method:
 * 
 * <code>List<FastaSequence> sequences = CmdParser
				.openInputStream("PATH_TO_INPUT_FILE");
   </code>
 * 
 * 2) The boolean parameters telling the system whether the results should be
 * normalized or not. Normalized results have values between 0 and 1. Please
 * note however, that some results cannot be normalized. In such a case, the
 * system returns not normalized values, and log the issue to the standard error
 * stream. The following formula is used for normalization n = (d - dmin)/(dmax
 * - dmin) Negative results first converted to positive by adding an absolute
 * value of the most negative result.
 * 
 * 3) The ExecutorService object which is used to parallel the calculations. The
 * ExecutorService initialized with the number of threads equals to the number
 * of cores on the executing machine can be obtained as follows:
 * 
 * <code>
    int corenum = Runtime.getRuntime().availableProcessors();
	ExecutorService executor = Executors.newFixedThreadPool(corenum);
	
	........DO CALCULATIONS.........
	
	// shutdown the Executor
	executor.shutdown();
	
   </code>
 * 
 * Please take care to initialize and pass only one executor to all the methods
 * to avoid the waist of resources. After use, the executor must be disposed of,
 * it can be done as follows:
 * 
 * <code>
 	executor.shutdown(); 
  </code>
 * 
 * An example use of this class for calculating conservation is below:
 * 
 * <code>
 	// Determine the number of CPU cores available on the system.  
 	int corenum = Runtime.getRuntime().availableProcessors();
 	
 	// Initialize the Executor instance with a number of cores
	ExecutorService executor = Executors.newFixedThreadPool(corenum);

	// Load the data from the file containing either Clustal formatted alignment 
	// or a list of FASTA formatted sequences
	List<FastaSequence> sequences = CmdParser
				.openInputStream("test/data/small.align");
	
	// Calculate conservation scores using all the methods available.
	Map<Method, double[]> result = getConservation(sequences, true,
	EnumSet.complementOf(EnumSet.of(Method.SMERFS)), executor);
	
	// Print the results to the console.  
	ConservationFormatter.outputScoreLine(result, System.out);
 
  </code>
 * 
 * @author Peter Troshin
 * @version 1.0 October 2010
 */
public class ConservationCalculator {

	/**
	 * Calculates the conservation by all the methods defined by {@link Method}
	 * enumeration apart from SMERFS.
	 * 
	 * @param alignment
	 *            the list of FastaSequence objects holding the alignment. All
	 *            sequences must be of the same length.
	 * @param normalize
	 *            true if the resulting scores should be normalized, false
	 *            otherwise.
	 * 
	 * @param methods
	 *            the methods to be used for the calculation, all {@link Method}
	 *            but SMERFS can be used. {@link EnumSet} provides a number of
	 *            convenience methods which can be used to provide a set of
	 *            methods for an input. For example, to use {@link Method#KABAT}
	 *            for conservation calculation one could construct a set in the
	 *            following way: {@link EnumSet#of(Method#KABAT)}) for all the
	 *            methods, but SMERFS use the following construct:
	 * 
	 *            <code>EnumSet.complementOf(EnumSet.of(Method.SMERFS))</code>
	 * 
	 *            for a set of methods including KABAT, JORES, SCHNEIDER,
	 *            SHENKIN, GERSTEIN use the following construct:
	 * 
	 *            <code>EnumSet.range(Method.KABAT, Method.GERSTEIN)</code>
	 * 
	 * @param executor
	 *            the {@link ExecutorService} to be used to paralelize the
	 *            calculations
	 * @return the Map of Method->double[] conservation scores
	 * @throws InterruptedException
	 *             if the calculating Thread was interrupted
	 * @see EnumSet for further help on preparing Method sets
	 */
	public static Map<Method, double[]> getConservation(
			List<FastaSequence> alignment, boolean normalize,
			Set<Method> methods, ExecutorService executor)
			throws InterruptedException {
		final Map<Method, double[]> results = new EnumMap<Method, double[]>(
				Method.class);

		// Make a matrix out of the alignment
		AminoAcidMatrix alignMatrix = new AminoAcidMatrix(alignment, null);

		Conservation scores = new Conservation(alignMatrix, normalize, executor);
		List<MethodWrapper> tasks = new ArrayList<MethodWrapper>();

		for (Method method : methods) {
			// Start other methods capable of running in multiple threads
			// from the main thread.
			if (method == Method.LANDGRAF || method == Method.SANDER
					|| method == Method.KARLIN || method == Method.VALDAR) {
				double[] result = scores.calculateScore(method);
				results.put(method, result);
				continue;
			}
			MethodWrapper wrapper = new MethodWrapper(method, scores, null);
			tasks.add(wrapper);
		}
		waitResults(executor, tasks, results);
		return results;
	}

	private static void waitResults(ExecutorService executor,
			List<MethodWrapper> tasks, Map<Method, double[]> results)
			throws InterruptedException {
		List<Future<MethodWrapper>> rawResults = executor.invokeAll(tasks);
		for (Future<MethodWrapper> rawResult : rawResults) {
			MethodWrapper entry = null;
			try {
				entry = rawResult.get();
			} catch (ExecutionException e) {
				System.err.println("Exception while executing method");
				throw new RuntimeException(e.getCause());
			}
			results.put(entry.method, entry.conservation);
		}
	}

	/**
	 * Calculating the SMERFS score with custom parameters
	 * 
	 * @param alignment
	 *            the List of FastaSequence objects holding each sequence from
	 *            the alignment
	 * @param windowWidth
	 *            the window size parameter for SMERFS algorithm
	 * @param scoringMethod
	 *            the {@link SMERFSColumnScore}
	 * @param gapTreshold
	 *            the gap threshold for SMERFS algorithm
	 * @param normalize
	 *            the boolean value indicating whether the resulting score
	 *            should be normalized, true if it does.
	 * @param service
	 *            the {@link ExecutorService} to be used to parallel
	 *            calculations
	 * @return the array of double values - SMERFS scores for each position in
	 *         the alignment
	 */
	public static double[] getSMERFSScore(Alignment alignment, int windowWidth,
			SMERFSColumnScore scoringMethod, float gapTreshold,
			boolean normalize, ExecutorService service) {

		// Make a matrix out of the alignment
		AminoAcidMatrix alignMatrix = new AminoAcidMatrix(
				alignment.getSequences(), null);
		Conservation scores = new Conservation(alignMatrix, normalize, service);

		// Start SMERFS from the main thread, as it has
		// its own means of dividing the tasks and executing in
		// parallel
		return scores.getSMERFS(windowWidth, scoringMethod, gapTreshold);
	}

	/**
	 * Calculating the SMERFS score with custom parameters. This method uses
	 * default parameters for SMERFS algorithm.
	 * 
	 * @param alignment
	 *            the List of FastaSequence objects holding each sequence from
	 *            the alignment the gap threshold for SMERFS algorithm
	 * @param normalize
	 *            the boolean value indicating whether the resulting score
	 *            should be normalized, true if it does.
	 * @param service
	 *            the {@link ExecutorService} to be used to parallel
	 *            calculations
	 * @return the array of double values - SMERFS scores for each position in
	 *         the alignment
	 * 
	 * @see SMERFSColumnScore
	 */
	public static double[] getSMERFSScore(List<FastaSequence> alignment,
			boolean normalize, ExecutorService service) {

		// Make a matrix out of the alignment
		AminoAcidMatrix alignMatrix = new AminoAcidMatrix(alignment, null);
		Conservation scores = new Conservation(alignMatrix, normalize, service);

		return scores.getSMERFS(SMERFSColumnScore.DEFAULT_WINDOW_SIZE,
				SMERFSColumnScore.MID_SCORE,
				SMERFSColumnScore.DEFAULT_GAP_THRESHOLD);
	}

}
