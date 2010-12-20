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
import compbio.data.sequence.ConservationMethod;
import compbio.data.sequence.FastaSequence;
import compbio.data.sequence.SMERFSConstraints;

/**
 * A public API for conservation calculation methods.
 * 
 * Input in all methods are <br/>
 * 1) The list of FastaSequence objects holding the aligned sequences. All
 * sequences must be of the same length (as this is an alignment). The input can
 * be conveniently loaded from either the Clustal formatted alignment file or
 * the file containing a list of FASTA formatted sequences using the following
 * method:
 * 
 * <pre>
 * {@code
 * 	List<FastaSequence> sequences = CmdParser
 * 			.openInputStream(<PATH_TO_INPUT_FILE>);
 * }
 * </pre>
 * 
 * 2) The boolean parameters telling the system whether the results should be
 * normalized or not. Normalized results have values between 0 and 1. Please
 * note however, that some results cannot be normalized. In such a case, the
 * system returns not normalized values, and log the issue to the standard error
 * stream. The following formula is used for normalization n = (d - dmin)/(dmax
 * - dmin) Negative results first converted to positive by adding an absolute
 * value of the most negative result. <br/>
 * 3) The ExecutorService object which is used to parallel the calculations. The
 * ExecutorService initialized with the number of threads equals to the number
 * of cores on the executing machine can be obtained as follows: <br/>
 * 
 * <pre>
 * {@code int corenum = Runtime.getRuntime().availableProcessors();
 *   ExecutorService executor = Executors.newFixedThreadPool(corenum);
 *   
 *   ........DO CALCULATIONS.........
 *   
 *   // shutdown the Executor 
 *   executor.shutdown(); 
 *   }
 * </pre>
 * 
 * 
 * Please take care to initialize and pass only one executor to all the methods
 * to avoid the waist of resources. After use, the executor must be disposed of,
 * it can be done as follows: {@code executor.shutdown(); }
 * 
 * An example use of this class for calculating conservation is below: <br/>
 * 
 * <pre>{@code
 
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
 } 
  </pre>
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
	 *            but SMERFS can be used. The {@link EnumSet} class provides a
	 *            number of convenience methods which can be used to prepare a
	 *            set of methods for the input. For example, to use
	 *            {@link Method#KABAT} for conservation calculation one could
	 *            construct a set in the following way:
	 *            {@code EnumSet.of(Method.KABAT)}.<br/>
	 * 
	 *            For all the methods, but SMERFS use the following construct:
	 * 
	 * <br/> {@code EnumSet.complementOf(EnumSet.of(Method.SMERFS))} <br/>
	 * 
	 *            For a set of methods including KABAT, JORES, SCHNEIDER,
	 *            SHENKIN, GERSTEIN use the following construct:
	 * 
	 * <br/> {@code EnumSet.range(Method.KABAT, Method.GERSTEIN)}<br/>
	 * 
	 *            For further details see {@link EnumSet}
	 * 
	 * @param executor
	 *            the {@link ExecutorService} to be used to parallelize the
	 *            calculations
	 * @return the Map of Method->double[] conservation scores
	 * @throws InterruptedException
	 *             if the calculating Thread was interrupted
	 */
	public static Map<ConservationMethod, double[]> getConservation(
			List<FastaSequence> alignment, boolean normalize,
			Set<ConservationMethod> methods, ExecutorService executor)
			throws InterruptedException {
		final Map<ConservationMethod, double[]> results = new EnumMap<ConservationMethod, double[]>(
				ConservationMethod.class);

		// Make a matrix out of the alignment
		AminoAcidMatrix alignMatrix = new AminoAcidMatrix(alignment, null);

		Conservation scores = new Conservation(alignMatrix, normalize, executor);
		List<MethodWrapper> tasks = new ArrayList<MethodWrapper>();

		for (ConservationMethod method : methods) {
			// Start SMERFS from the main thread, as it has
			// its own means of dividing the tasks and executing in
			// parallel
			if (method == ConservationMethod.SMERFS) {
				double[] conservation = scores.getSMERFS(
						SMERFSConstraints.DEFAULT_WINDOW_SIZE,
						SMERFSConstraints.MID_SCORE,
						SMERFSConstraints.DEFAULT_GAP_THRESHOLD);
				results.put(ConservationMethod.SMERFS, conservation);
				continue;
			}
			// Start other methods capable of running in multiple threads
			// from the main thread.
			if (method == ConservationMethod.LANDGRAF
					|| method == ConservationMethod.SANDER
					|| method == ConservationMethod.KARLIN
					|| method == ConservationMethod.VALDAR) {
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
			List<MethodWrapper> tasks, Map<ConservationMethod, double[]> results)
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
			SMERFSConstraints scoringMethod, float gapTreshold,
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

	/*
	 * TODO remove public static void main(String[] args) throws
	 * InterruptedException, IOException, UnknownFileFormatException { int
	 * corenum = Runtime.getRuntime().availableProcessors(); // Initialize the
	 * Executor instance with a number of cores ExecutorService executor =
	 * Executors.newFixedThreadPool(corenum); // Load the data from the file
	 * containing either Clustal formatted // alignment // or a list of FASTA
	 * formatted sequences. Assuming that small.align // file is // in the same
	 * directory as this program List<FastaSequence> sequences = CmdParser
	 * .openInputStream("test/data/small.align"); // Calculate conservation
	 * scores using all the methods available. Map<Method, double[]> result =
	 * getConservation(sequences, true,
	 * EnumSet.complementOf(EnumSet.of(Method.SMERFS)), executor); // Print the
	 * results to the console. FileOutputStream outfile = new
	 * FileOutputStream("results.txt");
	 * ConservationFormatter.formatResults(result, outfile); outfile.close();
	 * ConservationFormatter.formatResults(result, "test.txt",
	 * Format.RESULT_NO_ALIGNMENT, sequences); executor.shutdown(); }
	 */
}
