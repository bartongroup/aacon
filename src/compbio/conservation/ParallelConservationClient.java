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

import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import compbio.data.sequence.FastaSequence;
import compbio.util.Timer;

/**
 * Command line client for AAconservation methods.
 * 
 * @author Peter Troshin with input from A. Golicz
 */
final class ParallelConservationClient {

	private final Map<Method, double[]> results = new EnumMap<Method, double[]>(
			Method.class);

	private static class SMERFSParams {

		private SMERFSColumnScore colScoreSchema = SMERFSColumnScore.MID_SCORE;
		private double SMERFSGapTreshold = SMERFSColumnScore.DEFAULT_GAP_THRESHOLD;
		private int SMERFSWidth = SMERFSColumnScore.DEFAULT_WINDOW_SIZE;

		/**
		 * 
		 * @param SMERFSargs
		 * @throws IllegalArgumentException
		 */
		SMERFSParams(String[] SMERFSargs) {
			parseArguments(SMERFSargs);
		}

		private void parseArguments(String[] SMERFSargs) {
			try {
				if (SMERFSargs != null) {
					SMERFSWidth = Integer.parseInt(SMERFSargs[0]);
					colScoreSchema = SMERFSColumnScore
							.getSMERFSColumnScore(SMERFSargs[1]);
					SMERFSGapTreshold = Double.parseDouble(SMERFSargs[2]);
				}
			} catch (NumberFormatException e) {
				throwIllegalSMERFSParamException();
			} catch (ArrayIndexOutOfBoundsException e) {
				throwIllegalSMERFSParamException();
			}
		}

		private void throwIllegalSMERFSParamException() {
			throw new IllegalArgumentException(
					"To run SMERFS three arguments are"
							+ " needed, window width, how to give "
							+ "scores to columns and a gap treshold.");
		}

	}

	ParallelConservationClient(String[] cmd) throws IOException,
			InterruptedException {

		Timer timer = Timer.getMilliSecondsTimer();

		String statFile = CmdParser.getStatFilePath(cmd);
		if (statFile == null) {
			timer.setStatOutput(null);
		} else {
			timer.setStatOutput(new FileOutputStream(statFile));
		}

		Set<Method> methods = CmdParser.getMethodNames(cmd);
		// assume all methods are required
		if (methods.isEmpty()) {
			methods = EnumSet.allOf(Method.class);
			timer.println("No methods are request assuming all are required.");
		}

		String inFilePath = CmdParser.getInputFilePath(cmd);

		if (!methods.isEmpty() && inFilePath != null) {
			String format = CmdParser.getFormat(cmd);
			String outFilePath = CmdParser.getOutputFilePath(cmd);
			if (outFilePath == null) {
				timer
						.println("No output file is provided, writing results to the standard output.");
			}
			Format outFormat = Format.RESULT_NO_ALIGNMENT;
			if (format != null) {
				Format userFormat = Format.getFormat(format);
				if (userFormat == null) {
					timer
							.println("Cannot recognise format '" + format
									+ "' Assuming "
									+ Format.RESULT_NO_ALIGNMENT.toString()
									+ " format");
				} else {
					outFormat = userFormat;
					timer.println("Setting output format to " + userFormat);
				}
			} else {
				timer
						.println("No format is provided assuming RESULT_NO_ALIGNMENT is required");
			}

			String[] SMERFSDetails = CmdParser.getSMERFSDetails(cmd);
			// This will throw en exception if parameters supplied but not valid
			SMERFSParams smerfsPar = new SMERFSParams(SMERFSDetails);

			boolean normalize = CmdParser.getNormalize(cmd);
			String[] gap = CmdParser.getGapChars(cmd);
			Character[] gapChars = CmdParser.extractGapChars(gap);

			ExecutorFactory.initExecutor(CmdParser.getThreadNumber(cmd), timer
					.getStatWriter());
			ExecutorService executor = ExecutorFactory.getExecutor();

			List<FastaSequence> sequences = CmdParser
					.openInputStream(inFilePath);

			if (sequences != null) {
				AminoAcidMatrix alignment = new AminoAcidMatrix(sequences,
						gapChars);
				timer.println("Start time: " + CmdParser.getDateTime());
				timer.println("Alignment loaded in: " + timer.getStepTime()
						+ " ms");
				timer.println("Alignment has: " + alignment.numberOfRows()
						+ " sequences.");

				Conservation scores = new Conservation(alignment, normalize,
						ExecutorFactory.getExecutor());

				MethodWrapper wrapper = null;
				List<MethodWrapper> tasks = new ArrayList<MethodWrapper>();

				for (Method method : methods) {
					// Start SMERFS from the main thread, as it has
					// its own means of dividing the tasks and executing in
					// parallel
					if (method == Method.SMERFS) {
						double[] result = runSMERFS(scores, smerfsPar, timer);
						results.put(Method.SMERFS, result);
						continue;
					}
					// Start other methods capable of runing in multiple threads
					// from the main thread.
					if (method == Method.LANDGRAF || method == Method.SANDER
							|| method == Method.KARLIN
							|| method == Method.VALDAR) {
						double[] result = runParallelMethod(scores, method,
								timer);
						results.put(method, result);
						continue;
					}
					wrapper = new MethodWrapper(method, scores, timer);
					tasks.add(wrapper);
				}
				List<Future<MethodWrapper>> rawResults = executor
						.invokeAll(tasks);
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
				executor.shutdown();

				ConservationFormatter.formatResults(results, outFilePath,
						outFormat, alignment);

				timer.println("Total calculation time: "
						+ timer.getTotalTime(TimeUnit.SECONDS) + " s");
				timer.println("End time: " + CmdParser.getDateTime());
				timer.getStatWriter().close();
			} else {
				System.out.println("No input found in " + inFilePath
						+ " ! Exiting");
			}
		}

	}

	static void checkArguments(String[] args) {
		if (args == null) {
			System.out.println("No parameters were suppled");
			System.out.println();
			System.out.print(CmdParser.CONSERVATION_HELP);
			System.exit(0);
		}
		if (args.length < 2) {
			System.out
					.println("Method names, input file paths are required. Application will"
							+ " not run until these 2 arguments are provided.");
			System.out
					.println("If you want results printed, both format an input "
							+ "file path have to be provided");
			System.out.println();
			System.out.print(CmdParser.CONSERVATION_HELP);
			System.exit(0);
		}
	}

	private double[] runSMERFS(Conservation scores, SMERFSParams sparams,
			Timer timer) {
		timer.getStepTime();
		double[] conservation = scores.getSMERFS(sparams.SMERFSWidth,
				sparams.colScoreSchema, sparams.SMERFSGapTreshold);
		timer.println(Method.SMERFS.toString() + " " + timer.getStepTime()
				+ " ms");
		return conservation;
	}

	private double[] runParallelMethod(Conservation scores, Method method,
			Timer timer) {
		timer.getStepTime();
		double[] results = scores.calculateScore(method);
		timer.println(method.toString() + " " + timer.getStepTime() + " ms");
		return results;
	}

	public static void main(String[] args) {

		checkArguments(args);
		try {
			ParallelConservationClient cons = new ParallelConservationClient(
					args);
		} catch (IOException e) {
			System.err.println("Fail to write to the file system! "
					+ e.getLocalizedMessage());
			e.printStackTrace();
		} catch (InterruptedException e) {
			System.err.println("Interrupted!");
			e.printStackTrace();
		}
	}
}
