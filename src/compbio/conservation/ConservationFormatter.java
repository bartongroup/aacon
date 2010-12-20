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
import java.io.OutputStream;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import compbio.data.sequence.FastaSequence;
import compbio.data.sequence.Method;
import compbio.data.sequence.SequenceUtil;

/**
 * 
 * Conservation calculation results writer. A helper class for outputting the
 * conservation results.
 * 
 * @author Agnieszka Golicz & Peter Troshin
 * 
 */
public final class ConservationFormatter {

	/**
	 * Number if digits after the comma to print to results output
	 */
	public static final int PRECISION = 3;

	public static final NumberFormat NUMBER_FORMAT = NumberFormat
			.getNumberInstance(Locale.UK);
	static {
		NUMBER_FORMAT.setGroupingUsed(false);
	}

	/**
	 * Formats and prints results
	 * 
	 * @param scores
	 * @param outFilePath
	 * @param format
	 * @throws IOException
	 */
	static void formatResults(Map<Method, double[]> scores, String outFilePath,
			Format format, AminoAcidMatrix alignment) throws IOException {
		if (alignment == null) {
			throw new NullPointerException("Alignment must be provided!");
		}
		formatResults(scores, outFilePath, format, alignment.getAlignment());
	}

	/**
	 * Formats and prints results. If you do not need to output the alignment
	 * please use the other method
	 * {@link ConservationFormatter#formatResults(Map, OutputStream)}
	 * 
	 * @param scores
	 *            the results on the calculation, the Map with Method keys and
	 *            double[] values.
	 * @param outFilePath
	 *            the path to the output file. Optional. If null System.out is
	 *            used.
	 * @param format
	 *            the output format. Optional. Defaults to
	 *            {@value Format#RESULT_NO_ALIGNMENT}
	 * @param alignment
	 *            the alignment for which score was calculated
	 * @throws IOException
	 *             if the file cannot be written/created for whatever reason.
	 * @throws NullPointerException
	 *             if the alignment is null
	 */
	public static void formatResults(Map<Method, double[]> scores,
			String outFilePath, Format format, List<FastaSequence> alignment)
			throws IOException {

		assert format != null : "Format must not be null";
		assert scores != null : "Scores must not be null";
		if (scores == null) {
			return;
		}
		if (format == null) {
			// default to results with no alignment
			format = Format.RESULT_NO_ALIGNMENT;
		}
		if (alignment == null) {
			throw new NullPointerException("Alignment must be provided");
		}
		// If outFile is null default to System.out
		OutputStream out = openPrintWriter(outFilePath);
		switch (format) {
		case RESULT_NO_ALIGNMENT:
			formatResults(scores, out);
			break;
		case RESULT_WITH_ALIGNMENT:
			SequenceUtil.writeFastaKeepTheStream(out, alignment, 80);
			formatResults(scores, out);
			break;
		}
		out.flush();
		out.close();
	}

	/**
	 * Use this method to save the conservation results into a file. The call to
	 * this method produce the same results as
	 * {@link ConservationFormatter#formatResults(Map, String, Format, List)}
	 * where the first parameter is scores, the second is a name of the file,
	 * third if {@value Format#RESULT_NO_ALIGNMENT} and the fourth is null.
	 * 
	 * @param scores
	 *            the results of the calculation - the Map<Method,double[]>
	 * @param outStream
	 *            the stream to write the results to. Please note that the
	 *            method leaves this stream open. It is up to the caller to
	 *            close it!
	 */
	public static void formatResults(Map<Method, double[]> scores,
			OutputStream outStream) {
		if (outStream == null) {
			throw new NullPointerException("Output stream must be provided!");
		}
		PrintWriter print = new PrintWriter(outStream);
		Iterator<Method> itr = scores.keySet().iterator();
		while (itr.hasNext()) {
			Method key = itr.next();
			print.print("#" + key.toString() + " ");
			ConservationAccessory.printArrayOfDouble(scores.get(key), print,
					PRECISION);
		}
		print.println();
		print.flush();
	}

	static OutputStream openPrintWriter(String outFilePath) throws IOException {
		OutputStream outstream = null;
		if (outFilePath == null || outFilePath.isEmpty()) {
			outstream = System.out;
		} else {
			outstream = new FileOutputStream(outFilePath);
		}
		return outstream;
	}
}
