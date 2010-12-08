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
import java.util.Locale;
import java.util.Map;

import compbio.data.sequence.Method;

/**
 * 
 * @author Agnieszka Golicz & Peter Troshin
 * 
 */
final class ConservationFormatter {

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

		assert format != null : "Format must not be null";
		assert scores != null : "Scores must not be null";
		OutputStream out = openPrintWriter(outFilePath);

		switch (format) {
		case RESULT_NO_ALIGNMENT:
			outputScoreLine(scores, out);
			break;
		case RESULT_WITH_ALIGNMENT:
			if (alignment == null) {
				throw new NullPointerException("Alignment must be provided!");
			}
			alignment.printAlignment(out);
			outputScoreLine(scores, out);
			break;
		}
		out.flush();
		out.close();
	}

	static void outputScoreLine(Map<Method, double[]> scores,
			OutputStream outStream) {
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
