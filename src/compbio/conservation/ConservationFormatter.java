package compbio.conservation;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.Map;

public class ConservationFormatter {

	/**
	 * Number if digits after the comma to print to results output
	 */
	public static final int PRECISION = 3;

	/**
	 * Formats results
	 * 
	 * @param <T>
	 * @param tag
	 *            any object(usually used for members of enumeration)
	 * @param result
	 *            array of results to be printed
	 * @param resultPrecision
	 * @param print
	 *            reference to PrintWriter object
	 */
	private static <T> void formatResult(T tag, double[] result,
			int resultPrecision, PrintWriter print) {

		String tagFormat = "%s";
		String resultFormat = "%." + resultPrecision + "f";
		print.printf(tagFormat, "#" + tag.toString());
		for (int i = 0; i < result.length; i++) {
			print.printf(resultFormat, result[i]);
			print.print(" ");
		}
		print.println();
	}

	/**
	 * @param <T>
	 * @param tag
	 *            any object(usually used for members of enumeration)
	 * @param result
	 *            array of results to be printed
	 * @param tagWidth
	 * @param resultWidth
	 * @param resultPrecision
	 * @param print
	 *            reference to PrintWriter object
	 */
	private static <T> void formatResultWithAlignment(T tag, double[] result,
			int tagWidth, int resultWidth, int resultPrecision,
			PrintWriter print) {

		String tagFormat = "%-" + tagWidth + "s";
		String resultFormat = "%-" + resultWidth + "." + resultPrecision + "f";
		print.printf(tagFormat, "#" + tag.toString());
		for (int i = 0; i < result.length; i++) {
			print.printf(resultFormat, result[i]);
		}
		print.println();
	}

	/**
	 * @param <T>
	 * @param alignment
	 * @param tag
	 *            any object(usually used for members of enumeration)
	 * @param result
	 * @param tagWidth
	 * @param resultWidth
	 * @param resultPrecision
	 * @param outputFile
	 * @throws IOException
	 */
	static <T> void printResultWithAlignment(AminoAcidMatrix alignment, T tag,
			double[] result, int tagWidth, int resultWidth,
			int resultPrecision, String outputFile) throws IOException {

		PrintWriter print = openPrintWriter(outputFile, true);

		formatResultWithAlignment(tag, result, tagWidth, resultWidth,
				resultPrecision, print);
		print.close();
	}

	/**
	 * @param <T>
	 * @param alignment
	 * @param tag
	 *            any object(usually used for members of enumeration)
	 * @param result
	 * @param resultPrecision
	 * @param outputFile
	 * @param append
	 * @throws IOException
	 */
	static <T> void printResultNoAlignment(AminoAcidMatrix alignment, T tag,
			double[] result, int resultPrecision, String outputFile,
			boolean append) throws IOException {

		PrintWriter print = openPrintWriter(outputFile, append);
		formatResult(tag, result, resultPrecision, print);
		print.close();
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

		PrintWriter print = null;
		assert format!=null : "Format must not be null"; 
		assert scores!=null : "Scores must not be null"; 
		
		Iterator<Method> itr = scores.keySet().iterator();

		switch (format) {
		case RESULT_NO_ALIGNMENT:
			print = openPrintWriter(outFilePath, false);
			if (print != null) {
				while (itr.hasNext()) {
					Method key = itr.next();
					print.print("#" + key.toString() + " ");
					ConservationAccessory.printArrayOfDouble(scores.get(key),
							print, PRECISION);
				}
			}
			break;
		case RESULT_WITH_ALIGNMENT:
			if (alignment == null) {
				throw new NullPointerException("Alignment must be provided!");
			}
			print = openPrintWriter(outFilePath, true);
			int tagWidth = 30;
			int fieldWidth = 10;
			String tagFormat = "%-" + tagWidth + "s";
			alignment.printAlignment(tagWidth, fieldWidth, outFilePath);
			if (print != null) {
				print.println();
				while (itr.hasNext()) {
					Method key = itr.next();
					print.printf(tagFormat, "#" + key.toString() + " ");
					ConservationAccessory.printArrayOfDoubleFieldWidth(scores
							.get(key), print, PRECISION, fieldWidth);
				}
			}
			break;
		}
		print.flush();
		print.close();
	}

	static PrintWriter openPrintWriter(String outFilePath, boolean autoFlush)
			throws IOException {
		PrintWriter print = null;
		if (outFilePath == null || outFilePath.isEmpty()) {
			print = new PrintWriter(System.out);
		} else {
			print = new PrintWriter(new BufferedWriter(new FileWriter(
					outFilePath, autoFlush)));
		}
		return print;
	}

}
