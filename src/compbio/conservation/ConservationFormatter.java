package compbio.conservation;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.Iterator;
import java.util.Locale;
import java.util.Map;

class ConservationFormatter {

	/**
	 * Number if digits after the comma to print to results output
	 */
	public static final int PRECISION = 3;

	public static final NumberFormat NUMBER_FORMAT = NumberFormat
			.getNumberInstance(Locale.UK);

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
	 * 
	 *            private static <T> void formatResult(T tag, double[] result,
	 *            int resultPrecision, PrintWriter print) {
	 * 
	 *            String tagFormat = "%s"; String resultFormat = "%." +
	 *            resultPrecision + "f"; print.printf(tagFormat, "#" +
	 *            tag.toString()); for (int i = 0; i < result.length; i++) {
	 *            print.printf(resultFormat, result[i]); print.print(" "); }
	 *            print.println(); }
	 */

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
	 * 
	 *            private static <T> void formatResultWithAlignment(T tag,
	 *            double[] result, int tagWidth, int resultWidth, int
	 *            resultPrecision, PrintWriter print) {
	 * 
	 *            String tagFormat = "%-" + tagWidth + "s"; String resultFormat
	 *            = "%-" + resultWidth + "." + resultPrecision + "f";
	 *            print.printf(tagFormat, "#" + tag.toString()); for (int i = 0;
	 *            i < result.length; i++) { print.printf(resultFormat,
	 *            result[i]); } print.println(); }
	 */
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
	 * 
	 *             static <T> void printResultWithAlignment(AminoAcidMatrix
	 *             alignment, T tag, double[] result, int tagWidth, int
	 *             resultWidth, int resultPrecision, String outputFile) throws
	 *             IOException {
	 * 
	 *             PrintWriter print = openPrintWriter(outputFile, true);
	 * 
	 *             formatResultWithAlignment(tag, result, tagWidth, resultWidth,
	 *             resultPrecision, print); print.close(); }
	 */

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
	 * 
	 *             static <T> void printResultNoAlignment(AminoAcidMatrix
	 *             alignment, T tag, double[] result, int resultPrecision,
	 *             String outputFile, boolean append) throws IOException {
	 * 
	 *             PrintWriter print = openPrintWriter(outputFile, append);
	 *             formatResult(tag, result, resultPrecision, print);
	 *             print.close(); }
	 */
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

		Iterator<Method> itr = scores.keySet().iterator();

		switch (format) {
		case RESULT_NO_ALIGNMENT:
			outputScoreLine(scores, out, itr);
			break;
		case RESULT_WITH_ALIGNMENT:
			if (alignment == null) {
				throw new NullPointerException("Alignment must be provided!");
			}
			alignment.printAlignment(out);
			outputScoreLine(scores, out, itr);
			break;
		}
		out.flush();
		out.close();
	}

	private static void outputScoreLine(Map<Method, double[]> scores,
			OutputStream outStream, Iterator<Method> itr) {
		PrintWriter print = new PrintWriter(outStream);
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
