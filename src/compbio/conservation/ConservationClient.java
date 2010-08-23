package compbio.conservation;

import java.io.FileOutputStream;
import java.io.IOException;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import compbio.data.sequence.FastaSequence;
import compbio.util.Timer;

class ConservationClient {

	private final Map<Method, double[]> results = new EnumMap<Method, double[]>(
			Method.class);

	/**
	 * Returns results of SMERFS calculation or null, if parameters provided are
	 * not appropriate.
	 * 
	 * @param alignment
	 *            reference to alignment
	 * @param width
	 *            with of the window
	 * @param score
	 *            tells which score given to the column, either the highest
	 *            score of all the windows it belongs to, or the middle column
	 *            is given the score of the window.
	 * @param normalize
	 *            if true results will be normalized
	 * @return
	 */
	public static double[] getSMERFS(AminoAcidMatrix alignment, int width,
			SMERFSColumnScore score, double gapTreshold, boolean normalize) {

		if (alignment == null) {
			throw new IllegalArgumentException("Matrix must not be null.");
		}
		double[] result = null;
		if (width <= 0 || width % 2 != 1 || width > alignment.numberOfColumns()
				|| score == null || gapTreshold < 0 || gapTreshold > 1) {
			if (width <= 0 || width % 2 != 1) {
				System.out
						.println("Column width for SMERFS not provided or smaller or equal zero or not an odd number or not an integer.");
			}
			if (width > alignment.numberOfColumns()) {
				System.out
						.println("Column width greater than the length of the alignment");
			}
			if (score == null) {
				System.out
						.println("Column score not privided or the type provided is not supported.");
				SMERFSColumnScore.supportedSMERFSColumnSores();
			}
			if (gapTreshold < 0 || gapTreshold > 1) {
				System.out
						.println("Gap treshold could not have been parsed as a double or it was smaller than zero or it was greater than one.");
			}
			return result;
		}
		Correlation corr = new Correlation(alignment, width, gapTreshold);
		result = corr.getCorrelationScore(score, normalize);
		return result;
	}

	/**
	 * Constructor
	 * 
	 * @param cmd
	 *            command line arguments
	 * @throws IOException
	 */
	ConservationClient(String[] cmd) throws IOException {

		Timer timer = Timer.getMilliSecondsTimer();

		// dafault values
		int SMERFSWidth = SMERFSColumnScore.DEFAULT_WINDOW_SIZE;
		SMERFSColumnScore score = SMERFSColumnScore.MID_SCORE;
		double SMERFSGapTreshold = SMERFSColumnScore.DEFAULT_GAP_THRESHOLD;

		Set<Method> methods = CmdParser.getMethodNames(cmd);

		String inFilePath = CmdParser.getInputFilePath(cmd);

		if (methods != null && inFilePath != null) {
			String format = CmdParser.getFormat(cmd);
			String outFilePath = CmdParser.getOutputFilePath(cmd);

			Format outFormat = Format.getFormat(format);
			String[] SMERFSDetails = CmdParser.getSMERFSDetails(cmd);
			if (SMERFSDetails != null) {
				if (SMERFSDetails.length == 3) {
					try {
						SMERFSWidth = Integer.parseInt(SMERFSDetails[0]);
					} catch (NumberFormatException e) {
						SMERFSWidth = -1;
					}
					score = SMERFSColumnScore
							.getSMERFSColumnScore(SMERFSDetails[1]);
					try {
						SMERFSGapTreshold = Double
								.parseDouble(SMERFSDetails[2]);
					} catch (NumberFormatException e) {
						SMERFSGapTreshold = -1;
					}
				} else {
					System.err.println("To run SMERFS three arguments are"
							+ " needed, window width, how to give "
							+ "scores to columns and a gap treshold.");
					System.exit(1);
				}
			}
			boolean normalize = CmdParser.getNormalize(cmd);
			String[] gap = CmdParser.getGapChars(cmd);
			Character[] gapChars = CmdParser.extractGapChars(gap);

			String statFile = CmdParser.getStatFilePath(cmd);
			if (statFile == null) {
				timer.setStatOutput(null);
			} else {
				timer.setStatOutput(new FileOutputStream(statFile));
			}

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

				ConservationScores2 scores = new ConservationScores2(alignment);
				double[] result = null;

				for (Method method : methods) {

					if (method == Method.SMERFS && SMERFSDetails != null) {
						result = getSMERFS(alignment, SMERFSWidth, score,
								SMERFSGapTreshold, normalize);
					} else {
						result = scores.calculateScore(method, normalize);
					}
					results.put(method, result);
					timer.println(method.toString() + " done "
							+ timer.getStepTime() + " ms");
				}

				ConservationFormatter.formatResults(results, outFilePath,
						outFormat, alignment);

				timer.println("End time: " + CmdParser.getDateTime());
				timer.getStatWriter().close();
			} else {
				System.out.println("No input found in " + inFilePath
						+ " ! Exiting");
			}
		}

	}

	Map<Method, double[]> getResults() {
		return results;
	}

	/**
	 * Application entry point. Command line format looks like.
	 * -m=method1,method2,method3 - method names -n normalization factor
	 * -f=format -i=inputPath - inputfile path -o=outputPath - outputfile path
	 * 
	 * @param args
	 *            command line arguments
	 */
	public static void main(String[] args) {

		if (args == null) {
			System.out.println("No parameters were suppled");
			System.out.println();
			System.out.print(CmdParser.CONSERVATION_HELP);
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
		}
		try {
			ConservationClient cons = new ConservationClient(args);
		} catch (IOException e) {
			System.err.println("Fail to write to the file system! "
					+ e.getLocalizedMessage());
			e.printStackTrace();
		}
	}

}
