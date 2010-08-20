package compbio.conservation;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;

import compbio.common.IllegalGapCharacterException;
import compbio.data.sequence.FastaSequence;
import compbio.data.sequence.SequenceUtil;

class FamiliesConservationClient {

	// private final Map<Method, double[]> results = new EnumMap<Method,
	// double[]>(Method.class);
	ColumnInfo[][] info1;
	ColumnInfo[][] info2;
	final static String pseparator = "=";
	final static String methodKey = "-m";
	final static String normalizationKey = "-n";
	final static String inputKey = "-i";
	final static String outputKey = "-o";
	final static String SMERFSDetailsKey = "-s";
	final static String gapKey = "-g";
	final static String statKey = "-d";
	final static String groupsKey = "-r";

	/**
	 * Gets method name from the command line
	 * 
	 * @param cmd
	 *            array of cmd arguments
	 * @return method name or null if no method name provided
	 */
	private static String getMethodName(String[] cmd) {

		for (int i = 0; i < cmd.length; i++) {
			String meths = cmd[i];
			if (meths.trim().toLowerCase().startsWith(methodKey + pseparator)) {
				return meths.substring(meths.indexOf(pseparator) + 1);
			}
		}
		return null;
	}

	/**
	 * Gets method name from the command line
	 * 
	 * @param cmd
	 *            array of cmd arguments
	 * @return method name or null if no method name provided
	 */
	private static String[] getGroups(String[] cmd) {

		for (int i = 0; i < cmd.length; i++) {
			String meths = cmd[i];
			if (meths.trim().toLowerCase().startsWith(groupsKey + pseparator)) {
				return meths.substring(meths.indexOf(pseparator) + 1)
						.split(";");
			}
		}
		return null;
	}

	private static int[][] parseGroups(String[] groups, Boolean proceed) {

		int[][] groupDouble = new int[groups.length][];
		for (int i = 0; i < groups.length; i++) {
			String[] gr = groups[i].split(",");
			groupDouble[i] = new int[gr.length];
			for (int j = 0; j < groupDouble[i].length; j++) {
				int d = -1;
				try {
					d = Integer.parseInt(gr[j]);
					groupDouble[i][j] = d;
				} catch (NumberFormatException e) {
					System.out
							.println("Number could not have buun parsed as double");
					proceed = false;
				}
				if (d < 0) {
					System.out.println("Column out of range");
					proceed = false;
				}
			}
		}
		return groupDouble;
	}

	/**
	 * Gets method name from the command line
	 * 
	 * @param cmd
	 *            array of cmd arguments
	 * @return method name or null if no method name provided
	 */
	private static String[] getSMERFSDetails(String[] cmd) {

		for (int i = 0; i < cmd.length; i++) {
			String meths = cmd[i];
			if (meths.trim().toLowerCase().startsWith(
					SMERFSDetailsKey + pseparator)) {
				return meths.substring(meths.indexOf(pseparator) + 1)
						.split(",");
			}
		}
		return null;
	}

	/**
	 * Gets the normalization status from the command line.
	 * 
	 * @param cmd
	 *            array of cmd arguments
	 * @return true if results to be normalized false if else, returns false if
	 *         no normalization status provided
	 */
	private static boolean getNormalize(String[] cmd) {

		for (int i = 0; i < cmd.length; i++) {
			String norm = cmd[i];
			if (norm.trim().toLowerCase().equals(normalizationKey)) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Gets output file path
	 * 
	 * @param cmd
	 * @return null if not provided
	 */
	private static String getOutputFilePath(String[] cmd) {

		for (int i = 0; i < cmd.length; i++) {
			String name = cmd[i];
			if (name.trim().toLowerCase().startsWith(outputKey + pseparator)) {
				return name.substring(name.indexOf(pseparator) + 1);
			}
		}
		return null;
	}

	/**
	 * Input file path.
	 * 
	 * @param cmd
	 * @return null if not provided
	 */
	private static String getInputFilePath(String[] cmd) {

		for (int i = 0; i < cmd.length; i++) {
			String name = cmd[i];
			if (name.trim().toLowerCase().startsWith(inputKey + pseparator)) {
				return name.substring(name.indexOf(pseparator) + 1);
			}
		}
		return null;
	}

	/**
	 * Gets statistics file path
	 * 
	 * @param cmd
	 * @return null if not provided
	 */
	private static String getStatFilePath(String[] cmd) {

		for (int i = 0; i < cmd.length; i++) {
			String name = cmd[i];
			if (name.trim().toLowerCase().startsWith(statKey + pseparator)) {
				return name.substring(name.indexOf(pseparator) + 1);
			}
		}
		return null;
	}

	/**
	 * Gets the gap character .
	 * 
	 * @param cmd
	 *            array of cmd arguments
	 * @return gap character or null if gap character not provided
	 */
	private static String[] getGapChars(String[] cmd) {

		for (int i = 0; i < cmd.length; i++) {
			String form = cmd[i];
			if (form.trim().toLowerCase().startsWith(gapKey + pseparator)) {
				return form.substring(form.indexOf(pseparator) + 1).split(",");
			}
		}
		return null;
	}

	// don't know how to set up a path, so far path read from the command line;
	// but probably it is not how it is normally done
	// waiting for suggestions
	// formats not yet decided
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
	FamiliesConservationClient(String[] cmd) throws IOException {

		String startStr = this.getDateTime();
		long startTime = System.currentTimeMillis();
		int SMERFSWidth = 7;
		SMERFSColumnScore score = SMERFSColumnScore.MID_SCORE;
		double SMERFSGapTreshold = 0.1;
		Boolean proceed = true;
		String method = getMethodName(cmd);
		if (method == null) {
			System.out
					.println("Methods not provided. Please provide methods in format -m=method1,marthod2,method3");
			Method.supportedMethods();
		}
		String inFilePath = getInputFilePath(cmd);
		if (inFilePath == null) {
			System.out
					.println("Input file path not provided. Please provide input file path in format -i=inputFile - where inputFile is a full path to FASTA formatted file.");
		}
		String[] groups = getGroups(cmd);
		if (groups == null) {
			System.out.println("Groups not provided");
		}
		if (method != null && inFilePath != null && groups != null) {
			int[][] grInt = parseGroups(groups, proceed);
			if (proceed == true) {
				String outFilePath = getOutputFilePath(cmd);
				// boolean proceed = true;
				String[] SMERFSDetails = getSMERFSDetails(cmd);
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
					}
					// else {
					// System.out.println("To run SMERFS three arguments are needed, window width,how to give scores to columns and a gap treshold.");
					// runSMERFS = false;
					// }
				}
				boolean normalize = getNormalize(cmd);
				String[] gap = getGapChars(cmd);
				Character[] gapChars;
				if (gap == null) {
					gapChars = null;
				} else {
					gapChars = new Character[gap.length];
					for (int i = 0; i < gap.length; i++) {
						if (gap[i].length() == 1) {
							gapChars[i] = gap[i].charAt(0);
						} else {
							throw new IllegalGapCharacterException(
									"Argument provided as gap charcetr is more than one character long.");
						}
					}
				}
				String statFile = getStatFilePath(cmd);
				// InputStream inStr = null;
				List<FastaSequence> sequences = this
						.openInputStream(inFilePath);
				if (sequences != null) {
					PrintWriter print = null;
					if (statFile != null) {
						print = ConservationFormatter.openPrintWriter(statFile,
								false);
					}
					if ((statFile == null && print == null)
							|| (statFile != null && print != null)) {
						AminoAcidMatrix alignment = new AminoAcidMatrix(
								sequences, gapChars);
						long loadedTime = System.currentTimeMillis();
						long loadTime = loadedTime - startTime;
						if (print != null) {
							print.println("Start time: " + startStr);
							print.println("Alignment loaded in: " + loadTime
									+ "ms.");
							print.println("Alignment has: "
									+ alignment.numberOfRows() + " sequences.");
						}
						SubFamiliesConservation sub = new SubFamiliesConservation(
								alignment, grInt);
						// double[] result = null;
						long time1 = System.currentTimeMillis();
						if (Method.getMethod(method) == Method.SMERFS) {
							if (SMERFSDetails != null
									&& SMERFSDetails.length != 3) {
								System.out
										.println("To run SMERFS three arguments are needed, window width,how to give scores to columns and a gap treshold.");
							} else {
								// result = getSMERFS(alignment, SMERFSWidth,
								// score, SMERFSGapTreshold, normalize);
								info1 = sub.subgrupsConservation3(SMERFSWidth,
										SMERFSGapTreshold, score, normalize);
								info2 = sub.subFamilyPairsConservation3(
										SMERFSWidth, score, SMERFSGapTreshold,
										normalize);
							}
						} else {
							// result = getMethod(methods[i], scores,
							// normalize);
							info1 = sub.subgrupsConservation2(Method
									.getMethod(method), normalize);
							info2 = sub.subFamilyPairsConservation2(Method
									.getMethod(method), normalize);
						}
						long time2 = System.currentTimeMillis();
						long time = time2 - time1;
						if (print != null) {
							print.println(Method.getMethod(method).toString()
									+ " done " + time + "ms.");
						}
						PrintWriter print2 = null;
						try {
							print2 = new PrintWriter(new BufferedWriter(
									new FileWriter(outFilePath)));
						} catch (IOException ex) {
							System.out.println("Problem writing" + outFilePath);
						}
						int width1 = 20;
						int width2 = 20;
						String format1 = "%-" + width1 + "s";
						String format2 = "%-" + width2 + "s";
						print2.printf(format1, "");
						print2.printf(format1, "Group name");
						print2.printf(format2, "Conservation");
						print2.println("Properties");
						for (int i = 0; i < info1.length; i++) {
							print2.println("Column nr: " + i);
							for (int j = 0; j < info1[i].length; j++) {
								print2.printf(format1, "");
								info1[i][j].printInfo(Method.getMethod(method),
										print2);
							}
							for (int j = 0; j < info2[i].length; j++) {
								print2.printf(format1, "");
								info2[i][j].printInfo(Method.getMethod(method),
										print2);
							}
						}
						print2.close();
						if (print != null) {
							print.println("End time: " + getDateTime());
							print.close();
						}
					}
				}
				// ConservationFormatter.formatResults(scores);
			}
		}
	}

	// Map<Method, double[]> getResults() {
	// return results;
	// }
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
			// System.out.print(info);
		}
		if (args.length < 2) {
			System.out
					.println("Method names, input file paths are required. Application will not run until these 2 arguments are provided.");
			System.out
					.println("If you want results printed, both format an input file path have to be provided");
			System.out.println();
			// System.out.print(info);
		}
		try {
			FamiliesConservationClient cons = new FamiliesConservationClient(
					args);
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println("Fail to write to the filesystem! "
					+ "Cannot output results!  " + e.getLocalizedMessage());
		}
	}

	/**
	 * Opens input stream
	 * 
	 * @param inStr
	 * @param fastaSeqs
	 * @param inFilePath
	 * @return
	 */
	List<FastaSequence> openInputStream(String inFilePath) {

		InputStream inStr = null;
		List<FastaSequence> fastaSeqs = null;
		try {
			inStr = new FileInputStream(inFilePath);
			fastaSeqs = SequenceUtil.readFasta(inStr);
		} catch (FileNotFoundException e) {
			System.out
					.println("Can not find file. Please provide a valid file path.");
			// I'm using system exit here to avoid an exception I get, but I'm
			// not sure it;s a good programming practice.
			// System.exit(0)
		} catch (IOException e) {
			System.out.println("Sth wrong with reading the file.");
			// I'm using system exit here to avoid an exception I get, but I'm
			// not sure it;s a good programming practice.
			// System.exit(0);
			fastaSeqs = null;
		}
		return fastaSeqs;
	}

	private String getDateTime() {

		DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		Date date = new Date();
		return dateFormat.format(date);
	}
}
