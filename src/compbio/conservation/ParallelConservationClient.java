package compbio.conservation;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.SynchronousQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;

import compbio.common.IllegalGapCharacterException;
import compbio.data.sequence.FastaSequence;
import compbio.data.sequence.SequenceUtil;

/**
 * TODO to complete
 * 
 * @author pvtroshin
 */
class ParallelConservationClient {

	private final Map<Method, double[]> results = new EnumMap<Method, double[]>(
			Method.class);
	private static volatile ExecutorService executor;
	final static String pseparator = "=";
	final static String methodKey = "-m";
	final static String normalizationKey = "-n";
	final static String formatKey = "-f";
	final static String inputKey = "-i";
	final static String outputKey = "-o";
	final static String SMERFSDetailsKey = "-s";
	final static String gapKey = "-g";
	final static String statKey = "-d";
	final static String info = "AA Connservation 1.0\n"
			+ "\n"
			+ "This program allows calculation of conservation of amino acids in\n"
			+ "multiple sequence alignments.\n"
			+ "It implements 15 different conservation scores as described by Valdar in\n"
			+ "his paper (Scoring Residue Conservation, PROTEINS: Structure, Function\n"
			+ "and Genetics 48:227-241 (2002)) and SMERFS scoring algorithm as described\n"
			+ "by Manning, Jefferson and Barton (The contrasting properties of conservation\n"
			+ "and correlated phylogeny in protein functional residue prediction,\n"
			+ "BMC Bioinformatics (2008)).\n"
			+ "The scores supported are:\n"
			+ "KABAT, JORES, SCHNEIDER, SHENKIN, GERSTEIN, TAYLOR_GAPS, TAYLOR_NO_GAPS, \n"
			+ "ZVELIBIL, KARLIN, ARMON, THOMPSON, NOT_LANCET, MIRNY, WILLIAMSON, \n"
			+ "LANDGRAF, SANDER, VALDAR, SMERFS\n"
			+ "\n"
			+ "Input format is a FASTA file. By default program prints the results\n"
			+ "to the command window. If the output file is provided the results are \n"
			+ "printed to the file in two possible formats (with or without an alignment).\n"
			+ "The default format prints the results without alignment.  \n"
			+ "The scores can be normalized or not. By default the scores are not \n"
			+ "normalized. If SMERFS is called and no parameters for SMERFS are provided\n"
			+ "it is run with the default arguments (window width of 7, column score given to \n"
			+ "the middle column, gap% cutoff of 1.0). If new parameters are needed all\n"
			+ "three of them have to be provided. Details of the program execution can\n"
			+ "be listed to a separate file if an appropriate file path is provided.\n"
			+ "Program accepts command line arguments. All the arguments are preceded by a certain key:\n"
			+ "-m= - precedes a comma separated list of method names\n"
			+ "      EXAMPLE: -m=KABAT,JORES,GERSTEIN\n"
			+ "-i= - precedes a full path to the input FASTA file\n"
			+ "-o= - precedes a full path to the output file\n"
			+ "-f= - precedes the format  of the results in the output file\n"
			+ "      two different formats are possible:\n"
			+ "      RESULTS_WITH_ALIGNMENT\n"
			+ "      RESULTS_NO_ALIGNMENT\n"
			+ "-s= - precedes a list of three comma separated parameters for SMERFS\n"
			+ "      the order of parameters is as following:\n"
			+ "      1. window width - an integer and an odd number\n"
			+ "      2. how to allocate window scores to columns, two ways are possible:\n"
			+ "	     MID_SCORE - gives the window score to the middle column\n"
			+ "	     MAX_SCORE - gives the column the highest score of all the windows it belongs to\n"
			+ "      3. gap percentage cutoff - a float greater than 0 and smaller or equal 1\n"
			+ "      EXAMPLE: -s=5,MID_SCORE,0.1\n"
			+ "      \n"
			+ "-d= - precedes a full path to a file where program execution details are to be listed  \n"
			+ "-g= - precedes comma separated list of gap characters provided by the user, if you're using an unusual gap character\n"
			+ "      (not a -,., ,*) you have to provide it. If you you provide this list you have to list all the gaps accepted. Including those\n"
			+ "      that were previously treated as a default.      \n"
			+ "-n - using this key causes the results to be normalized\n"
			+ "\n"
			+ "EXAMPLE HOW TO RUN THE PROGRAM:\n"
			+ "-m=KABAT,SMERFS -i=/homes/agolicz/alignments/prot1 -o=/homes/agolicz/alignments/prot1_results \n"
			+ "-f=RESULTS_NO_ALIGNMENT -n\n"
			+ "\n"
			+ "KABAT and SMERFS scores will be calculated. Input will come form /homes/agolicz/alignments/prot1 \n"
			+ "file and will be printed to /homes/agolicz/alignments/prot1_results without alignment. \n"
			+ "Results will be normalized between 0 and 1.\n"
			+ "To normalize the results n = (d - dmin)/(dmax - dmin) formula is used.\n"
			+ "If the results are negative they are first shifted by adding the \n"
			+ "absolute value of the most negative number. \n" + "";

	private void initExecutor(int procNum) {

		int corenum = Runtime.getRuntime().availableProcessors();
		if (procNum <= 1 || procNum > corenum * 2) {
			System.err.println("Number of processors must be more than 1 and "
					+ "less than the number of cores*2"
					+ "However given number of processors is " + procNum);
			System.err.println("Changing number of processors to " + corenum
					+ " - the number of cores");
			procNum = corenum;
		}
		if (executor == null) {
			synchronized (ParallelConservationClient.class) {
				if (executor == null) {
					executor = new ThreadPoolExecutor(procNum, procNum, 0L,
							TimeUnit.MILLISECONDS,
							new SynchronousQueue<Runnable>(),
							new ThreadPoolExecutor.CallerRunsPolicy());
				}
			}
		}
	}

	/**
	 * Gets method name from the command line
	 * 
	 * @param cmd
	 *            array of cmd arguments
	 * @return method name or null if no method name provided
	 */
	private static String[] getMethodNames(String[] cmd) {

		for (int i = 0; i < cmd.length; i++) {
			String meths = cmd[i];
			if (meths.trim().toLowerCase().startsWith(methodKey + pseparator)) {
				return meths.substring(meths.indexOf(pseparator) + 1)
						.split(",");
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
	 * Gets the format of the file for the output.
	 * 
	 * @param cmd
	 *            array of cmd arguments
	 * @return format of null if format not provided
	 */
	private static String getFormat(String[] cmd) {

		for (int i = 0; i < cmd.length; i++) {
			String form = cmd[i];
			if (form.trim().toLowerCase().startsWith(formatKey + pseparator)) {
				return form.substring(form.indexOf(pseparator) + 1);
			}
		}
		return null;
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
	 * Gets the gap charcter .
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
	 * Returns the results of method calculation or null if method not
	 * supported.
	 */
	static double[] getMethod(String method, ConservationScores2 scores,
			boolean normalize) {

		double[] result = null;
		if (Method.getMethod(method) == null) {
			System.out.println("Method: " + method + " is not supported");
			Method.supportedMethods();
			return result;
		} else {
			result = scores.calculateScore(Method.getMethod(method), normalize);
			return result;
		}
	}

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
	ParallelConservationClient(String[] cmd) throws IOException {

		String startStr = this.getDateTime();
		long startTime = System.currentTimeMillis();
		int SMERFSWidth = 7;
		SMERFSColumnScore score = SMERFSColumnScore.MID_SCORE;
		double SMERFSGapTreshold = 0.1;
		// boolean runSMERFS = true;
		String[] methods = getMethodNames(cmd);
		// if(methods == null) {
		// System.out.println("Methods not provided. Please provide methods in format -m=method1,marthod2,method3");
		// Method.supportedMethods();
		// }
		String inFilePath = getInputFilePath(cmd);
		// if (inFilePath == null) {
		// System.out.println("Input file path not provided. Please provide input file path in format -i=inputFile - where inputFile is a full path to FASTA formatted file.");
		// }
		if (methods != null && inFilePath != null) {
			String format = getFormat(cmd);
			String outFilePath = getOutputFilePath(cmd);
			boolean proceed = true;
			if (format == null && outFilePath != null) {
				System.out.println("Format not provided.");
				Format.supportedFormats();
			}
			if (outFilePath == null && format != null) {
				System.out
						.println("Output file path not provided. Please provide output file path in format -o=outputFile - where outputFile is a full path to the file where the results are to be printed.");
				proceed = false;
			}
			Format outFormat = Format.getFormat(format);
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
			if (proceed == true) {
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
						ConservationScores2 scores = new ConservationScores2(
								alignment);
						double[] result = null;
						for (int i = 0; i < methods.length; i++) {
							long time1 = System.currentTimeMillis();
							MethodWrapper wrapper = null;
							if (Method.getMethod(methods[i]) == Method.SMERFS) {
								if (SMERFSDetails != null
										&& SMERFSDetails.length != 3) {
									System.out
											.println("To run SMERFS three arguments are needed, window width,how to give scores to columns and a gap treshold.");
								} else {
									wrapper = new MethodWrapper(scores,
											normalize, SMERFSWidth,
											SMERFSGapTreshold);
									// result = getSMERFS(alignment,
									// SMERFSWidth,
									// score, SMERFSGapTreshold, normalize);
								}
							} else {
								wrapper = new MethodWrapper(methods[i], scores,
										normalize);
								// result = getMethod(methods[i], scores,
								// normalize);
							}
							Future<double[]> wresult = executor.submit(wrapper);
							// TODO executor.awaitTermination(TimeUnit.SECONDS,
							// 3600);
							long time2 = System.currentTimeMillis();
							long time = time2 - time1;
							if (result != null) {
								results.put(Method.getMethod(methods[i]),
										result);
								if (print != null) {
									print.println(Method.getMethod(methods[i])
											.toString()
											+ " done " + time + "ms.");
								}
							}
						}
						if (results.size() != 0) {
							// alignment.printAlignment(30, 10, outFilePath);
							ConservationFormatter.formatResults(results,
									outFilePath, outFormat, alignment);
						}
						if (print != null) {
							print.println("End time: " + getDateTime());
							print.close();
						}
					}
				}
			}
		}
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
			System.out.print(info);
		}
		if (args.length < 2) {
			System.out
					.println("Method names, input file paths are required. Application will not run until these 2 arguments are provided.");
			System.out
					.println("If you want results printed, both format an input file path have to be provided");
			System.out.println();
			System.out.print(info);
		}
		try {
			ParallelConservationClient cons = new ParallelConservationClient(
					args);
		} catch (IOException e) {
			System.err.println("Fail to write to the file system! "
					+ e.getLocalizedMessage());
			e.printStackTrace();
		}
	}

	/**
	 * Opens input stream
	 * 
	 * @param inStr
	 * @param fastaSeqs
	 * @param inFilePath
	 * @returnif (statFile != null) { print =
	 *           ConservationFormatter.openPrintWriter(statFile, false); }
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

	/**
	 * Reads printed results from the file.
	 * 
	 * @param inStream
	 *            input stream
	 * @return map mapping a map of results(with method names as keys and result
	 *         arrays as values) to the list of FASTA sequences that constitute
	 *         the alignment.
	 * @throws IOException
	 */
	static Map<Map<Method, double[]>, List<FastaSequence>> readFile(
			InputStream inStream) throws IOException {

		Map<Map<Method, double[]>, List<FastaSequence>> result = new HashMap<Map<Method, double[]>, List<FastaSequence>>();
		Map<Method, double[]> resultMap = new EnumMap<Method, double[]>(
				Method.class);
		BufferedReader inResults = new BufferedReader(new InputStreamReader(
				inStream));
		List<FastaSequence> seqList = new ArrayList<FastaSequence>();
		Pattern pattern = Pattern.compile("\\s+");
		String line;
		String lineString = null;
		// line = inResults.readLine();
		// line = inResults.readLine();
		do {
			line = inResults.readLine();
			// System.out.println(line);
			if (line == null || line.startsWith("#") || line.startsWith(">")) {
				if (lineString != null) {
					if (lineString.startsWith("#")) {
						parseResults(lineString, resultMap, pattern);
					}
					if (lineString.startsWith(">")) {
						parseSequences(lineString, seqList, pattern);
					}
				}
				lineString = line;
			} else {
				lineString += line;
			}
		} while (line != null);
		result.put(resultMap, seqList);
		return result;
	}

	static void parseResults(String resultStr, Map<Method, double[]> resultMap,
			Pattern pattern) {

		String resultStrTemp = pattern.matcher(resultStr.trim())
				.replaceAll(" ");
		String[] results = resultStrTemp.split(" ");
		String name = results[0].substring(1);
		System.out.println(name);
		double[] resultsNum = new double[results.length - 1];
		for (int i = 0; i < resultsNum.length; i++) {
			resultsNum[i] = Double.parseDouble(results[i + 1]);
			System.out.println(resultsNum[i]);
		}
		resultMap.put(Method.getMethod(name), resultsNum);
	}

	static void parseSequences(String lineStr, List<FastaSequence> list,
			Pattern pattern) {

		StringTokenizer tokens = new StringTokenizer(lineStr, " ");
		String name = tokens.nextToken().trim();
		System.out.println(name);
		Pattern pattern2 = Pattern.compile(name);
		String seqStrMod = pattern.matcher(
				pattern2.matcher(lineStr).replaceFirst("")).replaceAll("");
		System.out.println(seqStrMod);
		list.add(new FastaSequence(name, seqStrMod));
	}

	private String getDateTime() {

		DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		Date date = new Date();
		return dateFormat.format(date);
	}
}
