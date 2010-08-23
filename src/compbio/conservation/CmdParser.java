package compbio.conservation;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.EnumSet;
import java.util.List;
import java.util.Set;

import compbio.common.IllegalGapCharacterException;
import compbio.data.sequence.FastaSequence;
import compbio.data.sequence.SequenceUtil;

public final class CmdParser {

	final static DateFormat DATE_FORMAT = new SimpleDateFormat(
			"yyyy/MM/dd HH:mm:ss");

	final static String pseparator = "=";
	final static String methodKey = "-m";
	final static String normalizationKey = "-n";
	final static String formatKey = "-f";
	final static String inputKey = "-i";
	final static String outputKey = "-o";
	final static String SMERFSDetailsKey = "-s";
	final static String gapKey = "-g";
	final static String statKey = "-d";
	final static String CONSERVATION_HELP = "AA Connservation 1.0\n"
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

	private CmdParser() {
		// This is a utility class no instantiation
	}

	/**
	 * Gets method name from the command line
	 * 
	 * @param cmd
	 *            array of cmd arguments
	 * @return method name or null if no method name provided
	 */
	static String[] getSMERFSDetails(String[] cmd) {

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
	static boolean getNormalize(String[] cmd) {

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
	static String getFormat(String[] cmd) {

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
	static String getOutputFilePath(String[] cmd) {

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
	static String getInputFilePath(String[] cmd) {

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
	static String getStatFilePath(String[] cmd) {

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
	static String[] getGapChars(String[] cmd) {

		for (int i = 0; i < cmd.length; i++) {
			String form = cmd[i];
			if (form.trim().toLowerCase().startsWith(gapKey + pseparator)) {
				return form.substring(form.indexOf(pseparator) + 1).split(",");
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
	static Set<Method> getMethodNames(String[] cmd) {
		Set<Method> methods = EnumSet.noneOf(Method.class);
		for (int i = 0; i < cmd.length; i++) {
			String meths = cmd[i];
			if (meths.trim().toLowerCase().startsWith(methodKey + pseparator)) {
				String[] mets = meths.substring(meths.indexOf(pseparator) + 1)
						.split(",");
				for (String method : mets) {
					methods.add(Method.getMethod(method));
				}
			}
		}
		return methods;
	}

	static Character[] extractGapChars(String[] gap) {
		Character[] gapChars = null;
		if (gap != null) {
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
		return gapChars;
	}

	/**
	 * Opens input stream
	 * 
	 * @param inStr
	 * @param fastaSeqs
	 * @param inFilePath
	 * @throws IOException
	 * @returnif (statFile != null) { print =
	 *           ConservationFormatter.openPrintWriter(statFile, false); }
	 */
	static List<FastaSequence> openInputStream(String inFilePath)
			throws IOException {

		InputStream inStr = new FileInputStream(inFilePath);
		List<FastaSequence> fastaSeqs = SequenceUtil.readFasta(inStr);

		return fastaSeqs;
	}

	static String getDateTime() {
		return DATE_FORMAT.format(new Date());
	}

}
