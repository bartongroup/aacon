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

/**
 * Makes sense of the command line arguments
 * 
 * @author Peter Troshin & Agnieszka Golicz
 * 
 */
final class CmdParser {

	final static DateFormat DATE_FORMAT = new SimpleDateFormat(
			"yyyy/MM/dd HH:mm:ss");

	final static String pseparator = "=";
	final static String methodKey = "-m";
	final static String normalizationKey = "-n";
	final static String threadNumberKey = "-t";
	final static String formatKey = "-f";
	final static String inputKey = "-i";
	final static String outputKey = "-o";
	final static String SMERFSDetailsKey = "-s";
	final static String gapKey = "-g";
	final static String statKey = "-d";
	final static String CONSERVATION_HELP = "\r\n"
			+ "AA Conservation version 1.0b (6 September 2010)\r\n"
			+ "\r\n"
			+ "This program allows calculation of conservation of amino acids in\r\n"
			+ "multiple sequence alignments.\r\n"
			+ "It implements 15 different conservation scores as described by Valdar in\r\n"
			+ "his paper (Scoring Residue Conservation, PROTEINS: Structure, Function\r\n"
			+ "and Genetics 48:227-241 (2002)) and SMERFS scoring algorithm as described\r\n"
			+ "by Manning, Jefferson and Barton (The contrasting properties of conservation\r\n"
			+ "and correlated phylogeny in protein functional residue prediction,\r\n"
			+ "BMC Bioinformatics (2008)).\r\n"
			+ "\r\n"
			+ "The conservation algorithms supported are:\r\n"
			+ "\r\n"
			+ "KABAT, JORES, SCHNEIDER, SHENKIN, GERSTEIN, TAYLOR_GAPS, TAYLOR_NO_GAPS, \r\n"
			+ "ZVELIBIL, KARLIN, ARMON, THOMPSON, NOT_LANCET, MIRNY, WILLIAMSON, \r\n"
			+ "LANDGRAF, SANDER, VALDAR, SMERFS\r\n"
			+ "\r\n"
			+ "Input format is either a FASTA formatted file containing aligned sequences with \r\n"
			+ "gaps or a Clustal alignment. The valid gap characters are *, -, space character,\r\n"
			+ "X and . (a dot). By default program prints the results to the command window. \r\n"
			+ "If the output file is provided the results are printed to the file in two \r\n"
			+ "possible formats with or without an alignment.\r\n"
			+ "If format is not specified, the program outputs conservation scores without \r\n"
			+ "alignment. The scores are not normalized by default but they can be (see below).\r\n"
			+ "SMERFS default parameters are window width of 7, column score is set to\r\n"
			+ "the middle column, gap% cutoff of 0.1. If different values for SMERFS parameters \r\n"
			+ "are required than all three parameters must be provided. Details of the program \r\n"
			+ "execution can be recorded to a separate file if an appropriate file path is \r\n"
			+ "provided.\r\n"
			+ "\r\n"
			+ "List of command line arguments:\r\n"
			+ "\r\n"
			+ "-m=  precedes a comma separated list of method names\r\n"
			+ "     EXAMPLE: -m=KABAT,JORES,GERSTEIN\r\n"
			+ "     Optional, if no method is specified request for all is assumed. \r\n"
			+ "\r\n"
			+ "-i=  precedes a full path to the input FASTA file, required\r\n"
			+ "\r\n"
			+ "-o=  precedes a full path to the output file, optional, if no output file is \r\n"
			+ "     provided the program will output to the standard out.  \r\n"
			+ "\r\n"
			+ "-t=  precedes the number of CPUs (CPU cores more precisely) to use. Optional, \r\n"
			+ "	 defaults to all processors available on the machine.  \r\n"
			+ "      \r\n"
			+ "-f=  precedes the format  of the results in the output file\r\n"
			+ "     two different formats are possible:\r\n"
			+ "      RESULT_WITH_ALIGNMENT\r\n"
			+ "      RESULT_NO_ALIGNMENT\r\n"
			+ "     Optional, if not specified RESULT_NO_ALIGNMENT is assumed \r\n"
			+ "\r\n"
			+ "-s=  precedes a list of three comma separated parameters for SMERFS\r\n"
			+ "     the order of parameters is as following:\r\n"
			+ "      1. window width - an integer and an odd number\r\n"
			+ "      2. how to allocate window scores to columns, two ways are possible:\r\n"
			+ "	     MID_SCORE - gives the window score to the middle column\r\n"
			+ "	     MAX_SCORE - gives the column the highest score of all the windows it \r\n"
			+ "	     belongs to\r\n"
			+ "      3. gap percentage cutoff - a float greater than 0 and smaller or equal 1\r\n"
			+ "     EXAMPLE: -s=5,MID_SCORE,0.1\r\n"
			+ "     Optional, default values are 7,MID_SCORE,0.1 \r\n"
			+ "      \r\n"
			+ "-d=  precedes a full path to a file where program execution details are to be \r\n"
			+ "     listed. Optional, if not provided, no execution statistics is produced.  \r\n"
			+ "      \r\n"
			+ "-g=  precedes comma separated list of gap characters provided by the user, if \r\n"
			+ "     you're using an unusual gap character (not a -,., ,*,X) you have to \r\n"
			+ "     provide it. If you you provide this list you have to list all the gaps \r\n"
			+ "     accepted. Including those that were previously treated as a default. \r\n"
			+ "     Optional.      \r\n"
			+ "      \r\n"
			+ "-n   using this key causes the results to be normalized. \r\n"
			+ "	 Normalized results have values between 0 and 1. Please note however, that \r\n"
			+ "	 some results cannot be normalized. In such a case, the system returns not \r\n"
			+ "	 normalized value, and log the issue to the standard error stream. \r\n"
			+ "	 The following formula is used for normalization \r\n"
			+ "			n = (d - dmin)/(dmax - dmin)\r\n"
			+ "	 Negative results first converted to positive by adding an absolute value of\r\n"
			+ "	 the most negative result. Optional. \r\n"
			+ "\r\n"
			+ "EXAMPLE HOW TO RUN THE PROGRAM:\r\n"
			+ "java -jar <jar name> -m=KABAT,SMERFS -i=prot1 -o=prot1_results -n\r\n"
			+ "\r\n"
			+ "As a result of the execution KABAT and SMERFS scores will be calculated. \r\n"
			+ "Input comes form prot1 file and an output without an alignment is recorded to \r\n"
			+ "prot1_results file. \r\n"
			+ "\r\n"
			+ "Authors: Agnieszka Golicz, Peter Troshin, David Martin and Geoff Barton.\r\n"
			+ "Please visit http://www.compbio.dundee.ac.uk for further information.\r\n"
			+ "\r\n ";

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
			if (meths.trim().toLowerCase()
					.startsWith(SMERFSDetailsKey + pseparator)) {
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
	 * Gets the number of threads to use from the command line.
	 * 
	 * By default the number of threads to be used is equals to the number of
	 * available processors
	 * 
	 * @param cmd
	 *            array of cmd arguments
	 * @return true if results to be normalized false if else, returns false if
	 *         no normalization status provided
	 */
	static int getThreadNumber(String[] cmd) {
		int threadNum = 0;
		for (int i = 0; i < cmd.length; i++) {
			String tnum = cmd[i].trim().toLowerCase();
			if (tnum.startsWith(threadNumberKey + pseparator)) {
				String num = tnum.substring(tnum.indexOf(pseparator) + 1);
				if (num != null) {
					threadNum = Integer.parseInt(num);
				}
			}
		}
		return threadNum;
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

	static char[] extractGapChars(String[] gap) {
		char[] gapChars = null;
		if (gap != null) {
			gapChars = new char[gap.length];
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
