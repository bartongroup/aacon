package compbio.conservation;

import java.io.*;
import java.util.*;

import compbio.util.FastaSequence;
import compbio.util.SequenceUtil;

class ConservationClient {
	
	private final Map<Method, double[]> results = new EnumMap<Method, double[]>(Method.class);
	
	final static String pseparator = "=";
	
	final static String methodKey = "-m";
	
	final static String normalizationKey = "-n";
	
	final static String formatKey = "-f";
	
	final static String inputKey = "-i";
	
	final static String outputKey = "-o";
	
	/**
	 * Gets method name from the command line
	 * 
	 * @param cmd array of cmd arguments
	 * @return method name or null if no method name provided
	 */
	
	private static String[] getMethodNames(String[] cmd) {
		
		for (int i = 0; i < cmd.length; i++) {
			
			String meths = cmd[i];
			
			if(meths.trim().toLowerCase().startsWith(methodKey + pseparator)) {
				
				return meths.substring(meths.indexOf(pseparator) + 1).split(",");
			}
		}
		
		return null;
		
	}
	
	/**
	 * Gets the normalization status from the command line.
	 * 
	 * @param cmd array of cmd arguments
	 * @return true if results to be normalized false if else, returns false if no normalization status provided
	 */
	
	private static boolean getNormalize(String[] cmd) {
		
		for (int i = 0; i < cmd.length; i++) {
			
			String norm = cmd[i];
			
			if(norm.trim().toLowerCase().equals(normalizationKey)) {
				
				return true;
			}
		}
		
		return false;
		
	}
	
	/**
	 * Gets the format of the file for the output.
	 * 
	 * @param cmd array of cmd arguments
	 * @return format of null if format not provided
	 */
	
	private static String getFormat(String[] cmd) {
		
		for (int i = 0; i < cmd.length; i++) {
			
			String form = cmd[i];
			
			if(form.trim().toLowerCase().startsWith(formatKey + pseparator)) {
				
				return form.substring(form.indexOf(pseparator) + 1);
			}
		}
		
		return null;
	}
	
	/**
	 * Gets output file path
	 * @param cmd
	 * @return null if not provided
	 */
	private static String getOutputFilePath(String[] cmd) {
		
		for (int i = 0; i < cmd.length; i++) {
			
			String name = cmd[i];
			
			if(name.trim().toLowerCase().startsWith(outputKey + pseparator)) {
				
				return name.substring(name.indexOf(pseparator) + 1);
			}
		}
		
		return null;
	}
	
	/**
	 * Input file path. 
	 * @param cmd
	 * @return null if not provided
	 */
	private static String getInputFilePath(String[] cmd) {
		
		for (int i = 0; i < cmd.length; i++) {
			
			String name = cmd[i];
			
			if(name.trim().toLowerCase().startsWith(inputKey + pseparator)) {
				
				return name.substring(name.indexOf(pseparator) + 1);
			}
		}
		
		return null;
	}
	
	
	
	// don't know how to set up a path, so far path read from the command line;
	// but probably it is not how it is normally done
	// waiting for suggestions
	// formats not yet decided
	
	/**
	 * Returns the results of method calculation.
	 */
	
	private double[] getMethod(Method method, AminoAcidMatrix matrix, boolean normalize) {
		
		double[] result = null;
		
		if (method == null) {
			
			System.out.println("Method not suppoted");
			
			return result;
		}
	
		else {
			
			ConservationScores2 scores = new ConservationScores2(matrix);
			
			result = scores.calculateScore(method, normalize);

			return result;
			
		}
		
	}
	/**
	 * Constructor
	 * @param cmd command line arguments
	 */
	
	ConservationClient(String[] cmd) {
		
		String[] methods = getMethodNames(cmd);
		
		if(methods == null) {
			
			System.out.println("Methods not provided. Please provide methods in format -m=method1,marthod2,method3");
			
			Method.supportedMethods();
			
		}

		String inFilePath = getInputFilePath(cmd);

		if (inFilePath == null) {
			
			System.out.println("Input file path not provided. Please provide input file path in format -i=inputFile - where inputFile is a full path to FASTA formatted file.");
			
		}
		
		String format = getFormat(cmd);
		
		String outFilePath = getOutputFilePath(cmd);
		
		if (format == null && outFilePath != null) {
			
			System.out.println("Format not provided. Please provide format in format -f=format.");
			
			Format.supportedFormats();
			
		}
		
		if (outFilePath == null && format != null) {
			
			System.out.println("Output file path not provided. Please provide output file path in format -o=outputFile - where outputFile is a full path to the file where the results are to be printed.");
			
		}
		
		boolean normalize = getNormalize(cmd);
		
		if (methods != null && inFilePath != null) {
		
		InputStream inStr = null;
		
		List<FastaSequence> fastaSeqs = null;
		
		try {
			
			inStr = new FileInputStream(inFilePath);
			
		}
		
		catch (IOException e) {
			
			System.out.println("Can not find file");
		
		}
		
		try {
			
			fastaSeqs = SequenceUtil.readFasta(inStr);
		}
		
		catch (IOException e) {
			
			System.out.println("Sth wrong with reading the file");
		}
			
		
		AminoAcidMatrix alignment = new AminoAcidMatrix(fastaSeqs);
		
		ConservationScores2 scores = new ConservationScores2(alignment);
		
		double[] result = null;
		
		for (int i = 0; i < methods.length; i++) {
			
			Method meth = Method.getMethod(methods[i]);
			
			result = scores.calculateScore(meth, normalize);
			
			results.put(meth, this.getMethod(meth, alignment, normalize));
			
			if (outFilePath != null && format != null) {
				
				if(Format.getFormat(format) == Format.RESULT_WITH_ALIGNMENT) {
			
					ConservationFormatter.printResultWithAlignment(alignment, meth, result, 20, 10, 3, outFilePath);
			
				}
				
				else {
					
					ConservationFormatter.printResultNoAlignment(alignment, meth, result, 20, 10, 3, outFilePath);
				}
			}
			
		}
		
		//ConservationFormatter.formatResults(scores);
		}
		
	}
	
	/**
	 * Application entry point. 
	 * Command line format looks like.
	 * -m=method1,method2,method3 - method names
	 * -n normalization factor
	 * -f=format
	 * -i=inputPath - inputfile path
	 * -o=outputPath - outputfile path
	 * 
	 * @param args command line arguments
	 */
	
	public static void main(String[] args) {
		
		if(args == null) {
			
			System.out.println ("No parameters were suppled");
		}
		
		if(args.length < 2) {
			
			System.out.println("Method names, input file paths are required. Application will not run until these 2 arguments are provided.");
			System.out.println("If you want results printed, both format an input file path have to be provided");
		}
		
		ConservationClient cons = new ConservationClient(args);
		
	}

}
