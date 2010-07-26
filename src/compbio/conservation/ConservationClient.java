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
	
	private double[] getMethod(String method, AminoAcidMatrix matrix, ConservationScores2 scores, boolean normalize) {
		
		double[] result = null;
		
		if (Method.getMethod(method) == null) {
			
			System.out.println("Method: " + method + "is not supported");
			
			Method.supportedMethods();
			
			return result;
		}
	
		else {
			
			result = scores.calculateScore(Method.getMethod(method), normalize);

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
		
		if(methods != null && inFilePath != null) {
		
		String format = getFormat(cmd);
		
		String outFilePath = getOutputFilePath(cmd);
		
		boolean proceed = true;
		
		if (format == null && outFilePath != null) {
			
			System.out.println("Format not provided. Default format will be used");
			
			Format.supportedFormats();
			
			
		}
		
		if (outFilePath == null && format != null) {
			
			System.out.println("Output file path not provided. Please provide output file path in format -o=outputFile - where outputFile is a full path to the file where the results are to be printed.");
			
			proceed = false;
			
		}
		
		if (outFilePath == null && format == null) {
			
			System.out.println("Output file path and format not provided. Results will be printed to the command window.");
			
		}
		
		boolean normalize = getNormalize(cmd);
		
		if (proceed == true) {
		
		//InputStream inStr = null;
		
			List<FastaSequence> sequences = this.openInputStream(inFilePath);
		
			if (sequences != null) {
		
				AminoAcidMatrix alignment = new AminoAcidMatrix(sequences);
		
				ConservationScores2 scores = new ConservationScores2(alignment);
		
				double[] result = null;
		
				for (int i = 0; i < methods.length; i++) {
			
					result = this.getMethod(methods[i], alignment, scores, normalize);
			
					if(result != null) {
						
						results.put(Method.getMethod(methods[i]), result);
			
					}
		
				}
				
				if (results.size() != 0) {
					
				//alignment.printAlignment(30, 10, outFilePath);
				
				ConservationFormatter.formatResults(results, outFilePath, format, alignment);
				
				}
		
		}
		//ConservationFormatter.formatResults(scores);
		}
		
		}
		
	}
	
	Map<Method, double[]> getResults() {
		
		return results;
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
	
	/**
	 * Opens input stream
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
			
			}
		
			catch (FileNotFoundException e) {
			
			System.out.println("Can not find file. Please provide a valid file path.");
			
			// I'm using system exit here to avoid an exception I get, but I'm not sure it;s a good programming practice.
			//System.exit(0)
			
			}
		
			catch (IOException e) {
			
			System.out.println("Sth wrong with reading the file.");
			
			// I'm using system exit here to avoid an exception I get, but I'm not sure it;s a good programming practice.
			//System.exit(0);
			
			fastaSeqs = null;
			
			}
			
			return fastaSeqs;
	
		}

}
