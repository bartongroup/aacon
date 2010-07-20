package compbio.conservation;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.List;

import compbio.util.FastaSequence;
import compbio.util.SequenceUtil;

public class SMERFSClient {
	
//private final Map<Method, double[]> scores = new EnumMap<Method, double[]>(Method.class);
	
	double[] correlationScores;
	
	final static String pseparator = "=";
	
	//final static String normalizationKey = "-n";
	
	final static String formatKey = "-f";
	
	final static String inputKey = "-i";
	
	final static String outputKey = "-o";
	
	final static String widthKey = "-w";
	
	
	/**
	 * Gets the normalization status from the command line.
	 * 
	 * @param cmd array of cmd arguments
	 * @return true if results to be normalized false if else, returns false if no normalization status provided
	 */
	
	//static boolean getNormalize(String[] cmd) {
		
		//for (int i = 0; i < cmd.length; i++) {
			
			//String norm = cmd[i];
			
			//if(norm.trim().toLowerCase().equals(normalizationKey)) {
				
				//return true;
			//}
		//}
		
		//return false;
		
	//}
	
	/**
	 * Gets the normalization status from the command line.
	 * 
	 * @param cmd array of cmd arguments
	 * @return true if results to be normalized false if else, returns false if no normalization status provided
	 */
	
	private static Integer getColumnWidth(String[] cmd) {
		
		for (int i = 0; i < cmd.length; i++) {
			
			String width = cmd[i];
			
			if(width.trim().toLowerCase().startsWith(widthKey + pseparator)) {
				
				int colWidth = Integer.parseInt(width.substring(width.indexOf(pseparator) + 1));
				
				if ( colWidth > 0 && colWidth%2 != 0) {
					
					return colWidth;
				}
				
				else {
					
					return null;
				}
			}
		}
		
		return null;
		
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
	 * Constructor
	 * @param cmd command line arguments
	 */
	
	SMERFSClient(String[] cmd) {

		String inFilePath = getInputFilePath(cmd);

		if (inFilePath == null) {
			
			System.out.println("Input file path not provided. Please provide input file path in format -i=inputFile - where inputFile is a full path to FASTA formatted file.");
			
		}
		
		Integer width = getColumnWidth(cmd);
		
		if (width == null) {
			
			System.out.println("No column width was provided or the column width provided was not appropriate");
			
			System.out.println("Column width must be an odd number greater than 0. Provide column width in format -w=integer");
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
		
		//boolean normalize = getNormalize(cmd);
		
		if (width != null && inFilePath != null) {
		
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
		
		double[] result = Correlation.getCorrelationScore(alignment, width);
		
		correlationScores = result;
		
			if (outFilePath != null && format != null) {
				
				if(Format.getFormat(format) == Format.RESULT_WITH_ALIGNMENT) {
			
					ConservationFormatter.printResultWithAlignment(alignment, Smerfs.SMERFS, result, 20, 10, 3, outFilePath);
			
				}
				
				else {
					
					ConservationFormatter.printResultNoAlignment(alignment, Smerfs.SMERFS, result, 20, 10, 3, outFilePath);
				}
			
			}
			
			//scores.put(meth, this.getMethod(meth, alignment, normalize));
			
		}
		
		//ConservationFormatter.formatResults(scores);
		
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
			
			System.out.println("Column width, input file path are required. Application will not run until these 2 arguments are provided.");
			System.out.println("If you want results printed, both format an input file path have to be provided");
		}
		
		SMERFSClient smerfs = new SMERFSClient(args);
	}


}
