package compbio.conservation;

import java.io.*;
import java.util.*;
import java.util.regex.Pattern;

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
	
	final static String SMERFSDetailsKey = "-s";
	
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
	 * Gets method name from the command line
	 * 
	 * @param cmd array of cmd arguments
	 * @return method name or null if no method name provided
	 */
	
	private static String[] getSMERFSDetails(String[] cmd) {
		
		for (int i = 0; i < cmd.length; i++) {
			
			String meths = cmd[i];
			
			if(meths.trim().toLowerCase().startsWith(SMERFSDetailsKey + pseparator)) {
				
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
	 * Returns the results of method calculation or null if method not supported.
	 */
	
	private static double[] getMethod(String method, ConservationScores2 scores, boolean normalize) {
		
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
	 * Returns results of SMERFS calculation or null, if parameters provided are not appropriate. 
	 * @param alignment reference to alignment 
	 * @param width with of the window
	 * @param score tells which score given to the column, either the highest score of all the windows it belongs to, or the middle column is given the score of the window. 
	 * @param normalize if true results will be normalized
	 * @return
	 */
	
	public static double[] getSMERFS(AminoAcidMatrix alignment, int width, SMERFSColumnScore score, boolean normalize){
		
		if(alignment == null) {
			
			throw new IllegalArgumentException("Matrix must not be null.");
		}
		
		double[] result = null;
		
		if (width <= 0 || width%2 != 1) {
			
			System.out.println("Column width for SMERFS smaller or equal zero or not an odd number.");
			
			return result;
		}
		
		if (width > alignment.numberOfColumns()) {
			
			System.out.println("Column width greater than the length of the alignment");
			
			return result;
		}
		
		if (score == null) {
			
			System.out.println("Column score type not supported.");
			
			SMERFSColumnScore.supportedSMERFSColumnSores();
			
			return result;
		}
		
		Correlation corr = new Correlation(alignment, width);
		
		result = corr.getCorrelationScore(score, normalize);
		
		return result;
	}
	/**
	 * Constructor
	 * @param cmd command line arguments
	 */
	
	ConservationClient(String[] cmd) {
		
		int SMERFSWidth = 7;
		
		SMERFSColumnScore score = SMERFSColumnScore.MID_SCORE;
		
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
			
			System.out.println("Output file path and format not provided.");
			
		}
		
		String[] SMERFSDetails = getSMERFSDetails(cmd);
		
		
		if (SMERFSDetails != null) {
			
			if (SMERFSDetails.length == 2) {
			
			SMERFSWidth = Integer.parseInt(SMERFSDetails[0]);
			
			score = SMERFSColumnScore.getSMERFSColumnScore(SMERFSDetails[1]);
			
			}
			
			else {
				
				SMERFSWidth = -1;
				
				score = null;
			}
			
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
					
					if (Method.getMethod(methods[i]) == Method.SMERFS) {
						
						result = getSMERFS(alignment, SMERFSWidth, score, normalize);
					}
					
					else {
			
						result = getMethod(methods[i], scores, normalize);
						
					}
			
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
		
		static Map<Method, double[]> readScores (InputStream inStream) throws IOException {
			
			
			Map<Method, double[]> resultMap = new EnumMap<Method, double[]>(Method.class);
			
			BufferedReader inResults = new BufferedReader(new InputStreamReader(inStream));
			
			String line;
			
			String resultsString = null;
			
			String name = "";
			
			boolean begin = false;
			
			do {
				
				line = inResults.readLine();
				
				if (line.startsWith("#") && begin == false) {
					
					begin = true;
					
				}
				
				if(begin == true) {
				
					if (line ==  null || line.startsWith("#")) {
						
						if(resultsString != null) {
							
							String[] results = resultsString.split(" ");
							
							double[] resultsNum = new double[results.length];
							
							for (int i = 0; i < resultsNum.length; i++) {
								
								resultsNum[i] = Double.parseDouble(results[i].trim());
							}
							
							resultMap.put(Method.getMethod(name), resultsNum);
							
							name = line.substring(1).trim();
							
							resultsString = "";
							
						}
						
					}
					
					else {
						
						resultsString += line;
					}
				}
				
				else {
					
					
				}
				
			} while (line != null);
			
				inResults.close();
				
				return resultMap;
		}
		
		static List<Object> readScores2 (InputStream inStream) throws IOException {
			
			
			Map<Method, double[]> resultMap = new EnumMap<Method, double[]>(Method.class);
			
			BufferedReader inResults = new BufferedReader(new InputStreamReader(inStream));
			
			List<String> seqNames = new ArrayList<String>();
			
			List<char[]> seqList = new ArrayList<char[]>();
			
			Pattern pattern = Pattern.compile("//s+");
			
			String line;
			
			String resultsString = null;
			
			String seqString = null;
			
			boolean beginResults = false;
			
			do {
				
				line = inResults.readLine();
				
				if (line.startsWith("#") && beginResults == false) {
					
					beginResults = true;
					
				}
				
				if(beginResults == true) {
				
					if (line ==  null || line.startsWith("#")) {
						
						if(resultsString != null) {
							
							String[] results = resultsString.split(" ");
							
							String name = results[0].trim().substring(1); 
							
							double[] resultsNum = new double[results.length - 1];
							
							for (int i = 0; i < resultsNum.length; i++) {
								
								resultsNum[i] = Double.parseDouble(results[i + 1].trim());
							}
							
							resultMap.put(Method.getMethod(name), resultsNum);
							
							name = line.substring(1).trim();
						
						}
						
						resultsString = line;
						
					}
					
					else {
						
						resultsString += line;
					}
				}
				
				else {
					
					if (line ==  null || line.startsWith(">")) {
						
						if(resultsString != null) {
							
							StringTokenizer tokens = new StringTokenizer(resultsString, " ");
							
							String name = tokens.nextToken().trim();
							
							Pattern pattern2 = Pattern.compile(name);
							
							String seqStringMod = pattern2.matcher(seqString).replaceFirst("");
							
							String seqStringMod2 = pattern.matcher(seqStringMod).replaceAll("");
							
							char[] seq = seqStringMod2.toCharArray();
							
							seqList.add(seq);
							
							seqNames.add(name);
						
						}
						
						seqString = line;
						
					}
					
					else {
						
						seqString += line;
					}
					
					
				}
				
			} while (line != null);
			
				inResults.close();
				
				char[][] alignment = new char[seqList.size()][];
				
				for (int i = 0; i < alignment.length; i++) {
					
					alignment[i] = seqList.get(i);
					
				}
				
				String[] names = new String[1];
				
				String[] names2 = seqNames.toArray(names);
				
				AminoAcidMatrix matrix = new AminoAcidMatrix(alignment, names2);
				
				List<Object> results = new ArrayList<Object>();
				
				results.add(matrix);
				
				results.add(resultMap);
				
				return results;
				
				
		}
		
		static Map<Method, double[]> readScores3 (InputStream inStream) throws IOException {
			
			
			Map<Method, double[]> resultMap = new EnumMap<Method, double[]>(Method.class);
			
			BufferedReader inResults = new BufferedReader(new InputStreamReader(inStream));
			
			List<String> seqNames = new ArrayList<String>();
			
			List<char[]> seqList = new ArrayList<char[]>();
			
			Pattern pattern = Pattern.compile("//s+");
			
			
			String line = null;
			
			String previousLine;
			
			String resultsString = null;
			
			String seqString = null;
			
			do {
				
				previousLine = line;
				
				line = inResults.readLine();
				
					if (line ==  null || line.startsWith("#") || line.startsWith(">")) {
						
						if(resultsString != null) {
							
							if (line.startsWith("#")) {
								
								ConservationClient.parseResults(resultsString, resultMap);
							}
							
							if (line.startsWith(">")) {
								
								ConservationClient.parseSequences(seqString, seqList, seqNames, pattern);
							}
							
							if(line == null) {
								
								if (previousLine.startsWith("#")) {
									
									ConservationClient.parseResults(resultsString, resultMap);
								}
								
								if (previousLine.startsWith(">")) {
									
									ConservationClient.parseSequences(seqString, seqList, seqNames, pattern);
								}
								
							}
						}
						
					}
					
					else {
						
						resultsString += line;
					}

					
			} while (line != null);
			
				inResults.close();
				
				return resultMap;
		}
		
		
		static void parseResults (String resultsString, Map<Method, double[]>resultMap) {
			
			String[] results = resultsString.split(" ");
			
			String name = results[0].trim().substring(1); 
			
			double[] resultsNum = new double[results.length - 1];
			
			for (int i = 0; i < resultsNum.length; i++) {
				
				resultsNum[i] = Double.parseDouble(results[i + 1].trim());
			}
			
			resultMap.put(Method.getMethod(name), resultsNum);
			
		}
		
		static void parseSequences (String seqString, List<char[]> seqList, List<String> seqNames, Pattern pattern) {
			
			StringTokenizer tokens = new StringTokenizer(seqString, " ");
			
			String name = tokens.nextToken().trim();
			
			Pattern pattern2 = Pattern.compile(name);
			
			String seqStringMod = pattern2.matcher(seqString).replaceFirst("");
			
			String seqStringMod2 = pattern.matcher(seqStringMod).replaceAll("");
			
			char[] seq = seqStringMod2.toCharArray();
			
			seqList.add(seq);
			
			seqNames.add(name);
		}


}
