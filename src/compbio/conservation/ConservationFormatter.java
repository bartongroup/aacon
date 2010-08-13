package compbio.conservation;

import java.util.*;
import java.io.*;

public class ConservationFormatter {
	
	/**
	 * Formats results
	 * @param <T> 
	 * @param tag any object(usually used for members of enumeration)
	 * @param result array of results to be printed
	 * @param resultPrecision
	 * @param print reference to PrintWriter object
	 */
	
	private static <T> void formatResult(T tag,  double[] result, int resultPrecision, PrintWriter print) {
		
		String tagFormat = "%s";
		
		String resultFormat = "%." + resultPrecision + "f";
		
		print.printf(tagFormat, "#" + tag.toString());
		
		for (int i = 0; i < result.length; i++) {
			
			print.printf(resultFormat, result[i] );
			print.print(" ");
		}
		
		print.println();
	}
	
	/**
	 * 
	 * @param <T>
	 * @param tag any object(usually used for members of enumeration)
	 * @param result array of results to be printed
	 * @param tagWidth
	 * @param resultWidth
	 * @param resultPrecision
	 * @param print reference to PrintWriter object
	 */
	
	private static <T> void formatResultWithAlignment(T tag,  double[] result, int tagWidth, int resultWidth, int resultPrecision, PrintWriter print) {
		
		String tagFormat = "%-" + tagWidth + "s";
		
		String resultFormat = "%-" + resultWidth + "." + resultPrecision + "f";
		
		print.printf(tagFormat, "#" + tag.toString());
		
		for (int i = 0; i < result.length; i++) {
			
			print.printf(resultFormat, result[i] );
		}
		
		print.println();
	}
	
	/**
	 * 
	 * @param <T>
	 * @param alignment 
	 * @param tag any object(usually used for members of enumeration)
	 * @param result
	 * @param tagWidth
	 * @param resultWidth
	 * @param resultPrecision
	 * @param outputFile
	 */
	static <T> void printResultWithAlignment(AminoAcidMatrix alignment, T tag,  double[] result, int tagWidth, int resultWidth, int resultPrecision, String outputFile) {
		
		PrintWriter print = null;
		
		try {
			
			print = new PrintWriter( new BufferedWriter (new FileWriter(outputFile, true)));
		}
		
		catch(IOException ex) {
			
			System.out.println("Problem writing" + outputFile);
			
			System.exit(0);
			
		}
		
		formatResultWithAlignment(tag, result, tagWidth, resultWidth, resultPrecision, print);
		
		print.close();
		
		
	}
	/**
	 * 
	 * @param <T>
	 * @param alignment
	 * @param tag any object(usually used for members of enumeration)
	 * @param result
	 * @param resultPrecision
	 * @param outputFile
	 * @param append
	 */
	
	static <T> void printResultNoAlignment(AminoAcidMatrix alignment, T tag,  double[] result, int resultPrecision, String outputFile, boolean append) {
		
		PrintWriter print = null;
		
		try {
			
			print = new PrintWriter( new BufferedWriter (new FileWriter(outputFile, append)));
		}
		
		catch(IOException ex) {
			
			System.out.println("Problem writing" + outputFile);
			
		}
		
		formatResult(tag, result, resultPrecision, print);
		
		print.close();
		
		
	}
	
	/**
	 * Formats and prints results
	 * @param scores
	 * @param outFilePath
	 * @param format
	 */
	
	static void formatResults (Map<Method, double[]> scores, String outFilePath, String format, AminoAcidMatrix alignment ) {
		
		int precision = 3;
		
		Iterator<Method> itr = scores.keySet().iterator();
		
		PrintWriter print = null;

		if (outFilePath == null) {
			
			while(itr.hasNext()) {
				
				Method key = itr.next();
				
				System.out.print("#" + key.toString() + " ");
				
				ConservationAccessory.printArrayOfDouble(scores.get(key), print, precision);
				
			}
			
		}
		
		else {
				
				if(format == null) { 
					
					print = openPrintWriter(outFilePath, false);
					
					if (print != null) {
					
						while(itr.hasNext()) {
						
							Method key = itr.next();
						
							print.print("#" + key.toString() + " ");
						
							ConservationAccessory.printArrayOfDouble(scores.get(key), print, precision);
						
						}
					
					}
				}
				
				else {
					
					if(Format.getFormat(format) == Format.RESULT_NO_ALIGNMENT) {
						
						print = openPrintWriter(outFilePath, false);
						
						if (print != null) {
								
							while(itr.hasNext()) {
									
								Method key = itr.next();
									
								print.print("#" + key.toString() + " ");
									
								ConservationAccessory.printArrayOfDouble(scores.get(key), print, precision);
									
							}
							
						}
						
					}
			
					if (Format.getFormat(format) == Format.RESULT_WITH_ALIGNMENT) {
						
						print = openPrintWriter(outFilePath, true);
						
						int tagWidth = 30;
						
						int fieldWidth = 10;
						
						String tagFormat = "%-" + tagWidth + "s";
						
						alignment.printAlignment(tagWidth, fieldWidth, outFilePath);
						
						if (print != null) {
							
							print.println();
								
							while(itr.hasNext()) {
									
								Method key = itr.next();
									
								print.printf(tagFormat, "#" + key.toString() + " ");
									
								ConservationAccessory.printArrayOfDoubleFieldWidth(scores.get(key), print, precision, fieldWidth);
									
							}
							
						}
						
					}
				
				}
				
				
				print.close();
		}
			
	}
		
	
	static PrintWriter openPrintWriter (String outFilePath, boolean append) {
		
		PrintWriter print = null;
		
		try {
			
			print = new PrintWriter( new BufferedWriter (new FileWriter(outFilePath, append)));
		}
		
		catch(IOException ex) {
			
			System.out.println("Problem writing" + outFilePath);
			
			print = null;
			
			//System.exit(0);
			
		}
		
		return print;
		
	}

}

	


