package compbio.conservation;

import java.util.*;
import java.io.*;

public class ConservationFormatter {
	
	static <T> void formatResultWithNumbers(T tag,  double[] result, int tagWidth, int resultWidth, int resultPrecision, PrintWriter print) {
		
		String tagFormat = "%-" + tagWidth + "s";
		
		String resultFormat = "%-" + resultWidth + "." + resultPrecision + "f";
		
		print.printf(tagFormat, tag.toString());
		
		print.println();
		
		for (int i = 0; i < result.length; i++) {
			
			print.print("Columnn number: " + i + ": ");
			print.printf(resultFormat, result[i] );
			print.println();
		}
		
		print.println();
	}
	
	static <T> void formatResultNoNumbers(T tag,  double[] result, int tagWidth, int resultWidth, int resultPrecision, PrintWriter print) {
		
		String tagFormat = "%-" + tagWidth + "s";
		
		String resultFormat = "%-" + resultWidth + "." + resultPrecision + "f";
		
		print.printf(tagFormat, tag.toString());
		
		for (int i = 0; i < result.length; i++) {
			
			print.printf(resultFormat, result[i] );
		}
		
		print.println();
	}
	
	static <T> void printResultWithAlignment(AminoAcidMatrix alignment, T tag,  double[] result, int tagWidth, int resultWidth, int resultPrecision, String outputFile) {
		
		boolean first = true;
		
		PrintWriter print = null;
		
		try {
			
			print = new PrintWriter( new BufferedWriter (new FileWriter(outputFile)));
		}
		
		catch(IOException ex) {
			
			System.out.println("Problem writing" + outputFile);
			
		}
		
		if (first) {
			
			alignment.printAlignment(tagWidth, resultWidth, print);
			
			first = false;
		}
		
		formatResultNoNumbers(tag, result,tagWidth, resultWidth, resultPrecision, print);
		
		print.close();
		
		
	}
	
	static <T> void printResultNoAlignment(AminoAcidMatrix alignment, T tag,  double[] result, int tagWidth, int resultWidth, int resultPrecision, String outputFile) {
		
		PrintWriter print = null;
		
		try {
			
			print = new PrintWriter( new BufferedWriter (new FileWriter(outputFile)));
		}
		
		catch(IOException ex) {
			
			System.out.println("Problem writing" + outputFile);
			
		}
		
		formatResultWithNumbers(tag, result, tagWidth, resultWidth, resultPrecision, print);
		
		print.close();
		
		
	}

}

	


