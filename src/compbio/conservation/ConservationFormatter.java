package compbio.conservation;

import java.util.*;
import java.io.*;

public class ConservationFormatter {
	
	static void formatResultWithNumbers(Method method,  double[] result, int tagWidth, int resultWidth, int resultPrecision, PrintWriter print) {
		
		String tagFormat = "%-" + tagWidth + "s";
		
		String resultFormat = "%-" + resultWidth + "." + resultPrecision + "f";
		
		print.printf(tagFormat, method.toString());
		
		print.println();
		
		for (int i = 0; i < result.length; i++) {
			
			print.print("Columnn number: " + i + ": ");
			print.printf(resultFormat, result[i] );
			print.println();
		}
		
		print.println();
	}
	
	static void formatResultNoNumbers(Method method,  double[] result, int tagWidth, int resultWidth, int resultPrecision, PrintWriter print) {
		
		String tagFormat = "%-" + tagWidth + "s";
		
		String resultFormat = "%-" + resultWidth + "." + resultPrecision + "f";
		
		print.printf(tagFormat, method.toString());
		
		for (int i = 0; i < result.length; i++) {
			
			print.printf(resultFormat, result[i] );
		}
		
		print.println();
	}
	
	static void printResultWithAlignment(AminoAcidMatrix alignment, Method method,  double[] result, int tagWidth, int resultWidth, int resultPrecision, String outputFile) {
		
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
		
		formatResultNoNumbers(method, result,tagWidth, resultWidth, resultPrecision, print);
		
		print.close();
		
		
	}
	
	static void printResultNoAlignment(AminoAcidMatrix alignment, Method method,  double[] result, int tagWidth, int resultWidth, int resultPrecision, String outputFile) {
		
		PrintWriter print = null;
		
		try {
			
			print = new PrintWriter( new BufferedWriter (new FileWriter(outputFile)));
		}
		
		catch(IOException ex) {
			
			System.out.println("Problem writing" + outputFile);
			
		}
		
		formatResultWithNumbers(method, result, tagWidth, resultWidth, resultPrecision, print);
		
		print.close();
		
		
	}

}

	


