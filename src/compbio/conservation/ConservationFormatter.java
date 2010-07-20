package compbio.conservation;

import java.util.*;
import java.io.*;

public class ConservationFormatter {
	
	static void formatResult(Method method,  double[] result, int tagWidth, int resultWidth, int resultPrecision, PrintWriter print) {
		
		String tagFormat = "%-" + tagWidth + "s";
		
		String resultFormat = "%-" + resultWidth + "." + resultPrecision + "f";
		
		print.printf(tagFormat, method.toString());
		
		for (int i = 0; i < result.length; i++) {
			
			print.printf(resultFormat, result[i] );
		}
		
		print.println("");
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
		
		formatResult(method, result,tagWidth, resultWidth, resultPrecision, print);
		
		
	}
	
	static void printResultNoAlignment(AminoAcidMatrix alignment, Method method,  double[] result, int tagWidth, int resultWidth, int resultPrecision, String outputFile) {
		
		PrintWriter print = null;
		
		try {
			
			print = new PrintWriter( new BufferedWriter (new FileWriter(outputFile)));
		}
		
		catch(IOException ex) {
			
			System.out.println("Problem writing" + outputFile);
			
		}
		
		formatResult(method, result, tagWidth, resultWidth, resultPrecision, print);
		
		
	}

}

	


