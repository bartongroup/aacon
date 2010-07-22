package compbio.conservation;

import java.util.*;
import java.io.*;

public class ConservationFormatter {
	
	static <T> void formatResultWithHashSign(T tag,  double[] result, int resultPrecision, PrintWriter print) {
		
		String tagFormat = "%s";
		
		String resultFormat = "%." + resultPrecision + "f";
		
		print.printf(tagFormat, "# " + tag.toString());
		
		print.println();
		
		for (int i = 0; i < result.length; i++) {
			
			print.printf(resultFormat, result[i] );
			print.print(" ");
		}
		
		print.println();
	}
	
	static <T> void formatResult(T tag,  double[] result, int tagWidth, int resultWidth, int resultPrecision, PrintWriter print) {
		
		String tagFormat = "%-" + tagWidth + "s";
		
		String resultFormat = "%-" + resultWidth + "." + resultPrecision + "f";
		
		print.printf(tagFormat, tag.toString());
		
		for (int i = 0; i < result.length; i++) {
			
			print.printf(resultFormat, result[i] );
		}
		
		print.println();
	}
	
	static <T> void printResultWithAlignment(AminoAcidMatrix alignment, T tag,  double[] result, int tagWidth, int resultWidth, int resultPrecision, String outputFile) {
		
		alignment.printAlignment(tagWidth, resultWidth, outputFile);
		
		PrintWriter print = null;
		
		try {
			
			print = new PrintWriter( new BufferedWriter (new FileWriter(outputFile, true)));
		}
		
		catch(IOException ex) {
			
			System.out.println("Problem writing" + outputFile);
			
			System.exit(0);
			
		}
		
		formatResult(tag, result,tagWidth, resultWidth, resultPrecision, print);
		
		print.close();
		
		
	}
	
	static <T> void printResultNoAlignment(AminoAcidMatrix alignment, T tag,  double[] result, int resultPrecision, String outputFile, boolean append) {
		
		PrintWriter print = null;
		
		try {
			
			print = new PrintWriter( new BufferedWriter (new FileWriter(outputFile, append)));
		}
		
		catch(IOException ex) {
			
			System.out.println("Problem writing" + outputFile);
			
		}
		
		formatResultWithHashSign(tag, result, resultPrecision, print);
		
		print.close();
		
		
	}

}

	


