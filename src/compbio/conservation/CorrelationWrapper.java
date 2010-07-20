package compbio.conservation;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.List;
import java.util.*;

import compbio.util.FastaSequence;
import compbio.util.SequenceUtil;

public class CorrelationWrapper {
	
	static void printCoeffs(double[] coeffs, double[] coeffs2) {
		
		//assert coeffs.length == coeffs2.length;
		
		for (int i = 0; i < coeffs.length; i++){
			
			System.out.println( i + ":" + coeffs[i] + "-" + coeffs[i]);
		}
		
	}

	public static void main (String[] args) {

		
		String filePath  = "/homes/agolicz/alignments/alignment1";
		
		InputStream inStr = null;
		
		List<FastaSequence> fastaSeqs = null;
		
		try {
			
			inStr = new FileInputStream(filePath);
			
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
			
			
		AminoAcidMatrix matrix = new AminoAcidMatrix(fastaSeqs);
		
		// FIXME the width should not be hard coded it will come form the command line
		
		double[] coeffs = Correlation.calcPearson4(matrix, 7);
		
		double[] coeffs2 = Correlation.calcPearsonCoeff3(matrix, 7);
		
		//int[][] localSim = Correlation.localSimilarity2(matrix);
		
		//int[] globalSim = Correlation.globalSimilarity(matrix);
		
		//printCoeffs(globalSim);
		
		//printCoeffs(localSim[0]);
		
		printCoeffs(coeffs, coeffs2);
		
		System.out.println(matrix.numberOfColumns());
	}		
}
