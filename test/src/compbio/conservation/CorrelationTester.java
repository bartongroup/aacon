package compbio.conservation;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.List;

import compbio.util.FastaSequence;
import compbio.util.SequenceUtil;

//import org.testng.Assert;
import org.testng.annotations.Test;


public class CorrelationTester {
	
	@Test
	void Corr1Tester() {
		
		String filePath  = "/homes/agolicz/alignments/prot2";
		
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
		
		double[] coeffs = Correlation.calcPearson4(matrix);
		
		//double[] coeffs2 = Correlation.calcPearsonCoeff3(matrix);
		
		//int[][] localSim = Correlation.localSimilarity2(matrix);
		
		//int[] globalSim = Correlation.globalSimilarity(matrix);
		
		//printCoeffs(globalSim);
		
		//printCoeffs(localSim[0]);
		
	}

}
