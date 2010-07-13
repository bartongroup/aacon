package compbio.conservation;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.List;

import compbio.util.FastaSequence;
import compbio.util.SequenceUtil;

public class ConservationScoresWrapper {
	
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
			
		ConservationScores scores = new ConservationScores(matrix);
		
		double[] kabat = scores.kabatScore(); 
		
		double[] jores = scores.joresScore();
		
		double[] schneider = scores.schneiderScore();
		
		double[] shenkin = scores.shenkinScore();
		
		double[] gerstein = scores.gersteinScore();
		
		double[] taylorNoGaps = scores.SmallestTaylorSetNoGaps();
		
		double[] taylorGaps = scores.SmallestTaylorSetGaps();
		
		double[] zvelibil = scores.zvelibilScore();
		
		double[] karlin = scores.karlinScore();
		
		double[] armon = scores.armonScore();
		
		double[] thompson = scores.thompsonScore();
		
		double[] lancet  = scores.notLancetScore();
		
		double[] mirny = scores.mirnyScore();
		
		double[] williamson = scores.williamsonScore();
		
		double[] landgraf = scores.landgrafScore();
		
		double[] sander = scores.landgrafScore();
		
		double[] valdar = scores.valdarScore();
		
		}
		
}
