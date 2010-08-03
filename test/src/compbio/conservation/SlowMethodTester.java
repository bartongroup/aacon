package compbio.conservation;

import static org.testng.Assert.fail;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.testng.annotations.Test;

import compbio.util.FastaSequence;
import compbio.util.SequenceUtil;
import compbio.util.SequenceUtilTester;
import compbio.util.Timer;

public class SlowMethodTester {

    static final String SMALL_AL = "TO1296.fasta.align";
    static final String AVG_AL = "avg.aln.fa";
    static final String LARGE_AL = "1000x3000Dna.aln.fa";

    @Test
    public void testSadler() {
	try {
	    Timer timer = new Timer(TimeUnit.MILLISECONDS);
	    List<FastaSequence> sequences = SequenceUtil
		    .readFasta(new FileInputStream(new File(
			    SequenceUtilTester.DATA_PATH + File.separator
				    + AVG_AL)));
	    System.out.println("Loading sequences: " + timer.getStepTime());

	    AminoAcidMatrix alignment = new AminoAcidMatrix(sequences);
	    System.out.println("Converting to Matrix: " + timer.getStepTime());

	    ConservationScores2 scores = new ConservationScores2(alignment);
	    System.out.println("Constructing conservation scores: "
		    + timer.getStepTime());

	    double[] result = scores.calculateScore(Method.SANDER, true);
	    System.out.println("Calculating sadler scores: "
		    + timer.getStepTime());

	    // ConservationScore2Tester.printScores(result, "Sander");
	    System.out.println("Total: " + timer.getTotalTime());

	} catch (FileNotFoundException e) {
	    e.printStackTrace();
	    fail(e.getLocalizedMessage());
	} catch (IOException e) {
	    e.printStackTrace();
	    fail(e.getLocalizedMessage());
	}
    }

}
