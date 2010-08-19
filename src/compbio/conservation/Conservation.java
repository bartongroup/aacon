package compbio.conservation;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;

import compbio.data.sequence.FastaSequence;
import compbio.data.sequence.SequenceUtil;

/**
 * TODO to complete
 * 
 * API for external users
 * 
 * @author pvtroshin
 * 
 */
public class Conservation {

    public static void calculate() throws FileNotFoundException, IOException {
	List<FastaSequence> sequences = SequenceUtil
		.readFasta(new FileInputStream(new File(
			SlowMethodTester.DATA_PATH + File.separator + "file")));

	AminoAcidMatrix alignment = new AminoAcidMatrix(sequences, null);
	ConservationScores2 scores = new ConservationScores2(alignment);
	double[] result = scores.calculateScore(Method.SANDER, true);
    }

}
