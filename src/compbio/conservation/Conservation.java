package compbio.conservation;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.EnumMap;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;

import compbio.data.sequence.Alignment;
import compbio.data.sequence.ClustalAlignmentUtil;
import compbio.data.sequence.FastaSequence;
import compbio.data.sequence.SequenceUtil;
import compbio.data.sequence.UnknownFileFormatException;

/**
 * TODO to complete API for external users
 * 
 * @author pvtroshin
 */
public class Conservation {

	/**
	 * Read either clustal formatted alignment or list of fasta formatted
	 * sequences, aligned sequences. Return the Map Method name->double[]
	 * conservation prediction results
	 * 
	 * @param file
	 * @param methods
	 * @param normilizeResults
	 * @return
	 * @throws FileNotFoundException
	 * @throws IOException
	 * @throws UnknownFileFormatException
	 */
	public static Map<Method, double[]> getConservation(File file, EnumSet<Method> methods,
			boolean normilizeResults) throws FileNotFoundException, IOException,
			UnknownFileFormatException {

		AminoAcidMatrix alignmatrix = null;
		if (file == null) {
			throw new NullPointerException("File must be provided!");
		}
		// there is no need to close input stream as the read method will close
		// it
		FileInputStream fis = new FileInputStream(file);
		if (ClustalAlignmentUtil.isValidClustalFile(fis)) {
			Alignment alignment = ClustalAlignmentUtil.readClustalFile(fis);
			assert alignment != null : "Fails to read the alignement!";
			alignmatrix = new AminoAcidMatrix(alignment);
		} else {
			// assume the file contain a list of fasta sequences then
			List<FastaSequence> sequences = SequenceUtil.readFasta(fis);
			alignmatrix = new AminoAcidMatrix(sequences, null);
		}
		return calculateConservation(alignmatrix, methods, normilizeResults);
	}

	/**
	 * @param alignment
	 * @param methods
	 * @param normilizeResults
	 * @return
	 */
	public static Map<Method, double[]> getConservation(Alignment alignment,
			EnumSet<Method> methods, boolean normilizeResults) {

		AminoAcidMatrix alignmatrix = null;
		if (alignment == null) {
			throw new NullPointerException("Alignment must be provided!");
		}
		// there is no need to close input stream as the read method will close
		// it
		alignmatrix = new AminoAcidMatrix(alignment);
		return calculateConservation(alignmatrix, methods, normilizeResults);
	}

	/**
	 * @param sequences
	 * @param methods
	 * @param normilizeResults
	 * @return
	 */
	public static Map<Method, double[]> getConservation(List<FastaSequence> sequences,
			EnumSet<Method> methods, boolean normilizeResults) {

		AminoAcidMatrix alignmatrix = null;
		if (sequences == null || sequences.isEmpty()) {
			throw new NullPointerException("Sequences must be provided!");
		}
		alignmatrix = new AminoAcidMatrix(sequences, null);
		return calculateConservation(alignmatrix, methods, normilizeResults);
	}

	private static Map<Method, double[]> calculateConservation(AminoAcidMatrix alignmatrix,
			EnumSet<Method> methods, boolean normilizeResults) {

		ConservationScores2 scores = new ConservationScores2(alignmatrix);
		Map<Method, double[]> result = new EnumMap<Method, double[]>(Method.class);
		for (Method method : methods) {
			double[] singleRes = scores.calculateScore(method, normilizeResults);
			assert singleRes != null && singleRes.length > 0;
			result.put(method, singleRes);
		}
		return result;
	}

	public static double[] getSMERFS(List<FastaSequence> sequences, EnumSet<Method> methods,
			boolean normilizeResults) {

		return null;
	}

}
