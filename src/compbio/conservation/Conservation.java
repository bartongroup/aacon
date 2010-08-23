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

	private AminoAcidMatrix alignMatrix;
	private Map<Method, double[]> results;

	public Conservation() {
	}

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
	public Map<Method, double[]> getConservation(File file,
			EnumSet<Method> methods, boolean normilizeResults)
			throws FileNotFoundException, IOException,
			UnknownFileFormatException {

		if (file == null) {
			throw new NullPointerException("File must be provided!");
		}
		// there is no need to close input stream as the read method will close
		// it
		FileInputStream fis = new FileInputStream(file);
		// the method closes the input stream
		boolean isClustalFile = ClustalAlignmentUtil.isValidClustalFile(fis);
		fis = new FileInputStream(file);
		if (isClustalFile) {
			Alignment alignment = ClustalAlignmentUtil.readClustalFile(fis);
			assert alignment != null : "Fails to read the alignement!";
			alignMatrix = new AminoAcidMatrix(alignment);
		} else {
			// assume the file contain a list of fasta sequences then
			List<FastaSequence> sequences = SequenceUtil.readFasta(fis);
			alignMatrix = new AminoAcidMatrix(sequences, null);
		}
		results = calculateConservation(methods, normilizeResults);
		return results;
	}

	/**
	 * @param alignment
	 * @param methods
	 * @param normilizeResults
	 * @return
	 */
	public Map<Method, double[]> getConservation(Alignment alignment,
			EnumSet<Method> methods, boolean normilizeResults) {

		if (alignment == null) {
			throw new NullPointerException("Alignment must be provided!");
		}
		// there is no need to close input stream as the read method will close
		// it
		alignMatrix = new AminoAcidMatrix(alignment);
		results = calculateConservation(methods, normilizeResults);
		return results;
	}

	/**
	 * @param sequences
	 * @param methods
	 * @param normilizeResults
	 * @return
	 */
	public Map<Method, double[]> getConservation(List<FastaSequence> sequences,
			EnumSet<Method> methods, boolean normilizeResults) {

		if (sequences == null || sequences.isEmpty()) {
			throw new NullPointerException("Sequences must be provided!");
		}
		alignMatrix = new AminoAcidMatrix(sequences, null);
		results = calculateConservation(methods, normilizeResults);
		return results;
	}

	private Map<Method, double[]> calculateConservation(
			EnumSet<Method> methods, boolean normilizeResults) {

		ConservationScores2 scores = new ConservationScores2(alignMatrix);
		Map<Method, double[]> result = new EnumMap<Method, double[]>(
				Method.class);
		for (Method method : methods) {
			double[] singleRes = scores
					.calculateScore(method, normilizeResults);
			assert singleRes != null && singleRes.length > 0;
			result.put(method, singleRes);
		}
		this.results = result;
		return result;
	}

	public void outputResults(File outFile, Format format) throws IOException {
		ConservationFormatter.formatResults(results, outFile.getAbsolutePath(),
				format, alignMatrix);
	}

	public void printResults() {
		try {
			ConservationFormatter.formatResults(results, null, null, null);
		} catch (IOException ignored) {
			// this will never happen as no writing to real file happens
			// in the call to the function above
			ignored.printStackTrace();
		}
	}

}
