package compbio.conservation;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.EnumMap;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import compbio.data.sequence.Alignment;
import compbio.data.sequence.ClustalAlignmentUtil;
import compbio.data.sequence.FastaSequence;
import compbio.data.sequence.SequenceUtil;
import compbio.data.sequence.UnknownFileFormatException;

/**
 * Calculates conservation
 * 
 * @author agolicz & pvtroshin
 */
public class Conservation {

	private AminoAcidMatrix alignMatrix;
	private final boolean normalize;
	private final Map<Method, double[]> results;

	private final ExecutorService executor;

	Conservation(AminoAcidMatrix alignment, boolean normalize,
			ExecutorService executor) {
		this.alignMatrix = alignment;
		this.normalize = normalize;
		this.executor = executor;
		this.results = new EnumMap<Method, double[]>(Method.class);
	}

	AminoAcidMatrix getAlignment() {
		return alignMatrix;
	}

	Map<Method, double[]> calculateScores(final EnumSet<Method> methods) {
		return calculateConservation(methods);
	}

	/**
	 * @param method
	 * @param normalize
	 * @return score for the given method
	 */
	synchronized double[] calculateScore(final Method method) {
		List<Callable<Object>> tasks = new ArrayList<Callable<Object>>();
		double[] result = new double[alignMatrix.numberOfColumns()];
		double[] normalized = null;
		switch (method) {
		case KABAT:
			for (int i = 0; i < alignMatrix.numberOfColumns(); i++) {
				result[i] = ColumnScores.kabatScore(alignMatrix, i);
			}
			if (normalize) {
				normalized = ConservationAccessory.inversedNormalize01(result,
						method);
			}
			break;
		case JORES:
			for (int i = 0; i < alignMatrix.numberOfColumns(); i++) {
				result[i] = ColumnScores.joresScore(alignMatrix, i);
			}
			if (normalize) {
				normalized = ConservationAccessory.inversedNormalize01(result,
						method);
			}
			break;
		case SCHNEIDER:
			for (int i = 0; i < alignMatrix.numberOfColumns(); i++) {
				result[i] = ColumnScores.schneiderScore(alignMatrix, i);
			}
			if (normalize) {
				normalized = ConservationAccessory.inversedNormalize01(result,
						method);
			}
			break;
		case SHENKIN:
			for (int i = 0; i < alignMatrix.numberOfColumns(); i++) {
				result[i] = ColumnScores.shenkinScore(alignMatrix, i);
			}
			if (normalize) {
				normalized = ConservationAccessory.inversedNormalize01(result,
						method);
			}
			break;
		case GERSTEIN:
			for (int i = 0; i < alignMatrix.numberOfColumns(); i++) {
				result[i] = ColumnScores.gersteinScore(alignMatrix, i);
			}
			if (normalize) {
				normalized = ConservationAccessory.inversedNormalize01(result,
						method);
			}
			break;
		case TAYLOR_GAPS:
			for (int i = 0; i < alignMatrix.numberOfColumns(); i++) {
				result[i] = ColumnScores.taylorScoreGaps(alignMatrix, i);
			}
			if (normalize) {
				normalized = ConservationAccessory.inversedNormalize01(result,
						method);
			}
			break;
		case TAYLOR_NO_GAPS:
			for (int i = 0; i < alignMatrix.numberOfColumns(); i++) {
				result[i] = ColumnScores.taylorScoreNoGaps(alignMatrix, i);
			}
			if (normalize) {
				normalized = ConservationAccessory.inversedNormalize01(result,
						method);
			}
			break;
		case ZVELIBIL:
			for (int i = 0; i < alignMatrix.numberOfColumns(); i++) {
				result[i] = ColumnScores.zvelibilScore(alignMatrix, i);
			}
			if (normalize) {
				normalized = ConservationAccessory.normalize01(result, method);
			}
			break;
		case KARLIN:
			for (int i = 0; i < alignMatrix.numberOfColumns(); i++) {
				tasks.add(Executors.callable(new TaskRunner(i, Method.KARLIN,
						result)));
			}
			executeAndWait(tasks);
			if (normalize) {
				normalized = ConservationAccessory.normalize01(result, method);
			}
			break;
		case ARMON:
			for (int i = 0; i < alignMatrix.numberOfColumns(); i++) {
				result[i] = ColumnScores.armonScore(alignMatrix, i);
			}
			if (normalize) {
				normalized = ConservationAccessory.inversedNormalize01(result,
						method);
			}
			break;
		case THOMPSON:
			for (int i = 0; i < alignMatrix.numberOfColumns(); i++) {
				result[i] = ColumnScores.thompsonScore(alignMatrix, i);
			}
			if (normalize) {
				normalized = ConservationAccessory.inversedNormalize01(result,
						method);
			}
			break;
		case NOT_LANCET:
			for (int i = 0; i < alignMatrix.numberOfColumns(); i++) {
				result[i] = ColumnScores.notLancetScore(alignMatrix, i);
			}
			if (normalize) {
				normalized = ConservationAccessory.normalize01(result, method);
			}
			break;
		case MIRNY:
			for (int i = 0; i < alignMatrix.numberOfColumns(); i++) {
				result[i] = ColumnScores.mirnyScore(alignMatrix, i);
			}
			if (normalize) {
				normalized = ConservationAccessory.normalize01(result, method);
			}
			break;
		case WILLIAMSON:
			for (int i = 0; i < alignMatrix.numberOfColumns(); i++) {
				result[i] = ColumnScores.williamsonScore(alignMatrix, i);
			}
			if (normalize) {
				normalized = ConservationAccessory.normalize01(result, method);
			}
			break;
		case LANDGRAF:
			for (int i = 0; i < alignMatrix.numberOfColumns(); i++) {
				tasks.add(Executors.callable(new TaskRunner(i, Method.LANDGRAF,
						result)));
			}
			executeAndWait(tasks);
			if (normalize) {
				normalized = ConservationAccessory.inversedNormalize01(result,
						method);
			}
			break;
		case SANDER:
			for (int i = 0; i < alignMatrix.numberOfColumns(); i++) {
				tasks.add(Executors.callable(new TaskRunner(i, Method.SANDER,
						result)));
			}
			executeAndWait(tasks);
			if (normalize) {
				normalized = ConservationAccessory.normalize01(result, method);
			}
			break;
		case VALDAR:
			for (int i = 0; i < alignMatrix.numberOfColumns(); i++) {
				tasks.add(Executors.callable(new TaskRunner(i, Method.VALDAR,
						result)));
			}
			executeAndWait(tasks);
			if (normalize) {
				normalized = ConservationAccessory.normalize01(result, method);
			}
			break;
		case SMERFS:
			result = getSMERFS(SMERFSColumnScore.DEFAULT_WINDOW_SIZE,
					SMERFSColumnScore.MID_SCORE,
					SMERFSColumnScore.DEFAULT_GAP_THRESHOLD);
			break;
		default:
			throw new RuntimeException("You should never ever get here");
		} // end of switch

		if (normalized != null) {
			results.put(method, normalized);
			return normalized;
		}
		// store results in this object
		results.put(method, result);
		return result;
	}

	private final void executeAndWait(List<Callable<Object>> tasks) {
		try {
			executor.invokeAll(tasks);
		} catch (InterruptedException e) {
			throw new RuntimeException(
					"The program was stopped in the middle of the calculation");
		}
	}

	private final class TaskRunner implements Runnable {

		private final int iteration;
		private final Method method;
		private final double[] result;

		public TaskRunner(final int i, Method method, double[] result) {
			this.iteration = i;
			this.method = method;
			this.result = result;
		}

		@Override
		public void run() {
			switch (method) {
			case KARLIN:
				result[iteration] = ColumnScores.karlinScore(alignMatrix,
						iteration);
				break;
			case VALDAR:
				result[iteration] = ColumnScores.valdarScore(alignMatrix,
						iteration);
				break;
			case LANDGRAF:
				result[iteration] = ColumnScores.landgrafScore(alignMatrix,
						iteration);
				break;
			case SANDER:
				result[iteration] = ColumnScores.sanderScore(alignMatrix,
						iteration);
				break;
			default:
				throw new IllegalArgumentException(
						"Only KARLIN, VALDAR, LANDGRAF and "
								+ "SANDER methods can be executed via TaskRunner!");
			}

		}
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
	public static Conservation getConservation(File file, boolean normalize,
			ExecutorService executor) throws FileNotFoundException,
			IOException, UnknownFileFormatException {

		if (file == null) {
			throw new NullPointerException("File must be provided!");
		}
		AminoAcidMatrix alignMatrix = null;
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

		return new Conservation(alignMatrix, normalize, executor);
	}

	/**
	 * @param alignment
	 * @param methods
	 * @param normilizeResults
	 * @return
	 */
	public Map<Method, double[]> getConservation(Alignment alignment,
			EnumSet<Method> methods) {

		if (alignment == null) {
			throw new NullPointerException("Alignment must be provided!");
		}
		// there is no need to close input stream as the read method will close
		// it
		alignMatrix = new AminoAcidMatrix(alignment);
		calculateConservation(methods);
		return Collections.unmodifiableMap(results);
	}

	/**
	 * @param sequences
	 * @param methods
	 * @param normilizeResults
	 * @return
	 */
	public Map<Method, double[]> getConservation(List<FastaSequence> sequences,
			EnumSet<Method> methods) {

		if (sequences == null || sequences.isEmpty()) {
			throw new NullPointerException("Sequences must be provided!");
		}
		alignMatrix = new AminoAcidMatrix(sequences, null);
		calculateConservation(methods);
		return results;
	}

	private synchronized Map<Method, double[]> calculateConservation(
			EnumSet<Method> methods) {
		for (Method method : methods) {
			double[] singleRes = calculateScore(method);
			assert singleRes != null && singleRes.length > 0;
			results.put(method, singleRes);
		}
		return Collections.unmodifiableMap(results);
	}

	void outputResults(File outFile, Format format) throws IOException {
		ConservationFormatter.formatResults(results, outFile.getAbsolutePath(),
				format, alignMatrix);
	}

	void printResults(Format format) throws IOException {
		ConservationFormatter.formatResults(results, null, format, alignMatrix);
	}

	public static void printResults(Map<Method, double[]> results) {
		try {
			ConservationFormatter.formatResults(results, null,
					Format.RESULT_NO_ALIGNMENT, null);
		} catch (IOException ignored) {
			// this will never happen as no writing to real file happens
			// in the call to the function above
			ignored.printStackTrace();
		}
	}

	public static void printResults(double[] result, Method method) {
		Map<Method, double[]> results = new EnumMap<Method, double[]>(
				Method.class);
		results.put(method, result);
		printResults(results);
	}

	/**
	 * Returns results of SMERFS calculation or null, if parameters provided are
	 * not appropriate.
	 * 
	 * @param alignment
	 *            reference to alignment
	 * @param width
	 *            with of the window
	 * @param score
	 *            tells which score given to the column, either the highest
	 *            score of all the windows it belongs to, or the middle column
	 *            is given the score of the window.
	 * @param normalize
	 *            if true results will be normalized
	 * @return TODO refactor?
	 */
	public double[] getSMERFS(int width, SMERFSColumnScore score,
			double gapTreshold) {

		if (alignMatrix == null) {
			throw new IllegalArgumentException("Matrix must not be null.");
		}
		double[] result = null;
		if (width <= 0 || width % 2 != 1
				|| width > alignMatrix.numberOfColumns() || score == null
				|| gapTreshold < 0 || gapTreshold > 1) {
			if (width <= 0 || width % 2 != 1) {
				System.out
						.println("Column width for SMERFS not provided or "
								+ "smaller or equal zero or not an odd number or not an integer.");
			}
			if (width > alignMatrix.numberOfColumns()) {
				System.out
						.println("Column width greater than the length of the alignment");
			}
			if (score == null) {
				System.out
						.println("Column score not privided or the type provided is not supported.");
				SMERFSColumnScore.supportedSMERFSColumnSores();
			}
			if (gapTreshold < 0 || gapTreshold > 1) {
				System.out
						.println("Gap treshold could not have been parsed as a double "
								+ "or it was smaller than zero or it was greater than one.");
			}
			return result;
		}
		Correlation corr = new Correlation(alignMatrix, width, gapTreshold,
				executor);

		try {
			result = corr.getCorrelationScore(score, normalize);
		} catch (InterruptedException e) {
			e.printStackTrace();
			throw new RuntimeException("Calculation was interrupted!");
		} catch (ExecutionException e) {
			throw new RuntimeException("Error during calculation!"
					+ e.getLocalizedMessage(), e.getCause());
		}

		return result;
	}

}
