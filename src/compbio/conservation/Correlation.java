package compbio.conservation;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import compbio.common.ColumnTooWideException;

/**
 * Class has static methods used to create similarity matrices, and correlation
 * matrices.
 * 
 * @author agolicz & Peter Troshin
 * 
 */
class Correlation {

	private final AminoAcidMatrix alignment;
	private final int winWidth;
	private final double gapTreshold;

	private double[] coeffs;
	private int[] global;
	// number of rows
	private final int numofSequences;
	private final ExecutorService executor;

	/**
	 * Calculates similarity between two sequences, similarity is calculated as
	 * a sum of blosum scores for pairs formed by corresponding amino acids in
	 * two sequences.
	 * 
	 * @param seq1
	 *            sequence 1
	 * @param seq2
	 *            sequence 2
	 * @return similarity score
	 */
	public Correlation(AminoAcidMatrix alignment, int winWidth,
			double gapTreshold) {
		// check preconditions
		if (alignment == null) {
			throw new IllegalArgumentException("Alignment must not be null.");
		}
		if (winWidth < 0 || winWidth % 2 != 1
				|| winWidth > alignment.numberOfColumns()) {
			throw new IllegalArgumentException(
					"ColWidth smaller than zero or an even number or largrt tah the number of columns");
		}
		if (winWidth > alignment.numberOfColumns()) {
			throw new ColumnTooWideException(
					"The width of the window is greater than the length of the allignment.");
		}
		this.alignment = alignment;
		this.winWidth = winWidth;
		this.gapTreshold = gapTreshold;
		this.numofSequences = alignment.numberOfRows();
		this.executor = ParallelConservationClient.getExecutor();
	}

	/**
	 * Calculates sequence similarity as a sum of blosum scores.
	 * 
	 * @param seq1
	 *            sequence 1
	 * @param seq2
	 *            sequence 2
	 * 
	 * @return similarity
	 */
	private static int sequenceSimilartyBlosum(char[] seq1, char[] seq2) {
		assert seq1.length == seq2.length;
		int similarity = 0;
		for (int i = 0; i < seq1.length; i++) {
			// it's the wrong blosum matrix, need the unmodulated one with
			// integers TODO FIXME!
			similarity += ConservationMatrices.BlosumPair(seq1[i], seq2[i]);
		}
		return similarity;
	}

	/**
	 * Creates a global similarity matrix. Always use this method along with
	 * local similarity matrix, the internal numbering scheme assures that
	 * Pearson correlation will be calculated properly
	 * 
	 * @param matrix
	 *            reference to the matrix containing sequences
	 * @return similarity matrix stored as a single array
	 * @throws ExecutionException
	 * @throws InterruptedException
	 * @throws ExecutionException
	 * @throws InterruptedException
	 */
	private void globalSimilarity() throws InterruptedException,
			ExecutionException {

		global = new int[numofSequences * (numofSequences - 1) / 2];
		int index = 0;

		List<Future<?>> tasks = new ArrayList<Future<?>>();

		for (int i = 0; i < numofSequences; i++) {
			tasks.add(executor.submit(new GlobalSimilarityWrapper(i, index)));
			index += numofSequences - i - 1;
		}
		// make sure all tasks are completed
		for (Future<?> future : tasks) {
			future.get();
		}

	}

	private class GlobalSimilarityWrapper implements Runnable {

		final int i;
		final int startIndex;
		final char[] alrow;

		public GlobalSimilarityWrapper(int i, int startIndex) {
			this.i = i;
			this.startIndex = startIndex;
			this.alrow = alignment.getRow(i);
		}

		@Override
		public void run() {
			int index = startIndex;
			for (int j = i + 1; j < numofSequences; j++) {
				global[index] = Correlation.sequenceSimilartyBlosum(alrow,
						alignment.getRow(j));
				index++;
			}
		}
	}

	/**
	 * Calculates correlation score for two vectors. If one of the vectors
	 * consists of single value only the correlation value is 0.
	 * 
	 * @param arr1
	 *            vector 1
	 * @param arr2
	 *            vector 2
	 * @return score
	 */
	private static double pearson(int[] arr1, int[] arr2) {
		assert arr1.length == arr2.length;
		int arr1Sum = 0;
		int arr2Sum = 0;
		for (int i = 0; i < arr1.length; i++) {
			arr1Sum += arr1[i];
			arr2Sum += arr2[i];
		}
		double arr1Ave = (double) arr1Sum / arr1.length;
		double arr2Ave = (double) arr2Sum / arr2.length;
		double sumProduct = 0;
		double a1Sum = 0;
		double a2Sum = 0;
		for (int i = 0; i < arr1.length; i++) {
			double s1 = arr1[i] - arr1Ave;
			double s2 = arr2[i] - arr2Ave;
			sumProduct += s1 * s2;
			a1Sum += s1 * s1;
			a2Sum += s2 * s2;
		}
		if (a1Sum == 0 || a2Sum == 0) {
			return 0;
		}
		assert Math.sqrt(a1Sum) * Math.sqrt(a2Sum) != 0;
		double result = sumProduct / (Math.sqrt(a1Sum) * Math.sqrt(a2Sum));
		return result;
	}

	/**
	 * Calculates correlation for a set number of windows.
	 * 
	 * @param nrOfWindows
	 * @param begin
	 * @param end
	 * @return
	 */
	private final int[][] localSimilarity2(final int nrOfWindows,
			final int begin, final int end) {

		final int[][] localSim = new int[nrOfWindows][numofSequences
				* (numofSequences - 1) / 2];
		// int res = (numofSequences * (numofSequences - 1) / 2);
		// System.out.println("!![]" + nrOfWindows + "[" + res + "]");

		int sum = 0;
		int globalIndex = 0;
		int windowNr = 0;

		for (int i = 0; i < numofSequences; i++) {
			final char[] rowI = alignment.getRow(i);
			for (int j = i + 1; j < numofSequences; j++) {
				final char[] rowJ = alignment.getRow(j);
				windowNr = 0;
				sum = 0;

				for (int z = begin; z < begin + winWidth; z++) {
					int index = 24 * ConservationMatrices.getIndex(rowI[z])
							+ ConservationMatrices.getIndex(rowJ[z]);
					sum += ConservationMatrices.blosum2[index];
				}
				localSim[windowNr][globalIndex] = sum;
				windowNr++;
				for (int k = begin + winWidth; k <= end; k++) {
					int index1 = 24 * ConservationMatrices.getIndex(rowI[k])
							+ ConservationMatrices.getIndex(rowJ[k]);
					int index2 = 24
							* ConservationMatrices.getIndex(rowI[k - winWidth])
							+ ConservationMatrices.getIndex(rowJ[k - winWidth]);
					sum = sum - ConservationMatrices.blosum2[index2]
							+ ConservationMatrices.blosum2[index1];
					localSim[windowNr][globalIndex] = sum;
					windowNr++;
				}
				globalIndex++;
			}
		}
		return localSim;
	}

	/**
	 * Also calculates correlation. Takes less memory than calcPearson3.
	 * 
	 * @return
	 * @throws InterruptedException
	 * @throws ExecutionException
	 */
	private double[] calcPearson() throws InterruptedException,
			ExecutionException {

		int nrOfWindows = ((alignment.numberOfColumns() - alignment
				.numberOfColumns()
				% winWidth)
				/ winWidth
				+ ((alignment.numberOfColumns() - alignment.numberOfColumns()
						% winWidth)
						/ winWidth - 1) * (winWidth - 1) + alignment
				.numberOfColumns()
				% winWidth);

		coeffs = new double[nrOfWindows];
		int coeffsIdx = 0;
		final int WIN_SIZE = 50;
		int tail = nrOfWindows % WIN_SIZE;
		final int turns = (nrOfWindows - tail) / WIN_SIZE;
		int start = 0;
		int end = (winWidth - 1) + WIN_SIZE - 1;

		// calculate global similarity
		globalSimilarity();

		List<Callable<Object>> locSimTasks = new ArrayList<Callable<Object>>();

		for (int i = 0; i < turns; i++) {
			locSimTasks.add(Executors.callable(new LocalSimilarityWrapper(
					WIN_SIZE, start, end, coeffsIdx)));
			coeffsIdx += WIN_SIZE;
			start = end - (winWidth - 1) + 1;
			end = start + (winWidth - 1) + WIN_SIZE - 1;
		}
		if (tail != 0) {
			locSimTasks.add(Executors.callable(new LocalSimilarityWrapper(tail,
					start, alignment.numberOfColumns() - 1, coeffsIdx)));
		}
		executor.invokeAll(locSimTasks);
		return coeffs;
	}

	private void calcLocalSimilarity(int chunkSize, int start, int end,
			int coeffsIdx) {

		int[][] result = localSimilarity2(chunkSize, start, end);

		for (int a = 0; a < result.length; a++) {
			coeffs[coeffsIdx] = pearson(result[a], global);
			coeffsIdx++;
		}
	}

	private class LocalSimilarityWrapper implements Runnable {

		final int chunkSize;
		final int start;
		final int end;
		final int coeffsIdx;

		public LocalSimilarityWrapper(int chunkSize, int start, int end,
				int coeffsIdx) {
			this.chunkSize = chunkSize;
			this.start = start;
			this.end = end;
			this.coeffsIdx = coeffsIdx;
		}

		@Override
		public void run() {
			Correlation.this.calcLocalSimilarity(chunkSize, start, end,
					coeffsIdx);
		}
	}

	/**
	 * Returns correlation score.
	 * 
	 * @param score
	 * @param normalize
	 * @return
	 */
	double[] getCorrelationScore(SMERFSColumnScore score, boolean normalize) {
		if (alignment == null) {
			throw new IllegalArgumentException("Alignment must not be null");
		}
		double[] results = null;

		try {
			results = calcPearson();
		} catch (InterruptedException e) {
			e.printStackTrace();
			throw new RuntimeException("Calculation was interrupted!");
		} catch (ExecutionException e) {
			throw new RuntimeException("Error during calculation!"
					+ e.getLocalizedMessage(), e.getCause());
		}

		double[] columnResults;
		if (score == SMERFSColumnScore.MAX_SCORE) {
			columnResults = giveMaxToColumn(results);
		} else {
			columnResults = giveMidToColumn(results);
		}
		rejectOverTreshold(columnResults);
		if (normalize) {
			double[] normalized = ConservationAccessory.normalize01(
					columnResults, Method.SMERFS);
			return normalized;
		} else {
			return columnResults;
		}
	}

	/**
	 * Finds max within a part of an array, both begin and end delimeters are
	 * also considered.
	 * 
	 * @param scores
	 *            array
	 * @param begin
	 * @param end
	 * @return max
	 */
	private double findMax(double[] scores, int begin, int end) {
		if (end < begin) {
			throw new IllegalArgumentException("End is smaller than the begin.");
		}
		if (begin == end) {
			return scores[begin];
		}
		double max = scores[begin];
		for (int i = begin; i < end + 1; i++) {
			if (scores[i] > max) {
				max = scores[i];
			}
		}
		return max;
	}

	/**
	 * Gives a score to the column. The score given is the max score of all the
	 * windows it belongs to.
	 * 
	 * @param windowScores
	 * @return
	 */
	private double[] giveMaxToColumn(double[] windowScores) {
		double[] scores = new double[alignment.numberOfColumns()];
		for (int i = 0; i < winWidth - 1; i++) {
			scores[i] = findMax(windowScores, 0, i);
		}
		for (int i = winWidth - 1; i < scores.length - (winWidth - 1); i++) {
			scores[i] = findMax(windowScores, i - (winWidth - 1), i);
		}
		int begin = windowScores.length - 1 - (winWidth - 2);
		int end = windowScores.length - 1;
		for (int i = scores.length - (winWidth - 1); i < scores.length; i++) {
			scores[i] = findMax(windowScores, begin, end);
			begin++;
		}
		return scores;
	}

	/**
	 * Gives scores to columns. The middle column gets the window score
	 * 
	 * @param windowScores
	 * @return
	 */
	private double[] giveMidToColumn(double[] windowScores) {
		double[] columnResults = new double[alignment.numberOfColumns()];
		int ends = (winWidth - 1) / 2;
		for (int i = 0; i < ends; i++) {
			columnResults[i] = windowScores[0];
			columnResults[(columnResults.length - 1) - i] = windowScores[windowScores.length - 1];
		}
		for (int i = 0; i < windowScores.length; i++) {
			columnResults[i + ends] = windowScores[i];
		}
		return columnResults;
	}

	private void rejectOverTreshold(double[] results) {
		for (int i = 0; i < alignment.numberOfColumns(); i++) {
			Map<Character, Integer> colMap = alignment.getTotalAcidsFreqByCol()
					.get(i);
			if (colMap.containsKey('-')) {
				if (colMap.get('-') / numofSequences > gapTreshold) {
					results[i] = 0.0;
				}
			}
		}
	}
}
