package compbio.conservation;

public class ConservationScores2 {

	private final AminoAcidMatrix alignment;

	ConservationScores2(AminoAcidMatrix alignment) {
		this.alignment = alignment;
	}

	AminoAcidMatrix getAlignment() {
		return alignment;
	}

	/**
	 * @param method
	 * @param normalize
	 * @return scoer for teh given method
	 */
	double[] calculateScore(final Method method, final boolean normalize) {
		double[] result = new double[alignment.numberOfColumns()];
		switch (method) {
		case KABAT:
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				result[i] = ColumnScores.kabatScore(alignment, i);
			}
			// scores.put(Method.KABAT_SCORE, result);
			if (normalize == true) {
				double[] normalized = ConservationAccessory
						.inversedNormalize01(result, method);
				return normalized;
			}
			return result;
		case JORES:
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				result[i] = ColumnScores.joresScore(alignment, i);
			}
			// scores.put(Method.JORES_SCORE, result);
			if (normalize == true) {
				double[] normalized = ConservationAccessory
						.inversedNormalize01(result, method);
				return normalized;
			}
			return result;
		case SCHNEIDER:
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				result[i] = ColumnScores.schneiderScore(alignment, i);
			}
			// scores.put(Method.SCHNEIDER_SCORE, result);
			if (normalize == true) {
				double[] normalized = ConservationAccessory
						.inversedNormalize01(result, method);
				return normalized;
			}
			return result;
		case SHENKIN:
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				result[i] = ColumnScores.shenkinScore(alignment, i);
			}
			// scores.put(Method.SHENKIN_SCORE, result);
			if (normalize == true) {
				double[] normalized = ConservationAccessory
						.inversedNormalize01(result, method);
				return normalized;
			}
			return result;
		case GERSTEIN:
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				result[i] = ColumnScores.gersteinScore(alignment, i);
			}
			// scores.put(Method.GERSTEIN_SCORE, result);
			if (normalize == true) {
				double[] normalized = ConservationAccessory
						.inversedNormalize01(result, method);
				return normalized;
			}
			return result;
		case TAYLOR_GAPS:
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				result[i] = ColumnScores.taylorScoreGaps(alignment, i);
			}
			// scores.put(Method.TAYLOR_SCORE_GAPS, result);
			if (normalize == true) {
				double[] normalized = ConservationAccessory
						.inversedNormalize01(result, method);
				return normalized;
			}
			return result;
		case TAYLOR_NO_GAPS:
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				result[i] = ColumnScores.taylorScoreNoGaps(alignment, i);
			}
			// scores.put(Method.TAYLOR_SCORE_NO_GAPS, result);
			if (normalize == true) {
				double[] normalized = ConservationAccessory
						.inversedNormalize01(result, method);
				return normalized;
			}
			return result;
		case ZVELIBIL:
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				result[i] = ColumnScores.zvelibilScore(alignment, i);
			}
			// scores.put(Method.ZVELIBIL_SCORE, result);
			if (normalize == true) {
				double[] normalized = ConservationAccessory.normalize01(result,
						method);
				return normalized;
			}
			return result;
		case KARLIN:
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				result[i] = ColumnScores.karlinScore(alignment, i);
			}
			// scores.put(Method.KARLIN_SCORE, result);
			if (normalize == true) {
				double[] normalized = ConservationAccessory.normalize01(result,
						method);
				return normalized;
			}
			return result;
		case ARMON:
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				result[i] = ColumnScores.armonScore(alignment, i);
			}
			// scores.put(Method.ARMON_SCORE, result);
			if (normalize == true) {
				double[] normalized = ConservationAccessory
						.inversedNormalize01(result, method);
				return normalized;
			}
			return result;
		case THOMPSON:
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				result[i] = ColumnScores.thompsonScore(alignment, i);
			}
			// scores.put(Method.THOMPSON_SCORE, result);
			if (normalize == true) {
				double[] normalized = ConservationAccessory
						.inversedNormalize01(result, method);
				return normalized;
			}
			return result;
		case NOT_LANCET:
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				result[i] = ColumnScores.notLancetScore(alignment, i);
			}
			// scores.put(Method.NOT_LANCET_SCORE, result);
			if (normalize == true) {
				double[] normalized = ConservationAccessory.normalize01(result,
						method);
				return normalized;
			}
			return result;
		case MIRNY:
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				result[i] = ColumnScores.mirnyScore(alignment, i);
			}
			// scores.put(Method.MIRNY_SCORE, result);
			if (normalize == true) {
				double[] normalized = ConservationAccessory.normalize01(result,
						method);
				return normalized;
			}
			return result;
		case WILLIAMSON:
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				result[i] = ColumnScores.williamsonScore(alignment, i);
			}
			// scores.put(Method.WILLIAMSON_SCORE, result);
			if (normalize == true) {
				double[] normalized = ConservationAccessory.normalize01(result,
						method);
				return normalized;
			}
			return result;
		case LANDGRAF:
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				result[i] = ColumnScores.landgrafScore(alignment, i);
			}
			// scores.put(Method.LANDGRAF_SCORE, result);
			if (normalize == true) {
				double[] normalized = ConservationAccessory
						.inversedNormalize01(result, method);
				return normalized;
			}
			return result;
		case SANDER:
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				result[i] = ColumnScores.sanderScore(alignment, i);
			}
			// scores.put(Method.SANDER_SCORE, result);
			if (normalize == true) {
				double[] normalized = ConservationAccessory.normalize01(result,
						method);
				return normalized;
			}
			return result;
		case VALDAR:
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				result[i] = ColumnScores.valdarScore(alignment, i);
			}
			// scores.put(Method.VALDAR_SCORE, result);
			if (normalize == true) {
				double[] normalized = ConservationAccessory.normalize01(result,
						method);
				return normalized;
			}
			return result;
		case SMERFS:
			result = ConservationClient.getSMERFS(alignment,
					SMERFSColumnScore.DEFAULT_WINDOW_SIZE,
					SMERFSColumnScore.MID_SCORE,
					SMERFSColumnScore.DEFAULT_GAP_THRESHOLD, normalize);
			return result;
		default:
			throw new RuntimeException("You should never ever get here");
		}
	}
}
