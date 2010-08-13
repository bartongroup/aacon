package compbio.conservation;

import java.util.EnumMap;
import java.util.Map;

public class ConservationScores2 {
	
	//private final Map<Method, double[]> scores = new EnumMap<Method, double[]>(Method.class);
	
	private final AminoAcidMatrix alignment;
	
	ConservationScores2 (AminoAcidMatrix alignment) {
		
		this.alignment = alignment;
		
	}

	/**
	 * 
	 * @param method
	 * @param normalize
	 * @return scoer for teh given method
	 */
double[] calculateScore(Method method, boolean normalize) {
		
		if (method == Method.KABAT) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.kabatScore(alignment, i);
				
			}
			
			//scores.put(Method.KABAT_SCORE, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.inversedNormalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
		
		}
	
		if (method == Method.JORES) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.joresScore(alignment, i);
				
			}
			
			//scores.put(Method.JORES_SCORE, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.inversedNormalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
		
		}
		
		if (method == Method.SCHNEIDER) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.schneiderScore(alignment, i);
				
			}
			
			//scores.put(Method.SCHNEIDER_SCORE, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.inversedNormalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.SHENKIN) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.shenkinScore(alignment, i);
				
			}
			
			//scores.put(Method.SHENKIN_SCORE, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.inversedNormalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.GERSTEIN) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.gersteinScore(alignment, i);
				
			}
			
			//scores.put(Method.GERSTEIN_SCORE, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.inversedNormalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.TAYLOR_GAPS) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.taylorScoreGaps(alignment, i);
				
			}
			
			//scores.put(Method.TAYLOR_SCORE_GAPS, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.inversedNormalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.TAYLOR_NO_GAPS) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.taylorScoreNoGaps(alignment, i);
				
			}
			
			//scores.put(Method.TAYLOR_SCORE_NO_GAPS, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.inversedNormalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.ZVELIBIL) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.zvelibilScore(alignment, i);
				
			}
			
			//scores.put(Method.ZVELIBIL_SCORE, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.normalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.ZVELIBIL) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.zvelibilScore(alignment, i);
				
			}
			
			//scores.put(Method.ZVELIBIL_SCORE, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.inversedNormalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.KARLIN) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.karlinScore(alignment, i);
				
			}
			
			//scores.put(Method.KARLIN_SCORE, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.normalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.ARMON) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.armonScore(alignment, i);
				
			}
			
			//scores.put(Method.ARMON_SCORE, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.inversedNormalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.THOMPSON) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.thompsonScore(alignment, i);
				
			}
			
			//scores.put(Method.THOMPSON_SCORE, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.inversedNormalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.NOT_LANCET) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.notLancetScore(alignment, i);
				
			}
			
			//scores.put(Method.NOT_LANCET_SCORE, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.normalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.MIRNY) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.mirnyScore(alignment, i);
				
			}
			
			//scores.put(Method.MIRNY_SCORE, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.normalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.WILLIAMSON) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.williamsonScore(alignment, i);
				
			}
			
			//scores.put(Method.WILLIAMSON_SCORE, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.normalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.LANDGRAF) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.landgrafScore(alignment, i);
				
			}
			
			//scores.put(Method.LANDGRAF_SCORE, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.inversedNormalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.SANDER) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.sanderScore(alignment, i);
				
			}
			
			//scores.put(Method.SANDER_SCORE, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.normalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.VALDAR) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.valdarScore(alignment, i);
				
			}
			
			//scores.put(Method.VALDAR_SCORE, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.normalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		System.out.println("You should never ever get here");
		
		return null;
		
		}

	//Map<Method, double[]> getScores() {
		
		//return this.scores;
	//}

	}
