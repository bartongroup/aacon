package compbio.conservation;

import java.util.EnumMap;
import java.util.Map;

public class ConservationScores2 {
	
	private final Map<Method, double[]> scores = new EnumMap<Method, double[]>(Method.class);
	
	private final AminoAcidMatrix alignment;
	
	ConservationScores2 (AminoAcidMatrix alignment) {
		
		this.alignment = alignment;
		
	}

double[] calculateScore(Method method, boolean normalize) {
		
		if (method == Method.kabatScore) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.kabatScore(alignment, i);
				
			}
			
			scores.put(Method.kabatScore, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.inversedNormalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
		
		}
	
		if (method == Method.joresScore) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.joresScore(alignment, i);
				
			}
			
			scores.put(Method.joresScore, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.inversedNormalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
		
		}
		
		if (method == Method.schneiderScore) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.schneiderScore(alignment, i);
				
			}
			
			scores.put(Method.schneiderScore, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.inversedNormalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.shenkinScore) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.shenkinScore(alignment, i);
				
			}
			
			scores.put(Method.shenkinScore, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.inversedNormalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.gersteinScore) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.gersteinScore(alignment, i);
				
			}
			
			scores.put(Method.gersteinScore, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.inversedNormalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.SmallestTaylorSetGaps) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.SmallestTaylorSetGaps(alignment, i);
				
			}
			
			scores.put(Method.SmallestTaylorSetGaps, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.inversedNormalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.SmallestTaylorSetNoGaps) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.SmallestTaylorSetNoGaps(alignment, i);
				
			}
			
			scores.put(Method.SmallestTaylorSetNoGaps, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.inversedNormalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.zvelibilScore) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.zvelibilScore(alignment, i);
				
			}
			
			scores.put(Method.zvelibilScore, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.normalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.zvelibilScore) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.zvelibilScore(alignment, i);
				
			}
			
			scores.put(Method.zvelibilScore, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.inversedNormalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.karlinScore) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.karlinScore(alignment, i);
				
			}
			
			scores.put(Method.karlinScore, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.normalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.armonScore) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.armonScore(alignment, i);
				
			}
			
			scores.put(Method.armonScore, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.inversedNormalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.thompsonScore) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.thompsonScore(alignment, i);
				
			}
			
			scores.put(Method.thompsonScore, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.inversedNormalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.notLancetScore) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.notLancetScore(alignment, i);
				
			}
			
			scores.put(Method.notLancetScore, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.inversedNormalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.mirnyScore) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.mirnyScore(alignment, i);
				
			}
			
			scores.put(Method.mirnyScore, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.normalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.mirnyScore) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.mirnyScore(alignment, i);
				
			}
			
			scores.put(Method.mirnyScore, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.normalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.williamsonScore) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.williamsonScore(alignment, i);
				
			}
			
			scores.put(Method.williamsonScore, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.normalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.landgrafScore) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.landgrafScore(alignment, i);
				
			}
			
			scores.put(Method.landgrafScore, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.inversedNormalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.sanderScore) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.sanderScore(alignment, i);
				
			}
			
			scores.put(Method.sanderScore, result);
			
			if (normalize == true) {
				
				double[] normalized = ConservationAccessory.normalize01(result);
				
				return normalized;
			}
			
			else {
				
				return result;
			}
			
		}
		
		if (method == Method.valdarScore) {
			
			double[] result = new double[alignment.numberOfColumns()];
			
			for (int i = 0; i < alignment.numberOfColumns(); i++) {
				
				result[i] = ColumnScores.valdarScore(alignment, i);
				
			}
			
			scores.put(Method.valdarScore, result);
			
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

	}
