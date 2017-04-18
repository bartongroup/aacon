/*
 * Copyright (c) 2010 Agnieszka Golicz & Peter Troshin 
 * 
 * Amino Acid Conservation @version: 1.0 
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the Apache License version 2 as published by the
 * Apache Software Foundation This library is distributed in the hope that it
 * will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Apache
 * License for more details. A copy of the license is in apache_license.txt. It
 * is also available here: http://www.apache.org/licenses/LICENSE-2.0.txt 
 * Any republication or derived work distributed in source code form must 
 * include this copyright and license notice.
 * 
 */
package compbio.conservation;

import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import compbio.data.sequence.ConservationMethod;

public class ConservationScore3Tester {
		
	// row1
	char p1 = 'D';
	char p2 = 'A';
	char p3 = 'S';
	char p4 = 'P';
	char p5 = 'F';
	
	// row2
	char p6 = 'D';
	char p7 = 'T';
	char p8 = 'S';
	char p9 = 'K';
	char p10 = '-';
	
	// row3
	char p11 = 'D';
	char p12 = 'P';
	char p13 = 'T';
	char p14 = 'Y';
	char p15 = '-';
	
	static ExecutorFactory efactory;
	AminoAcidMatrix alignment;
	Conservation scores;
	Conservation scores_norm;

	@BeforeClass
	public void init() {
		ExecutorFactory.initExecutor();
		alignment = new AminoAcidMatrix(p1, p2, p3, p4, p5, p6, p7, p8, p9, 
										p10, p11, p12, p13, p14, p15);
		
		// try on normalized = true and default = false
		scores = new Conservation(alignment, false,
				ExecutorFactory.getExecutor());
		scores_norm = new Conservation(alignment, true,
				ExecutorFactory.getExecutor());
	}

	static void printScores(double[] results, double[] results_norm, String tag) {

		System.out.println(tag);

		for (int i = 0; i < results.length; i++) {

			System.out.printf("%.4f, %.4f\n", results[i], results_norm[i]);
		}
	}

	@Test
	public void kabatTester() {

		double[] result = scores.calculateScore(ConservationMethod.KABAT);
		double[] result_norm = scores_norm.calculateScore(ConservationMethod.KABAT);
		
		printScores(result, result_norm, "Kabat");

	}



	@Test
	public void joresTester() {

		double[] result = scores.calculateScore(ConservationMethod.JORES);
		double[] result_norm = scores_norm.calculateScore(ConservationMethod.JORES);
		
		printScores(result, result_norm, "Jores");

	}

	@Test
	public void schneiderTester() {

		double[] result = scores.calculateScore(ConservationMethod.SCHNEIDER);
		double[] result_norm = scores_norm.calculateScore(ConservationMethod.SCHNEIDER);
		
		printScores(result, result_norm, "Schneider");

	}

	@Test
	public void shenkinTester() {

		double[] result = scores.calculateScore(ConservationMethod.SHENKIN);
		double[] result_norm = scores_norm.calculateScore(ConservationMethod.SHENKIN);
		
		printScores(result, result_norm, "Shenkin");

	}

	@Test
	public void gersteinTester() {

		double[] result = scores.calculateScore(ConservationMethod.GERSTEIN);
		double[] result_norm = scores_norm.calculateScore(ConservationMethod.GERSTEIN);
		
		printScores(result, result_norm, "Gerstein");

	}

	@Test
	public void smallestTaylosetGapsTester() {

		double[] result = scores.calculateScore(ConservationMethod.TAYLOR_GAPS);
		double[] result_norm = scores_norm.calculateScore(ConservationMethod.TAYLOR_GAPS);
		
		printScores(result, result_norm, "SmallestTaylorSetGaps");

	}

	@Test
	public void smallestTaylorSetNoGapsTester() {

		double[] result = scores.calculateScore(ConservationMethod.TAYLOR_NO_GAPS);
		double[] result_norm = scores_norm.calculateScore(ConservationMethod.TAYLOR_NO_GAPS);
		
		printScores(result, result_norm, "SmallestTaylorSetNoGaps");

	}

	@Test
	public void zvelibilTester() {

		double[] result = scores.calculateScore(ConservationMethod.ZVELIBIL);
		double[] result_norm = scores_norm.calculateScore(ConservationMethod.ZVELIBIL);
		
		printScores(result, result_norm, "Zvelibil");

	}

	@Test
	public void karlinTester() {

		double[] result = scores.calculateScore(ConservationMethod.KARLIN);
		double[] result_norm = scores_norm.calculateScore(ConservationMethod.KARLIN);
		
		printScores(result, result_norm, "Karlin");

	}

	@Test
	public void armonTester() {

		double[] result = scores.calculateScore(ConservationMethod.ARMON);
		double[] result_norm = scores_norm.calculateScore(ConservationMethod.ARMON);
		
		printScores(result, result_norm, "Armon");

	}

	@Test
	public void thompsonTester() {

		double[] result = scores.calculateScore(ConservationMethod.THOMPSON);
		double[] result_norm = scores_norm.calculateScore(ConservationMethod.THOMPSON);
		
		printScores(result, result_norm, "Thompson");

	}

	@Test
	public void notLancetTester() {

		double[] result = scores.calculateScore(ConservationMethod.NOT_LANCET);
		double[] result_norm = scores_norm.calculateScore(ConservationMethod.NOT_LANCET);
		
		printScores(result, result_norm, "NotLancet");

	}

	@Test
	public void mirnyTester() {

		double[] result = scores.calculateScore(ConservationMethod.MIRNY);
		double[] result_norm = scores_norm.calculateScore(ConservationMethod.MIRNY);
		
		printScores(result, result_norm, "Mirny");

	}

	@Test
	public void williamsonTester() {

		double[] result = scores.calculateScore(ConservationMethod.WILLIAMSON);
		double[] result_norm = scores_norm.calculateScore(ConservationMethod.WILLIAMSON);
		
		printScores(result, result_norm, "Williamson");

	}

	@Test
	public void landgrafTester() {

		double[] result = scores.calculateScore(ConservationMethod.LANDGRAF);
		double[] result_norm = scores_norm.calculateScore(ConservationMethod.LANDGRAF);
		
		printScores(result, result_norm, "Landgraf");

	}

	@Test
	public void sanderTester() {

		double[] result = scores.calculateScore(ConservationMethod.SANDER);
		double[] result_norm = scores_norm.calculateScore(ConservationMethod.SANDER);
		
		printScores(result, result_norm, "Sander");

	}

	@Test
	public void valdarTester() {

		double[] result = scores.calculateScore(ConservationMethod.VALDAR);
		double[] result_norm = scores_norm.calculateScore(ConservationMethod.VALDAR);
		
		printScores(result, result_norm, "Valdar");

	}
	
//	@Test
//	public void totalAcidsWillSetsTester() {
//
//		Map<String, Integer> sets = alignment.totalAcidsWillSets();
//
//		Iterator<String> itr = sets.keySet().iterator();
//
//		while (itr.hasNext()) {
//
//			String key = itr.next();
//
//			// System.out.println("Key" + key + "Value" + sets.get(key));
//		}
//
//	}
//	
//	@Test
//	public void tester() {
//
//		List<Map<Character, Integer>> map = alignment.getTotalAcidsFreqByCol();
//
//		System.out.println("freqs");
//
//		for (int i = 0; i < map.size(); i++) {
//
//			Iterator<Character> itr = map.get(i).keySet().iterator();
//
//			while (itr.hasNext()) {
//
//				Character key = itr.next();
//
//				System.out.println(i + " key " + key + " value "
//						+ map.get(i).get(key));
//			}
//
//		}
//	}
}
