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

public class ConservationScore2Tester {

	char a = 'D';
	char b = 'F';
	char c = 'P';

	char d = 'D';
	char e = 'K';
	char f = 'P';

	char g = 'D';
	char h = '-';
	char i = 'K';

	static ExecutorFactory efactory;
	AminoAcidMatrix alignment;
	Conservation scores;

	@BeforeClass
	public void init() {
		ExecutorFactory.initExecutor();
		alignment = new AminoAcidMatrix(a, b, c, d, e, f, g, h, i);

		scores = new Conservation(alignment, true,
				ExecutorFactory.getExecutor());
	}

	static void printScores(double[] results, String tag) {

		System.out.println(tag);

		for (int i = 0; i < results.length; i++) {

			System.out.println(results[i]);
		}
	}

	@Test
	public void kabatTester() {

		double[] result = scores.calculateScore(ConservationMethod.KABAT);

		printScores(result, "Kabat");

	}

	@Test
	public void totalAcidsWillSetsTester() {

		Map<String, Integer> sets = alignment.totalAcidsWillSets();

		Iterator<String> itr = sets.keySet().iterator();

		while (itr.hasNext()) {

			String key = itr.next();

			// System.out.println("Key" + key + "Value" + sets.get(key));
		}

	}

	@Test
	public void joresTester() {

		double[] result = scores.calculateScore(ConservationMethod.JORES);

		printScores(result, "Jores");

	}

	@Test
	public void schneiderTester() {

		double[] result = scores.calculateScore(ConservationMethod.SCHNEIDER);

		printScores(result, "Schneider");

	}

	@Test
	public void shenkinTester() {

		double[] result = scores.calculateScore(ConservationMethod.SHENKIN);

		printScores(result, "Shenkin");

	}

	@Test
	public void gersteinTester() {

		double[] result = scores.calculateScore(ConservationMethod.GERSTEIN);

		printScores(result, "Gerstein");

	}

	@Test
	public void smallestTaylosetGapsTester() {

		double[] result = scores.calculateScore(ConservationMethod.TAYLOR_GAPS);

		printScores(result, "SmallestTaylorSetGaps");

	}

	@Test
	public void smallestTaylorSetNoGapsTester() {

		double[] result = scores
				.calculateScore(ConservationMethod.TAYLOR_NO_GAPS);

		printScores(result, "SmallestTaylorSetNoGaps");

	}

	@Test
	public void zvelibilTester() {

		double[] result = scores.calculateScore(ConservationMethod.ZVELIBIL);

		printScores(result, "Zvelibil");

	}

	@Test
	public void karlinTester() {

		double[] result = scores.calculateScore(ConservationMethod.KARLIN);

		printScores(result, "Karlin");

	}

	@Test
	public void armonTester() {

		double[] result = scores.calculateScore(ConservationMethod.ARMON);

		printScores(result, "Armon");

	}

	@Test
	public void thompsonTester() {

		double[] result = scores.calculateScore(ConservationMethod.THOMPSON);

		printScores(result, "Thompson");

	}

	@Test
	public void notLancetTester() {

		double[] result = scores.calculateScore(ConservationMethod.NOT_LANCET);

		printScores(result, "NotLancet");

	}

	@Test
	public void mirnyTester() {

		double[] result = scores.calculateScore(ConservationMethod.MIRNY);

		printScores(result, "Mirny");

	}

	@Test
	public void williamsonTester() {

		double[] result = scores.calculateScore(ConservationMethod.WILLIAMSON);

		printScores(result, "Williamson");

	}

	@Test
	public void landgrafTester() {

		double[] result = scores.calculateScore(ConservationMethod.LANDGRAF);

		printScores(result, "Landgraf");

	}

	@Test
	public void sanderTester() {

		double[] result = scores.calculateScore(ConservationMethod.SANDER);

		printScores(result, "Sander");

	}

	@Test
	public void valdarTester() {

		double[] result = scores.calculateScore(ConservationMethod.VALDAR);

		printScores(result, "Valdar");

	}

	@Test
	public void tester() {

		List<Map<Character, Integer>> map = alignment.getTotalAcidsFreqByCol();

		System.out.println("freqs");

		for (int i = 0; i < map.size(); i++) {

			Iterator<Character> itr = map.get(i).keySet().iterator();

			while (itr.hasNext()) {

				Character key = itr.next();

				System.out.println(i + " key " + key + " value "
						+ map.get(i).get(key));
			}

		}
	}
}
