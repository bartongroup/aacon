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

import org.testng.annotations.Test;

import compbio.data.sequence.ConservationMethod;

public class AccessoryTester {

	// @Test
	// public void VoronoiTester() {

	// char a = 'D';
	// char b = 'D';
	// char c = 'D';
	// char d = 'D';
	// char e = 'D';
	// char f = 'E';
	// /char g = 'E';
	// char h = 'E';
	// char i = 'F';
	// char j = 'F';

	// AminoAcidMatrix matrix = new AminoAcidMatrix(a, b, c, d, e, f, g, h, i,
	// j);

	// double result = ConservationAccessory.voronoiWeights(matrix, 1000)[0];

	// Assert.assertEquals(result, 1.0, 0.1);
	// }

	@Test
	public void normalizationTester() {

		double[] scores = { 15, 13, 13, -32, 0 };

		double[] normalized = ConservationAccessory.normalize01(scores,
				ConservationMethod.KABAT);

		System.out.println("normal");

		for (int i = 0; i < normalized.length; i++) {

			System.out.println(normalized[i]);
		}
	}

	@Test
	public void inversedNormalizationTester() {

		double[] scores = { 15, 13, 13, -32, 0 };

		double[] normalized = ConservationAccessory.inversedNormalize01(scores,
				ConservationMethod.KABAT);

		System.out.println("inversed");

		for (int i = 0; i < normalized.length; i++) {

			System.out.println(normalized[i]);
		}
	}

}
