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

/**
 * Enumeration defining two constraints for SMERFS columns score calculation.
 * MAX_SCORE gives the highest core of all the windows the column belongs to.
 * MID_SCORE gives the window score to the column in the middle.
 * 
 * @author Agnieszka Golicz & Peter Troshin
 */
public enum SMERFSColumnScore {
	MAX_SCORE, MID_SCORE;

	/**
	 * Default window size value for SMERFS algorithm
	 */
	public static final int DEFAULT_WINDOW_SIZE = 7;

	/**
	 * Default gap threshold value for SMERFS algorithm
	 */
	public static final double DEFAULT_GAP_THRESHOLD = 0.1;

	static SMERFSColumnScore getSMERFSColumnScore(String score) {

		score = score.trim().toLowerCase();
		if (score.equalsIgnoreCase(SMERFSColumnScore.MAX_SCORE.toString())) {
			return SMERFSColumnScore.MAX_SCORE;
		}
		if (score.equalsIgnoreCase(SMERFSColumnScore.MID_SCORE.toString())) {
			return SMERFSColumnScore.MID_SCORE;
		}
		return null;
	}

}
