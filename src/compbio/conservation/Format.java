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
 * Conservation output Format.
 * 
 * @author Agnieszka Golicz & Peter Troshin
 * 
 */
public enum Format {

	/**
	 * Output just the conservation scores, but no alignment in the following
	 * format: <br/>
	 * #Method name scores <br/>
	 * For example:
	 * 
	 * <pre>{@code
	 * #SANDER 0 0.345 0.543 0.667 1 0.2
	 * #SMERFS 1 1 1 1 1 1 1 1 1 1 
	 * }</pre>
	 */
	RESULT_WITH_ALIGNMENT,

	/**
	 * Output the conservation scores and the alignment
	 * 
	 * The alignment is saved in FASTA format prior to the conservation scores.
	 * The scores themselves are saved in the same format as with
	 * RESULT_WITH_ALIGNMENT
	 * 
	 */
	RESULT_NO_ALIGNMENT;

	/**
	 * Converts from String format name to enum Format
	 * 
	 * @param format
	 *            the format to parse
	 * @return the Format Enum value
	 */
	public static Format getFormat(String format) {

		format = format.trim().toLowerCase();
		if (format.equalsIgnoreCase(Format.RESULT_WITH_ALIGNMENT.toString())) {
			return Format.RESULT_WITH_ALIGNMENT;
		}
		if (format.equalsIgnoreCase(Format.RESULT_NO_ALIGNMENT.toString())) {
			return Format.RESULT_NO_ALIGNMENT;
		}
		return null;
	}

}
