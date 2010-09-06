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
 * 
 * @author Agnieszka Golicz
 * 
 */
public enum Format {
	RESULT_WITH_ALIGNMENT, RESULT_NO_ALIGNMENT;

	static Format getFormat(String format) {

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
