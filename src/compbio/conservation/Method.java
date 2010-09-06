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
 * Enumeration listing all the supported methods.
 * 
 * @author Agnieszka Golicz
 */
public enum Method {
	KABAT, JORES, SCHNEIDER, SHENKIN, GERSTEIN, TAYLOR_GAPS, TAYLOR_NO_GAPS, ZVELIBIL, KARLIN, ARMON, THOMPSON, NOT_LANCET, MIRNY, WILLIAMSON, LANDGRAF, SANDER, VALDAR, SMERFS;

	static Method getMethod(String meth) {
		meth = meth.trim().toLowerCase();
		if (meth.equalsIgnoreCase(KABAT.toString())) {
			return KABAT;
		}
		if (meth.equalsIgnoreCase(JORES.toString())) {
			return JORES;
		}
		if (meth.equalsIgnoreCase(SCHNEIDER.toString())) {
			return SCHNEIDER;
		}
		if (meth.equalsIgnoreCase(SHENKIN.toString())) {
			return SHENKIN;
		}
		if (meth.equalsIgnoreCase(GERSTEIN.toString())) {
			return GERSTEIN;
		}
		if (meth.equalsIgnoreCase(TAYLOR_GAPS.toString())) {
			return TAYLOR_GAPS;
		}
		if (meth.equalsIgnoreCase(TAYLOR_NO_GAPS.toString())) {
			return TAYLOR_NO_GAPS;
		}
		if (meth.equalsIgnoreCase(ZVELIBIL.toString())) {
			return ZVELIBIL;
		}
		if (meth.equalsIgnoreCase(KARLIN.toString())) {
			return KARLIN;
		}
		if (meth.equalsIgnoreCase(ARMON.toString())) {
			return ARMON;
		}
		if (meth.equalsIgnoreCase(THOMPSON.toString())) {
			return THOMPSON;
		}
		if (meth.equalsIgnoreCase(NOT_LANCET.toString())) {
			return NOT_LANCET;
		}
		if (meth.equalsIgnoreCase(MIRNY.toString())) {
			return MIRNY;
		}
		if (meth.equalsIgnoreCase(WILLIAMSON.toString())) {
			return WILLIAMSON;
		}
		if (meth.equalsIgnoreCase(LANDGRAF.toString())) {
			return LANDGRAF;
		}
		if (meth.equalsIgnoreCase(SANDER.toString())) {
			return SANDER;
		}
		if (meth.equalsIgnoreCase(VALDAR.toString())) {
			return VALDAR;
		}
		if (meth.equalsIgnoreCase(SMERFS.toString())) {
			return SMERFS;
		}
		return null;
	}

}
