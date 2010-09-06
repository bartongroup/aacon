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
 * @author Peter Troshin
 * 
 */
final class OverflowChecker {

	static final boolean preAdd(double left, double right)
			throws ArithmeticException {
		if (right > 0 ? left > Double.MAX_VALUE - right
				: left < Double.MIN_VALUE - right) {
			return false;
		}
		return true;
	}

	static final boolean preSubtract(double left, double right)
			throws ArithmeticException {
		if (right > 0 ? left < Double.MIN_VALUE + right
				: left > Double.MAX_VALUE + right) {
			return false;
		}
		return true;
	}

	static final boolean preMultiply(double left, double right)
			throws ArithmeticException {
		if (right > 0 ? left > Double.MAX_VALUE / right
				|| left < Double.MIN_VALUE / right
				: (right < -1 ? left > Double.MIN_VALUE / right
						|| left < Double.MAX_VALUE / right : right == -1
						&& left == Double.MIN_VALUE)) {
			return false;
		}
		return true;
	}

	static final boolean preDivide(double left, double right)
			throws ArithmeticException {
		if ((left == Double.MIN_VALUE) && (right == -1)) {
			return false;
		}
		return true;
	}

	static final boolean preAbs(double a) throws ArithmeticException {
		if (a == Double.MIN_VALUE) {
			return false;
		}
		return true;
	}

	static final boolean preAdd(int left, int right) throws ArithmeticException {
		if (right > 0 ? left > Integer.MAX_VALUE - right
				: left < Integer.MIN_VALUE - right) {
			return false;
		}
		return true;
	}

	static final boolean preSubtract(int left, int right)
			throws ArithmeticException {
		if (right > 0 ? left < Integer.MIN_VALUE + right
				: left > Integer.MAX_VALUE + right) {
			return false;
		}
		return true;
	}

	static final boolean preMultiply(int left, int right)
			throws ArithmeticException {
		if (right > 0 ? left > Integer.MAX_VALUE / right
				|| left < Integer.MIN_VALUE / right
				: (right < -1 ? left > Integer.MIN_VALUE / right
						|| left < Integer.MAX_VALUE / right : right == -1
						&& left == Integer.MIN_VALUE)) {
			return false;
		}
		return true;
	}

	static final boolean preDivide(int left, int right)
			throws ArithmeticException {
		if ((left == Integer.MIN_VALUE) && (right == -1)) {
			return false;
		}
		return true;
	}

	static final boolean preAbs(int a) throws ArithmeticException {
		if (a == Integer.MIN_VALUE) {
			return false;
		}
		return true;
	}

}
