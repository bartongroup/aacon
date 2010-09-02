package compbio.conservation;

public class OverflowChecker {

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
