package compbio.util;

import java.text.SimpleDateFormat;
import java.util.HashMap;
import java.util.Random;

public class Util {

	public static final boolean isEmpty(String string) {
		return string == null || string.trim().length() == 0;
	}

	public static final SimpleDateFormat datef = new SimpleDateFormat(
			"HH:mm-dd/MM/yyyy");

	public static <K, V> HashMap<K, V> getNewHashMap() {
		return new HashMap<K, V>();
	}

	final static Random rand = new Random();

	/*
	 * Returns random integers with value in range from min to max
	 * 
	 * @param min
	 * 
	 * @param max
	 * 
	 * @return
	 */
	public static int getRandomNumber(int start, int end) {
		if (start >= end) {
			throw new IllegalArgumentException("Start cannot exceed End.");
		}
		// get the range, casting to long to avoid overflow problems
		long range = (long) end - (long) start + 1;
		// compute a fraction of the range, 0 <= frac < range
		long fraction = (long) (range * rand.nextDouble());
		return (int) (fraction + start);
	}

	public static double getRandomNumber(double min, double max) {
		return (max - min) * rand.nextDouble() + min;
	}
}
