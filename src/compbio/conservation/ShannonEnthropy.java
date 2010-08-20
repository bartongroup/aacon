package compbio.conservation;

import java.util.Collection;
import java.util.Iterator;
import java.util.Map;

class ShannonEnthropy {

	/**
	 * Calculates Shanon enthropy. Uses logarithm with base 2.
	 * 
	 * @param map
	 * @param nrSequences
	 * @return
	 */
	static double ShannonLog2(final Map<Character, Integer> map, double nrSequences) {

		if (map == null) {
			throw new IllegalArgumentException("Map must not be null");
		}
		assert !map.isEmpty() : "Shannon has been fed an empty map";
		assert nrSequences > 0;
		double sum = 0;
		Collection<Character> keys = map.keySet();
		Iterator<Character> itr = keys.iterator();
		while (itr.hasNext()) {
			Character key = itr.next();
			if (key != '-') {
				double value = map.get(key) / nrSequences;
				sum = sum + (value * (Math.log(value) / Math.log(2.0)));
			}
		}
		// assert sum > 0 : "Shannon has been fed an empty column";
		return -sum;
	}

	/**
	 * Calculates Shannon enthropy. Uses natural logarithm
	 * 
	 * @param map
	 * @param nrSequences
	 * @return
	 */
	static double ShannonLn(final Map<Character, Integer> map, int nrSequences) {

		if (map == null) {
			throw new IllegalArgumentException("Map must not be null");
		}
		assert !map.isEmpty() : "Shannon has been fed an empty map";
		assert nrSequences > 0;
		double sum = 0;
		Collection<Character> keys = map.keySet();
		Iterator<Character> itr = keys.iterator();
		while (itr.hasNext()) {
			Character key = itr.next();
			if (key != '-') {
				double value = map.get(key) / (double) nrSequences;
				sum = sum + (value * Math.log(value));
			}
		}
		// assert sum > 0 : "Shannon enthropy has been fed an empty column";
		return -sum;
	}
}
