package compbio.alt;

import java.util.EnumMap;

/**
 * Example of a conservation counter
 * 
 * @author Petr Troshin
 * 
 */
public class ConservationCounter {
	/**
	 * AA alphabet, non AA residue must be represented by u.
	 * 
	 * @author Petr Troshin
	 * 
	 */
	static enum ALPHABET {
		R, H, K, D, E, S, T, N, Q, C, G, P, A, I, L, M, F, W, Y, V, u;

		/**
		 * Counts occurrence of the characters in the char array
		 * 
		 * @param column
		 *            a list of characters. Characters must all be in upper case
		 *            to match the alphabet, whitespace is not allowed
		 * @return a EnumMap with ALPHABET chars as keys and the frequency of
		 *         occurence are the values
		 * @throws IllegalArgumentException
		 *             if character is not defined in the alphabet
		 * @throws NullPointerException
		 *             is character is null
		 */
		static EnumMap<ALPHABET, Integer> occurenceCounter(final char[] column) {
			final EnumMap<ALPHABET, Integer> occ = new EnumMap<ALPHABET, Integer>(
					ALPHABET.class);
			for (final char ch : column) {
				final ALPHABET key = ALPHABET.valueOf(Character.toString(ch));
				Integer count = occ.get(key);
				if (count == null) {
					count = new Integer(0);
				}
				occ.put(key, ++count);
			}
			return occ;
		}
	}

	/**
	 * Usage example
	 * 
	 * @param args
	 *            not used
	 */
	public static void main(final String[] args) {
		final StringBuilder build = new StringBuilder();

		for (int i = 0; i < 10000; i++) {
			for (final ALPHABET A : ALPHABET.values()) {
				build.append(A.toString());
			}
		}
		final char[] chars = build.toString().toCharArray();
		System.out.println(chars.length);
		assert chars.length > 10000;
		final long st = System.nanoTime();
		EnumMap<ALPHABET, Integer> map = null;
		for (int i = 0; i < 100; i++) {
			map = ALPHABET.occurenceCounter(chars);
		}
		System.out.println(map);
		System.out.println((System.nanoTime() - st) / 1000000000 + "s");
	}
}
