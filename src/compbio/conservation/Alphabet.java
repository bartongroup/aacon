package compbio.conservation;

import java.util.*;

/** 
 * This class provides an amino acid alphabet
*  There is no special character for unknown amino acid, it is treated as a gap
*
*   @Author agolicz
*/

	class Alphabet {
	
	/** 
	 * A set containing the amino acid alphabet
	 */
	
	private static final Set<Character> ALPHABET;
	
	/** 
	 * An array containing 20 basic amino acids, no gap character
	 */
	
	private static final char[] alpArray = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'};

	
	static {
		
		Set<Character> alph = new HashSet<Character>();
			
			alph.add('R');
			alph.add('H');
			alph.add('K');
			alph.add('D');
			alph.add('E');
			alph.add('S');
			alph.add('T');
			alph.add('Q');
			alph.add('C');
			alph.add('G');
			alph.add('P');
			alph.add('A');
			alph.add('I');
			alph.add('L');
			alph.add('M');
			alph.add('F');
			alph.add('W');
			alph.add('Y');
			alph.add('V');
			alph.add('N');
			alph.add('-');
			
			ALPHABET = Collections.unmodifiableSet(alph);
			
			}
	/** 
	 * Returns the array of 20 basic amino acids
	 * 
	 * @return alphabet array
	 */
	static char[] alphabetArray() {

		return alpArray;

		}
	/** 
	 * Returns a set of amino acids in the alphabet.
	 * Gap represented as a dash
	 * 
	 * @return the alphabet set
	 */

	static Set<Character> alphabet() {
		
		return ALPHABET;
	
		}
	
	/** 
	 * Calculates how many times a particular amino acid is present in a given array of characters representing amino acids.
	 * 
	 * @param column
	 * @return A Map, amino acid characters are keys, the frequencies are values
	 */
	
	static Map<Character, Integer> calculateOccurance( final char[] column) {
    	
    	if (column == null) {
    		
    		throw new IllegalArgumentException("Column must not be  null");
    	}
    		
        Set<Character> alph = ALPHABET;
    	
    	Map<Character,Integer> charCount = new HashMap<Character,Integer>();
    	
    	for (char ch : column) {
    		
    	if(ch == '.' || ch == '*' || ch == ' ' || ch =='X') {
    		
    		ch = '-';
    	}
    	
    	assert alph.contains(ch) : "Illegal character in the column";
    	
        Integer count = charCount.get(ch);
    	
    	if (count == null) {
    	    charCount.put(ch, new Integer(1));
    	    
    	} else {
    	    charCount.put(ch, count + 1);
    	}
    	
    			
    	}
    	
    	//addToOccurance(ch);
    	    
        assert !charCount.isEmpty();
        
    	return charCount;
    	
        }
	
	
}
	
	


