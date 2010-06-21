package compbio.conservation;

import java.util.*;

// This class provides an amino acid alphabet
// Alphabet comes in two varieties as a set, containing a gap character and as an array that contain only 20 basic amino acids
// There is no special character for unknown amino acid, it is treated as a gap


class Alphabet {
	
	private static final Set<Character> ALPHABET;
	
	// takes into consideration only basic 20 amino acids and no gaps
	
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
			alph.add('-');
			
			ALPHABET = Collections.unmodifiableSet(alph);
			
			}
	
	static char[] alphabetArray() {

		return alpArray;

		}
	

	static Set<Character> alphabet() {
		
		return ALPHABET;
	
		}
	
	// calculates occurance of aa, aa not defined treats as gap
	
// calculates occurance of aa, aa not defined treats as gap
	
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
	
	


