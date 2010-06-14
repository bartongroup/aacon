package compbio.conservation;

import java.util.*;


class Alphabet {
	
	//private final Set<Character> alph;
	
	///private final Map<Character,Integer> charCount;
	
	//private final Map<Character,Integer> charCountNoGaps;
	
	static Set<Character> alphabet() {
		
		Set<Character> alph = new HashSet<Character>();
		
	//	charCount = new HashMap<Character,Integer>();
		
	//	charCountNoGaps = new HashMap<Character,Integer>();
		
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
		
		return alph;
		
		}
	
   static Map<Character, Integer> calculateOccurance( final char[] column) {
    	
    	if (column == null) {
    		
    		throw new IllegalArgumentException("Column must not be  null");
    	}
    		
        Set<Character> alph = Alphabet.alphabet();
    	
    	Map<Character,Integer> charCount = new HashMap<Character,Integer>();
    	
    	for (char ch : column) {
    		
    	if(ch == '.' || ch == '*' || ch == ' ') {
    		
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

       // void addToOccurance(char ch) {
        	
    	//Integer count = charCount.get(ch);
    	
    	//if (count == null) {
    	   // charCount.put(ch, new Integer(1));
    	    
    	//} else {
    	   // charCount.put(ch, count + 1);
    	//}
    	
        //}
        
  //      Map<Character, Integer> calculateOccuranceNoGaps( final char[] column) {
  //    	
  //  	if (column == null) {
        		
  //      		throw new IllegalArgumentException("Column must not be  null");
  //      	}
  //     		
  //      	for (char ch : column) {
        	
  //      	assert alph.contains(ch) : "Illegal character in the column";
        			
  //      	    addToOccurance(ch);
        	    
  //      		}
        	    
        	
  //          assert !charCountNoGaps.isEmpty();
            
  //      	return charCountNoGaps;
  //          }

  //       void addToOccuranceNoGaps(char ch) {
        	
  //      	Integer count = charCountNoGaps.get(ch);
        	
  //      	if (count == null) {
  //      	    charCountNoGaps.put(ch, new Integer(1));
        	    
  //      	} else {
  //      	    charCountNoGaps.put(ch, count + 1);
  //      	}
        	
  //          }

}

	
	


