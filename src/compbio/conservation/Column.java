package compbio.conservation;
import java.util.*;
import org.testng.Assert;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;




public class Column {

	private final char[] columnArr;
    private final Map<Character,Integer> acidsIntMap;

    // TO BE DELETED, THIS CONSTRUCTOR USED FOR TESTS ONLY

    //public Column(char a, char b, char c, char d, char e, char f,
	 //  char g, char h, char i, char j) {

	//columnArr = new char[10];

	//columnArr[0] = a;
	//columnArr[1] = b;
	//columnArr[2] = c;
	//columnArr[3] = d;
	//columnArr[4] = e;
	//columnArr[5] = f;
	//columnArr[6] = g;
	//columnArr[7] = h;
	//columnArr[8] = i;
	//columnArr[9] = j;


    //}

    // constructor, 

    public Column(AminoAcidMatrix m, int column) {

	if (m == null) {
	    throw new IllegalArgumentException("Matrix must not be null");
	}
	assert column >= 0 && column <= m.numberOfColumns();

	int colLength = m.numberOfRows();

	StringBuilder s = new StringBuilder();

	// creates a string 

	for (int i = 0; i < colLength; i++) {

	    char aacid = m.getMatrixPosition(i, column);
	    s = s.append(aacid);
	}

	String columnStr = s.toString().toUpperCase();

	columnArr = columnStr.toCharArray();

// calls mathod to create a map
	
	Alphabet alp = new Alphabet();
	
	acidsIntMap = alp.calculateOccurance(columnArr);
	
	}
	
    // check whether column is empty 
    
    boolean isEmpty() {
    	
    	if (acidsIntMap.containsKey('-') && acidsIntMap.get('-') == columnArr.length) {
    			
    			return true;
    		}
    	
    	else {
    		
    		return false;
    		
    	}
    	
    }
	// returns true if all but one residues in the column are gaps
    boolean allButOneGaps() {
    	
    	if (acidsIntMap.containsKey('-') && acidsIntMap.get('-') == columnArr.length - 1) {
			
			return true;
		}
	
	else {
		
		return false;
		
	}
	
    }
    // returns true if there is only one residue type in the column
    // this residue is not a gap
    boolean oneResidueType() {
    	
    	if (acidsIntMap.size() == 1 && acidsIntMap.containsKey('-') == false) {
    		
    		return true;
    		
    	}
    	
    	else {
    		
    		return false;
    		
    	}
    	
    }
    
    // returns true if column contains gaps 
    
    boolean containsGaps() {
    	
    	if (acidsIntMap.containsKey('-')) {
    		
    		return true;
    	}
    	
    	else {
    		
    		return false;
    		
    	}
    	
    }
 
    // returns the nr of aa types in the column
    // counts gap as 21 aminoacid
    int numberOfAcidsWithGap() {
    	
    	return acidsIntMap.size(); 
    	
    }
    
    //don't want gap to count gaps as 21 amino acids
    //does not coutn gap
    
    int numberOfAcidsNoGap() {
    
    if(this.containsGaps() == true) {
    	
    	return acidsIntMap.size() - 1;
    }
    
    else {
    	
    	return acidsIntMap.size();
    
    }
    
    }
    // look out fo this mathod if column empty return most common nr as 0
    
    int mostCommonNumber() {
    	
    	int max = 0;
    	
    	Set<Character> keys = acidsIntMap.keySet();
    	
    	Iterator<Character> itr = keys.iterator();
    	
    	while(itr.hasNext()) {
    		
    	Character key = itr.next();
    		
    		if (key != '-' && acidsIntMap.get(key) > max) {
    			
    			max = acidsIntMap.get(key);
    			
    		}
    	
    	}
    	
    	assert max != 0 : "Zero in the most Common Number";
    	
    	return max;
    	
    }
    
    Map<Character,Integer> getAcidsIntMap() {
    
    return acidsIntMap;
    
    
    }
    			
 int length() {
	 
	 return columnArr.length;
	 
 }
 
 // this method gets the number of elements in the smallest set  
 // among the Taylor sets that covers all the aa in the column
 
 int SmallestTaylorSetGaps(Map<String, HashSet<Character>> setMap) {
 
	if (setMap == null) {
		
		throw new IllegalArgumentException("setMap must not be null");
	}
	
	Map<String,Integer> repSets = new HashMap<String,Integer>();
	
	Set<String> setMapKeys = setMap.keySet();
	
	Iterator<String> itr = setMapKeys.iterator();
	
	while(itr.hasNext()) {
		
		Object key = itr.next();
		
		if (setMap.get(key).contains(acidsIntMap.keySet())) {
			
			repSets.put((String) key, new Integer(setMap.get(key).size()));
				
		}
	}
		
	int smallestSetSize = Collections.min(repSets.values());
	
	assert smallestSetSize > 0 && smallestSetSize < 21;
	
	return smallestSetSize;
	
	}
 
 int SmallestTaylorSetNoGaps(Map<String, HashSet<Character>> setMap) {
	 
		if (setMap == null) {
			
			throw new IllegalArgumentException("setMap must not be null");
		}
		Alphabet alp = new Alphabet();
		
		Map<Character, Integer> acidsMapNoGaps = alp.calculateOccuranceNoGaps(columnArr);
		
		Map<String,Integer> repSets = new HashMap<String,Integer>();
		
		Set<String> setMapKeys = setMap.keySet();
		
		Iterator<String> itr = setMapKeys.iterator();
		
		while(itr.hasNext()) {
			
			Object key = itr.next();
			
			if (setMap.get(key).contains(acidsMapNoGaps.keySet())) {
				
				repSets.put((String) key, new Integer(setMap.get(key).size()));
					
			}
		}
			
		int smallestSetSize = Collections.min(repSets.values());
		
		assert smallestSetSize > 0 && smallestSetSize < 21;
		
		return smallestSetSize;
		
		}
	
 	int zvelibilScore(Map<String, HashSet<Character>> setMap){
 		
 		int result = 0;
 		
 		if (setMap == null) {
 			
 			throw new IllegalArgumentException("setMap must not be null");
 		}
 		
 		Set<String> keys = setMap.keySet();
 		
 		Iterator<String> itr = keys.iterator();
 		
 		while(itr.hasNext()) {
 			
 			if(setMap.get(itr.next()).contains(acidsIntMap.keySet())) {
 				
 				result++;
 			}
 			
 		}
 		
 		assert result > 0 && result < 11;
 		
 		return result;
 	}




}
 
  
    	
    	
 

