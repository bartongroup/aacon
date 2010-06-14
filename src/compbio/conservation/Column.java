package compbio.conservation;
import java.util.*;





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
	
	//Alphabet alp = new Alphabet();
	
	acidsIntMap = Alphabet.calculateOccurance(columnArr);
	
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
		
		String key = itr.next();
		
		if (setMap.get(key).containsAll(acidsIntMap.keySet())) {
			
			repSets.put(key, new Integer(setMap.get(key).size()));
				
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
		//Alphabet alp = new Alphabet();
		
		Map<Character,Integer> acidsMapNoGaps = new HashMap<Character,Integer>(acidsIntMap);
		
		if(acidsMapNoGaps.containsKey('-')) {
			
			acidsMapNoGaps.remove('-');
			
		}
		
		Map<String,Integer> repSets = new HashMap<String,Integer>();
		
		Set<String> setMapKeys = setMap.keySet();
		
		Iterator<String> itr = setMapKeys.iterator();
		
		while(itr.hasNext()) {
			
			String key = itr.next();
			
			if (setMap.get(key).containsAll(acidsMapNoGaps.keySet())) {
				
				repSets.put(key, new Integer(setMap.get(key).size()));
					
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
 			
 			if(setMap.get(itr.next()).containsAll(acidsIntMap.keySet())) {
 				
 				result++;
 			}
 			
 		}
 		
 		assert result > 0 && result < 11;
 		
 		return result;
 	}

 	// done exactly like in jon's code: gaps not counted, if all but 1 gaps conservation score 0
 	// although it looks to me that the formula suggests that if one doesn't count gaps
 	// than the score of all but one gaps would be undefined
 	// will read the original paper
 	
 	double karlinScore() {
 		
 		double blosumSum = 0;
 		
 		double finalSum = 0;
 		
 		for( int a = 0; a < columnArr.length; a++) {
 			
 			if(columnArr[a] != '-') {
 				
 				for(int b = a + 1; b < columnArr.length; b++) {
 					
 					if(columnArr[b] != '-') {
 						
 						int pairScore = ConservationAccessory.BlosumPair(columnArr[a], columnArr[b]);
 						
 						int aSelf = ConservationAccessory.BlosumPair(columnArr[a], columnArr[a]);
 						
 						int bSelf = ConservationAccessory.BlosumPair(columnArr[a], columnArr[b]);
 						
 						blosumSum = blosumSum + (pairScore / (Math.sqrt(aSelf * bSelf)));
 					
 					}
 					
 				}
 			}
 		}
 		
 		finalSum = blosumSum * (2/ (columnArr.length * (columnArr.length - 1)));
 		
 		assert finalSum >= -1 && finalSum <= 1;
 		
 		return finalSum;
 	}
 	
 // creates an array containg all the aminoacids and gaps present in the colum
 // iterates through that array twice(nested loops), finds all the possible pairs 
 // that can be formed by aa present
 // gap is considered the 21 aminoacid
 	double armonScore() {
 		
 		double scoreSum = 0;
 		
 		//Object[] acidsPresent = acidsIntMap.keySet().toArray();
 		
 		int arrayLength = acidsIntMap.keySet().size();
 		
 		Character[] acidsPresent = new Character[arrayLength];
 		
 		int arrayIndex = 0;
 		
 		Set<Character> keys = acidsIntMap.keySet();
 		
 		Iterator<Character> itr = keys.iterator();
 		
 		while(itr.hasNext()) {
 			
 			acidsPresent[arrayIndex] = itr.next();
 		}
 		
 		for ( int a = 0; a < acidsPresent.length; a++ ) {
 			
 			char charA = acidsPresent[a];
 			
 			for (int b = a + 1; b < acidsPresent.length; b++) {
 				
 				char charB = acidsPresent[b];
 				
 				scoreSum = scoreSum + ConservationAccessory.BlosumPair(charA, charB);
 				
 			}
 		}
 		
 		return scoreSum;
 		
 	}
// amino acids are viewed as points in k dimensional space
 	
 	double thompsonScore(){
 		
 		int[] sum = null;
 		
 		double[] meanPoint = null;
 		
 		double distance = 0;
 		
 		char[] alp = ConservationAccessory.alphabetArray();
 		
 		int[][] points = new int[columnArr.length][alp.length];
 		
 		for (int i = 0; i < columnArr.length; i++) {
 			
 			for (int j = 0; j < alp.length; j++) {
 				
 				points[i][j] = ConservationAccessory.BlosumPair(columnArr[i], alp[j]);
 				
 			}
 			
 			if (sum == null) {
 				
 				sum = points[i];
 				
 			}
 				
 				else {
 					
 					sum = ConservationAccessory.addPoints(sum, points[i]);
 				}
 			
 		}
 		
 		assert sum != null;
 		
 		meanPoint = ConservationAccessory.multPointByScalar(sum, 1.0/(double) columnArr.length);
 		
 		for (int i = 0; i < columnArr.length; i++) {
 			
 			distance = distance + ConservationAccessory.pointDistance(meanPoint, points[i]);
 			
 		}
 		
 		double nonGapsFraction = (columnArr.length - acidsIntMap.get('-')) / columnArr.length;
 		
 		double result = nonGapsFraction * 1/columnArr.length * distance;
 		
 		return result;
 			
 		}
 		
 		
	double lancetScore() {
		
		double result = 0.0;
		
		Set<Character> keys = acidsIntMap.keySet();
		
		Iterator<Character> itr1 = keys.iterator();
		
		Iterator<Character> itr2 = keys.iterator();
		
		while (itr1.hasNext()) {
			
			char key1 = itr1.next();
			
			while (itr2.hasNext()) {
				
				char key2 = itr2.next();
				
				double blosum = ConservationAccessory.BlosumPair(key1, key2);
				
				if (blosum == 0) { blosum = 0.000000000000001;}
				
				result = result + ((acidsIntMap.get(key1)/columnArr.length * acidsIntMap.get(key1)/columnArr.length)/blosum);
			}
		}
	
	return result;
	
	}
	
	
}

 

    	
    	
 

