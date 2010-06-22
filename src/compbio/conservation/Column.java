package compbio.conservation;
import java.util.*;

// This class provides all the basic calculations on column level
// Has all the methods for calculating conservation scores





public class Column {

	private final char[] columnArr;
    private final Map<Character,Integer> acidsIntMap;
    // a reference to a matrix a column is based on is stored as instance variable 
    // not sure I like this idea, have to think about some better one, but don't have it yet
    private final AminoAcidMatrix matrix;

    // TO BE DELETED, THIS CONSTRUCTOR USED FOR TESTS ONLY

    public Column(char a, char b, char c, char d, char e, char f,
	  char g, char h, char i, char j, AminoAcidMatrix m) {

	columnArr = new char[10];

	columnArr[0] = a;
	columnArr[1] = b;
	columnArr[2] = c;
	columnArr[3] = d;
	columnArr[4] = e;
	columnArr[5] = f;
	columnArr[6] = g;
	columnArr[7] = h;
	columnArr[8] = i;
	columnArr[9] = j;


    acidsIntMap = Alphabet.calculateOccurance(columnArr);
    
    assert !acidsIntMap.isEmpty() && acidsIntMap != null;
    
    matrix = m;
    
    }

    // constructor, 

    public Column(AminoAcidMatrix m, int columnNr) {
    	
    // 

	if (m == null) {
	    throw new IllegalArgumentException("Matrix must not be null");
	}
	
	matrix  = m;
	
	assert columnNr >= 0 && columnNr <= m.numberOfColumns();

	String columnStr = new String(matrix.getColumn(columnNr));

	String upColumnString = columnStr.toUpperCase();

	columnArr = upColumnString.toCharArray();

// calls method to create a map
	
	//Alphabet alp = new Alphabet();
	
	Map<Character, Integer> map = Alphabet.calculateOccurance(columnArr);
	
	acidsIntMap = Collections.unmodifiableMap(map);
	
	assert acidsIntMap != null && !acidsIntMap.isEmpty();
	
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
    boolean oneResidueTypeNoGaps() {
    	
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
    // look out for this method if column empty return most common nr as 0
    
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
 
 double kabatScore() {
	 
        assert this.mostCommonNumber() > 0 && this.mostCommonNumber() < columnArr.length + 1;	 
         		
		double result = this.length() * this.numberOfAcidsNoGap()/ this.mostCommonNumber(); 
	
		return result;
		
	}
 
 double joresScore() {
	 
	 double result = 0.0;
	 
	 if (this.oneResidueTypeNoGaps() == true) {
		 
		 result = 1.0;
		 
		 return result;
	 
	 }
	 
	 if (this.allButOneGaps() == true) {
		 
		 return result;
		 
	 }
	 
	 int samePairs = 0;
	 
	 int differentPairs = 0;
	 
	 Map<Character,Integer> acidsIntMapCopy = new HashMap<Character,Integer>(acidsIntMap);
	 
	 Set<Character> keys = acidsIntMapCopy.keySet();
	 
	 keys.remove('-');
	 
	 int types = keys.size();
	 
	 differentPairs = types * (types - 1) / 2;
	 
	 Iterator<Integer> itr = acidsIntMapCopy.values().iterator();
	 
	 while (itr.hasNext()) {
		 
		 if (itr.next() > 1) {
			 
			 samePairs++;
			 
		 }
		 
	 }
	 
	 int totalPairs = samePairs + differentPairs;
	 
	 Iterator<Character> itr2 = acidsIntMapCopy.keySet().iterator();
	 
	 int max1 = 0;
	 
	 Character maxKey = null;
	 
	 while(itr2.hasNext()) {
		 
		 Character key = itr2.next();
		 
		 if (acidsIntMapCopy.get(key) > max1) {
			 
			 maxKey = key;
			 
			 max1 = acidsIntMapCopy.get(key);
			 
		 }
		 
	 }
	 
	 acidsIntMapCopy.remove(maxKey);
	 
	 Iterator<Integer> itr3 = acidsIntMapCopy.values().iterator();
	 
	 int max2 = 0;
	 
	 while(itr3.hasNext()) {
		 
		 int value = itr3.next();
		 
		 if(value > max2) {
			 
			 max2 = value;
			 
		 }
		 
	 }
	 
	 int mostFreqNr = 0;
	 
	 if (max2 == 0) {
		 
		 mostFreqNr = (max1) * (max1 - 1) / 2;
		 
	 }
	 
	 else {
		 
		 if (max1 == max2) {
			 
			 mostFreqNr = max1 * max2;
			 
		 }
		 
		 else {
			 
			 int same = (max1) * (max1 - 1) / 2;
			 
			 int diff = max1 * max2;
			 
			 if (same > diff) {
				 
				 mostFreqNr = same;
				 
			 }
				 
			 else {
					 
					 mostFreqNr = diff;
				 }
			 
		 	}
		 
	 	}
		 
 	result = ((double) totalPairs / (double) mostFreqNr) * ((columnArr.length) * (columnArr.length -1 ) / 2); 
 
 	return result;
 	
 	}
	 
	 
	// Symbol Enthropy Scores

	double schneiderScore() {
		
		double normal = 1.0 / 20.0;
		
		double result = ShannonEnthropy.ShannonLn(acidsIntMap, columnArr.length) * normal;
			
		return result;
	}


	double shenkinScore() {
		
		double result = Math.pow( 2.0, ShannonEnthropy.ShannonLog2(acidsIntMap, columnArr.length)) * 6.0;
			
		return result;
	}

	double gersteinScore() {
		
		double result = ShannonEnthropy.ShannonLn(matrix.totalAcidsFrequency(), matrix.numberOfColumns() * matrix.numberOfRows()) - ShannonEnthropy.ShannonLn(acidsIntMap, columnArr.length);
		
		return result;
		
	}

 
 // this method gets the number of elements in the smallest set  
 // among the Taylor sets that covers all the aa in the column
 // the set that includes a gap symbol, is the biggest set of all(all aa + gap symbol) 
 int SmallestTaylorSetGaps() {
 
	//if (setMap == null) {
		
		//throw new IllegalArgumentException("setMap must not be null");
	//}
	
	Map<String, HashSet<Character>> setMap = ConservationSets.taylorSets();
	
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
	
	assert smallestSetSize > 0 && smallestSetSize < 25;
	
	return smallestSetSize;
	
	}
 
 // does a very similar thing as a previous one, but does not take gaps into account
 int SmallestTaylorSetNoGaps() {
	 
		//if (setMap == null) {
			
			//throw new IllegalArgumentException("setMap must not be null");
		//}
		//Alphabet alp = new Alphabet();
	 
	    Map<String, HashSet<Character>> setMap = ConservationSets.taylorSets();
		
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
		
		assert smallestSetSize > 0 && smallestSetSize < 25;
		
		return smallestSetSize;
		
		}
	// gives score 0-10 based on the fact that all the characters belong to certain sets
 	// there are 20 sets(they consist of 10 sets based on characteristics and 10 complements
 	int zvelibilScore(){
 		
 		int result = 0;
 		
 		Map<String, HashSet<Character>> setMap = ConservationSets.zvelibilSets();
 		
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
 		
 		for( int a = 0; a < columnArr.length; a++) {
 			
 			if(columnArr[a] != '-') {
 			
 				for(int b = a + 1; b < columnArr.length; b++) {
 					
 					if(columnArr[b] != '-') {
 						
 						double pairScore =  ConservationMatrices.BlosumPair(columnArr[a], columnArr[b]);
 						
 						double aSelf =  ConservationMatrices.BlosumPair(columnArr[a], columnArr[a]);
 						
 						assert aSelf > 0;
 						
 						double bSelf = ConservationMatrices.BlosumPair(columnArr[b], columnArr[b]);
 						
 						assert bSelf > 0;
 						
 						blosumSum = blosumSum + ((pairScore) / (Math.sqrt(aSelf * bSelf)));
 					
 					}
 					
 				}
 		
 			}
 			
 		}
 		
 		
 		double finalSum = blosumSum * (2.0 / ( columnArr.length * (columnArr.length - 1)));
 		
 	    assert finalSum >= -1 && finalSum <= 1;
 		
 		return finalSum;
 	}
 	
 // creates an array containing all the amino acids and gaps present in the column
 // iterates through that array twice(nested loops), finds all the possible pairs 
 // that can be formed by aa present
 // gap is considered the 21 aminoacid
 	double armonScore() {
 		
 		double scoreSum = 0;
 		
 		int arrayLength = acidsIntMap.keySet().size();
 		
 		Character[] acidsPresent = new Character[arrayLength];
 		
 		int arrayIndex = 0;
 		
 		Set<Character> keys = acidsIntMap.keySet();
 		
 		Iterator<Character> itr = keys.iterator();
 		
 		while(itr.hasNext()) {
 			
 			acidsPresent[arrayIndex] = itr.next();
 			
 			arrayIndex++;
 			
 		}
 		
 		for ( int a = 0; a < acidsPresent.length; a++ ) {
 			
 			char charA = acidsPresent[a];
 			
 			for (int b = a + 1; b < acidsPresent.length; b++) {
 				
 				char charB = acidsPresent[b];
 				
 				scoreSum = scoreSum + ConservationMatrices.miyataArmonPair(charA, charB);
 				
 			}
 		}
 		
 		return scoreSum;
 		
 	}
// amino acids are viewed as points in k dimensional space
// an average point is calculated
// score is the distance between the av point and the actual point
 	
 	double thompsonScore(){
 		
 		double[] sum = null;
 		
 		double[] meanPoint = null;
 		
 		double distance = 0.0;
 		
 		double nonGapsFraction = 0.0;
 		
 		char[] alp = Alphabet.alphabetArray();
 		
 		assert alp != null && alp.length != 0;
 		
 		double[][] points = new double[columnArr.length][alp.length];
 		
 		for (int i = 0; i < columnArr.length; i++) {
 			
 			for (int j = 0; j < alp.length; j++) {
 				
 				points[i][j] = ConservationMatrices.BlosumPair(columnArr[i], alp[j]);
 				
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
 			
 			distance = distance + ConservationAccessory.pointDistance(points[i], meanPoint);
 			
 		}
 		
 		if (acidsIntMap.keySet().contains('-')) {
 			
 			nonGapsFraction = (double) (columnArr.length - acidsIntMap.get('-')) / (double) columnArr.length;
 			
 		}
 		
 		else {
 			
 			nonGapsFraction = 1.0;
 			
 		}
 		double result = nonGapsFraction * 1.0/columnArr.length * distance;
 		
 		return result;
 			
 		}
 		
 	// causes some math problem because denominator can be 0, that's a formula flaw
 	// nothing can be done about it
 	double lancetScore() {
		
		double result = 0.0;
		
		Set<Character> keys = acidsIntMap.keySet();
		
		Iterator<Character> itr1 = keys.iterator();
		
		while (itr1.hasNext()) {
			
			char key1 = itr1.next();
			
			Iterator<Character> itr2 = keys.iterator();
			
			while (itr2.hasNext()) {
				
				char key2 = itr2.next();
				
				double blosum = ConservationMatrices.BlosumPair(key1, key2);
				
				if (blosum == 0.0) { blosum = 0.00000000000000000000000001;}
				
				result = result + ((((double) acidsIntMap.get(key1)/ (double) columnArr.length * (double) acidsIntMap.get(key1)/ (double) columnArr.length)) / blosum);
			}
		}
	
	return result;
	
	}
	
 	//calculates shannon enthropy but based on the sets, has to calculate number of amino acids tat belong 
 	//a particular set, stores them in a hashmap
 	//reads in a hashmap with the sets needed and creates a hashmap with set names as keys and number of aa belonging to set as value
 	
	double mirnyScore() {
		
		double mirnySum = 0.0;
		
		Map<String, HashSet<Character>> mirnySets= ConservationSets.mirnySets();
		
		assert mirnySets != null && !mirnySets.isEmpty();
		
		Set<String> mirnyKeys = mirnySets.keySet();
		
		assert !mirnyKeys.isEmpty();
		
		Iterator<String> mirnyKeysItr = mirnyKeys.iterator();
		
		Set<Character> acInKeys = acidsIntMap.keySet(); 
		
		assert !acInKeys.isEmpty();
		
		Map<String,Integer> setsFreq = new HashMap<String,Integer>();  
		
		while (mirnyKeysItr.hasNext()) {
			
			String mirnyKey = mirnyKeysItr.next();
			
			Iterator<Character> acInKeysItr = acInKeys.iterator();
			
			while (acInKeysItr.hasNext()) {
				
				Character acInKey = acInKeysItr.next();
				
				if (mirnySets.get(mirnyKey).contains(acInKey)) {
					
					Integer count = setsFreq.get(mirnyKey);
					
					if (count == null) {
						
						setsFreq.put(mirnyKey, acidsIntMap.get(acInKey));
					}
					
					else {
						
						setsFreq.put(mirnyKey, count + acidsIntMap.get(acInKey));
						
					}
				}
			}
		}
		
		// this assertion will not work if we feed an empty column
		
		assert !setsFreq.isEmpty();
		
		Set<String> setsFreqKeys = setsFreq.keySet();
		
		Iterator<String> setsFreqKeysItr = setsFreqKeys.iterator();
		
		while(setsFreqKeysItr.hasNext()) {
			
			String setFreqKey = setsFreqKeysItr.next();
			
			double pI = (double) setsFreq.get(setFreqKey) / (double) columnArr.length; 
			
			mirnySum = mirnySum + (pI * Math.log(pI));
			
		}
		
		return mirnySum;
		
	}
	

	double williamsonScore() {
		  
		double willSum = 0.0;
		
		Map<String, HashSet<Character>> willSets= ConservationSets.williamsonSets();
		
		assert willSets != null && ! willSets.isEmpty();
		
		Set<String> willKeys = willSets.keySet();
		
		assert !willKeys.isEmpty();
		
		Iterator<String> willKeysItr = willKeys.iterator();
		
		Set<Character> acInKeys = acidsIntMap.keySet(); 
		
		assert !acInKeys.isEmpty();
		
		Map<String,Integer> setsFreq = new HashMap<String,Integer>();  
		
		while (willKeysItr.hasNext()) {
			
			String willKey = willKeysItr.next();
			
			Iterator<Character> acInKeysItr = acInKeys.iterator();
			
			while (acInKeysItr.hasNext()) {
				
				Character acInKey = acInKeysItr.next();
				
				if (willSets.get(willKey).contains(acInKey)) {
					
					Integer count = setsFreq.get(willKey);
					
					if (count == null) {
						
						setsFreq.put(willKey, acidsIntMap.get(acInKey));
					}
					
					else {
						
						setsFreq.put(willKey, count + acidsIntMap.get(acInKey));
						
					}
				}
			}
		}
		
		// this assertion will not work if we feed an empty column
		
		assert !setsFreq.isEmpty();
		
		Set<String> setsFreqKeys = setsFreq.keySet();
		
		Iterator<String> setsFreqKeysItr = setsFreqKeys.iterator();
		
		while(setsFreqKeysItr.hasNext()) {
			
			String setFreqKey = setsFreqKeysItr.next();
			
			// FIXME Pi in the logarithm needs to be divided by average pi, do it once u get classes to be nested
			
			double pI = (double) setsFreq.get(setFreqKey) / (double) columnArr.length; 
			
			double piAve = (double) matrix.TotalAcidsWillSets().get(setFreqKey);
			
			willSum = willSum + (pI * Math.log(pI/piAve));
			
		}
		
		return willSum;
		
	}
	
	double landgrafScore() {
		
		double sum = 0.0;
		
		for (int i = 0; i < columnArr.length; i++) {
			
			for (int j = i + 1; j < columnArr.length; j++) {
				
				double disIJ = ConservationMatrices.dissimilarity(columnArr[i], columnArr[j]);
				
				double disJI = ConservationMatrices.dissimilarity(columnArr[j], columnArr[i]);
				
				sum = sum + (ConservationAccessory.voronoiWeights(matrix, 1000)[i] * disIJ) + (ConservationAccessory.voronoiWeights(matrix, 1000)[j] * disJI);
			}
			
		}
		
		double result = sum / columnArr.length;
		
		return result;
	}
	
	double sanderScore() {
		
		double sum = 0.0;
		
		double moderator = 0.0;
		
		for (int i = 0; i < columnArr.length; i++) {
			
			for (int j = i + 1; j < columnArr.length; j++) {
				
				sum = sum + (1 - ConservationAccessory.percentIdentity(matrix.getRow(i), matrix.getRow(j)) * ConservationMatrices.pam250Pair(columnArr[i], columnArr[j]));
			}
			
		}
		
		for (int i = 0; i < columnArr.length; i++) {
			
			for (int j = i + 1; j < columnArr.length; j++) {
			
			moderator = moderator + (1 - ConservationAccessory.percentIdentity(matrix.getRow(i), matrix.getRow(j)));
			
			}
		}

		double result = sum * moderator;
		
		return result;
		
		}

	
	double valdarScore() {
		
		double sum = 0.0;
		
		double moderator = 0.0;
		
		for (int i = 0; i < columnArr.length; i++) {
			
			for (int j = i + 1; j < columnArr.length; j++) {
				
				sum = sum + ConservationAccessory.weightOfSequenceVingronArgos(i, matrix) * ConservationAccessory.weightOfSequenceVingronArgos(i, matrix) * ConservationMatrices.pet91Pair(columnArr[i], columnArr[j]);
			}
			
		}
		
		for (int i = 0; i < columnArr.length; i++) {
			
			for (int j = i + 1; j < columnArr.length; j++) {
			
				moderator = moderator + ConservationAccessory.weightOfSequenceVingronArgos(i, matrix) * ConservationAccessory.weightOfSequenceVingronArgos(i, matrix);
			}
			
		}
		
		double result = sum * moderator;
		
		return result;
		
	}
	
} 

    	
    	
 

