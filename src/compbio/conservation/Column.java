package compbio.conservation;
import java.util.*;

/**
 * This class provides all the basic calculations on column level.
 * Has all the methods for calculating conservation scores.
 * The conservation scorse are based on Valdar's paper
 * 
 * @author agolicz
 */

public class Column {

	/**
	 * Array consisting of character representation of amino acids in the column
	 */
	private final char[] columnArr;
	
	/**
	 * Maps a characters representing amino acids with their occurrence.
	 */
    private final Map<Character,Integer> acidsIntMap;

    /** 
     * Reference to the matrix object the column was taken from. Indices of the amino acids in the column correspond to the indices of the the sequences they belong to in this matrix.
     */
    private final AminoAcidMatrix matrix;

    /**
     * Constructor that allows manual creation of the column of the length of 10.
     * @param a
     * @param b
     * @param c
     * @param d
     * @param e
     * @param f
     * @param g
     * @param h
     * @param i
     * @param j
     * @param m
     */

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


    /**
     * Constructor. Creates a column object by copying the appropriate column form the matrix.
     * Takes a reference to the alignment matrix and the index of the corresponding column in the matrix. 
     * @param mmatrix
     * @param columnNr
     */
    public Column(AminoAcidMatrix matrix, int columnNr) {
    	 

	if (matrix == null) {
	    throw new IllegalArgumentException("Matrix must not be null");
	}
	
	assert columnNr >= 0 && columnNr <  matrix.numberOfColumns();
	
	this.matrix = matrix;

	String columnStr = new String(matrix.getColumn(columnNr));

	String upColumnString = columnStr.toUpperCase();

	columnArr = upColumnString.toCharArray();
	
	Map<Character, Integer> map = Alphabet.calculateOccurance(columnArr);
	
	acidsIntMap = Collections.unmodifiableMap(map);
	
	assert acidsIntMap != null && !acidsIntMap.isEmpty();
	
	}
    
    /**
     * Return array of characters representing amino acids in the column.
     */
    
    char[] getColumnAcids() {
    	
    	return this.columnArr;
    }
    
    
	/**
	 * Checks if all but one residues in the column are gaps.
	 * 
	 * @return true if all but one residues are gaps, false if not 
	 */
	
    boolean allButOneGaps() {
    	
    	if (acidsIntMap.containsKey('-') && acidsIntMap.get('-') == columnArr.length - 1) {
			
			return true;
		}
	
	else {
		
		return false;
		
	}
	
    }
    /**
     * Checks whether there is only one residue type in the column. Gap is not counted as a residue type.
     * 
     * @return returns true if there is only one residue type in the column, false if not
     */
    boolean oneResidueTypeNoGaps() {
    	
    	if (acidsIntMap.size() == 1 && acidsIntMap.containsKey('-') == false) {
    		
    		return true;
    		
    	}
    	
    	else {
    		
    		return false;
    		
    	}
    	
    }
    
    
    /** 
     * Checks whether column contains gaps.
     * 
     * @return true if column contains gaps
     */
    
    boolean containsGaps() {
    	
    	if (acidsIntMap.containsKey('-')) {
    		
    		return true;
    	}
    	
    	else {
    		
    		return false;
    		
    	}
    	
    }
 
    /**
     * Counts the number of different amino acids in the column.
     * Gap is counted as a 21 amino acid.
     * 
     * @return the number of different amino acids
     */
    
    int numberOfAcidsWithGap() {
    	
    	return acidsIntMap.size(); 
    	
    }
    
    /**
     * Counts the number of different amino acids in the column.
     * Gap not counted as a 21 amino acid.
     * 
     * @return the number of different amino acids
     */
    
    int numberOfAcidsNoGap() {
    
    if(this.containsGaps() == true) {
    	
    	return acidsIntMap.size() - 1;
    }
    
    else {
    	
    	return acidsIntMap.size();
    
    }
    
    }
    /** 
     * Calculates the number of times the most common amino acid in the column occurs.
     * Does not count gap as 21 amino acid.
     * 
     * @return times the most common amino acid occurs
     */
    
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
    /**
     * Gets length of the column.
     * 
     * @return length of the column
     */
    			
    int length() {
	 
	 return columnArr.length;
	 
    }
 
    /**
     * Calculates Kabat score for the column.
     *  
     * @return Kabat score
     */
    double kabatScore() {
	 
        assert this.mostCommonNumber() > 0 && this.mostCommonNumber() < columnArr.length + 1;	 
         		
		double result = this.length() * (double) this.numberOfAcidsNoGap()/ (double) this.mostCommonNumber(); 
	
		return result;
		
	}
 
    /**
     * Calculates Jores score for the column.
     * When calculating the number of distinct pair formed calculates a combination of two out of the aminoacid types
     * present in the column. Than adds the self pairs. However when there is only one of a type
     * in the column ne self pair is formed. 
     * 
     * @return Jores score 
     */
    double joresScore() {
	 
	 double result = 0.0;
	 
	 // special case #1 one residue type only
	 
	 if (this.oneResidueTypeNoGaps() == true) {
		 
		 result = 1.0;
		 
		 return result;
	 
	 }
	
	 // special case #2 all but one are gaps
	 
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

    /**
     * Calculates Schneider score.
     */
	double schneiderScore() {
		
		double normal = 1.0 / Math.log(20.0);
		
		double result = ShannonEnthropy.ShannonLn(acidsIntMap, columnArr.length) * normal;
		
		assert result >= 0 && result <= 1;
			
		return result;
	}
	/**
	 * Calculates Shenkin score.
	 * 
	 * @return Shenkin score
	 */

	double shenkinScore() {
		
		double result = Math.pow( 2.0, ShannonEnthropy.ShannonLog2(acidsIntMap, columnArr.length)) * 6.0;
			
		assert result >= 6 && result <= 120;
		
		return result;
	}

	/**
	 * Calculates Gerstein score.
	 * 
	 * @return Gerstein score
	 */
	
	double gersteinScore() {
		
		double result = ShannonEnthropy.ShannonLn(matrix.totalAcidsFrequency(), matrix.numberOfColumns() * matrix.numberOfRows()) - ShannonEnthropy.ShannonLn(acidsIntMap, columnArr.length);
		
		return result;
		
	}

 
	/**
	 * Calculates Taylor score. Returns the number of elements in the smallest set covering the whole
	 * column. Gaps sore is included in the largest set. Therefore a column containing gaps is automatically given the lowest possibel score.  
	 *
	 * @return Taylor score
	 */
	
	int SmallestTaylorSetGaps() {
	
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
 
	/**
	 * Does a very similar thing to SmallestTaylorSetGaps but does not take gaps into account at all.
	 *  
	 * @return Taylor score
	 */
	int SmallestTaylorSetNoGaps() {
	 
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
	/**
	 * Gives a score of one to ten based on whether all the amino acids though the column
	 * maintain or fail to maintain a certain trait.
	 *  
	 * @return
	 */
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
 		
 		assert result >= 0 && result < 11;
 		
 		return result;
 	}

 	/**
 	 * Calculates Karlin score.
 	 *  done exactly like in jon's code: gaps not counted, if all but 1 gaps conservation score 0
 	 *  although it looks to me that the formula suggests that if one doesn't count gaps
 	 *  than the score of all but one gaps would be undefined.
 	 *  
 	 * @return Karlin score
 	 */
 	
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
 	
 	/**
 	 * Calculates Armon score.
 	 * 
 	 * @return Armon score.
 	 */
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
 	
 	/** 
 	 * Calculates Thompson score.
 	 * Gaps are accounted for. 
 	 * 
 	 * @return
 	 */
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
 	double notLancetScore() {
		
		double result = 0.0;
		
		Set<Character> keys = acidsIntMap.keySet();
		
		Iterator<Character> itr1 = keys.iterator();
		
		while (itr1.hasNext()) {
			
			char key1 = itr1.next();
			
			Iterator<Character> itr2 = keys.iterator();
			
			while (itr2.hasNext()) {
				
				char key2 = itr2.next();
				
				double blosum = ConservationMatrices.BlosumPair(key1, key2);
				
				result = result + ((((double) acidsIntMap.get(key1)/ (double) columnArr.length * (double) acidsIntMap.get(key1)/ (double) columnArr.length)) * blosum);
			}
		}
	
	return result;
	
	}
	/**
	 * Calculates Mirny Score.
	 * 
	 * @return mirny score
	 */
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
	/**
	 * Calculates Williamson score.
	 * 
	 * 
	 * @return Williamson score.
	 */

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
			
			double piAve = (double) matrix.totalAcidsWillSets().get(setFreqKey);
			
			willSum = willSum + (pI * Math.log(pI/piAve));
			
		}
		
		return willSum;
		
	}
	
	/**
	 * Calculates Landgraf score.
	 * 
	 * @return landgraf score
	 */
	double landgrafScore() {
		
		double sum = 0.0;
		
		for (int i = 0; i < columnArr.length; i++) {
			
			for (int j = i + 1; j < columnArr.length; j++) {
				
				double disIJ = ConservationMatrices.dissimilarity(columnArr[i], columnArr[j]);
				
				double disJI = ConservationMatrices.dissimilarity(columnArr[j], columnArr[i]);
				
				sum = sum + (matrix.getVoronoiWeights(1000)[i] * disIJ) + (matrix.getVoronoiWeights(1000)[j] * disJI);
			}
			
		}
		
		double result = sum / columnArr.length;
		
		return result;
	}
	
	/** 
	 * Calculates sander score.
	 * 
	 * @return sander score
	 */ 
	
	double sanderScore() {
		
		double sum = 0.0;
		
		double moderator = 0.0;
		
		for (int i = 0; i < columnArr.length; i++) {
			
			for (int j = i + 1; j < columnArr.length; j++) {
				
				sum = sum + (1 - matrix.getPercentIdentity()[i][j]) * ConservationMatrices.pam250Pair(columnArr[i], columnArr[j]);
			}
			
		}
		
		for (int i = 0; i < columnArr.length; i++) {
			
			for (int j = i + 1; j < columnArr.length; j++) {
			
			moderator = moderator + (1 - matrix.getPercentIdentity()[i][j]);
			
			}
		}

		double result = sum * moderator;
		
		return result;
		
		}

	/**
	 * Calculates Valdar score.
	 * 
	 * @return Valdar score
	 */
	double valdarScore() {
		
		double sum = 0.0;
		
		double moderator = 0.0;
		
		for (int i = 0; i < columnArr.length; i++) {
			
			for (int j = i + 1; j < columnArr.length; j++) {
				
				sum = sum + matrix.vingronArgosWeights()[i] * matrix.vingronArgosWeights()[j] * ConservationMatrices.pet91Pair(columnArr[i], columnArr[j]);
			}
			
		}
		
		for (int i = 0; i < columnArr.length; i++) {
			
			for (int j = i + 1; j < columnArr.length; j++) {
			
				moderator = moderator + matrix.vingronArgosWeights()[i] * matrix.vingronArgosWeights()[j];
			}
			
		}
		
		double result = sum * moderator;
		
		return result;
		
	}
	
} 

    	
    	
 

